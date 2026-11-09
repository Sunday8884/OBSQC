#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
ppp.py — 精密单点定位（PPP）主引擎。

模型概述
--------
- 观测值：双频无电离层组合伪距 P_IF 与载波 L_IF（GPS/BDS/Galileo 多系统）；
- 卫星星历：SP3 精密星历 + 精密钟差（可选 RINEX CLK，缺省用 SP3 内嵌钟差）；
- 误差改正：卫星钟差、对流层（Collins 先验 + 湿延迟随机游走估计）、
  固体潮（度-2 主项）、相位缠绕、Sagnac（geodist 内置）；
- 参数估计：逐历元扩展卡尔曼滤波（EKF），
  状态 = [坐标(3), 各系统接收机钟差(白噪声), 天顶湿延迟(随机游走), 各弧段浮点模糊度(常量)]；
- 周跳探测：MW + GF 组合（ppp/slip.py），周跳或长中断后模糊度重置；
- 粗差处理：逐量测标准化新息 > outlier_nsigma 剔除。

接口
----
ppp(obsHead, obsData, sp3Data, clkData, bandChoose, startDatetime, endDatetime,
    config=None, posFile=None, plot_show=False, plot_save=None,
    snx_file=None, station=None)
    -> (posData, XMEAN, YMEAN, ZMEAN)   # 与 spos.spp() 接口保持一致
"""
import datetime
import os
import numpy as np

from com.gnssParameter import CLIGHT, bandDictGetFreqPos
from com.xyz2blh import xyz2blh, geodist, bl2enuRot
from com.sat2siteAngle import getEle
from com.gnssTime import datetime2doys
from com.pub import printPanda
from com.readNav import readNav

from spos.trop import tropColins
from plot.plotppp import plot_ppp
from spos.spos import read_snx_coord

from ppp.pppConfig import DEFAULT_CONFIG, load_config
from ppp.satposSp3 import satPosSp3
from ppp.slip import detectCycleSlip
from ppp.qcweight import qc_grey_weights, qc_abs_weights, abs_calibrate
from ppp.corrections import (sunPosECEF, solidEarthTide, phaseWindUp, tropMapping,
                             saastamoinenZenith, nmfMapf, poleTide)
from com.readAtx import readAntex, satAntCorr, recAntCorr, bandToAtxTag
from com.readBlq import readBlq, oceanLoadingDisp
from com.readDcb import readDcb, applyDcbSat, dcb_stat
from com.readErp import readErp, getErp
from ppp.pppfilter import KalmanFilter
from ppp.pppresult import writePosData


# ======================================================== 辅助：IGGⅢ 等价权
def _iggiFactor(vhist, k0, k1, wmin):
    """
    IGGⅢ 三段式等价权（抗差估计，杨元喜），自适应标准化版。
    vhist : 该星该类型最近若干历元的标准化新息 v = innov/√R（含待评估的
            最后一历元；R 为该行量测噪声方差）。
    两级归一化：
      ① v 用量测噪声 σ_row 归一（消除码/载波量纲差）；
      ② σ̂ = max(1.4826·MAD(|v|历史窗, 不含 v_last), 1.0)——该星历史鲁棒
        噪声水平，floor=1.0 表示名义噪声。语义为"当前新息相对自身历史
        是否突跳"：孤立粗差持续触发降权，星噪声水平整体缓变时 σ̂ 随窗
        滚动自适应（实为 n 倍名义噪声 → σ̂≈n，不会持续误压权）。
    u = |v_last|/σ̂ 后套 IGGⅢ：
      |u| ≤ k0 ：w = 1（正常，不降权）
      k0<|u|≤k1：w = (k0/|u|)·((k1−|u|)/(k1−k0))²（平滑降权）
      |u| > k1 ：w = wmin（软化剔除，避免删行破坏批结构）
    """
    v = vhist[-1] if vhist else None
    if v is None or not np.isfinite(v):
        return 1.0
    arr = np.abs(np.asarray(vhist[:-1], dtype=float))
    mad = float(np.median(arr)) if arr.size else 0.0
    sigma = max(1.4826 * mad, 1.0)
    u = abs(v) / sigma
    if u <= k0:
        return 1.0
    if u >= k1:
        return wmin
    return max((k0 / u) * ((k1 - u) / (k1 - k0)) ** 2, wmin)


# ======================================================== 辅助：码定位初值
def _codePosition(satObs, satClkInEpoch, gSysList, x0=None, maxIter=12, rejThresh=100.0):
    """
    仅用无电离层组合伪距做最小二乘，给出 PPP 初始坐标（带粗差剔除）。
    satObs : {prn: {'P': P_IF, 'sys': 'G', 'rs': np.array(6)}}
    rejThresh : 残差（去中位数后）剔除阈值（米），逐迭代剔除最差卫星。
    返回 (x, y, z) 或 None（卫星不足/不收敛）。
    """
    gSysIndex = {s: i for i, s in enumerate(gSysList)}
    nSys = len(gSysList)
    X = np.zeros((3 + nSys, 1))
    if x0 is not None:
        X[:3, 0] = x0

    active = dict(satObs)
    for _ in range(maxIter):
        nUse = 0
        G = np.zeros((len(active), 3 + nSys))
        y = np.zeros((len(active), 1))
        for prn, so in active.items():
            rs = so['rs'][:3]
            rr = X[:3, 0]
            rho, e = geodist(rs, rr)
            satClk = satClkInEpoch.get(prn)
            if satClk is None or not np.isfinite(rho) or not np.isfinite(so['P']):
                continue
            G[nUse, :3] = e
            G[nUse, 3 + gSysIndex[so['sys']]] = 1.0
            y[nUse, 0] = so['P'] - (rho - satClk * CLIGHT + X[3 + gSysIndex[so['sys']], 0])
            nUse += 1
        if nUse < 4:
            return None
        G = G[:nUse]
        y = y[:nUse]
        try:
            dx = np.linalg.pinv(G.T @ G) @ G.T @ y
        except Exception:
            return None
        X = X + dx

        # --- 残差检查与粗差剔除（去中位数后超出阈值者，逐迭代剔除一个）---
        if not np.isfinite(X[:3]).all() or np.linalg.norm(dx[:3]) > 1e7:
            return None
        resid = []
        for prn, so in active.items():
            rho, _ = geodist(so['rs'][:3], X[:3, 0])
            satClk = satClkInEpoch.get(prn)
            if satClk is None or not np.isfinite(rho):
                continue
            resid.append((prn, so['P'] - (rho - satClk * CLIGHT + X[3 + gSysIndex[so['sys']], 0])))
        if len(resid) < 4:
            return None
        med = float(np.median([r for _, r in resid]))
        bad = [(abs(r - med), prn) for prn, r in resid if abs(r - med) > rejThresh]
        if bad:
            bad.sort(reverse=True)
            del active[bad[0][1]]
            continue
        if np.max(np.abs(dx[:3])) < 1e-4:
            break

    if len(active) < 4 or not np.isfinite(X[:3]).all():
        return None
    # 顺便返回各系统接收机钟差（米）——用于滤波钟差状态初始化
    # （RTKLIB udclk_ppp 同源：initx(CLIGHT*dtr, VAR_CLK)，避免首历元
    #   钟差 0 初值与真实钟差相差数十公里导致载波新息超限被拒）
    return (float(X[0, 0]), float(X[1, 0]), float(X[2, 0]),
            [float(X[3 + i, 0]) for i in range(nSys)])


# ======================================================== 主引擎
def ppp(obsHead, obsData, sp3Data, clkData, bandChoose, startDatetime, endDatetime,
        config=None, posFile=None, plot_show=False, plot_save=None,
        snx_file=None, station=None):
    """
    精密单点定位主函数。参数与返回值见模块 docstring。
    """
    start_time = datetime.datetime.now()

    # ---------------- 配置 ----------------
    cfg = dict(DEFAULT_CONFIG)
    if isinstance(config, dict):
        cfg.update(config)
    elif isinstance(config, str):
        cfg.update(load_config(config))

    cutoff = np.radians(float(cfg['elev_cutoff_deg']))
    mode = str(cfg.get('mode', 'static')).lower()
    gSysList = list(bandChoose.keys())
    gSysIndex = {s: i for i, s in enumerate(gSysList)}
    nSys = len(gSysList)

    # 对流层梯度估计（MacMillan 模型，RTKLIB TROPOPT_ESTG 同源）；
    # 梯度作用于湿延迟，必须同时估计 ZWD
    estZwd = bool(cfg.get('est_zwd', True))
    estGrad = bool(cfg.get('est_trop_grad', False))
    if estGrad and not estZwd:
        print('警告：est_trop_grad 需要开启 est_zwd，对流层梯度估计未启用。')
        estGrad = False
    # 对流层模型：saas_nmf（RTKLIB TROPOPT_EST 同源：ZTD 总天顶延迟状态 +
    # Saastamoinen 干延迟先验 + NMF 干/湿分离映射）| collins（旧：Collins 全延迟先验 + 残余 ZWD）
    tropSaas = str(cfg.get('trop_prior', 'collins')).lower() == 'saas_nmf'

    # 波段频率/波长
    bandChoose = bandDictGetFreqPos(bandChoose)

    # ---------------- 精密产品（ANTEX / BLQ / DCB / ERP）----------------
    atxData = None
    if cfg.get('atx_file'):
        atxData = readAntex(cfg['atx_file'])
        cfg['sp3_velocity'] = True   # 卫星姿态（沿轨单位矢量）需要速度
        print('ANTEX：卫星 %d 条记录，接收机天线 %d 类 -> %s'
              % (atxData['n_sat'], atxData['n_rec'], cfg['atx_file']))
    blqData = None
    if cfg.get('blq_file'):
        blqAll = readBlq(cfg['blq_file'])
        blqStation = None
        if station and station.upper() in blqAll:
            blqStation = station.upper()
        else:
            # IGS 9 字符站名约定：MARKER NAME 前 4 字符为 IGS 站名
            # （如 WUH200CHN -> WUH2），用于匹配 BLQ 测站
            for cand in sorted(blqAll.keys()):
                if station and station.upper()[:4] == cand[:4]:
                    blqStation = cand
                    break
        if blqStation is None:
            print('警告：BLQ 文件中无测站 %s，海潮载荷改正未启用。' % (station or '?'))
        else:
            station = blqStation
            blqData = blqAll
            print('BLQ：测站 %s 海潮载荷已启用 -> %s' % (station, cfg['blq_file']))
    dcbData = None
    if cfg.get('dcb_file'):
        dcbData = readDcb(cfg['dcb_file'])
        nDcb, dcbType = dcb_stat(dcbData)
        print('DCB：条目 %d（%s）-> %s' % (nDcb, dcbType, cfg['dcb_file']))
    erpData = None
    if cfg.get('erp_file'):
        erpData = readErp(cfg['erp_file'])
        print('ERP：%d 天（MJD %.1f ~ %.1f）-> %s'
              % (len(erpData['mjd']), erpData['mjd'][0], erpData['mjd'][-1],
                 cfg['erp_file']))

    # 每系统两个波段的 ANTEX 频率标签（卫星/接收机 PCO·PCV 用）
    atxTags = {}
    if atxData is not None:
        for s, bands in bandChoose.items():
            bl = list(bands.keys())[:2]
            if len(bl) < 2:
                continue
            tag1 = bandToAtxTag(s, bl[0])
            tag2 = bandToAtxTag(s, bl[1])
            if tag1 and tag2:
                atxTags[s] = (bl[0], tag1, bl[1], tag2)

    recAntType = cfg.get('rec_antenna') or obsHead.get('Antenna Type') or None
    antDel = obsHead.get('Antenna Delta') or None
    if antDel is not None and antDel[0] is None:
        antDel = None   # 头文件未给出 DELTA H/E/N
    if atxData is not None and recAntType is None:
        print('警告：RINEX 头未含天线类型且未配置 rec_antenna，接收机 PCO/PCV 未启用。')

    # ---------------- 周跳探测 ----------------
    slipMap = detectCycleSlip(obsData, bandChoose, startDatetime, endDatetime, cfg)

    # ---------------- 状态初始化 ----------------
    sigmaPosInit = float(cfg['sigma_pos_init'])
    sigmaClk = float(cfg['sigma_clk'])
    sigmaZwdInit = float(cfg['sigma_zwd_init'])
    sigmaAmb = float(cfg['sigma_amb'])

    # 初始坐标：RINEX 头近似坐标，无效则用首历元码定位
    x0 = None
    approx = obsHead.get('Approx Position')
    if approx is not None and approx[0] is not None and abs(approx[0]) > 1e5:
        x0 = [float(approx[0]), float(approx[1]), float(approx[2])]

    kf = None
    tracks = {}   # prn -> {'idx': 状态索引, 'lastEpoch': datetime, 'windup': 累计缠绕(周)}

    # 状态布局：[pos(3), clk(nSys), d(nSys), zwd, gradN, gradE, amb...]
    # d = 码-相钟差差（载波量测含 d、码不含）：接收机钟漂使码观测漂移而
    # 载波不同步，钟差状态由码驱动漂移后，载波新息含全星一致共模 ±数百
    # 米；d 状态（白噪声每历元重置）正规吸收该共模，EKF 协方差自洽
    # （等效 RTKLIB 由白噪声钟差吸收的行为，且不污染 amb/ZWD 协方差）
    dIdx = 3 + nSys                    # 码-相钟差差状态起始索引
    zwdIdx = dIdx + nSys               # ZWD 状态固定索引
    gradNIdx = zwdIdx + 1 if estGrad else None   # 梯度北分量（MacMillan Gn）
    gradEIdx = zwdIdx + 2 if estGrad else None   # 梯度东分量（MacMillan Ge）
    nBase = zwdIdx + (1 if estZwd else 0) + (2 if estGrad else 0)

    # ---------------- 结果 ----------------
    posData = {}
    Xs, Ys, Zs = [], [], []
    dbgRes = None   # DEBUG: 残差转储缓冲（环境变量 PPP_DEBUG_RES 指定输出路径时启用）
    if os.environ.get('PPP_DEBUG_RES'):
        dbgRes = []
    qcState = {}    # QC 定权滑窗状态：{prn: {gf, gfep, l1m, mp, cmc, p, ep}}
    absHist = {}    # qc_method='abs' 跨历元指标历史：{metric: [值...]}（标定缓冲）
    absN = 0        # abs 热身历元计数
    resState = {}   # IGGⅢ 抗差状态：{prn: {'P': [v...], 'L': [v...], 'ep'}}
                    # v = 标准化新息 innov/√R（上一历元，RAE 单历元延迟反馈）
    qcwDbg = None   # DEBUG: QC 权值转储（环境变量 PPP_DEBUG_QCW）
    if os.environ.get('PPP_DEBUG_QCW'):
        qcwDbg = []

    # ---------------- 历元循环 ----------------
    for epoch in sorted(obsData.keys()):
        if startDatetime is not None and epoch < startDatetime:
            continue
        if endDatetime is not None and epoch > endDatetime:
            continue
        sodDbg = (epoch - epoch.replace(hour=0, minute=0, second=0,
                                        microsecond=0)).total_seconds()

        # --- 卫星位置/钟差（精密星历）---
        satPosInEpoch, satClkInEpoch, satVarInEpoch = satPosSp3(
            bandChoose, obsData, sp3Data, clkData, epoch, cfg)
        if len(satPosInEpoch) < 4:
            continue

        # --- 组装本历元观测 ---
        exSats = set(s for s in (cfg.get('exclude_sats') or []) if s)
        satObs = {}
        for prn, rs in satPosInEpoch.items():
            if prn in exSats:
                continue
            gnssSys = prn[0]
            if gnssSys not in bandChoose:
                continue
            bands = list(bandChoose[gnssSys].keys())[:2]
            if len(bands) < 2:
                continue
            rec = obsData[epoch].get(prn)
            if not rec:
                continue
            f1 = bandChoose[gnssSys][bands[0]]['freq']
            f2 = bandChoose[gnssSys][bands[1]]['freq']
            lam1 = bandChoose[gnssSys][bands[0]]['lambda']
            lam2 = bandChoose[gnssSys][bands[1]]['lambda']
            if not f1 or not f2 or abs(f1 - f2) < 1.0:
                continue

            P1 = rec.get(bands[0])
            P2 = rec.get(bands[1])
            if P1 is None or P2 is None:
                continue
            # DCB：观测码归算到精密钟差码基准（P1/P2、B1I/B3I、E1/E5a 等）
            dcb1 = applyDcbSat(dcbData, prn, bands[0])
            dcb2 = applyDcbSat(dcbData, prn, bands[1])
            if dcb1:
                P1 += dcb1
            if dcb2:
                P2 += dcb2
            PIF = (f1 ** 2 * P1 - f2 ** 2 * P2) / (f1 ** 2 - f2 ** 2)

            # 载波（周 → 无电离层组合，米）；LLI 非零按周跳处理（由 slipMap/弧段重置负责）
            L1 = rec.get('L' + bands[0][1:])
            L2 = rec.get('L' + bands[1][1:])
            LIF = None
            if L1 is not None and L2 is not None:
                lli1 = rec.get('L' + bands[0][1:] + 'LLI', 0)
                lli2 = rec.get('L' + bands[1][1:] + 'LLI', 0)
                if lli1 == 0 and lli2 == 0:
                    LIF = (f1 ** 2 * lam1 * L1 - f2 ** 2 * lam2 * L2) / (f1 ** 2 - f2 ** 2)

            lamIF = CLIGHT / (f1 + f2)
            varSat = float(satVarInEpoch.get(prn, 0.0)) \
                if cfg.get('use_sat_var', True) else 0.0

            # QC 指标原值（qcweight 灰色关联定权用，差分/窗口量在装配后
            # 统一计算，状态在周跳/中断后重置）：
            #   MP 组合：P1 − (f1²+f2²)/(f1²−f2²)·λ1L1 + 2f2²/(f1²−f2²)·λ2L2
            #   CMC    ：P1 − λ1L1
            #   GF     ：λ1L1 − λ2L2（PN/IOD 用）
            #   L1 米  ：λ1L1（HOD 用）
            mpRaw = cmcRaw = gfRaw = l1mRaw = None
            if L1 is not None and L2 is not None:
                mpRaw = P1 - (f1 ** 2 + f2 ** 2) / (f1 ** 2 - f2 ** 2) * lam1 * L1 \
                    + 2 * f2 ** 2 / (f1 ** 2 - f2 ** 2) * lam2 * L2
                cmcRaw = P1 - lam1 * L1
                gfRaw = lam1 * L1 - lam2 * L2
                l1mRaw = lam1 * L1
            s1 = rec.get('S' + bands[0][1:])
            s2 = rec.get('S' + bands[1][1:])
            snrVals = [float(v) for v in (s1, s2) if v is not None]
            snr = float(np.mean(snrVals)) if snrVals else None

            satObs[prn] = {'P': PIF, 'L': LIF, 'lamIF': lamIF,
                           'sys': gnssSys, 'rs': rs, 'varSat': varSat,
                           'mp': mpRaw, 'cmc': cmcRaw, 'gf': gfRaw,
                           'l1m': l1mRaw, 'snr': snr}

        # 本历元最大 SNR（snr_weight 相对权参考，缺失时退化高度角模型）
        epochMaxSnr = max((float(so['snr']) for so in satObs.values()
                           if so.get('snr')), default=48.0)

        # --- QC 灰色关联定权（每历元跨可见卫星，替换高度角随机模型）---
        # 指标（qc/ 模块全部逐星动态指标，周跳/中断重置状态）：
        #   mpd=|ΔMP|  cmcd=|ΔCMC|  snr=SNR均值
        #   pn=相位噪声（GF 三阶差分/√40）
        #   hod=高阶差分噪声（L1 四阶差分/√70）
        #   iod=电离层变化率 |ΔGF/Δt|
        # min-max 正向化 → 熵权 → 灰色关联系数(ρ) → 加权关联度 p_i → w=p_i
        if cfg.get('qc_weight', False) and satObs:
            satInd = {}
            newQC = {}
            for prn, so in satObs.items():
                st = qcState.get(prn)
                slipped = epoch in slipMap.get(prn, ())
                if st is None or slipped \
                        or (epoch - st['ep']).total_seconds() > 90.0:
                    st = {'gf': [], 'gfep': None, 'l1m': [],
                          'mp': None, 'cmc': None, 'p': []}
                ind = {}
                mp, cmc, gf, l1m = so.get('mp'), so.get('cmc'), \
                    so.get('gf'), so.get('l1m')
                # 码类差分（消除模糊度常数）
                if mp is not None and st['mp'] is not None:
                    ind['mpd'] = abs(mp - st['mp'])
                if cmc is not None and st['cmc'] is not None:
                    ind['cmcd'] = abs(cmc - st['cmc'])
                # IOD：GF 电离层变化率 (m/s)
                if gf is not None and st['gf'] and st['gfep'] is not None:
                    ind['iod'] = abs(
                        (gf - st['gf'][-1]) /
                        (epoch - st['gfep']).total_seconds())
                # PN：GF 三阶差分 / √40（4 点滑窗）
                if gf is not None and len(st['gf']) >= 3:
                    seq = st['gf'][-3:] + [gf]
                    d = seq[3] - 3 * seq[2] + 3 * seq[1] - seq[0]
                    if abs(d) <= 100.0:
                        ind['pn'] = abs(d) / np.sqrt(40.0)
                # HOD：L1(米) 四阶差分 / √70（5 点滑窗）
                if l1m is not None and len(st['l1m']) >= 4:
                    seq = st['l1m'][-4:] + [l1m]
                    d = seq[4] - 4 * seq[3] + 6 * seq[2] - 4 * seq[1] + seq[0]
                    if abs(d) <= 100.0:
                        ind['hod'] = abs(d) / np.sqrt(70.0)
                if so.get('snr') is not None:
                    ind['snr'] = so['snr']
                if ind:
                    satInd[prn] = ind
                # 更新滑窗状态（GF/L1 缓冲、MP/CMC/IOD 前值）
                if gf is not None:
                    st['gf'].append(gf)
                    if len(st['gf']) > 3:
                        st['gf'].pop(0)
                    st['gfep'] = epoch
                if l1m is not None:
                    st['l1m'].append(l1m)
                    if len(st['l1m']) > 4:
                        st['l1m'].pop(0)
                st['mp'] = mp
                st['cmc'] = cmc
                st['ep'] = epoch
                newQC[prn] = st
            qcState = newQC
            # 指标子集过滤（qc_indicators：逗号分隔，空 = 全部 6 指标；
            # pn/hod/iod 高度相关，可只留 pn+iod 降低熵权冗余）
            _keep = {s.strip() for s in
                     str(cfg.get('qc_indicators', '') or '').split(',') if s.strip()}
            if _keep:
                satInd = {p: {k: v for k, v in ind.items() if k in _keep}
                          for p, ind in satInd.items()}
            # winsorize：列向截尾到 [q, 1−q] 分位（单星异常值不再压扁
            # 跨星 min-max 归一化域），p 时间平滑（降低逐历元抖动）
            _meth = str(cfg.get('qc_method', 'gra'))
            if _meth == 'abs':
                # 绝对阈值标定法：指标对跨历元分位数 [q,1−q] 打绝对分，
                # p 跨历元可比（整体质量差时相对评分集体失真的根治）；
                # 热身（qc_abs_warmup 历元）前回退 gra 相对评分
                absN += 1
                _warm = max(int(cfg.get('qc_abs_warmup', 300)), 10)
                for ind in satInd.values():
                    for k, v in ind.items():
                        if v is not None and np.isfinite(v):
                            absHist.setdefault(k, []).append(float(v))
                _cap = _warm * 6
                for k in absHist:
                    if len(absHist[k]) > _cap:
                        del absHist[k][:len(absHist[k]) - _cap]
                if absN >= _warm:
                    pmap = qc_abs_weights(satInd, abs_calibrate(absHist))
                else:
                    pmap = qc_grey_weights(
                        satInd, rho=float(cfg.get('qc_rho', 0.5)),
                        winsor=float(cfg.get('qc_winsor', 0.1)), method='gra')
            else:
                pmap = qc_grey_weights(
                    satInd, rho=float(cfg.get('qc_rho', 0.5)),
                    winsor=float(cfg.get('qc_winsor', 0.1)), method=_meth)
            _nS = max(int(cfg.get('qc_smooth', 1)), 1)
            for prn, so in satObs.items():
                p = pmap.get(prn)
                if p is not None:
                    st = newQC.get(prn)
                    if st is not None:
                        st['p'].append(p)
                        if len(st['p']) > _nS:
                            st['p'].pop(0)
                        p = float(np.mean(st['p']))
                so['qcW'] = p if p is not None else 1.0
                so['qcPbar'] = 1.0
            # p̄ 取实际施加值的均值（含未评分星的 1.0），mul 归一化基准
            _pbar = float(np.mean([so['qcW'] for so in satObs.values()])) \
                if satObs else 1.0
            for so in satObs.values():
                so['qcPbar'] = _pbar
            if qcwDbg is not None:
                for prn, so in satObs.items():
                    si = satInd.get(prn, {})
                    qcwDbg.append((sodDbg, prn, so['qcW'],
                                   si.get('mpd', float('nan')),
                                   si.get('cmcd', float('nan')),
                                   si.get('pn', float('nan')),
                                   si.get('hod', float('nan')),
                                   si.get('iod', float('nan')),
                                   si.get('snr', float('nan'))))

        # --- 首历元：无有效初值时用码定位确定 ---
        if kf is None:
            x0Code = _codePosition(satObs, satClkInEpoch, gSysList, x0)
            if x0Code is None:
                continue
            x0 = list(x0Code[:3])
            clk0 = x0Code[3]
            xInit = np.zeros(nBase)
            xInit[:3] = x0
            # 接收机钟差状态用码解初值（RTKLIB udclk_ppp 同源）
            for i in range(3, 3 + nSys):
                xInit[i] = clk0[i - 3]
            if estZwd and tropSaas:
                # RTKLIB udtrop 同源：ZTD 状态 = 总天顶延迟，初值 Saastamoinen（含标准湿度湿延迟）
                _blh0 = xyz2blh(x0[0], x0[1], x0[2])
                _zhd0, _zwd0 = saastamoinenZenith(_blh0, humi=0.7)
                xInit[zwdIdx] = _zhd0 + _zwd0
            PInit = np.zeros((nBase, nBase))
            PInit[:3, :3] = np.diag([sigmaPosInit ** 2] * 3)
            for i in range(3, 3 + nSys):
                PInit[i, i] = sigmaClk ** 2
            # d（码-相钟差差）：初值 0、白噪声大方差（覆盖钟漂累积量）
            for i in range(dIdx, dIdx + nSys):
                PInit[i, i] = sigmaClk ** 2
            if estZwd:
                PInit[zwdIdx, zwdIdx] = sigmaZwdInit ** 2
            if estGrad:
                # RTKLIB 同源：梯度状态初值 1e-6，初始方差 VAR_GRA = SQR(0.01)
                for gi in (gradNIdx, gradEIdx):
                    xInit[gi] = 1e-6
                    PInit[gi, gi] = float(cfg.get('sigma_grad_init', 0.01)) ** 2
            kf = KalmanFilter(xInit, PInit)

        # --- 当前状态 ---
        rr = kf.x[:3, 0]
        lat, lon, h = xyz2blh(rr[0], rr[1], rr[2])
        enuRot = bl2enuRot(lat, lon)
        doys = datetime2doys(epoch)
        # saas_nmf 模式：天顶干延迟（历元级，与卫星无关）
        zhdEpoch = 0.0
        if tropSaas and estZwd and h < 200000.0:
            zhdEpoch = saastamoinenZenith((lat, lon, h), humi=0.0)[0]
        dt = 0.0
        if tracks:
            lastList = [t['lastEpoch'] for t in tracks.values() if t['lastEpoch'] is not None]
            if lastList:
                dt = (epoch - max(lastList)).total_seconds()

        # --- 过期弧段清理（长时间未见的卫星，其模糊度状态在下次出现时重建）---
        for prn in list(tracks.keys()):
            tr = tracks[prn]
            if tr['lastEpoch'] is not None and \
                    (epoch - tr['lastEpoch']).total_seconds() > 2.0 * float(cfg['gap_reset_sec']):
                del tracks[prn]

        # --- 过程噪声 ---
        # 接收机钟差与 d（码-相钟差差）：白噪声，每历元重置。RTKLIB
        # udclk_ppp 经 initx() 重置时会清零该状态与其他状态的全部协方差
        # （P 行列置 0 后设对角 VAR_CLK）；若只累加 Q 而保留协方差，巨幅
        # 载波更新建立的强负相关会在 S 收缩后产生病态增益 → 发散。
        for i in list(range(3, 3 + nSys)) + list(range(dIdx, dIdx + nSys)):
            kf.P[i, :] = 0.0
            kf.P[:, i] = 0.0
            kf.P[i, i] = sigmaClk ** 2
        Q = np.zeros((kf.n, kf.n))
        if mode == 'kinematic' and dt > 0:
            qPos = float(cfg['q_pos_kinematic']) * dt
            Q[np.ix_(range(3), range(3))] += np.eye(3) * qPos
        if estZwd and dt > 0:
            Q[zwdIdx, zwdIdx] = float(cfg['q_zwd']) * dt
        if estGrad and dt > 0:
            # RTKLIB udtrop：梯度过程噪声 SQR(prn_trop*0.1)·dt（prn_trop 缺省 1e-3）
            qGrad = float(cfg.get('q_grad', 1e-8)) * dt
            Q[gradNIdx, gradNIdx] = qGrad
            Q[gradEIdx, gradEIdx] = qGrad

        kf.predict(Q)

        # --- 每卫星几何与量测构建 ---
        # 行结构：保存几何分量，更新时用当前状态重算非线性模型 h(x)，
        # 避免 H@x 与 |rs-x| 的线性化偏差（e·rs 常数项）污染 innovation
        rows = []
        codeRows = []  # 供 DOP 计算

        # ERP：UT1−UTC（GMST / 太阳月亮方向 / 固体潮 / 相位缠绕统一使用）
        dut1 = getErp(erpData, epoch)['dut1'] if erpData is not None else 0.0
        sunPos = sunPosECEF(epoch, dut1) \
            if (cfg.get('windup', True) or cfg.get('solid_tide', True)
                or atxData is not None) else None

        # 固体潮：改正到卫星视线上等价于改正测站坐标（此处直接修正几何距离）
        tideCorr = None
        if cfg.get('solid_tide', True):
            tideCorr = solidEarthTide(epoch, rr, dut1)

        # 海潮载荷：本历元 ENU 位移（mm~cm 级），投影到视线
        oceanCorr = None
        if blqData is not None:
            oceanDisp = oceanLoadingDisp(blqData, station, epoch)
            if oceanDisp is not None:
                oceanCorr = enuRot.T @ oceanDisp   # ENU → ECEF

        # 极潮：极移激发的测站位移（~2cm，准常值），RTKLIB tide_pole 同源
        poleCorr = None
        if erpData is not None and cfg.get('pole_tide', True):
            erpNow = getErp(erpData, epoch)
            poleDispENU = poleTide(epoch, (lat, lon, h),
                                   erpNow['xp'], erpNow['yp'])
            poleCorr = enuRot.T @ poleDispENU      # ENU → ECEF

        for prn, so in satObs.items():
            gnssSys = prn[0]
            rs = so['rs'][:3]
            rho, e = geodist(rs, rr)
            if not np.isfinite(rho):
                continue

            # 高度角/方位角
            _, az, elev = getEle(rs, rr, enuRot)
            if elev < cutoff:
                continue

            sinE = max(np.sin(elev), 1e-3)
            elevFactor = float(cfg['sigma_elev_a']) + float(cfg['sigma_elev_b']) / sinE
            # 观测定权（σ² = σ0²·obsFactor²），分解为基础因子 f × 质量因子 g：
            #   基础因子 f：snr_weight 开启时用相对 SNR 权代替高度角权
            #              f = SNR_max/SNR（逐历元相对、均值保持，限幅 [0.25,4]），
            #              卫星无 SNR 或未开启时回落高度角模型 f = a+b/sinE
            #   质量因子 g：qc_weight 开启时按 qc_mode 组合灰色关联评分 p
            #       replace ：σ = σ0/√p           完全替换（f 不参与）
            #       hybrid  ：σ = σ0·f/√p         非均值保持（平均噪声被 1/√p̄ 抬升）
            #       mul     ：σ = σ0·f·√(p̄/p)    乘性混合，g=p̄/p 逐历元均值保持
            #       mul_pure：σ = σ0·√(p̄/p)      无 f 项（消融对照）
            #   即 snr_weight 只替换 f 位置的高度角项，与 qc_mode 正交组合；
            #   二者均关闭时 σ = σ0·elevFactor（RTKLIB 同源默认）
            if cfg.get('snr_weight', False) and so.get('snr'):
                fBase = min(max(float(epochMaxSnr) / max(float(so['snr']), 1.0),
                                0.25), 4.0)
            else:
                fBase = elevFactor
            if cfg.get('qc_weight', False):
                _mode = cfg.get('qc_mode', 'replace')
                _p = max(float(so.get('qcW', 1.0)), 1e-3)
                if _mode == 'replace':
                    obsFactor = 1.0 / np.sqrt(_p)
                elif _mode == 'hybrid':
                    obsFactor = fBase / np.sqrt(_p)
                elif _mode == 'mul':
                    obsFactor = fBase * np.sqrt(
                        max(float(so.get('qcPbar', 1.0)), 1e-3) / _p)
                elif _mode == 'mul_pure':
                    obsFactor = np.sqrt(max(float(so.get('qcPbar', 1.0)),
                                            1e-3) / _p)
                else:
                    obsFactor = fBase
            else:
                obsFactor = fBase
            # --- IGGⅢ 残差抗差层（RAE：上一历元标准化新息反馈，可选）---
            # 码/载波分别取等价权；新弧段或无历史时 w=1 不干预
            iggiWP = iggiWL = 1.0
            if cfg.get('robust_iggi', False):
                stR = resState.get(prn)
                if stR is not None \
                        and (epoch - stR['ep']).total_seconds() <= 90.0:
                    _k0 = float(cfg.get('robust_k0', 2.0))
                    _k1 = float(cfg.get('robust_k1', 3.5))
                    _wm = float(cfg.get('robust_wmin', 0.1))
                    if stR['P']:
                        iggiWP = _iggiFactor(stR['P'], _k0, _k1, _wm)
                    if stR['L']:
                        iggiWL = _iggiFactor(stR['L'], _k0, _k1, _wm)
            mfw = tropMapping(elev)
            cotz = 1.0 / np.tan(elev)   # 梯度映射用（MacMillan）

            # 对流层先验：
            #   saas_nmf：NMF 干映射 × 天顶干延迟（RTKLIB trop_model_prec 同源）
            #   collins ：Collins 经验模型全延迟
            tropPrior = 0.0
            if tropSaas:
                mfwH, mfw = nmfMapf(doys, (lat, lon, h), elev)
                tropPrior = mfwH * zhdEpoch
            elif h < 200000.0:
                tropPrior = tropColins(elev, lat, lon, h, doys)

            # 固体潮 + 海潮载荷：测站位移 d 使视距变化 Δρ = e·d
            # （e 为卫星→接收机方向；RTKLIB res_ppp 将 disp 加到 rr 后再算
            #   geodist，等价于 corr = +e·disp）
            # 固体潮 + 海潮载荷：测站位移 d 使视距变化 Δρ = e·d
            # （e 为卫星→接收机方向；RTKLIB res_ppp 将 disp 加到 rr 后再算
            #   geodist，等价于 corr = +e·disp）
            corr = 0.0
            if tideCorr is not None:
                corr += float(np.dot(tideCorr, e))
            if oceanCorr is not None:
                corr += float(np.dot(oceanCorr, e))
            if poleCorr is not None:
                corr += float(np.dot(poleCorr, e))

            satClk = satClkInEpoch.get(prn)
            if satClk is None:
                continue

            # 卫星天线 PCO/PCV（相位中心相对质心沿视线偏移 + 天底角相位变化）
            if atxData is not None and gnssSys in atxTags:
                b1, tag1, b2, tag2 = atxTags[gnssSys]
                f1v = bandChoose[gnssSys][b1]['freq']
                f2v = bandChoose[gnssSys][b2]['freq']
                # PCO 轴系（RTKLIB satantoff 同源）：Z 指向地心，Y = ez×es(sat→太阳)，
                # X = ey×ez —— ANTEX 卫星 PCO X 轴指向太阳（太阳指向姿态），非速度沿轨！
                ez = -rs / np.linalg.norm(rs)                 # Z 轴：指向地心
                es = (sunPos - rs) / np.linalg.norm(sunPos - rs)
                ey = np.cross(ez, es)
                ey /= np.linalg.norm(ey)
                ex = np.cross(ey, ez)
                ex /= np.linalg.norm(ex)
                # 真天底角（RTKLIB satantpcv 同源）：nadir = acos(dot(e, ez))，
                # e(星→接收机) 与 ez(星→地心) 均为卫星发出方向，天顶时≈0°
                nadir = np.degrees(np.arccos(np.clip(float(np.dot(ez, e)), -1.0, 1.0)))
                pcoSat, pcvSat = satAntCorr(atxData, prn, epoch, gnssSys,
                                            tag1, f1v, tag2, f2v, nadir)
                if pcoSat is not None and cfg.get('sat_pco', True):
                    dEC = ex * pcoSat[0] + ey * pcoSat[1] + ez * pcoSat[2]
                    corr += -float(np.dot(e, dEC)) + pcvSat
                # 接收机天线 PCO/PCV + RINEX 头 DELTA H/E/N（RTKLIB antdel 同源）
                if recAntType is not None or antDel is not None:
                    pcoRec, pcvRec = None, None
                    if atxData is not None and recAntType is not None \
                            and cfg.get('rec_pco', True):
                        zen = np.degrees(np.pi / 2.0 - elev)      # 天顶角
                        pcoRec, pcvRec = recAntCorr(atxData, recAntType, gnssSys,
                                                    tag1, f1v, tag2, f2v, zen)
                    dNEU = None
                    if pcoRec is not None:
                        dNEU = np.array([pcoRec[0], pcoRec[1], pcoRec[2]])
                    if antDel is not None and cfg.get('rec_ant_delta', False):
                        # DELTA H/E/N：ARP 相对标石的位移，仅当解算基准为
                        # 标石坐标（或参考坐标为标石）时才应施加。IGS SINEX
                        # 参考坐标指 ARP，施加 antdel 会使估出位置 = ARP +
                        # antdel（U 偏高 antdel_U）→ 默认关闭（RTKLIB GUI
                        # 默认 antdel=0，同源行为）。
                        dNEU = (dNEU if dNEU is not None else np.zeros(3)) \
                            + np.array([antDel[2], antDel[1], antDel[0]])
                    if dNEU is not None:
                        # 接收机相位中心位移 dRec 使视距变化 Δρ = +e·dRec
                        # （e=卫星→接收机；抬高 H 时贡献 −H·sin(el)，视距变短；
                        #   RTKLIB 等价于 dantr=e·dRec 从观测中扣除）
                        dRecEC = enuRot.T @ np.array([dNEU[1], dNEU[0], dNEU[2]])
                        corr += float(np.dot(e, dRecEC)) + (pcvRec or 0.0)

            rowBase = {'prn': prn, 'sys': gnssSys, 'e': e,
                       'rho': rho, 'corr': corr, 'dtsC': satClk * CLIGHT,
                       'tropPrior': tropPrior, 'mfw': mfw,
                       'cotz': cotz, 'az': az, 'varSat': so['varSat'],
                       'zhd': zhdEpoch}

            # 伪距量测（R 叠加卫星轨道/钟差方差 varSat，RTKLIB var_rs 同源；
            # IGGⅢ 等价权作用于量测噪声项，产品方差不膨胀）
            rows.append(dict(rowBase, type='P', z=so['P'],
                             R=(float(cfg['sigma_p']) * obsFactor) ** 2
                             / iggiWP + rowBase['varSat']))
            codeRows.append((e, gnssSys))

            # 载波量测（弧段管理延迟到相位-码跳变检测之后统一执行）
            if so['L'] is None:
                continue
            tr = tracks.get(prn)
            if cfg.get('windup', True):
                w = phaseWindUp(prn, rs, rr, sunPos, tr['windup'] if tr else 0.0)
            else:
                w = tr['windup'] if tr else 0.0

            isNewSat = tr is None
            gapNew = (not isNewSat and tr['lastEpoch'] is not None and
                      (epoch - tr['lastEpoch']).total_seconds() > float(cfg['gap_reset_sec']))
            slipped = epoch in slipMap.get(prn, ())

            # 新弧段模糊度候选初值：当前码-相残差（标准做法），
            # 避免首历元巨幅新息（载波整周计数含锁定时刻距离）反复触发剔除
            clockNow = float(kf.x[3 + gSysIndex[gnssSys], 0])
            zwdNow = float(kf.x[zwdIdx, 0]) if estZwd else 0.0
            if tropSaas:
                zwdNow = zwdNow - zhdEpoch   # ZTD 状态 → 湿延迟残余
            tropWet = zwdNow * mfw
            if estGrad:
                gnNow = float(kf.x[gradNIdx, 0])
                geNow = float(kf.x[gradEIdx, 0])
                tropWet *= 1.0 + cotz * (gnNow * np.cos(az) + geNow * np.sin(az))
            hNoAmb = rho + corr - satClk * CLIGHT + clockNow + tropPrior + tropWet
            bInit = (so['L'] - so['lamIF'] * w) - hNoAmb

            rows.append(dict(rowBase, type='L',
                             z=so['L'] - so['lamIF'] * w,
                             R=(float(cfg['sigma_l']) * obsFactor) ** 2
                             / iggiWL
                             + rowBase['varSat'],
                             w=w, isNewSat=isNewSat, gapNew=gapNew,
                             slipped=slipped, bInit=bInit))

        # --- 相位-码跳变改正（RTKLIB udbias 同源机制）---
        # 每历元对连续弧段卫星计算码新息，均值超过 0.0005·c 判定为接收机
        # 毫秒级钟跳（MW/GF 对码钟跳不敏感或全星同跳，无法区分）：
        # 连续弧段模糊度整体平移 offset（保持相位-码相干性，避免全部重置），
        # 并撤销本历元周跳标记使弧段保留。
        jumpShift = 0.0
        if cfg.get('phase_code_jump', True) and tracks:
            pRows = {row['prn']: row for row in rows if row['type'] == 'P'}
            offs = []
            contIdxs = []
            for prn, tr in tracks.items():
                row = pRows.get(prn)
                if row is None or tr['lastEpoch'] is None:
                    continue
                if (epoch - tr['lastEpoch']).total_seconds() > float(cfg['gap_reset_sec']):
                    continue
                clock = float(kf.x[3 + gSysIndex[row['sys']], 0])
                zwdS0 = float(kf.x[zwdIdx, 0]) if estZwd else 0.0
                if tropSaas:
                    zwdS0 = zwdS0 - row['zhd']
                tropWet0 = zwdS0 * row['mfw']
                if estGrad:
                    gn0 = float(kf.x[gradNIdx, 0])
                    ge0 = float(kf.x[gradEIdx, 0])
                    tropWet0 *= 1.0 + row['cotz'] * (gn0 * np.cos(row['az'])
                                                     + ge0 * np.sin(row['az']))
                hHat0 = row['rho'] + row['corr'] - row['dtsC'] + clock \
                    + row['tropPrior'] + tropWet0
                offs.append(row['z'] - hHat0)
                contIdxs.append(tr['idx'])
            if len(offs) >= 2:
                offMean = float(np.mean(offs))
                if abs(offMean) > 0.0005 * CLIGHT:
                    jumpShift = offMean
                    for idxA in contIdxs:
                        kf.x[idxA, 0] -= jumpShift   # amb −= 码新息均值（保持 L-P 相干）
                    for prn in slipMap:
                        slipMap[prn].discard(epoch)
                    print('相位-码跳变改正：offset=%.2f m（%d 颗连续卫星），'
                          '模糊度整体平移，弧段保留。' % (jumpShift, len(offs)))

        # --- 载波弧段管理（跳变改正后统一执行）---
        for row in rows:
            if row['type'] != 'L':
                continue
            prn = row['prn']
            tr = tracks.get(prn)
            if tr is None or row['gapNew'] or (row['slipped'] and jumpShift == 0.0):
                # 新弧段：用候选初值新建模糊度状态
                idx = kf.add_state(row['bInit'], sigmaAmb ** 2)
                tracks[prn] = {'idx': idx, 'lastEpoch': epoch, 'windup': row['w']}
                row['contArc'] = False
            else:
                # 连续弧段（含被跳变改正挽救的弧段）
                tr['windup'] = row['w']
                tr['lastEpoch'] = epoch
                row['contArc'] = True
            row['ambIdx'] = tracks[prn]['idx']

        # ---- RTKLIB 同源批量量测更新（整历元单次联合 filter）----
        # RTKLIB pppos 将全部卫星的码+载波残差堆为一个向量，用同一状态
        # 快照一次 filter：最小二乘自动把共模钟漂分给钟差/d 状态、码-相
        # 残差之差分给各星模糊度，位置只吸收几何差分信息。
        # 任何"逐星/逐行"序贯处理都有序贯污染：白噪声钟差方差在第一颗星
        # 更新后即塌缩到 ~R，后续星的个体差新息以增益≈1 直接冲击位置
        # （实测 G19 amb 被甩到 -3888、G20 新息 +11279 → 爆炸，fix17）。
        # 批量更新无逐行拒绝门（RTKLIB ppp 同源），仅保留 guard 防护。
        nMeas = 0
        nRej = 0

        def _measModel(row):
            """非线性量测模型 h(x) 与设计阵 H（用当前状态）。"""
            idxClk = 3 + gSysIndex[row['sys']]
            idxD = dIdx + gSysIndex[row['sys']]
            clock = float(kf.x[idxClk, 0])
            zwdS = float(kf.x[zwdIdx, 0]) if estZwd else 0.0
            if tropSaas:
                zwdS = zwdS - row['zhd']
            gradFac = 1.0
            if estGrad:
                gnS = float(kf.x[gradNIdx, 0])
                geS = float(kf.x[gradEIdx, 0])
                # MacMillan 模型：m_w = m_0·(1 + cot(el)·(Gn·cosA + Ge·sinA))
                gradFac = 1.0 + row['cotz'] * (gnS * np.cos(row['az'])
                                               + geS * np.sin(row['az']))
            tropWet = zwdS * row['mfw'] * gradFac if estZwd else 0.0
            hHat = row['rho'] + row['corr'] - row['dtsC'] + clock \
                + row['tropPrior'] + tropWet
            H = np.zeros((1, kf.n))
            H[0, :3] = row['e']
            H[0, idxClk] = 1.0
            if row['type'] == 'L':
                # 载波量测含 d（码-相钟差差）：钟差由码驱动吸收钟漂，载波
                # 与码钟差之差由 d 白噪声状态吸收（RTKLIB 钟差吸收共模同效）
                H[0, idxD] = 1.0
                hHat += float(kf.x[idxD, 0])
            if estZwd:
                H[0, zwdIdx] = row['mfw'] * gradFac
            if estGrad:
                # ∂trop/∂G = m_0·cot(el)·zwd·{cosA, sinA}（RTKLIB dtdx 同式）
                cg = row['mfw'] * row['cotz'] * zwdS
                H[0, gradNIdx] = cg * np.cos(row['az'])
                H[0, gradEIdx] = cg * np.sin(row['az'])
            return hHat, H

        # 全部行（码+载波、全部卫星）堆叠为一次联合批量 filter
        H = np.zeros((len(rows), kf.n))
        z = np.zeros(len(rows))
        R = np.zeros((len(rows), len(rows)))
        hh = np.zeros(len(rows))
        for i, row in enumerate(rows):
            hHat, Hi = _measModel(row)
            if row['type'] == 'L':
                Hi[0, row['ambIdx']] = 1.0
                hHat += float(kf.x[row['ambIdx'], 0])
            H[i, :] = Hi[0, :]
            z[i] = row['z']
            R[i, i] = row['R']
            hh[i] = hHat
        innov, S, rejected = kf.update_multi(H, z, R, hHat=hh)
        if innov is not None:
            nMeas = len(rows)
            if dbgRes is not None:
                for i, row in enumerate(rows):
                    _sinel = float(np.clip(-np.dot(row['e'], enuRot[2]),
                                           -1.0, 1.0))
                    _ambDbg = float(kf.x[row['ambIdx'], 0]) \
                        if row['type'] == 'L' else 0.0
                    dbgRes.append((sodDbg, row['prn'], row['type'],
                                   np.degrees(np.arcsin(_sinel)),
                                   np.degrees(row['az']), float(innov[i]),
                                   float(S[i, i]), int(rejected), _ambDbg))
            # --- IGGⅢ 抗差反馈存储：标准化新息 v = innov/√R（纯量测噪声
            # 归一，码/载波同尺度）；新弧段载波历史作废。u 的自适应归一
            # （MAD）在下一历元应用时进行（见 _iggiFactor）---
            if cfg.get('robust_iggi', False):
                _rWin = max(int(cfg.get('robust_win', 8)), 1)
                for i, row in enumerate(rows):
                    stR = resState.setdefault(row['prn'],
                                              {'P': [], 'L': [], 'ep': epoch})
                    v = float(innov[i]) / np.sqrt(max(float(row['R']), 1e-12))
                    if row['type'] == 'L' and not row.get('contArc', True):
                        stR['L'] = []
                    buf = stR[row['type']]
                    buf.append(v)
                    if len(buf) > _rWin:
                        buf.pop(0)
                    stR['ep'] = epoch
        if nMeas < 4:
            continue

        # --- 失效模糊度状态剪枝（无任何弧段引用的旧状态）---
        usedIdx = set(tr['idx'] for tr in tracks.values())
        deadIdx = [i for i in range(nBase, kf.n) if i not in usedIdx]
        if deadIdx:
            stateMap = kf.remove_states(deadIdx)
            for prn in tracks:
                tracks[prn]['idx'] = stateMap.get(tracks[prn]['idx'], tracks[prn]['idx'])

        # --- 载波侧码-相一致性平移（已由量测更新内的 cL 共模扣除取代）---
        # 载波更新前已按连续弧段新息中位数 cL 扣除码-相钟差共模（见量测
        # 更新块），共模不再进入任何状态；此处若再把 offMed（≈cL）写进
        # amb，将与 cL 扣除双重处理（amb 每历元累积漂移 → 发散），故删除。
        # 毫秒级钟跳仍由量测更新前的码侧跳变改正（jumpShift）处理。

        # --- DOP（伪距几何）---
        if len(codeRows) >= 4:
            Gdop = np.zeros((len(codeRows), 3 + nSys))
            for i, (e, s) in enumerate(codeRows):
                Gdop[i, :3] = e
                Gdop[i, 3 + gSysIndex[s]] = 1.0
            try:
                Qm = np.linalg.pinv(Gdop.T @ Gdop)
                posData[epoch] = {
                    'HDOP': np.sqrt(Qm[0, 0] + Qm[1, 1]),
                    'VDOP': np.sqrt(Qm[2, 2]),
                    'PDOP': np.sqrt(Qm[0, 0] + Qm[1, 1] + Qm[2, 2]),
                    'GDOP': np.sqrt(Qm[0, 0] + Qm[1, 1] + Qm[2, 2] + Qm[3, 3]),
                }
            except Exception:
                posData[epoch] = {'HDOP': 0, 'VDOP': 0, 'PDOP': 0, 'GDOP': 0}
        else:
            continue

        # --- 记录结果 ---
        x, y, z = (float(kf.x[0, 0]), float(kf.x[1, 0]), float(kf.x[2, 0]))
        posData[epoch]['x'] = x
        posData[epoch]['y'] = y
        posData[epoch]['z'] = z
        # 坐标改正数：相对初始坐标（RINEX 近似坐标或首历元码定位值），ECEF 米
        posData[epoch]['dx'] = x - x0[0]
        posData[epoch]['dy'] = y - x0[1]
        posData[epoch]['dz'] = z - x0[2]
        posData[epoch]['clk'] = float(kf.x[3, 0])
        posData[epoch]['NSAT'] = len(codeRows)
        posData[epoch]['NREJ'] = nRej
        posData[epoch]['ZWD'] = (float(kf.x[zwdIdx, 0]) - zhdEpoch) if estZwd else 0.0
        posData[epoch]['SIG3D'] = float(np.sqrt(max(kf.P[0, 0] + kf.P[1, 1] + kf.P[2, 2], 0.0)))
        Xs.append(x)
        Ys.append(y)
        Zs.append(z)

    if len(Xs) == 0:
        printPanda('PPP 失败：无有效解算历元（请检查精密星历是否覆盖观测时段）。')
        return None, None, None, None

    # ---------------- ENU 参考坐标 ----------------
    refX, refY, refZ = _reference_xyz(posData, Xs, Ys, Zs, cfg, snx_file, station)

    print('PPP POS -> ' + '%15.3f' % refX + '%15.3f' % refY + '%15.3f' % refZ)
    print('PPP done : %d epochs, %.1f s' %
          (len(Xs), (datetime.datetime.now() - start_time).total_seconds()))

    # ---------------- 输出 ----------------
    if posFile is not None:
        writePosData(posData, posFile, refX, refY, refZ)
    if os.environ.get('PPP_DEBUG_RES') and dbgRes:
        with open(os.environ['PPP_DEBUG_RES'], 'w', encoding='utf-8') as f:
            f.write('sod prn type elev az innov S rej amb\n')
            for r in dbgRes:
                f.write('%.1f %s %s %.3f %.3f %.6f %.6f %d %.3f\n' % r)
    if os.environ.get('PPP_DEBUG_QCW') and qcwDbg:
        with open(os.environ['PPP_DEBUG_QCW'], 'w', encoding='utf-8') as f:
            f.write('sod prn qcW mpd cmcd pn hod iod snr\n')
            for r in qcwDbg:
                f.write('%.1f %s %.4f %.4f %.4f %.4f %.4f %.4g %.2f\n' % r)
    if plot_show or plot_save:
        plot_ppp(posData, refX, refY, refZ, plot_show=plot_show,
                 plot_save=plot_save,
                 conv_thresh=float(cfg.get('ppp_conv_thresh', 0.3)),
                 conv_hold_min=float(cfg.get('ppp_conv_hold', 5)))

    return posData, refX, refY, refZ


# ======================================================== 内部工具
def _reference_xyz(posData, Xs, Ys, Zs, cfg, snx_file, station):
    """
    ENU 参考：优先 SNX 周解；否则可用 ref_mean_skip_min 跳过滤波收敛段取均值。
    """
    if snx_file is not None and station is not None:
        try:
            ref = read_snx_coord(snx_file, station)
            print('SNX REF -> ' + '%15.3f' % ref[0] + '%15.3f' % ref[1] + '%15.3f' % ref[2])
            return float(ref[0]), float(ref[1]), float(ref[2])
        except Exception as e:
            print(f'读取 SNX 失败（{e}），改用坐标均值作为 ENU 参考。')

    skipMin = float(cfg.get('ref_mean_skip_min', 0))
    if skipMin > 0 and len(posData) > 10:
        epochs = sorted([ep for ep in posData if 'x' in posData[ep]])
        t0 = epochs[0]
        ex = [posData[ep]['x'] for ep in epochs
              if (ep - t0).total_seconds() >= skipMin * 60]
        ey = [posData[ep]['y'] for ep in epochs
              if (ep - t0).total_seconds() >= skipMin * 60]
        ez = [posData[ep]['z'] for ep in epochs
              if (ep - t0).total_seconds() >= skipMin * 60]
        if len(ex) >= 10:
            return float(np.mean(ex)), float(np.mean(ey)), float(np.mean(ez))

    return float(np.mean(Xs)), float(np.mean(Ys)), float(np.mean(Zs))


if __name__ == '__main__':
    print('本模块为 PPP 引擎，请使用 ppp/run_ppp.py 作为命令行入口。')
