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
from ppp.corrections import sunPosECEF, solidEarthTide, phaseWindUp, tropMapping
from com.readAtx import readAntex, satAntCorr, recAntCorr, bandToAtxTag
from com.readBlq import readBlq, oceanLoadingDisp
from com.readDcb import readDcb, applyDcbSat, dcb_stat
from com.readErp import readErp, getErp
from ppp.pppfilter import KalmanFilter
from ppp.pppresult import writePosData


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
    return float(X[0, 0]), float(X[1, 0]), float(X[2, 0])


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
        if not station or station.upper() not in (blqData := readBlq(cfg['blq_file'])):
            print('警告：BLQ 文件中无测站 %s，海潮载荷改正未启用。' % (station or '?'))
            blqData = None
        else:
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

    zwdIdx = 3 + nSys   # ZWD 状态固定索引
    gradNIdx = zwdIdx + 1 if estGrad else None   # 梯度北分量（MacMillan Gn）
    gradEIdx = zwdIdx + 2 if estGrad else None   # 梯度东分量（MacMillan Ge）
    nBase = zwdIdx + (1 if estZwd else 0) + (2 if estGrad else 0)

    # ---------------- 结果 ----------------
    posData = {}
    Xs, Ys, Zs = [], [], []

    # ---------------- 历元循环 ----------------
    for epoch in sorted(obsData.keys()):
        if startDatetime is not None and epoch < startDatetime:
            continue
        if endDatetime is not None and epoch > endDatetime:
            continue

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
            satObs[prn] = {'P': PIF, 'L': LIF, 'lamIF': lamIF,
                           'sys': gnssSys, 'rs': rs, 'varSat': varSat}

        # --- 首历元：无有效初值时用码定位确定 ---
        if kf is None:
            x0Code = _codePosition(satObs, satClkInEpoch, gSysList, x0)
            if x0Code is None:
                continue
            x0 = list(x0Code)
            xInit = np.zeros(nBase)
            xInit[:3] = x0
            PInit = np.zeros((nBase, nBase))
            PInit[:3, :3] = np.diag([sigmaPosInit ** 2] * 3)
            for i in range(3, 3 + nSys):
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
        Q = np.zeros((kf.n, kf.n))
        if mode == 'kinematic' and dt > 0:
            qPos = float(cfg['q_pos_kinematic']) * dt
            Q[np.ix_(range(3), range(3))] += np.eye(3) * qPos
        for i in range(3, 3 + nSys):
            Q[i, i] = sigmaClk ** 2          # 接收机钟差：白噪声（每历元重置方差）
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
            if (cfg.get('windup', True) or cfg.get('solid_tide', True)) else None

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
            mfw = tropMapping(elev)
            cotz = 1.0 / np.tan(elev)   # 梯度映射用（MacMillan）

            # 对流层：Collins 先验（全延迟作干分量先验）+ 估计的残余湿延迟
            tropPrior = 0.0
            if h < 200000.0:
                tropPrior = tropColins(elev, lat, lon, h, doys)

            # 固体潮 + 海潮载荷：投影到视线方向（Δρ = −e·d）
            corr = 0.0
            if tideCorr is not None:
                corr -= float(np.dot(tideCorr, e))
            if oceanCorr is not None:
                corr -= float(np.dot(oceanCorr, e))

            satClk = satClkInEpoch.get(prn)
            if satClk is None:
                continue

            # 卫星天线 PCO/PCV（相位中心相对质心沿视线偏移 + 天底角相位变化）
            if atxData is not None and gnssSys in atxTags:
                b1, tag1, b2, tag2 = atxTags[gnssSys]
                f1v = bandChoose[gnssSys][b1]['freq']
                f2v = bandChoose[gnssSys][b2]['freq']
                ez = -rs / np.linalg.norm(rs)                 # 指向地心
                v = np.asarray(so['rs'][3:6], dtype=np.float64)
                vT = v - ez * float(np.dot(v, ez))            # 沿轨（正交化）
                nv = np.linalg.norm(vT)
                ex = vT / nv if nv > 1e-6 else np.cross(ez, np.array([0.0, 0.0, 1.0]))
                ex /= np.linalg.norm(ex)
                ey = np.cross(ez, ex)
                nadir = np.degrees(np.arccos(np.clip(float(np.dot(ez, e)), -1.0, 1.0)))
                pcoSat, pcvSat = satAntCorr(atxData, prn, epoch, gnssSys,
                                            tag1, f1v, tag2, f2v, nadir)
                if pcoSat is not None:
                    dEC = ex * pcoSat[0] + ey * pcoSat[1] + ez * pcoSat[2]
                    corr += -float(np.dot(e, dEC)) + pcvSat
                # 接收机天线 PCO/PCV
                if recAntType is not None:
                    zen = np.degrees(np.pi / 2.0 - elev)      # 天顶角
                    pcoRec, pcvRec = recAntCorr(atxData, recAntType, gnssSys,
                                                tag1, f1v, tag2, f2v, zen)
                    if pcoRec is not None:
                        dRecEC = enuRot.T @ np.array([pcoRec[1], pcoRec[0], pcoRec[2]])
                        corr += -float(np.dot(e, dRecEC)) + (pcvRec or 0.0)

            rowBase = {'prn': prn, 'sys': gnssSys, 'e': e,
                       'rho': rho, 'corr': corr, 'dtsC': satClk * CLIGHT,
                       'tropPrior': tropPrior, 'mfw': mfw,
                       'cotz': cotz, 'az': az, 'varSat': so['varSat']}

            # 伪距量测（R 叠加卫星轨道/钟差方差 varSat，RTKLIB var_rs 同源）
            rows.append(dict(rowBase, type='P', z=so['P'],
                             R=(float(cfg['sigma_p']) * elevFactor) ** 2
                             + rowBase['varSat']))
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
            tropWet = zwdNow * mfw
            if estGrad:
                gnNow = float(kf.x[gradNIdx, 0])
                geNow = float(kf.x[gradEIdx, 0])
                tropWet *= 1.0 + cotz * (gnNow * np.cos(az) + geNow * np.sin(az))
            hNoAmb = rho + corr - satClk * CLIGHT + clockNow + tropPrior + tropWet
            bInit = (so['L'] - so['lamIF'] * w) - hNoAmb

            rows.append(dict(rowBase, type='L',
                             z=so['L'] - so['lamIF'] * w,
                             R=(float(cfg['sigma_l']) * elevFactor) ** 2
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
            else:
                # 连续弧段（含被跳变改正挽救的弧段）
                tr['windup'] = row['w']
                tr['lastEpoch'] = epoch
            row['ambIdx'] = tracks[prn]['idx']

        # --- 量测更新（逐量测 + 粗差剔除）---
        nMeas = 0
        nRej = 0
        for row in rows:
            idxClk = 3 + gSysIndex[row['sys']]
            clock = float(kf.x[idxClk, 0])
            zwdS = float(kf.x[zwdIdx, 0]) if estZwd else 0.0
            gradFac = 1.0
            if estGrad:
                gnS = float(kf.x[gradNIdx, 0])
                geS = float(kf.x[gradEIdx, 0])
                # MacMillan 模型：m_w = m_0·(1 + cot(el)·(Gn·cosA + Ge·sinA))
                gradFac = 1.0 + row['cotz'] * (gnS * np.cos(row['az'])
                                               + geS * np.sin(row['az']))
            tropWet = zwdS * row['mfw'] * gradFac if estZwd else 0.0
            # 非线性量测模型 h(x)
            hHat = row['rho'] + row['corr'] - row['dtsC'] + clock \
                + row['tropPrior'] + tropWet
            H = np.zeros((1, kf.n))
            H[0, :3] = row['e']
            H[0, idxClk] = 1.0
            if estZwd:
                H[0, zwdIdx] = row['mfw'] * gradFac
            if estGrad:
                # ∂trop/∂G = m_0·cot(el)·zwd·{cosA, sinA}（RTKLIB dtdx 同式）
                cg = row['mfw'] * row['cotz'] * zwdS
                H[0, gradNIdx] = cg * np.cos(row['az'])
                H[0, gradEIdx] = cg * np.sin(row['az'])
            if row['type'] == 'L':
                H[0, row['ambIdx']] = 1.0
                hHat += float(kf.x[row['ambIdx'], 0])
            innov, S, rejected = kf.update(H, row['z'], row['R'], hHat=hHat,
                                           nsigma=float(cfg['outlier_nsigma']))
            if innov is None:
                continue
            nMeas += 1
            if rejected:
                nRej += 1
        if nMeas < 4:
            continue

        # --- 失效模糊度状态剪枝（无任何弧段引用的旧状态）---
        usedIdx = set(tr['idx'] for tr in tracks.values())
        deadIdx = [i for i in range(nBase, kf.n) if i not in usedIdx]
        if deadIdx:
            stateMap = kf.remove_states(deadIdx)
            for prn in tracks:
                tracks[prn]['idx'] = stateMap.get(tracks[prn]['idx'], tracks[prn]['idx'])

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
        posData[epoch]['ZWD'] = float(kf.x[zwdIdx, 0]) if cfg.get('est_zwd', True) else 0.0
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
