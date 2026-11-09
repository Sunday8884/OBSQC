#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
slip.py — PPP 周跳探测（MW 宽巷组合 + GF 几何无关组合）。

原理
----
- MW 组合：Nw = L1 - L2 - (f1-f2)/(f1+f2) * (P1+P2)，历元间差分，
  超过 max(滑动窗标准差 × mw_sigma, 最小阈值) 判定周跳；
- GF 组合：L1 - L2（几何无关，含电离层变化），历元间差分超过
  slip_gf_thresh（米）判定周跳；
- LLI 非零的历元直接标记周跳。

返回 {prn: set(epochs)}，供 ppp 主引擎做模糊度重置。
"""
import numpy as np

# MW 滑动统计窗口（历元数）
MW_WINDOW = 50
# MW 绝对最小阈值（周），防止小样本下阈值过紧
MW_MIN_THRESH = 1.0


def _get_lc(obsRec, bandC1, bandC2):
    """从单星观测记录中取 L1/L2（周）、P1/P2（m）及 LLI；缺数返回 None。"""
    bandL1 = 'L' + bandC1[1:]
    bandL2 = 'L' + bandC2[1:]
    L1 = obsRec.get(bandL1)
    L2 = obsRec.get(bandL2)
    P1 = obsRec.get(bandC1)
    P2 = obsRec.get(bandC2)
    lli1 = obsRec.get(bandL1 + 'LLI', 0)
    lli2 = obsRec.get(bandL2 + 'LLI', 0)
    if L1 is None or L2 is None or P1 is None or P2 is None:
        return None
    return L1, L2, P1, P2, lli1, lli2


def detectCycleSlip(obsData, bandChooseWithFreq, startDatetime=None, endDatetime=None, config=None):
    """
    obsData            : {epoch: {prn: {band: value}}}
    bandChooseWithFreq : 已填入 freq 的双频配置 {'G': {'C1C': {'freq':...}, 'C2W': {...}}}
    config             : slip_mw_sigma / slip_gf_thresh

    返回 slipMap = {prn: set(epochs)}（发生周跳的历元集合）。
    """
    if config is None:
        config = {}
    mwSigmaFactor = float(config.get('slip_mw_sigma', 4.0))
    gfThresh = float(config.get('slip_gf_thresh', 0.10))

    slipMap = {}
    epochList = sorted(obsData.keys())

    # 第一遍：收集全星 MW 值（星间共模消除用）
    mwAll = {}   # {epoch: {prn: mwVal}}
    lcMeta = {}  # {prn: [(epoch, L1, L2, lli1, lli2)]} 有完整四量的记录
    for prn in _all_prn(obsData, startDatetime, endDatetime):
        gnssSys = prn[0]
        if gnssSys not in bandChooseWithFreq:
            continue
        bands = list(bandChooseWithFreq[gnssSys].keys())[:2]
        if len(bands) < 2:
            continue
        f1 = bandChooseWithFreq[gnssSys][bands[0]]['freq']
        f2 = bandChooseWithFreq[gnssSys][bands[1]]['freq']
        if f1 is None or f2 is None or abs(f1 - f2) < 1.0:
            continue
        lambdaW = 299792458.0 / (f1 - f2)
        coef = (f1 - f2) / (f1 + f2)

        for epoch in epochList:
            if startDatetime is not None and epoch < startDatetime:
                continue
            if endDatetime is not None and epoch > endDatetime:
                continue
            rec = obsData[epoch].get(prn)
            if not rec:
                continue
            lc = _get_lc(rec, bands[0], bands[1])
            if lc is None:
                continue
            L1, L2, P1, P2, lli1, lli2 = lc
            # MW 组合（周）：Nw = φ1 − φ2 − (f1−f2)/(f1+f2)·(P1·f1 + P2·f2)/c
            # 注意码项必须按频率加权（P·f/c），几何项才能与 φ1−φ2 的
            # ρ·(f1−f2)/c 抵消；若写成 (P1+P2)/λw，几何残留 = 几何变化率
            # （实测 ~15000 周/30s），将每历元误报周跳、模糊度弧段无法维持
            mwAll.setdefault(epoch, {})[prn] = \
                L1 - L2 - coef * (P1 * f1 + P2 * f2) / 299792458.0
            # GF 组合（米）：φ1·λ1 − φ2·λ2，几何完全抵消，历元间差分
            # 只剩电离层变化（mm~cm 级）；若写成 (L1−L2)·λw，几何残留
            # (f1−f2)/(f1+f2)·Δρ（实测 ~2km/30s）同样每历元误报周跳
            lcMeta.setdefault(prn, []).append(
                (epoch, lambdaW, lli1, lli2,
                 L1 * 299792458.0 / f1 - L2 * 299792458.0 / f2))

    # 全星 MW 中位数（接收机钟漂/钟跳为全星共模，中位数差分可消除；
    # 单星周跳在星间残差中保留）
    medMap = {ep: float(np.median(list(v.values())))
              for ep, v in mwAll.items() if len(v) >= 4}

    for prn, metas in lcMeta.items():
        slipSet = set()

        # LLI 直接标记
        for epoch, lambdaW, lli1, lli2, gfVal in metas:
            if lli1 != 0 or lli2 != 0:
                slipSet.add(epoch)

        # ---- MW 历元间差分 + 滑动窗标准差（星间共模消除）----
        # 接收机钟漂/钟跳使 P1、P2 历元间同量跳变，MW 码项产生数百周
        # 共模变化（实测本数据 ~1400 周/30s），直接差分将每历元误报周跳、
        # 模糊度弧段无法维持；改用 MW − 全星中位数 的星间残差做历元差分
        mwList = [(ep, mwAll[ep][prn]) for ep in sorted(mwAll) if prn in mwAll[ep]]
        mwEpochs = [e for e, _ in mwList]
        mwVals = np.array([v for _, v in mwList])
        medVals = np.array([medMap.get(e, np.nan) for e in mwEpochs])
        resMW = mwVals - medVals            # 星间残差（消除共模钟漂）
        for i in range(1, len(mwVals)):
            diff = resMW[i] - resMW[i - 1]
            if (mwEpochs[i] - mwEpochs[i - 1]).total_seconds() > 120.0:
                slipSet.add(mwEpochs[i])   # 长中断后重置
                continue
            lo = max(0, i - MW_WINDOW)
            window = np.abs(np.diff(resMW[lo:i]))
            sigma = float(np.std(window)) if len(window) >= 5 else 1.0
            thresh = max(mwSigmaFactor * sigma, MW_MIN_THRESH)
            if abs(diff) > thresh:
                slipSet.add(mwEpochs[i])

        # ---- GF 历元间差分（纯相位组合，天然免疫钟漂）----
        gfEpochs = [e for e, _, _, _, gf in metas]
        gfVals = np.array([gf for _, _, _, _, gf in metas])
        for i in range(1, len(gfVals)):
            diff = gfVals[i] - gfVals[i - 1]
            if (gfEpochs[i] - gfEpochs[i - 1]).total_seconds() > 120.0:
                continue                   # 长中断已由 MW 标记
            if abs(diff) > gfThresh:
                slipSet.add(gfEpochs[i])

        if slipSet:
            slipMap[prn] = slipSet

    return slipMap


def _all_prn(obsData, startDatetime=None, endDatetime=None):
    """收集时间范围内的所有卫星号。"""
    prnSet = set()
    for epoch in sorted(obsData.keys()):
        if startDatetime is not None and epoch < startDatetime:
            continue
        if endDatetime is not None and epoch > endDatetime:
            continue
        for prn in obsData[epoch]:
            if isinstance(prn, str) and len(prn) >= 2:
                prnSet.add(prn)
    return sorted(prnSet)
