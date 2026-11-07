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

        mwList = []       # (epoch, Nw)
        gfList = []       # (epoch, GF)
        slipSet = set()

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
            if lli1 != 0 or lli2 != 0:
                slipSet.add(epoch)
            mwList.append((epoch, L1 - L2 - coef * (P1 + P2) / lambdaW))
            gfList.append((epoch, (L1 - L2) * lambdaW))

        # ---- MW 历元间差分 + 滑动窗标准差 ----
        mwEpochs = [e for e, _ in mwList]
        mwVals = np.array([v for _, v in mwList])
        for i in range(1, len(mwVals)):
            diff = mwVals[i] - mwVals[i - 1]
            if (mwEpochs[i] - mwEpochs[i - 1]).total_seconds() > 120.0:
                slipSet.add(mwEpochs[i])   # 长中断后重置
                continue
            lo = max(0, i - MW_WINDOW)
            window = np.abs(np.diff(mwVals[lo:i]))
            sigma = float(np.std(window)) if len(window) >= 5 else 1.0
            thresh = max(mwSigmaFactor * sigma, MW_MIN_THRESH)
            if abs(diff) > thresh:
                slipSet.add(mwEpochs[i])

        # ---- GF 历元间差分 ----
        gfEpochs = [e for e, _ in gfList]
        gfVals = np.array([v for _, v in gfList])
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
