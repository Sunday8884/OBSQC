#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
readBlq.py — BLQ 海潮载荷位移文件（如 FES2004）读取与计算。

能力
----
- readBlq()            : 解析测站的 11 分潮（M2 S2 N2 K2 K1 O1 P1 Q1 MF MM SSA）
                         振幅（m）与相位（度），行序：径向 / EW / NS
- oceanLoadingDisp()   : 给定测站与历元，计算海潮载荷 ENU 位移（米）

说明
----
- BLQ 约定：位移正向为向上、向南、向西，相位滞后为正（相对格林尼治）；
- 幅角采用 Doodson 变量低精度展开（s/h/p 取 Meeus 主项），
  未加交点（nodal）改正，分潮振幅 mm~cm 级时误差可忽略；
- 返回 ENU 位移可直接与固体潮叠加：视线改正 Δρ = −e·(Rᵀ·dENU)。
"""
import math
import numpy as np

# 分潮固定顺序（BLQ 文件列顺序）
CONSTITUENTS = ['M2', 'S2', 'N2', 'K2', 'K1', 'O1', 'P1', 'Q1', 'MF', 'MM', 'SSA']


def readBlq(blqFile):
    """
    解析 BLQ 文件。
    返回 {station: {'amp': (3,11) m [radial,ew,ns], 'phase': (3,11) 度, 'lonlat': (lon,lat,h) or None}}
    """
    stations = {}
    pendingName = None
    rows = []

    with open(blqFile, 'r', encoding='utf-8', errors='ignore') as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith('$$'):
                continue
            toks = s.split()
            # 测站名行：单个 token（或名称+注释），后随 6 行数据
            if len(toks) == 1 and not _isNumber(toks[0]):
                if pendingName and len(rows) >= 6:
                    _store(stations, pendingName, rows)
                pendingName = toks[0].upper()
                rows = []
                continue
            if _isNumber(toks[0]) and len(toks) >= 10:
                rows.append([float(v) for v in toks])
            # 其它内容（表头等）忽略

    if pendingName and len(rows) >= 6:
        _store(stations, pendingName, rows)
    return stations


def _isNumber(t):
    try:
        float(t)
        return True
    except Exception:
        return False


def _store(stations, name, rows):
    try:
        amp = np.array(rows[0:3], dtype=np.float64)      # (3,11) m
        phase = np.array(rows[3:6], dtype=np.float64)    # (3,11) deg
    except Exception:
        return
    stations[name] = {'amp': amp, 'phase': phase, 'lonlat': None}


# ---------------------------------------------------------------- 天文幅角
def _astroArgs(epoch):
    """
    计算 Doodson 幅角基础量（弧度）：
        T = 格林尼治平恒星时角, s = 月亮平黄经, h = 太阳平黄经, p = 月亮近地点平黄经
    """
    jd = epoch.toordinal() + 0.5 + \
        (epoch.hour * 3600 + epoch.minute * 60 + epoch.second + epoch.microsecond * 1e-6) / 86400.0
    tut = (jd - 2451545.0) / 36525.0

    # 平恒星时（度）
    Tdeg = 280.46061837 + 360.98564736629 * (jd - 2451545.0) \
        + 0.000387933 * tut * tut - tut ** 3 / 38710000.0
    # 月亮平黄经（度）
    s = 218.3164477 + 481267.88123421 * tut - 0.0015786 * tut * tut
    # 太阳平黄经（度）
    h = 280.46646 + 36000.76983 * tut
    # 月亮近地点平黄经（度）
    p = 83.3532465 + 4069.0137287 * tut

    T = math.radians(Tdeg % 360.0)
    return T, math.radians(s % 360.0), math.radians(h % 360.0), math.radians(p % 360.0)


def _constituentArg(name, T, s, h, p):
    """各分潮幅角（弧度），基于 τ = T − s + h 的 Doodson 展开。"""
    if name == 'M2':
        return 2.0 * (T - s + h)
    if name == 'S2':
        return 2.0 * T
    if name == 'N2':
        return 2.0 * T - 3.0 * s + 2.0 * h + p
    if name == 'K2':
        return 2.0 * T + 2.0 * h
    if name == 'K1':
        return T + h
    if name == 'O1':
        return T - 2.0 * s + h
    if name == 'P1':
        return T - h
    if name == 'Q1':
        return T - 3.0 * s + 2.0 * h + p
    if name == 'MF':
        return 2.0 * s
    if name == 'MM':
        return s - p
    if name == 'SSA':
        return 2.0 * h
    return None


# ---------------------------------------------------------------- 载荷位移
def oceanLoadingDisp(blqData, station, epoch):
    """
    计算测站海潮载荷位移（ENU，米）。
    blqData : readBlq 返回
    station : 测站名（4 字符，自动大写；无记录返回 None）
    返回 np.array([dE, dN, dU])（米）或 None。
    """
    st = blqData.get(station.upper())
    if st is None:
        return None

    T, s, h, p = _astroArgs(epoch)
    amp = st['amp']       # (3,11) radial, ew, ns
    phase = st['phase']   # (3,11) deg

    dU = 0.0
    dE = 0.0
    dN = 0.0
    for i, cname in enumerate(CONSTITUENTS):
        chi = _constituentArg(cname, T, s, h, p)
        if chi is None:
            continue
        c = math.cos(chi - math.radians(phase[0, i]))
        dU += amp[0, i] * c
        dE -= amp[1, i] * math.cos(chi - math.radians(phase[1, i]))  # 正向西
        dN -= amp[2, i] * math.cos(chi - math.radians(phase[2, i]))  # 正向南

    return np.array([dE, dN, dU])
