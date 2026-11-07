#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
readErp.py — IGS ERP（地球定向参数）文件读取与应用。

能力
----
- readErp() : 解析 IGS ERP（version 1/2），支持 -- 所有数据列
- getErp()  : 线性内插任意历元的极移 xp/yp 与 UT1-UTC

应用
----
- UT1-UTC 用于 GMST 计算（太阳/月亮 ECEF 方向、相位缠绕参考系），
  当前 UTC 替代 UT1 引入的恒星时误差 ~0.1 s ≈ 1.5″；
- xp/yp（角秒）可用于极移改正（本引擎暂未用于测站坐标变换，保留接口）。
"""
import numpy as np


def readErp(erpFile):
    """
    解析 ERP 文件。
    返回 {
      'mjd':  np.array float,   # MJD（v2 文件带 .50 等小数）
      'xp':   np.array arcsec,  # 极移 X（原值 10^-6 角秒）
      'yp':   np.array arcsec,
      'dut1': np.array s,       # UT1-UTC（原值 0.1 µs）
      'lod':  np.array s,       # 日长变化（原值 0.1 µs/day）
    }
    """
    mjds, xps, yps, dut1s, lods = [], [], [], [], []
    with open(erpFile, 'r', encoding='utf-8', errors='ignore') as f:
        lines = f.readlines()

    started = False
    for line in lines:
        s = line.strip()
        if not s:
            continue
        low = s.lower()
        if low.startswith('version'):
            continue
        # 数据行以 MJD 开头（如 59700.50）
        try:
            v0 = float(s.split()[0])
        except Exception:
            continue
        if v0 < 13000 or v0 > 90000:   # 不是合理的 MJD
            continue
        toks = s.split()
        try:
            mjd = float(toks[0])
            xp = float(toks[1]) * 1e-6   # 10^-6" → "
            yp = float(toks[2]) * 1e-6
            dut1 = float(toks[3]) * 0.1e-6   # 0.1 µs → s
            lod = float(toks[4]) * 0.1e-6 if len(toks) > 4 else 0.0
        except (ValueError, IndexError):
            continue
        mjds.append(mjd)
        xps.append(xp)
        yps.append(yp)
        dut1s.append(dut1)
        lods.append(lod)
        started = True

    if not started:
        raise ValueError('ERP 文件无有效数据行: %s' % erpFile)
    order = np.argsort(mjds)
    return {
        'mjd': np.array(mjds)[order],
        'xp': np.array(xps)[order],
        'yp': np.array(yps)[order],
        'dut1': np.array(dut1s)[order],
        'lod': np.array(lods)[order],
    }


def getErp(erpData, epoch):
    """
    线性内插历元 epoch（datetime）的 ERP 值。
    返回 {'xp': ", 'yp': ", 'dut1': s, 'lod': s}；超出范围就近外推一格。
    """
    mjd = epoch.toordinal() - 678576 + \
        (epoch.hour * 3600 + epoch.minute * 60 + epoch.second + epoch.microsecond * 1e-6) / 86400.0
    m = erpData['mjd']

    i = np.searchsorted(m, mjd)
    i0 = min(max(i - 1, 0), len(m) - 2)
    i1 = i0 + 1
    if m[i1] == m[i0]:
        a = 0.0
    else:
        a = (mjd - m[i0]) / (m[i1] - m[i0])
    a = min(max(a, 0.0), 1.0)

    return {
        'xp': float(erpData['xp'][i0] * (1 - a) + erpData['xp'][i1] * a),
        'yp': float(erpData['yp'][i0] * (1 - a) + erpData['yp'][i1] * a),
        'dut1': float(erpData['dut1'][i0] * (1 - a) + erpData['dut1'][i1] * a),
        'lod': float(erpData['lod'][i0] * (1 - a) + erpData['lod'][i1] * a),
    }
