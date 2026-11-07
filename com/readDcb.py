#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
readDcb.py — 码偏差（DCB/DSB）文件读取与应用。

支持两类文件
------------
1. CODE 月度 P1-C1 DCB 文件（旧格式，如 P1C12205.DCB）
   - readDcb() 自动识别（无 %=BIA 头）；
   - applyDcbSat(): GPS C/A 码伪距归算到 P1 基准。

2. Bias SINEX 文件（MGEX 长文件名格式，如 CAS0MGXRAP_*_DCB.BSX）
   - readDcb() 按 '%=BIA' 头自动识别；
   - 解析 +BIASES/SOLUTION 块中的 DSB（差分码偏差）行；
   - applyDcbSat(): 将任意观测码归算到该系统精密钟差的 IF 组合码基准。

钟差码基准（MGEX/WUM/CODE 约定）：
    GPS     P1/P2     (C1W/C2W)
    GLONASS P1/P2     (C1P/C2P)
    BDS     B1I/B3I   (C2I/C6I)
    Galileo E1/E5a    (C1X/C5X)
    QZSS    L1/L2     (C1C/C2W)

符号约定：DSB(s1,s2) = b(s1) − b(s2)，b 为卫星硬件延迟（直接体现在伪距上）。
把码 s1 的伪距归算到基准码 datum：P_corr = P_obs − [b(s1) − b(datum)]·c
                                        = P_obs − DSB(s1,datum)·c
CODE P1C1 文件的值为 b(P1) − b(C1) = DSB(P1,C1)，故 C1 归算改正 = +值（与旧实现一致）。
"""
import re

CLIGHT_NS = 0.299792458  # ns → m

DCB_LINE_SAT = re.compile(r'^\s*([GRECJI]\d{2})\s+(-?\d+\.\d+)\s+(\d+\.\d+)')
DCB_LINE_REC = re.compile(r'^\s*([A-Z0-9]{4})\s+(-?\d+\.\d+)\s+(\d+\.\d+)')
BSX_SATID = re.compile(r'^[GRECJI]\d{2}$')

# 各系统精密钟差的无电离层组合码基准（按 RINEX 波段号分组）
# 结构: {系统: {波段号: 基准码}}；未列出的波段号回退到最近的已列条目
CLK_DATUM = {
    'G': {1: 'C1W', 2: 'C2W'},            # GPS: P1/P2
    'R': {1: 'C1P', 2: 'C2P'},            # GLONASS: P1/P2
    'C': {2: 'C2I', 6: 'C6I'},            # BDS: B1I/B3I
    'E': {1: 'C1X', 5: 'C5X'},            # Galileo: E1/E5a
    'J': {1: 'C1C', 2: 'C2W'},            # QZSS: L1/L2
}
# 波段号分组回退（同频点组内共用基准槽位）
_DATUM_GROUP = {
    'G': {1: 1, 2: 2, 3: 2, 5: 2},        # GPS L1→P1；L2/L3/L5→P2
    'R': {1: 1, 2: 2, 3: 2},
    'C': {1: 2, 2: 2, 5: 6, 6: 6, 7: 6, 8: 6},   # BDS B1*→B1I；B1C/B2*/B3*→B3I 槽
    'E': {1: 1, 5: 5, 6: 5, 7: 5, 8: 5},  # GAL E1→E1；E5/E6→E5a 槽
    'J': {1: 1, 2: 2, 5: 2},
}


def _band_num(band):
    """C1C → 1，C6I → 6。"""
    try:
        return int(band[1])
    except (IndexError, ValueError):
        return 0


def _datum_code(sys, band):
    """返回该系统/波段对应的钟差基准码（无定义返回 None）。"""
    tab = CLK_DATUM.get(sys)
    grp = _DATUM_GROUP.get(sys)
    if not tab or not grp:
        return None
    n = _band_num(band)
    if n in grp:
        return tab.get(grp[n])
    return None


def _readBiasSinex(dcbFile):
    """
    解析 Bias SINEX（BSX）文件中的 DSB（差分码偏差）行。
    注：部分产品（如 CAS MGEX）DSB 行不带 +BIASES/SOLUTION 块标记，
    故直接全文扫描行首为 DSB 的数据行。
    返回 {'type': 'DSB-BSX', 'sat': {}, 'rec': {}, 'dsb': {(prn, s1, s2): value_ns}}
    """
    out = {'type': 'DSB-BSX', 'sat': {}, 'rec': {}, 'dsb': {}}
    with open(dcbFile, 'r', encoding='utf-8', errors='ignore') as f:
        for line in f:
            fields = line.split()
            if not fields or fields[0] not in ('DSB', 'OSB'):
                continue
            # DSB <SVN|-> <SATID> <sig1> <sig2> <start> <end> <unit> <value> [std]
            if len(fields) < 9:
                continue
            svn, satid, s1, s2, unit, val = \
                fields[1], fields[2], fields[3], fields[4], fields[7], fields[8]
            if not BSX_SATID.match(satid):
                continue  # 测站偏差（接收机端被钟差参数吸收）
            try:
                v = float(val)
            except ValueError:
                continue
            if unit == 'm':
                v /= CLIGHT_NS  # 统一为 ns
            elif unit != 'ns':
                continue
            # OSB 为绝对偏差、DSB 为相对偏差，仅使用 DSB（IF 组合归算需要相对量）
            if fields[0] == 'OSB':
                continue
            out['dsb'][(satid, s1, s2)] = v
            out['sat'][satid] = 1.0  # 兼容旧统计字段
    return out


def _isBiasSinex(dcbFile):
    up = dcbFile.upper()
    if up.endswith('.BSX'):
        return True
    try:
        with open(dcbFile, 'r', encoding='utf-8', errors='ignore') as f:
            for _ in range(20):
                line = f.readline()
                if not line:
                    break
                if line.startswith('%=BIA'):
                    return True
    except OSError:
        return False
    return False


def readDcb(dcbFile):
    """
    解析 DCB / Bias SINEX 文件（按文件头自动识别格式）。

    返回：
        CODE  格式: {'type': 'P1-C1', 'sat': {prn: ns}, 'rec': {station: ns}}
        BSX   格式: {'type': 'DSB-BSX', 'sat': {prn: 1.0},
                     'rec': {}, 'dsb': {(prn, sig1, sig2): ns}}
    """
    if _isBiasSinex(dcbFile):
        return _readBiasSinex(dcbFile)

    out = {'type': 'P1-C1', 'sat': {}, 'rec': {}}
    inRecSection = False
    with open(dcbFile, 'r', encoding='utf-8', errors='ignore') as f:
        for line in f:
            up = line.upper()
            if 'FOR SATELLITES AND RECEIVERS' in up:
                inRecSection = False
                continue
            if 'BIASES FOR STATION' in up or 'FOR RECEIVERS' in up:
                inRecSection = True
                continue
            if line.strip().startswith('***'):
                continue

            m = DCB_LINE_SAT.match(line)
            if m and not inRecSection:
                out['sat'][m.group(1)] = float(m.group(2))
                continue
            if inRecSection:
                m2 = DCB_LINE_REC.match(line)
                if m2 and not re.match(r'^\s*[GRECJI]\d{2}\s', line):
                    out['rec'][m2.group(1)] = float(m2.group(2))
    return out


def _dsb_lookup(dsb, prn, s1, s2):
    """查 DSB(s1,s2)；若键序相反则取负。返回 ns 或 None。"""
    v = dsb.get((prn, s1, s2))
    if v is not None:
        return v
    v = dsb.get((prn, s2, s1))
    if v is not None:
        return -v
    return None


def _dsb_chain(dsb, prn, s1, s2):
    """经中间码 X 链式归算：DSB(s1,X) − DSB(s2,X)。返回 ns 或 None。"""
    xs = set()
    for (sat, a, b) in dsb:
        if sat != prn:
            continue
        if a not in (s1, s2):
            xs.add(a)
        if b not in (s1, s2):
            xs.add(b)
    best = None
    for x in xs:
        v1 = _dsb_lookup(dsb, prn, s1, x)
        v2 = _dsb_lookup(dsb, prn, s2, x)
        if v1 is None or v2 is None:
            continue
        diff = v1 - v2
        if best is None or abs(diff) < abs(best):
            best = diff
    return best


def applyDcbSat(dcbData, prn, band):
    """
    返回把该观测码伪距归算到系统钟差基准所需的改正量（米，直接加到伪距上）。

    - dcbData 为 None              → 0
    - CODE P1-C1 格式              → 旧逻辑（仅 GPS C/A 码，+DCB）
    - DSX-BSX 格式                 → −DSB(band, datum)·c，找不到时尝试链式归算
    - 波段本身即基准码（如 GPS C1W）→ 0
    """
    if dcbData is None:
        return 0.0
    sys = prn[0] if prn else ''
    band = band.upper()

    if dcbData.get('type') == 'DSB-BSX':
        dsb = dcbData.get('dsb', {})
        if not dsb:
            return 0.0
        datum = _datum_code(sys, band)
        if datum is None or band == datum:
            return 0.0
        v = _dsb_lookup(dsb, prn, band, datum)
        if v is None:
            v = _dsb_chain(dsb, prn, band, datum)
        if v is None:
            return 0.0
        return -v * CLIGHT_NS

    # 旧格式：CODE P1-C1（GPS C/A 码归算到 P1 基准）
    if not prn.startswith('G'):
        return 0.0
    if not band.startswith('C1') or band in ('C1W', 'C1P'):
        return 0.0
    d = dcbData['sat'].get(prn)
    if d is None:
        return 0.0
    return d * CLIGHT_NS


# 兼容旧接口名
def applyDcbP1C1(dcbData, prn, band):
    """旧接口：GPS C/A 码归算 P1 基准（现自动路由到通用实现）。"""
    return applyDcbSat(dcbData, prn, band)


def dcb_stat(dcbData):
    """返回 (卫星条目数, 格式类型) 用于加载日志。"""
    if dcbData is None:
        return 0, 'None'
    if dcbData.get('type') == 'DSB-BSX':
        return len(dcbData.get('dsb', {})), dcbData['type']
    return len(dcbData.get('sat', {})), dcbData.get('type', '?')
