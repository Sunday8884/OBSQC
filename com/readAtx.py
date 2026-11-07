#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
readAtx.py — ANTEX 相位中心改正文件（igsXX.atx）读取与应用。

能力
----
- readAntex()      : 解析卫星与接收机天线的 PCO（mm）与 PCV（NOAZI / 方位角网格）
- satAntCorr()     : 卫星 PCO/PCV → 无电离层组合的距离改正量（米）
- recAntCorr()     : 接收机 PCO（NEU）/PCV → 距离改正量（米）

约定
----
- ANTEX 卫星 PCO 给在卫星固定坐标系：NORTH=沿轨 x、EAST=横轨 y、UP=指向地心 z（IGS 约定）；
- 卫星姿态近似：ez 指向地心，ex 沿速度方向（正交化），ey = ez×ex（无偏航机动建模）；
- PCV 仅用 NOAZI（无方位角依赖），卫星按天底角、接收机按天顶角线性内插；
- IF 组合改正量 = (f1²·corr1 − f2²·corr2) / (f1² − f2²)。
"""
import datetime
import numpy as np


# ---------------------------------------------------------------- 波段 → ANTEX 频率标签
def bandToAtxTag(sys, band):
    """
    RINEX 观测码 → ANTEX 频率标签（如 'G01' / 'C06' / 'E05'）。
    未识别返回 None。
    """
    band = band.upper()
    m = {
        'G': {'C1': 'G01', 'C2': 'G02', 'C5': 'G05'},
        'C': {'C2': 'C01', 'C6': 'C06', 'C1': 'C02',
              'C7': 'C07', 'C5': 'C08'},
        'E': {'C1': 'E01', 'C5': 'E05', 'C7': 'E07', 'C8': 'E08', 'C6': 'E06'},
    }
    tab = m.get(sys)
    if tab is None:
        return None
    return tab.get(band[:2])


# ---------------------------------------------------------------- 解析
def _parse_valid(tokens):
    try:
        return datetime.datetime(int(tokens[0]), int(tokens[1]), int(tokens[2]),
                                 int(tokens[3]), int(tokens[4]), int(float(tokens[5])))
    except Exception:
        return None


def _parse_pcv_line(line):
    """NOAZI 或方位角行 → [值(mm)...]；失败返回 None。"""
    parts = line.split()
    try:
        if parts[0].upper() == 'NOAZI':
            return [float(v) for v in parts[1:-1]]  # 末尾字段是注释名
        return [float(v) for v in parts[:-1]]
    except Exception:
        return None


def readAntex(atxFile):
    """
    解析 ANTEX 文件。

    返回
    ----
    {
      'sat': {prn: [rec, ...]},          # rec 按VALID FROM 升序
      'rec': {antType: [rec, ...]},
      'n_sat': int, 'n_rec': int,
    }
    rec = {
      'type': 天线类型, 'serial': 序列号/'NONE',
      'validFrom': datetime|None, 'validUntil': datetime|None,
      'zen1': float, 'dzen': float,      # PCV 天顶角网格（度）
      'freq': {tag: {'NEU': (n,e,u) m, 'pcv': [mm,...] or None}},
    }
    """
    out = {'sat': {}, 'rec': {}, 'n_sat': 0, 'n_rec': 0}
    cur = None          # 当前 antenna 记录
    inFreq = False      # 是否在频率块内

    with open(atxFile, 'r', encoding='utf-8', errors='ignore') as f:
        for line in f:
            label = line[60:].strip()

            if label == 'START OF ANTENNA':
                cur = {'type': None, 'serial': None,
                       'validFrom': None, 'validUntil': None,
                       'zen1': 0.0, 'dzen': 1.0, 'freq': {}}
                inFreq = False
                continue

            if cur is None:
                continue

            if label == 'TYPE / SERIAL NO':
                cur['type'] = line[0:20].strip()
                cur['serial'] = line[20:40].strip() or 'NONE'
                # 兼容部分接收机记录 serial 未对齐第 21 列的写法：
                # 若 serial 为 NONE 但 type 字段内含尾随 token，则拆出
                if cur['serial'] == 'NONE' and ' ' in cur['type']:
                    toks = cur['type'].split()
                    if len(toks) > 1:
                        cur['serial'] = toks[-1]
                        cur['type'] = ' '.join(toks[:-1])
            elif label == 'ZEN1 / ZEN2 / DZEN':
                try:
                    cur['zen1'] = float(line[2:8])
                    cur['dzen'] = float(line[13:19]) if line[13:19].strip() else 1.0
                except Exception:
                    pass
            elif label == 'VALID FROM':
                cur['validFrom'] = _parse_valid(line.split()[:6])
            elif label == 'VALID UNTIL':
                cur['validUntil'] = _parse_valid(line.split()[:6])
            elif label == 'START OF FREQUENCY':
                inFreq = True
                tag = line[3:9].strip()
                cur['freq'][tag] = {'NEU': None, 'pcv': None}
            elif label == 'END OF FREQUENCY':
                inFreq = False
            elif inFreq and label == 'NORTH / EAST / UP':
                try:
                    n, e, u = (float(line[0:10]), float(line[10:20]), float(line[20:30]))
                    # 找到最近一个 START OF FREQUENCY 的 tag
                    tag = list(cur['freq'].keys())[-1]
                    cur['freq'][tag]['NEU'] = (n / 1000.0, e / 1000.0, u / 1000.0)
                except Exception:
                    pass
            elif inFreq and (line.strip().startswith('NOAZI') or
                             (line.strip()[:1].isdigit() and cur.get('_dazi', 0))):
                vals = _parse_pcv_line(line)
                if vals:
                    tag = list(cur['freq'].keys())[-1]
                    cur['freq'][tag]['pcv'] = vals
            elif label == 'END OF ANTENNA':
                isSat = (len(cur['type'].split()) >= 2 and cur['type'].split()[0] in
                         ('BLOCK', 'GPS', 'GLO', 'GALILEO', 'BDS', 'BDS-3', 'QZS', 'QZSS', 'IRNSS', 'SBAS')) \
                    or (len(cur['serial']) == 3 and cur['serial'][1:].isdigit() and cur['serial'][0] in 'GRECJIS')
                if isSat:
                    out['sat'].setdefault(cur['serial'], []).append(cur)
                    out['n_sat'] += 1
                else:
                    out['rec'].setdefault(cur['type'], []).append(cur)
                    out['n_rec'] += 1
                cur = None

    # 卫星记录按有效期排序
    for prn in out['sat']:
        out['sat'][prn].sort(key=lambda r: r['validFrom'] or datetime.datetime(1980, 1, 1))
    return out


def _pick_sat_record(atxData, prn, epoch):
    recs = atxData['sat'].get(prn)
    if not recs:
        return None
    for rec in recs:  # 已排序，取最后一个满足时间窗的
        vf, vu = rec['validFrom'], rec['validUntil']
        if vf and epoch < vf:
            continue
        if vu and epoch > vu:
            continue
        return rec
    return recs[-1]  # 无时间信息时兜底取最后一条


def _pick_rec_record(atxData, antType):
    recs = atxData['rec'].get(antType)
    if not recs:
        return None
    for rec in recs:  # 优先 serial=NONE 的通用记录
        if rec['serial'] == 'NONE':
            return rec
    return recs[0]


def _interpPcv(pcvList, zen1, dzen, angleDeg):
    """PCV 网格线性内插（angleDeg 为天底角/天顶角，度）。"""
    if pcvList is None or len(pcvList) == 0:
        return 0.0
    idxF = (angleDeg - zen1) / dzen
    i0 = int(np.floor(idxF))
    if i0 < 0:
        return float(pcvList[0])
    if i0 >= len(pcvList) - 1:
        return float(pcvList[-1])
    a = idxF - i0
    return float(pcvList[i0] * (1.0 - a) + pcvList[i0 + 1] * a)


# 频率标签回退：igs20 中 BDS-2（GEO/IGSO 及部分 MEO）的 B1I 标为 C02 而非 C01
_TAG_FALLBACK = {
    ('C', 'C01'): 'C02',   # BDS-2 B1I(1561.098MHz) → igs20 记录标签
}


def _freq_of(rec, sys, tag):
    """取记录中指定标签的频率数据，必要时按回退表替换标签。"""
    fr = rec['freq'].get(tag)
    if fr is not None:
        return fr
    alt = _TAG_FALLBACK.get((sys, tag))
    if alt:
        return rec['freq'].get(alt)
    return None


def satAntCorr(atxData, prn, epoch, sys, tag1, f1, tag2, f2, nadirDeg):
    """
    卫星天线相位中心改正（无电离层组合，米）。

    tag1/tag2 : 两个波段对应的 ANTEX 频率标签
    f1/f2     : 两波段频率（Hz）
    nadirDeg  : 天底角（度）

    返回 (pcoVec, pcv)：
        pcoVec : PCO 向量（卫星固定系 x/y/z，米；x=沿轨, y=横轨, z=指向地心）
        pcv    : IF 组合 PCV（米，该天底角处）
    卫星无记录时返回 (None, None)。
    """
    rec = _pick_sat_record(atxData, prn, epoch)
    if rec is None:
        return None, None
    fr1 = _freq_of(rec, sys, tag1)
    fr2 = _freq_of(rec, sys, tag2)
    if fr1 is None or fr2 is None or fr1['NEU'] is None or fr2['NEU'] is None:
        return None, None

    def _comb(attr, valFn):
        a1, a2 = valFn(fr1), valFn(fr2)
        return (f1 ** 2 * a1 - f2 ** 2 * a2) / (f1 ** 2 - f2 ** 2)

    pco = np.array([_comb(0, lambda fr: fr['NEU'][i]) for i in range(3)])
    pcv = 0.0
    if fr1['pcv'] and fr2['pcv']:
        pcv = (f1 ** 2 * _interpPcv(fr1['pcv'], rec['zen1'], rec['dzen'], nadirDeg) -
               f2 ** 2 * _interpPcv(fr2['pcv'], rec['zen1'], rec['dzen'], nadirDeg)) \
            / (f1 ** 2 - f2 ** 2)
    # NEU 解析时已换算为米；PCV 仍为 mm 需 /1000
    return pco, pcv / 1000.0


def recAntCorr(atxData, antType, sys, tag1, f1, tag2, f2, zenithDeg):
    """
    接收机天线相位中心改正（无电离层组合，米）。

    antType   : RINEX 头 ANT # / TYPE 中的天线类型
    zenithDeg : 观测天顶角（度，= 90° − 高度角）

    返回 (pcoNEU, pcv)：pcoNEU = (N, E, U)（米）；无记录返回 (None, None)。
    """
    rec = _pick_rec_record(atxData, antType)
    if rec is None:
        return None, None
    fr1 = _freq_of(rec, sys, tag1)
    fr2 = _freq_of(rec, sys, tag2)
    if fr1 is None or fr2 is None or fr1['NEU'] is None:
        return None, None

    pcoNEU = np.array([(f1 ** 2 * fr1['NEU'][i] - f2 ** 2 * (fr2['NEU'][i] or 0.0))
                       / (f1 ** 2 - f2 ** 2) for i in range(3)])
    pcv = 0.0
    if fr1['pcv'] and fr2['pcv']:
        pcv = (f1 ** 2 * _interpPcv(fr1['pcv'], rec['zen1'], rec['dzen'], zenithDeg) -
               f2 ** 2 * _interpPcv(fr2['pcv'], rec['zen1'], rec['dzen'], zenithDeg)) \
            / (f1 ** 2 - f2 ** 2)
    return pcoNEU, pcv / 1000.0
