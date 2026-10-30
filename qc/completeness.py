# -*- coding: utf-8 -*-
"""
卫星数据完整率（Data Completeness）统计
================================================

【为什么需要这个指标】
现有 freq 模块给出的是「某频段观测次数 / 文件总历元数」。这个值对移动卫星
会被"卫星本就不在该时段可见"严重稀释——例如 G01 的 L1C 只有 670/2880 = 23%，
但这 23% 里绝大部分是卫星本就在地平线以下，并无任何数据质量问题。
本模块以**卫星可见弧段（presence arc）**为基准计算完整率，才反映真实的数据缺口。

【指标定义】
对每颗卫星 PRN 及其每个选定频段 band（**载波观测码全称**，如 'L2I'）：

  presence arc（可见弧段）
      = 该 PRN **实际有值**的历元集合，按连续段切分后各段历元数之和。
      PRN「有值」定义：该历元中该 PRN 至少有一个观测码的数值非 None。
      弧段之外的历元视为「卫星不可见」，不计入缺测。
      （这是本模块与 freq 模块的关键区别）

      **为什么按连续段切分**：卫星可能在时段中途长时间消失（接收机换星座、
      遮挡、信号中断等）后又重新出现。若仍按「首个有值 ~ 末个有值」取弧段，
      中间的长时间不可见期会被计入分母，把完整率严重拉低——实测 C19 实际
      有值 1165 历元，却因跨越 878 历元的中断而只算出 56.8%。
      切分后其完整率为 99% 左右，才反映真实的数据质量。

  频段分组粒度 = 频点号 + 跟踪模式（band key，如 '2I'），**不跨跟踪模式合并**。
  例如 GPS 的 L1C / L1W / L1X 分列为三行，各自独立统计完整率。
  同组内含该 band key 的全部观测码（伪距 C / 伪距 P / 载波 L / 载噪比 S /
  多普勒 D），判定「有值」时按 C > P > L > S > D 优先级取第一个非 None 值。

  N_arc      = 弧段内历元总数                              （分母）
  N_band     = 弧段内该频段有非 None 值的历元数
  R_arc      = N_band / N_arc                             → 【数据完整率】主指标
  N_present  = 弧段内该 PRN 至少一个观测码有值的历元数
  R_present = N_band / N_present                         → 频段层面缺测率
  R_total    = N_band / total_epochs                     → 相对整个观测时段的数据覆盖率

  N_span     = 该频段自身首个有值历元 ~ 末个有值历元之间的历元总数
  R_span     = N_band / N_span                           → 该频段数据自身的连续完整性
               （与 R_arc 互补：R_arc 衡量「相对整星可见期」，R_span 衡量
                「该频段一旦开始后是否连续」）

  gap        = 弧段内连续缺测段（连续 ≥ gap_min_epochs 个历元该频段无值）
  gap_count  = gap 段数
  max_gap    = 最长 gap 的历元数 / 时长

【可选：理论可见性（需广播星历）】
  N_visible  = 理论可见历元数（卫星高度角 ≥ cutoff 且星历健康可用）
  R_visible  = N_band / N_visible
  n_vis_miss = 理论可见但实际无数据的历元数

  注意两点：
  1. 星历不可用（无对应星历 / SVhealth 异常）的历元排除在分母之外，
     否则会把「星历问题」误判成「数据缺失」；
  2. 广播星历间隔大（常为 2 小时）且部分卫星 SVhealth 异常，
     插值精度有限，R_visible 允许略微超过 100%（此时说明星历口径偏保守）。

【数据结构】
completenessData = {
    'G': {
        'G01': {
            'L1C': {
                'obs': 670,            # 该频段有值历元数
                'present': 698,        # 该 PRN 有值历元数
                'arc': 698,            # 可见弧段历元数
                'total': 2880,         # 文件总历元数
                'ratio_arc': 0.9599,   # 数据完整率（主指标）
                'ratio_present': 0.9599,
                'ratio_total': 0.2326,
                'gap_count': 3,
                'max_gap': 28,
                'max_gap_sec': 840.0,
                'visible': None,       # 有星历时为理论可见历元数
                'ratio_visible': None,
                'vis_miss': None,
                'nav_skipped': None,   # 因星历不可用被跳过的历元数
            },
            ...
        },
        ...
    },
    ...
}
"""

import datetime


# ========================= 内部工具 =========================
def _band_key(code):
    """观测码 -> 频段键（去掉观测类型前缀）。'L1C' -> '1C'，'C1C' -> '1C'"""
    if len(code) >= 3 and code[0] in ('L', 'C', 'P', 'D', 'S'):
        return code[1:]
    return code


# GLONASS 的中心频率（obsFreq/obsCodes 静态表中无 R 条目，此处单独补充）
# L1C/L1P = L1 (1602 MHz)，L2C/L2P = L2 (1246 MHz)，L3X = L3 (1202.025 MHz)
_GLONASS_FREQ = {
    'L1C': 1602.0, 'L1P': 1602.0,
    'L2C': 1246.0, 'L2P': 1246.0,
    'L3X': 1202.025,
}


def _band_center_freq(gSys, band_code):
    """
    频段键 -> 中心频率 (MHz)，取不到返回 None。

    band_code 为完整载波观测码，如 'L2I'。按 obsCodes/obsFreq 静态表反查，
    例如 BDS：L1X -> B1C 1575.42、L2I -> B1 1561.098、L6I -> B3 1268.52、
    L7I -> B2 1207.14、L7D/L7Z -> B2b 1207.14、L8X -> B2a+B2b 1191.795。
    GLONASS 静态表无条目，用内置映射兜底。
    """
    try:
        from com.gnssParameter import obsCodes, obsFreq
    except Exception:
        obsCodes, obsFreq = {}, {}
    for band, codes in (obsCodes.get(gSys) or {}).items():
        if band_code in codes:
            return (obsFreq.get(gSys) or {}).get(band)
    if gSys == 'R':
        return _GLONASS_FREQ.get(band_code)
    return None


def _collect_bands(obsType, gSys, band_choose):
    """
    确定某系统需要统计的频段，返回 {bandName: [观测码列表]}。

    bandName 统一为**载波观测码全称**（如 'L1C' / 'L2I' / 'L7D'），与 RINEX
    观测码书写一致，便于直接与观测头对照。

    分组粒度 = 频点号 + 跟踪模式（即 band key，如 '2I'），**不跨跟踪模式合并**。
    同一 band key 下的其他观测码（伪距 C / 伪距 P / 载波 L / 载噪比 S /
    多普勒 D）归入同一组，判定「有值」时按 C > P > L > S > D 的优先级
    取第一个非 None 的值。

    两种指定方式：
    1. 用户指定 band_choose（如 {'G': ['1C', '2W']}）时，严格按用户给定的
       band key 分组，bandName 取 'L' + band key（'2I' -> 'L2I'）。
    2. 自动模式（band_choose 为 None/空）时，遍历观测头中所有载波观测码，
       每个不同的 band key 生成一个条目（'L1C' / 'L1W' / 'L2W' ...）。
       不做跨跟踪模式合并，因此同一频点的不同跟踪模式会各自成行，
       例如 GPS 的 L1C / L1W / L1X 分列为三行。
    """
    if gSys not in obsType:
        return {}

    codes = obsType[gSys]

    if band_choose:
        wanted = band_choose.get(gSys) or []
        keys = []
        for code in codes:
            k = _band_key(code)
            if k in wanted and k not in keys:
                keys.append(k)
        # 用户指定的频段键可能观测头中不存在，保留以示"完全缺失"
        for k in wanted:
            if k not in keys:
                keys.append(k)
    else:
        # 自动模式：以载波观测码为锚点，逐个 band key 生成条目
        keys = []
        for code in codes:
            if not code.startswith('L'):
                continue
            k = _band_key(code)
            if k not in keys:
                keys.append(k)

    prio = {'C': 0, 'P': 1, 'L': 2, 'S': 3, 'D': 4}
    result = {}
    for k in keys:
        group = [c for c in codes if _band_key(c) == k]
        group.sort(key=lambda c: (prio.get(c[0], 9), codes.index(c)))
        # bandName 用载波观测码全称；无载波码时退化为原 key
        result['L' + k] = group
    return result


def _has_value(obs_rec, code_group):
    """该频段在某个历元是否有有效值：按优先级取首个非 None 的码。"""
    for code in code_group:
        if obs_rec.get(code) is not None:
            return True
    return False


def _prn_has_any(obs_rec):
    """该历元该 PRN 是否有任一非 None 观测值"""
    for code, val in obs_rec.items():
        if code.endswith('LLI'):
            continue
        if val is not None:
            return True
    return False


def _visible_epochs(nav_data, obs_head, prn, epochs, cutoff_deg):
    """
    计算某卫星的理论可见历元集合。

    返回 (visible_set, skipped)：
      visible_set —— 高度角 ≥ cutoff_deg 且星历可用的历元集合
      skipped     —— 因星历不可用而跳过的历元数（不计入分母）
    """
    visible = set()
    skipped = 0

    if not nav_data or prn not in nav_data:
        return visible, len(epochs)

    try:
        import math
        import numpy as np
        from com.nav2posclk import getPosInNav
        from com.xyz2blh import xyz2blh, bl2enuRot
    except Exception:
        return visible, len(epochs)

    approx = obs_head.get('Approx Position') or [None, None, None]
    if approx[0] is None:
        return visible, len(epochs)

    try:
        lat, lon, _h = xyz2blh(approx[0], approx[1], approx[2])
        enuRot = bl2enuRot(lat, lon)
    except Exception:
        return visible, len(epochs)

    rr = np.array([approx[0], approx[1], approx[2]])
    cut_rad = math.radians(cutoff_deg)

    for ep in epochs:
        try:
            rs = getPosInNav(nav_data, prn, ep, ep)
        except Exception:
            rs = (None,) * 6
        if rs[0] is None:
            skipped += 1
            continue
        try:
            d = np.array(rs[:3]) - rr
            east, north, up = enuRot @ d
            el = math.atan(up / math.sqrt(north * north + east * east))
        except Exception:
            skipped += 1
            continue
        if el >= cut_rad:
            visible.add(ep)

    return visible, skipped


# ========================= 主计算 =========================
def _split_present_segments(epochs, present_set, max_gap_epochs=1):
    """
    将卫星的「有值历元」按连续段切分。

    卫星可能在时段中途长时间消失（接收机换星座、遮挡、信号中断等），
    之后又重新出现。此时若仍按「首个有值 ~ 末个有值」取弧段，
    中间的长时间不可见期会被计入分母，把完整率严重拉低
    （例如 C19 可见 1165 历元、却因跨越 878 历元的中断而只算出 57%）。

    参数
    ----
    epochs          : 全部历元（时间序）
    present_set     : 该星有值的历元集合
    max_gap_epochs  : 允许的最大间断历元数，超过则切分为两段

    返回
    ----
    [(段内历元列表), ...]，按时间顺序。
    """
    if not present_set:
        return []
    present_eps = sorted(present_set)
    segments = []
    cur = [present_eps[0]]
    for ep in present_eps[1:]:
        # 两历元在 epochs 中的索引间隔
        try:
            gap = epochs.index(ep) - epochs.index(cur[-1]) - 1
        except ValueError:
            gap = 0
        if gap > max_gap_epochs:
            segments.append(cur)
            cur = [ep]
        else:
            cur.append(ep)
    segments.append(cur)
    return segments


def computeCompleteness(obsHead, obsData, start_datetime=None, end_datetime=None,
                        sat_list=None, sys_list=None, band_choose=None,
                        nav_data=None, cutoff=10.0, gap_min_epochs=2,
                        keep_empty_bands=False):
    """
    计算各卫星各频段的数据完整率。

    Parameters
    ----------
    obsHead, obsData : dict
        观测文件头与观测数据（readObsHead / readObs 返回值）
    start_datetime, end_datetime : datetime, optional
        时间范围过滤（闭区间）
    sat_list : list of str, optional
        卫星 PRN 列表，如 ['G01', 'C06']；None 表示全部
    sys_list : list of str, optional
        GNSS 系统列表，如 ['G', 'C']；None 表示全部
    band_choose : dict, optional
        频段键字典，如 {'G': ['1C', '2W']}；None 表示自动取全部载波频段。
        与 qc_common.band_choose_keys 同一格式，故 --bands 可直接约束本模块。
    nav_data : dict, optional
        广播星历数据；提供时可额外计算「理论可见」口径的完整率
    cutoff : float
        理论可见性判定的截止高度角（度），默认 10
    gap_min_epochs : int
        判定为「数据缺口」所需的最小连续缺测历元数，默认 2
    keep_empty_bands : bool
        是否保留「整个时段无任何数据」的频段（默认 False，即自动剔除）。
        RINEX3 MIXED 的观测头会声明全部可能类型，但接收机实际只提供其中一部分，
        例如 wuh20050.24o 中 BDS 各星的 L1X/C5X 全部无数据，
        若保留会产生大量 0% 的假告警。

    Returns
    -------
    completenessData : dict
        {gSys: {prn: {band_key: {...}}}}
    """
    # ---------- 时间范围内的历元序列 ----------
    all_epochs = sorted(obsData.keys())
    epochs = []
    for ep in all_epochs:
        if start_datetime and ep < start_datetime:
            continue
        if end_datetime and ep > end_datetime:
            continue
        epochs.append(ep)
    total_epochs = len(epochs)
    if total_epochs == 0:
        return {}

    # 采样间隔（用于 gap 时长换算）
    if len(epochs) >= 2:
        interval = (epochs[1] - epochs[0]).total_seconds()
    else:
        interval = 30.0

    # ---------- 卫星列表 ----------
    if sat_list is None:
        prn_list = list(obsHead.get('prn') or [])
    else:
        prn_list = list(sat_list)
    if sys_list:
        prn_list = [p for p in prn_list if p[0] in sys_list]

    obsType = obsHead.get('OBS TYPES') or {}

    # 每个 PRN 的「有值历元」序列（任一观测码非 None）
    prn_present = {}
    for ep in epochs:
        rec = obsData[ep]
        for prn in rec:
            if _prn_has_any(rec[prn]):
                prn_present.setdefault(prn, []).append(ep)

    completenessData = {}

    for prn in prn_list:
        present = prn_present.get(prn)
        if not present:
            continue
        gSys = prn[0]

        bands = _collect_bands(obsType, gSys, band_choose)
        if not bands:
            continue

        # ---- 可见弧段 ----
        # 卫星可能中途长时间消失后重新出现，此时按「首个有值 ~ 末个有值」
        # 取弧段会把中间的不可见期计入分母，严重拉低完整率。
        # 故按连续可见段切分，分母 = 各可见段内的历元数之和。
        segments = _split_present_segments(epochs, set(present))
        if not segments:
            continue
        seg_epoch_sets = [set(s) for s in segments]
        arc_epochs = [ep for ep in epochs if any(ep in s for s in seg_epoch_sets)]
        arc_n = len(arc_epochs)
        if arc_n == 0:
            continue

        present_set = set(present)

        # ---- 理论可见历元（可选）----
        visible = None
        nav_skipped = None
        if nav_data:
            visible, nav_skipped = _visible_epochs(nav_data, obsHead, prn,
                                                   arc_epochs, cutoff)

        sys_dict = completenessData.setdefault(gSys, {}).setdefault(prn, {})

        for band_key, code_group in bands.items():
            # 该频段在观测头中确实存在的观测码
            exist_codes = [c for c in code_group if c in (obsType.get(gSys) or [])]
            if not exist_codes:
                continue

            # 逐历元判定该频段是否有值
            band_mark = {}
            for ep in arc_epochs:
                obs_rec = obsData[ep].get(prn)
                band_mark[ep] = (obs_rec is not None
                                 and _has_value(obs_rec, exist_codes))

            n_band = sum(1 for ok in band_mark.values() if ok)

            # 全时段无数据的频段：默认剔除（观测头声明 ≠ 实际提供）
            if n_band == 0 and not keep_empty_bands:
                continue

            # ---- 连续缺测段（gap）----
            gap_count = 0
            cur_gap = 0
            max_gap = 0
            for ok in band_mark.values():
                if ok:
                    cur_gap = 0
                else:
                    cur_gap += 1
                    if cur_gap > max_gap:
                        max_gap = cur_gap
                    if cur_gap == gap_min_epochs:
                        gap_count += 1

            # 理论可见但无数据的历元数
            n_vis_miss = None
            ratio_visible = None
            if visible is not None:
                n_vis_miss = sum(1 for ep, ok in band_mark.items()
                                  if (not ok) and ep in visible)
                if len(visible) > 0:
                    ratio_visible = n_band / len(visible)

            n_present = sum(1 for ep in arc_epochs if ep in present_set)

            # 该频段自身的跨度（首个有值历元 ~ 末个有值历元）
            # 同样按可见段切分，避免跨越中断期
            band_eps = [ep for ep in arc_epochs if band_mark[ep]]
            if band_eps:
                band_seg_sets = _split_present_segments(
                    arc_epochs, set(band_eps))
                span_n = sum(len(s) for s in band_seg_sets)
                ratio_span = n_band / span_n if span_n else 0.0
            else:
                span_n = 0
                ratio_span = 0.0

            sys_dict[band_key] = {
                'obs': n_band,
                'present': n_present,
                'arc': arc_n,
                'span': span_n,
                'total': total_epochs,
                'ratio_arc': n_band / arc_n,
                'ratio_present': (n_band / n_present) if n_present else 0.0,
                'ratio_span': ratio_span,
                'ratio_total': n_band / total_epochs,
                'gap_count': gap_count,
                'max_gap': max_gap,
                'max_gap_sec': max_gap * interval,
                'visible': len(visible) if visible is not None else None,
                'ratio_visible': ratio_visible,
                'vis_miss': n_vis_miss,
                'nav_skipped': nav_skipped,
            }

    return completenessData


# ========================= 输出 =========================
def writeCompleteness(compData, outFile, total_epochs=0):
    """
    写出数据完整率统计文本，并返回参与 QCValue 汇总的字典。

    输出结构：
        +SYS
         gSys  BandNum   ObsNum   ArcNum   Ratio(%)   GapCnt   MaxGap
        +SAT
          PRN  Band   Freq(MHz)  ObsNum  ArcNum  SpanNum  Arc(%)  Span(%)  GapCnt  MaxGap  MaxGapSec

    PRN	         卫星号
    Band	     载波观测码全称，如 L2I、L6I、L7D
    Freq(MHz)	该观测码的中心频率
    ObsNum	弧段内该频段有值的历元数	分子
    ArcNum	该卫星的可见弧段历元数	首个有值历元~末个有值历元之间，分母
    SpanNum	该频段自身首末有值历元之间的历元数	频段级跨度
    Arc(%)	ObsNum / ArcNum	    衡量"相对整星可见期"的数据完整率
    Span(%)	ObsNum / SpanNum	衡量"该频段一旦开始后是否连续"
    GapCnt	连续缺测段数	连续 ≥2 历元无值记 1 段
    MaxGap	最长缺测段的历元数	—
    MaxGapSec	最长缺测段的秒数	MaxGap × 采样间隔
    """
    summary_dict = {}
    with open(outFile, 'w+') as f:
        f.write('# Inf       : Satellite Data Completeness [%]\n')
        f.write('# Time      : ' + datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S') + '\n')
        f.write('# Ratio     : ObsNum / ArcNum  (弧段内该频段有值历元数 / 可见弧段历元数)\n')
        f.write('# ArcNum    : 该卫星首个有值历元 ~ 末个有值历元之间的历元总数\n')
        f.write('# Band      : 载波观测码全称（如 L2I），同频点不同跟踪模式不合并\n')
        f.write('# Freq      : 该观测码对应的中心频率 (MHz)，取自 obsFreq 静态表\n')
        f.write('#             END OF HEADER\n\n')

        # ---- 测站级：按系统汇总 ----
        sys_stat = {}
        for gSys in compData:
            agg = {'obs': 0, 'arc': 0, 'ratio_sum': 0.0, 'n': 0,
                   'gap_count': 0, 'max_gap': 0}
            for prn in compData[gSys]:
                for band in compData[gSys][prn]:
                    r = compData[gSys][prn][band]
                    agg['obs'] += r['obs']
                    agg['arc'] += r['arc']
                    agg['ratio_sum'] += r['ratio_arc']
                    agg['n'] += 1
                    agg['gap_count'] += r['gap_count']
                    agg['max_gap'] = max(agg['max_gap'], r['max_gap'])
            if agg['n']:
                agg['ratio_mean'] = agg['ratio_sum'] / agg['n']
                sys_stat[gSys] = agg

        f.write('+SYS\n')
        f.write(' gSys  BandNum   ObsNum   ArcNum   Ratio(%)   GapCnt   MaxGap\n')
        f.write('--------------------------------------------------------------\n')
        for gSys in sorted(sys_stat):
            a = sys_stat[gSys]
            f.write(f'   {gSys:<4}{a["n"]:>8d}{a["obs"]:>9d}{a["arc"]:>9d}'
                    f'{a["ratio_mean"] * 100:>11.2f}{a["gap_count"]:>9d}{a["max_gap"]:>9d}\n')
        f.write('--------------------------------------------------------------\n')
        f.write('-SYS\n\n')

        # ---- 卫星级 ----
        f.write('+SAT\n')
        f.write('   PRN     Band  Freq(MHz)   ObsNum   ArcNum  SpanNum   Arc(%)  Span(%)'
                '   GapCnt   MaxGap   MaxGapSec\n')
        f.write('---------------------------------------------------------------------------------\n')
        for gSys in sorted(compData):
            for prn in sorted(compData[gSys]):
                for band in sorted(compData[gSys][prn]):
                    r = compData[gSys][prn][band]
                    freq = _band_center_freq(gSys, band)
                    freq_str = f'{freq:.3f}' if freq else '-'
                    f.write(f'  {prn:<5}{band:>7}{freq_str:>12}'
                            f'{r["obs"]:>9d}{r["arc"]:>9d}{r["span"]:>9d}'
                            f'{r["ratio_arc"] * 100:>8.2f}{r["ratio_span"] * 100:>9.2f}'
                            f'{r["gap_count"]:>9d}'
                            f'{r["max_gap"]:>9d}{r["max_gap_sec"]:>12.1f}\n')
                    # 参与 QCValue 汇总：主指标为弧段完整率
                    summary_dict.setdefault('completeness', {}) \
                                 .setdefault(prn, {})[band] = {'val': r['ratio_arc']}
        f.write('---------------------------------------------------------------------------------\n')
        f.write('-SAT\n')

    return summary_dict


def writeCompletenessDetailed(compData, obsData, obsHead, detailFile,
                              band_choose=None, gap_min_epochs=2,
                              start_datetime=None, end_datetime=None):
    """
    写出逐历元的完整率标记 CSV（宽表），用于缺口可视化。

        行 = 历元
        列 = PRN_BAND，值 = 1（有值）/ 0（缺测）/ 空（该历元该星不可见）

    同时返回各卫星各频段的连续缺测段明细：
        {gSys: {prn: {band: [{'start','end','n','sec'}, ...]}}}

    band_choose 必须与调用 computeCompleteness 时一致，否则频段分组会不匹配。
    """
    import csv

    obsType = obsHead.get('OBS TYPES') or {}

    epochs = []
    for ep in sorted(obsData.keys()):
        if start_datetime and ep < start_datetime:
            continue
        if end_datetime and ep > end_datetime:
            continue
        epochs.append(ep)

    if len(epochs) >= 2:
        interval = (epochs[1] - epochs[0]).total_seconds()
    else:
        interval = 30.0

    # 收集 (prn, band) 与对应观测码组（与 computeCompleteness 用同一套分组逻辑）
    entries = []   # [(prn, band, exist_codes)]
    for gSys in compData:
        bands = _collect_bands(obsType, gSys, band_choose)
        for prn in compData[gSys]:
            for band in compData[gSys][prn]:
                exist = [c for c in bands.get(band, [])]
                if exist:
                    entries.append((prn, band, exist))
    entries.sort(key=lambda x: (x[0], x[1]))

    # 逐历元逐条目计算标记
    marks = {}    # (prn,band) -> {epoch: 0/1}
    gaps = {}     # 缺测段明细
    for prn, band, exist in entries:
        gSys = prn[0]
        present_eps = []
        for ep in epochs:
            rec = obsData[ep].get(prn)
            if rec is not None and _prn_has_any(rec):
                present_eps.append(ep)
        if not present_eps:
            continue
        # 弧段按连续可见段切分（与 computeCompleteness 一致）：
        # 卫星中途长时间消失后再出现时，中间不可见期不计入分母
        segments = _split_present_segments(epochs, set(present_eps))
        seg_sets = [set(s) for s in segments]
        arc = [ep for ep in epochs if any(ep in s for s in seg_sets)]

        m = {}
        cur_start = None
        cur_n = 0
        segs = []
        for ep in arc:
            rec = obsData[ep].get(prn)
            ok = rec is not None and _has_value(rec, exist)
            m[ep] = 1 if ok else 0
            if not ok:
                if cur_start is None:
                    cur_start = ep
                    cur_n = 1
                else:
                    cur_n += 1
            else:
                if cur_start is not None and cur_n >= gap_min_epochs:
                    segs.append((cur_start, ep, cur_n))
                cur_start = None
                cur_n = 0
        if cur_start is not None and cur_n >= gap_min_epochs:
            segs.append((cur_start, arc[-1], cur_n))

        marks[(prn, band)] = m
        if segs:
            gl = gaps.setdefault(gSys, {}).setdefault(prn, {}).setdefault(band, [])
            for s, e, n in segs:
                gl.append({'start': s, 'end': e, 'n': n, 'sec': n * interval})

    # 写 CSV
    with open(detailFile, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['epoch'] + [f'{p}_{b}' for p, b, _ in entries])
        for ep in epochs:
            row = [ep.strftime('%Y-%m-%d %H:%M:%S')]
            for prn, band, _ in entries:
                m = marks.get((prn, band))
                if m is None:
                    row.append('')
                elif ep in m:
                    row.append(str(m[ep]))
                else:
                    row.append('')      # 该历元该星不可见
            writer.writerow(row)

    return gaps
