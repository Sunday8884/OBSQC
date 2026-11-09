#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
spp_gate_main —— 带观测质量门控的 SPP 独立入口
================================================================
以 spos.py 为基础(多文件通配/日期配对/波段选择/SNX参考全部保留),
在解算前增加"极端剔除门控": QC 指标超阈值的卫星不进入该历元最小二乘。

门控层次(三层, 可叠加):
  [1] 阈值极端剔除(默认开启, gate_spp.py V1, 实验证实有效)
      MP>mp_thresh 或 CNR<cnr_thresh 的卫星剔除, 剔后保底 min_keep 星
      口径: 优先合并单星宽表 QCValue_sats.csv, MP=弧段累积RMS
      WUH2/BRUX: 3D RMS +3.0% ~ +5.4%, 最差历元 8.0->3.0m
  [2] 卫星级 S 逐历元剔除(--gate_sat_s worst1|fixed, 慎用)
      每(历元,卫星) 8 指标批内秩归一化评分(sppgate.score.sat_score),
      worst1=每历元剔 S_sat 最差1颗(可见>=sat_s_min_sat);
      fixed=剔 S_sat>=固定阈值(--sat_s_thresh, 默认0.75)的星, 自动保底
      min_keep
      ⚠ 方向 bug 已修复(2026-11-06, 原实现把 MP 大的星判成优), 且修复了
        sat_score 权重双重归一化 bug(S_sat 曾被放大 1/0.898~1/0.433 倍,
        阈值语义失真; 修复后 S_sat∈[0,1])。全天 WUH2 重测仍全负:
        critic T=0.75 -0.20% / T=0.70 -3.51% / T=0.60 -16.96%,
        fixed T=0.70 -4.92% / T=0.60 -22.43%, worst1 -1.70% ——
        秩归一化"相对最差"剔除在几何上得不偿失, 仅作诊断/对照, 默认关闭
  [3] 历元级 S 筛选(--epoch_s_pct 10 20, 解后筛选)
      综合评分 S(sppgate.score.epoch_score, 9指标加权)最差的 p% 历元
      从统计中剔除(解不变, 不丢观测), 对基线/门控解均可叠加
      ⚠ gate_epoch_spp.py 实测站内信号弱(rho≈-0.11), 保留率换精度
        收益有限, 同样默认关闭

方案编排(--compare 时并列对比, 各方案独立解算不组合):
  基线 | 门控(阈值) | S剔星(仅卫星级S) ; 历元级 S 筛选只对基线解
  独立成行(S历元筛<p>%), 输出对应 csv

QC 指标计算(sppgate.indicators):
  仅阈值门控 -> 只算 MP/CNR(快速档);
  开启任一 S 门控 -> 全模块(CNR/SLIP/MP/CMC/IOD/PN/HOD/COMP), 供评分。

流程(每观测文件):
  1. readObs 读入
  2. sppgate.indicators 进程内计算指标宽表
  3. 构建门控表(+S表), 运行时包裹 spos.spp.initObs 注入剔除
  4. spp() 解算 -> *_gate_spp.txt + 剔除明细 CSV
  5. (可选)历元级 S 筛选 -> *_Ssieve<pct>_enu.csv
  --compare 时基线(不剔除)与门控各跑一遍, 输出 RMS 对比(需 --snx)

独立运行示例:
  python sppgate.py --obs ./test/data/wuh20050.24o \
      --nav ./test/data/brdc0050.24p --snx ./test/data/IGS0OPSSNX*.SNX \
      --out_dir ./spp_gate --compare
  # 叠加 S 门控:
  python sppgate.py --obs ... --nav ... --snx ... \
      --gate_sat_s worst1 --epoch_s_pct 10
================================================================
"""
import os
import glob
import argparse
import numpy as np
import pandas as pd

from com.readObs import readObsHead, readObs
from com.readNav import readNav
from com.xyz2neu import xyz2neu
from spos.spos import spp, read_snx_coord
# 复用 spp_main 的骨架辅助函数(多文件/配对/波段/测站名)
from spp import (extract_station_from_filename, resolve_snx_file,
                 parse_time, get_bandchoose, extract_date_key,
                 ensure_dir)
# 门控支持包(自包含): QC 指标 + 综合评分 S
from sppgate.indicators import compute_indicators, load_wide
from sppgate import score as qcscore

_DEFAULT_WEIGHT_CSV = os.path.join(
    os.path.dirname(os.path.abspath(__file__)),
    'sppgate', 'weights_QC9.csv')

_ORIG_INITOBS = None


# ========================= 门控注入 =========================
def install_gate(gate_table, mp_thresh, cnr_thresh, min_keep,
                 s_table=None, s_rule=None, s_thresh=0.75,
                 s_min_sat=7):
    """包裹 spos.spp.initObs: 阈值命中 + 卫星级S规则 即剔星。

    s_table: {epoch: {prn: S_sat}} 卫星级评分(可见星内相对, 1=最差);
    s_rule : 'worst1' | 'quantile' | None
    返回 state(可统计/关闭)。
    """
    global _ORIG_INITOBS
    import spos.spos as sppmod
    if _ORIG_INITOBS is None:
        _ORIG_INITOBS = sppmod.initObs

    state = {'on': True, 'n_gated': 0, 'n_dropped': 0, 'drop_detail': {},
             'epoch_flags': {}}

    def _s_flags(sat_use, epoch, flags):
        """卫星级 S 规则 -> 追加到 flags(不超过保底卫星数)"""
        row = s_table.get(epoch)
        if not row:
            return flags
        sub = {p: row[p] for p in sat_use if p in row}
        if not sub:
            return flags
        extra = set()
        if s_rule == 'worst1':
            if len(sat_use) >= s_min_sat and len(sub) >= 5:
                extra = {max(sub, key=sub.get)}
        elif s_rule == 'fixed':
            thr = s_thresh
            extra = {p for p, v in sub.items() if v >= thr}
            n_keep = len(sat_use) - len(flags) - len(extra)
            if n_keep < min_keep:   # 保底: 只剔评分最高的若干颗
                allowed = max(0, len(sat_use) - len(flags) - min_keep)
                extra = set(sorted(sub, key=sub.get, reverse=True)[:allowed])
        for p in extra:
            if p not in flags:
                flags[p] = f'S={sub[p]:.3f}'
        return flags

    def gated_initObs(bandChoose, obsData, epoch, satPosInEpoch, satClkInEpoch):
        o, isLess = _ORIG_INITOBS(bandChoose, obsData, epoch,
                                  satPosInEpoch, satClkInEpoch)
        if state['on'] and not isLess and (gate_table or s_table):
            flags = {}
            for prn in o['satUse']:
                d = gate_table.get(prn) if gate_table else None
                if d is None or epoch not in d.index:
                    continue
                row = d.loc[epoch]
                mp, cnr = row['mp'], row['cnr']
                if pd.notna(mp) and float(mp) > mp_thresh:
                    flags[prn] = f'MP={mp:.2f}'
                if pd.notna(cnr) and float(cnr) < cnr_thresh:
                    flags[prn] = (flags.get(prn, '') + ' ' if prn in flags
                                  else '') + f'CNR={cnr:.1f}'
            if s_table and s_rule:
                flags = _s_flags(o['satUse'], epoch, flags)
            if flags:
                keep = [p for p in o['satUse'] if p not in flags]
                if len(keep) >= min_keep:      # 安全保护: 保底卫星数
                    o['satUse'] = keep
                    o['obs'] = {p: o['obs'][p] for p in keep}
                    o['prnIndex'] = {p: i for i, p in enumerate(keep)}
                    state['n_gated'] += 1
                    state['n_dropped'] += len(flags)
                    for p in flags:
                        state['drop_detail'][p] = \
                            state['drop_detail'].get(p, 0) + 1
                    state['epoch_flags'][epoch] = flags
        return o, len(o['satUse']) < 4

    sppmod.initObs = gated_initObs
    return state


def restore_initObs():
    """恢复原始 initObs(多文件循环防重复包裹)"""
    global _ORIG_INITOBS
    if _ORIG_INITOBS is not None:
        import spos.spos as sppmod
        sppmod.initObs = _ORIG_INITOBS
        _ORIG_INITOBS = None


# ========================= 门控表构建 =========================
def build_gate_table(obs_file, nav, out_dir, bands=None, systems=None,
                     modules=('mp', 'cnr')):
    """QC 指标计算 -> {prn: DataFrame(mp, cnr)}

    先查缓存: qc_ind/QCValue_sats.csv 存在且口径一致(systems/bands/modules
    与本次一致, 见 _cache_meta.json)时直接复用, 不重复计算 QC。
    口径说明: 优先读合并单星宽表 QCValue_sats.csv(列 '<prn>_<指标>_<频段>',
    其 MP 为弧段累积RMS —— 门控实验验证有效的口径; 由 sppgate.indicators
    自动合并产出); 兼容旧的单星文件 QCValue_<sat>_*.csv; 都缺失时退回
    逐历元宽表(瞬时值)。modules: ('mp','cnr') 快速档; None=全模块(供 S)。
    """
    import json
    import re
    qc_out = os.path.join(out_dir, 'qc_ind')
    cache = os.path.join(qc_out, 'QCValue_sats.csv')
    meta_p = os.path.join(qc_out, '_cache_meta.json')
    meta = {'systems': systems or None, 'bands': bands or None,
            'modules': sorted(modules) if modules else None}
    found = None
    if os.path.exists(cache) and os.path.exists(meta_p):
        try:
            with open(meta_p, encoding='utf-8') as f:
                if json.load(f).get('meta') == meta:
                    found = {'QCValue_sats': cache}
                    print(f'[门控] 已有合并单星宽表且口径一致, '
                          f'跳过 QC 重复计算: {cache}')
        except (OSError, ValueError):
            pass
    if found is None:
        print(f'[门控] QC 指标计算 '
              f'({"全模块" if modules is None else "MP/CNR"} 口径)...')
        found = compute_indicators(obs_file, nav=nav,
                                   out_dir=qc_out,
                                   bands=bands, systems=systems,
                                   modules=list(modules) if modules else None)
        try:
            with open(meta_p, 'w', encoding='utf-8') as f:
                json.dump({'meta': meta}, f)
        except OSError:
            pass
    gate = {}

    # ---- 优先: 合并单星宽表 QCValue_sats.csv (累积RMS口径) ----
    merged_path = found.get('QCValue_sats')
    if merged_path and os.path.exists(merged_path):
        try:
            w = load_wide(merged_path)
            for prn in sorted({c.split('_')[0] for c in w.columns}):
                cols = [c for c in w.columns if c.startswith(prn + '_')]
                mp_cols = [c for c in cols if c.startswith(prn + '_MP_')]
                cnr_cols = [c for c in cols if c.startswith(prn + '_CNR_')]
                mp_s = w[mp_cols].max(axis=1) if mp_cols else \
                    pd.Series(np.nan, index=w.index)   # MP取频段最大
                cnr_s = w[cnr_cols].min(axis=1) if cnr_cols else \
                    pd.Series(np.nan, index=w.index)   # CNR取频段最小
                d = pd.DataFrame({'mp': mp_s, 'cnr': cnr_s}).dropna(how='all')
                if len(d):
                    gate[prn] = d
            print(f'[门控] 使用合并单星宽表 QCValue_sats.csv '
                  f'({len(gate)} 颗卫星, MP=弧段累积RMS口径)')
            return gate
        except Exception as e:
            print(f'  警告: 合并单星宽表读取失败({e}), 退回单星文件模式')

    # ---- 兼容: 旧的单星时间序列文件(累积RMS口径) ----
    for f in glob.glob(os.path.join(out_dir, '**', 'QCValue_*.csv'),
                       recursive=True):
        m = re.search(r'QCValue_([GCE]\d{2})_', os.path.basename(f))
        if not m:
            continue  # 汇总表(频段标识命名)被正则过滤
        prn = m.group(1)
        try:
            w = load_wide(f)
        except Exception:
            continue
        mp_cols = [c for c in w.columns if c.startswith('MP_')]
        cnr_cols = [c for c in w.columns if c.startswith('CNR_')]
        if not mp_cols and not cnr_cols:
            continue
        mp_s = w[mp_cols].max(axis=1) if mp_cols else \
            pd.Series(np.nan, index=w.index)
        cnr_s = w[cnr_cols].min(axis=1) if cnr_cols else \
            pd.Series(np.nan, index=w.index)
        d = pd.DataFrame({'mp': mp_s, 'cnr': cnr_s})
        gate[prn] = d
    if gate:
        print(f'[门控] 使用单星时间序列表 ({len(gate)} 颗卫星, '
              f'MP=弧段累积RMS口径)')
        return gate

    # ---- 后备: 逐历元宽表(瞬时值口径) ----
    print('[门控] 单星表缺失, 使用逐历元宽表(瞬时值口径)')
    for suffix, col in (('MP', 'mp'), ('CNR', 'cnr')):
        p = found.get(suffix)
        if p is None:
            print(f'  警告: 缺少 {suffix} 宽表, 该指标不参与门控')
            continue
        w = load_wide(p)
        for prn in {c.split('_')[0] for c in w.columns}:
            cols = [c for c in w.columns if c.startswith(prn + '_')]
            if not cols:
                continue
            s = w[cols].max(axis=1) if suffix == 'MP' \
                else w[cols].min(axis=1)   # MP取频段最大, CNR取频段最小
            s = s.dropna()
            if len(s):
                gate.setdefault(prn, pd.DataFrame(index=s.index)) \
                    .__setitem__(col, s)
    # 对齐 mp/cnr 两列
    for prn, d in gate.items():
        if 'mp' not in d.columns:
            d['mp'] = np.nan
        if 'cnr' not in d.columns:
            d['cnr'] = np.nan
    return gate


def build_s_tables(qc_dir, weight_csv, weight_mode='fixed',
                   score_method='wavg'):
    """S 评分表: (S_epoch Series, S_sat {epoch:{prn:S_sat}}, 权重表)

    weight_mode: fixed=外部权重表 CSV(默认); entropy/critic/equal/ahp/game=
    自包含赋权(由本文件指标矩阵确定权重, 不读 weight_csv; ahp=主观先验,
    game=博弈论组合熵权/CRITIC/AHP)。
    score_method: 历元级 S 评分器(wavg/topsis/gra/waspas/ecdf),
    详见 sppgate.score.SCORE_METHODS。
    """
    ind = qcscore.epoch_indicators(qc_dir)
    if weight_mode == 'fixed':
        wtab = qcscore.load_weight_table(weight_csv)
        print(f'[S评分] 权重表: {weight_csv}')
    else:
        wtab = qcscore.build_weight_table_local(ind, weight_mode)
        print(f'[S评分] 自包含赋权({weight_mode}): ' +
              ', '.join(f'{k} {v:.3f}' for k, v in
                        wtab['最终权重'].sort_values(ascending=False).items()))
    S_epoch, contrib = None, {}
    if len(ind):
        S_epoch, contrib = qcscore.epoch_score(ind, wtab,
                                               method=score_method)
        print(f'[S评分] 历元级 S({score_method}): {len(S_epoch)} 历元, '
              f'指标 {list(contrib)}')
    S_sat = qcscore.sat_score(qc_dir, wtab)
    n_rec = sum(len(v) for v in S_sat.values())
    print(f'[S评分] 卫星级 S_sat: {len(S_sat)} 历元 / {n_rec} 条(历元,卫星)')
    return S_epoch, S_sat, wtab


def filter_by_epoch_s(df, S_epoch, pct):
    """历元级 S 筛选(解后筛选): 剔 S 最差 pct% 历元, 返回筛选后 DataFrame"""
    common = df.index.intersection(S_epoch.index)
    S_c = S_epoch.loc[common]
    thr = S_c.quantile(1 - pct / 100)
    return df.loc[S_c[S_c <= thr].index]


# ========================= 结果处理 =========================
def pos_to_enu(posData, ref):
    """SPP posData -> DataFrame(index=epoch, E/N/U/PDOP/3D)"""
    refX, refY, refZ = ref
    rows = []
    for ep, v in posData.items():
        if 'x' not in v:
            continue
        N, E, U = xyz2neu(refX, refY, refZ, v['x'], v['y'], v['z'], 'WGS84')
        rows.append((ep, E, N, U, v.get('PDOP', np.nan)))
    df = pd.DataFrame(rows, columns=['epoch', 'E', 'N', 'U', 'PDOP'])
    df = df.set_index('epoch').sort_index()
    df['3D'] = np.sqrt(df['E'] ** 2 + df['N'] ** 2 + df['U'] ** 2)
    return df


def rms_table(df, name):
    r = {'方案': name, '历元数': len(df)}
    for c in ('E', 'N', 'U', '3D'):
        r[c] = float(np.sqrt((df[c] ** 2).mean()))
    r['3D中位'] = float(df['3D'].median())
    r['PDOP'] = float(df['PDOP'].mean())
    return r


# ========================= 单文件处理 =========================
def process_one_file(obs_file, navData, args, base_output_dir):
    file_base = os.path.splitext(os.path.basename(obs_file))[0]
    out_dir = os.path.join(base_output_dir, file_base)
    ensure_dir(out_dir)

    start_dt = parse_time(args.start) if args.start else None
    end_dt = parse_time(args.end) if args.end else None
    band_choose = get_bandchoose(args.systems, args.bands)
    if not band_choose:
        print('错误: 未为任何系统指定波段')
        return
    station = extract_station_from_filename(obs_file)

    sat_s_on = args.gate_sat_s != 'off'
    s_any_on = sat_s_on or bool(args.epoch_s_pct)

    print(f"{'=' * 60}\n处理观测文件: {obs_file}")
    obsHead = readObsHead(obs_file, needSatList=True)
    obsData = readObs(obs_file, obsHead=obsHead)
    print(f'历元数量: {len(obsData)}')

    # ---- 门控表 + S 表(QC 指标, 进程内) ----
    # 开启任一 S 门控时需全模块指标(9类), 否则只算 MP/CNR 快速档
    modules = None if s_any_on else ['mp', 'cnr']
    gate = None
    S_epoch, S_sat = None, None
    state = None
    if not args.no_gate or s_any_on:
        print(f'[门控] 阈值口径: MP>{args.mp_thresh}m 或 '
              f'CNR<{args.cnr_thresh}dB-Hz, 保底{args.min_keep}星')
        gate = build_gate_table(obs_file, args.nav, out_dir,
                                bands=args.bands, systems=args.systems,
                                modules=modules)
        n_over = sum(int(((d['mp'] > args.mp_thresh) |
                          (d['cnr'] < args.cnr_thresh)).sum())
                     for d in gate.values())
        print(f'[门控] {len(gate)} 颗卫星, 超阈值记录 {n_over} 条')
    if s_any_on:
        qc_dir = os.path.join(out_dir, 'qc_ind')
        S_epoch, S_sat, _ = build_s_tables(qc_dir, args.s_weight_csv,
                                           weight_mode=args.s_weight_mode,
                                           score_method=args.s_score_method)
        if sat_s_on and not S_sat:
            print('  警告: 卫星级 S 表为空, --gate_sat_s 不生效')
        if args.epoch_s_pct and S_epoch is None:
            print('  警告: 历元级 S 为空, --epoch_s_pct 不生效')

    # ---- SNX 参考 ----
    ref = None
    if args.snx:
        try:
            ref = read_snx_coord(args.snx, station)
        except Exception as e:
            print(f'警告: SNX 读取失败({e}), 不输出 RMS 摘要')
            ref = None

    # ---- 解算 ----
    # 方案编排: 各方案独立解算, 不组合 --
    #   基线 | 门控(阈值) | S剔星(仅卫星级S) ;
    #   S历元筛为解后处理, 只对基线解独立成行(不与门控/S剔星叠加)
    summaries, posFrames, df_schemes = [], {}, {}
    if args.compare:
        runs = [('基线', False, False)]
        if not args.no_gate:
            runs.append(('门控', True, False))
        if sat_s_on:
            runs.append(('S剔星', False, True))
    elif args.no_gate:
        runs = [('基线', False, False)]
    elif sat_s_on:
        runs = [('S剔星', False, True)]
    else:
        runs = [('门控', True, False)]
    for name, gate_on, sat_on in runs:
        use_gate = (gate is not None and gate_on)
        use_sat = (sat_s_on and S_sat and sat_on)
        if use_gate or use_sat:
            if use_sat:
                print(f'[门控] 叠加卫星级 S 剔除: {args.gate_sat_s}'
                      f'(fixed阈值={args.sat_s_thresh}, '
                      f'worst1需可见>={args.sat_s_min_sat})')
            state = install_gate(gate if use_gate else {},
                                 args.mp_thresh, args.cnr_thresh,
                                 args.min_keep,
                                 s_table=S_sat if use_sat else None,
                                 s_rule=args.gate_sat_s if use_sat else None,
                                 s_thresh=args.sat_s_thresh,
                                 s_min_sat=args.sat_s_min_sat)
            state['on'] = True
        else:
            state = install_gate({}, 1e9, -1e9, args.min_keep)
            state['on'] = False
        tag = name
        pos_file = os.path.join(out_dir, f'{file_base}_{tag}_spp.txt')
        print(f'[SPP] {tag} 解算中...')
        posData = spp(obsHead, obsData, navData, band_choose,
                      startDatetime=start_dt, endDatetime=end_dt,
                      posFile=pos_file)[0]
        restore_initObs()
        if posData is None:
            print(f'  {tag} 解算失败: 有效历元不足')
            continue
        if state['n_gated']:
            print(f"  触发 {state['n_gated']} 历元, 剔除 {state['n_dropped']} 颗次; "
                  f"剔除Top5: {dict(sorted(state['drop_detail'].items(), key=lambda kv: -kv[1])[:5])}")
        posFrames[tag] = posData
        if ref is not None:
            df = pos_to_enu(posData, ref)
            df['剔除'] = [','.join(f'{p}({r})' for p, r in
                                   state['epoch_flags'].get(ep, {}).items())
                          for ep in df.index]
            df.to_csv(os.path.join(out_dir, f'{file_base}_{tag}_enu.csv'),
                      encoding='utf-8-sig')
            df_schemes[tag] = df
            summaries.append(rms_table(df, tag))
            print(f'  已保存: {file_base}_{tag}_enu.csv')
            # ---- 历元级 S 筛选(解后筛选, 解不变; 仅基线, 独立成行) ----
            if (args.epoch_s_pct and S_epoch is not None
                    and tag == '基线'):
                for pct in args.epoch_s_pct:
                    dfk = filter_by_epoch_s(df, S_epoch, pct)
                    if dfk.empty:
                        print(f'  警告: S筛选{pct}% 后无历元, 跳过')
                        continue
                    summaries.append(rms_table(dfk, f'S历元筛{pct}%'))
                    df_schemes[f'S历元筛{pct}%'] = dfk
                    dfk.to_csv(os.path.join(
                        out_dir, f'{file_base}_Ssieve{pct}_enu.csv'),
                        encoding='utf-8-sig')
                    print(f'  S历元筛{pct}%: 保留 {len(dfk)}/{len(df)} 历元, '
                          f'已保存: {file_base}_Ssieve{pct}_enu.csv')

    # ---- 方案对比图 (E/N/U/3D 四子图, 逐历元折线) ----
    if df_schemes:
        from sppgate.plot import plot_schemes
        png = os.path.join(out_dir, f'{file_base}_对比.png')
        try:
            plot_schemes(df_schemes, png,
                         title=f'{file_base} 各方案 SPP 误差时序对比')
            print(f'-->方案对比图: {png}')
        except Exception as e:
            print(f'  警告: 对比图绘制失败({e})')

    # ---- 对比摘要 ----
    if summaries:
        st = pd.DataFrame(summaries).set_index('方案')
        if '基线' in st.index:
            b = st.loc['基线', '3D']
            st['3D提升%'] = (b - st['3D']) / b * 100
        print('\n[RMS 摘要]')
        print(st.to_string(float_format=lambda v: f'{v:.3f}'))
        st.to_csv(os.path.join(out_dir, f'{file_base}_gate_summary.csv'),
                  encoding='utf-8-sig')
    print(f'完成: 结果保存在 {out_dir}')


# ========================= 主流程(骨架同 spp_main) =========================
def run_spp(args):
    args.snx = resolve_snx_file(args.snx)
    if args.snx:
        print(f'使用SINEX文件: {args.snx}')

    obs_files = glob.glob(args.obs)
    if not obs_files:
        print(f"错误: 未找到匹配 '{args.obs}' 的观测文件")
        return
    nav_files = glob.glob(args.nav)
    if not nav_files:
        print(f"错误: 未找到匹配 '{args.nav}' 的导航文件")
        return
    print(f'找到 {len(obs_files)} 个观测文件, {len(nav_files)} 个导航文件')

    base_output_dir = args.out_dir
    ensure_dir(base_output_dir)

    if len(nav_files) == 1:
        navData = readNav(nav_files[0])
        if not navData:
            print('错误: 读取导航文件失败')
            return
        for obs_file in obs_files:
            process_one_file(obs_file, navData, args, base_output_dir)
    else:
        nav_dict = {}
        for nav_file in nav_files:
            key = extract_date_key(nav_file)
            if key is None:
                continue
            navData = readNav(nav_file)
            if navData:
                nav_dict[key] = navData
        for obs_file in obs_files:
            key = extract_date_key(obs_file)
            if key in nav_dict:
                process_one_file(obs_file, nav_dict[key], args, base_output_dir)
            else:
                print(f'警告: 观测 {obs_file} 无匹配导航文件, 跳过')

    print(f'\n所有处理完成, 结果保存在 {base_output_dir}')


def create_parser():
    parser = argparse.ArgumentParser(
        description='带观测质量门控的SPP(阈值极端剔除 + 卫星级S剔除 + '
                    '历元级S筛选) —— 多文件通配/日期配对/SNX参考同 spp_main')
    parser.add_argument('--obs', required=True, help='RINEX 观测文件(通配符)')
    parser.add_argument('--nav', required=True, help='广播星历(通配符)')
    parser.add_argument('--out_dir', default='./spp_gate', help='输出根目录')
    parser.add_argument('--start', help='起始时间(YYYY-MM-DD HH:MM:SS)')
    parser.add_argument('--end', help='结束时间(YYYY-MM-DD HH:MM:SS)')
    parser.add_argument('--systems', nargs='+', help='GNSS 系统, 如 G C E')
    parser.add_argument('--bands', help="波段, 如 'G:C1C,C2W;C:C2I,C6I'")
    parser.add_argument('--snx', default=None, help='SINEX 文件(ENU 参考)')
    # [1] 阈值极端剔除
    parser.add_argument('--no-gate', dest='no_gate', action='store_true',
                        help='关闭阈值门控(等价 spp_main)')
    parser.add_argument('--mp_thresh', type=float, default=1.0,
                        help='MP 极端阈值(m), 默认 1.0')
    parser.add_argument('--cnr_thresh', type=float, default=20.0,
                        help='CNR 极端下限(dB-Hz), 默认 20')
    parser.add_argument('--min_keep', type=int, default=5,
                        help='每历元最少保留卫星数, 默认 5')
    # [2] 卫星级 S 逐历元剔除
    parser.add_argument('--gate_sat_s', default='off',
                        choices=['off', 'worst1', 'fixed'],
                        help='卫星级S剔除规则(默认 off; worst1=每历元剔S最差'
                             '1颗; fixed=剔S>=固定阈值的星)')
    parser.add_argument('--sat_s_thresh', type=float, default=0.75,
                        help='fixed 规则的 S_sat 剔除阈值(0~1), 默认 0.75')
    parser.add_argument('--sat_s_min_sat', type=int, default=7,
                        help='worst1 规则要求的可见星数下限, 默认 7')
    # [3] 历元级 S 筛选
    parser.add_argument('--epoch_s_pct', type=int, nargs='+', default=None,
                        help='历元级S筛选剔除比例(%%), 可多个, 如 10 20')
    parser.add_argument('--s_weight_csv', default=_DEFAULT_WEIGHT_CSV,
                        help='S 评分权重表 CSV(默认 sppgate/weights_QC9.csv)')
    parser.add_argument('--s_weight_mode', default='fixed',
                        choices=['fixed', 'entropy', 'critic', 'equal',
                                 'ahp', 'game'],
                        help='S 赋权模式: fixed=权重表(默认); '
                             'entropy=熵权 critic=CRITIC equal=等权 '
                             'ahp=AHP先验等级 game=博弈论组合赋权 '
                             '(后五种自包含, 不读权重 CSV, S 仅文件内可比)')
    parser.add_argument('--s_score_method', default='wavg',
                        choices=['wavg', 'topsis', 'gra', 'waspas', 'ecdf'],
                        help='历元级 S 评分器: wavg=加权平均(默认); '
                             'topsis=理想解相对接近度; gra=灰色关联; '
                             'waspas=加权和+积; ecdf=经验分位秩')
    parser.add_argument('--compare', action='store_true',
                        help='基线与门控双跑并输出 RMS 对比(需 --snx)')
    return parser


if __name__ == '__main__':
    run_spp(create_parser().parse_args())
