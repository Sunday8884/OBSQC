#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
time 模式：逐历元时间序列质量汇总
==================================================

与其他模式的关系
----------------
  prn      : 单文件 → PRN 级汇总（_QCValue.csv）
  station  : 多文件 → 测站级汇总 + 多站对比表
  Sstation : 多文件（单站多日期）→ 行=日期 的时间序列表
  time     : 多文件 → **行=历元** 的时间序列表（QCValue.csv）

time 模式流程
-------------
  1. 逐文件执行 QC 各模块（与 station 模式完全相同的 10 个模块 + SPP）；
  2. 读取各模块的**逐历元明细 CSV**：
        _cnr.csv / _slip.csv / _multipath.csv / _cmc.csv / _iod.csv /
        _phasenoise.csv / _highOrderDiff.csv / _completeness.csv
  3. **只保留「计算 SPP 时实际选用的频段」**——按每个文件实际选出的 SPP
     伪距码对过滤列（非 SPP 频点的列一律丢弃）；
  4. 逐历元把各模块的明细汇总成一行（均值 / RMS / 计数 / 有效样本数），
     写入 <output_dir>/QCValue_<标签>.csv；
  5. 另出 <output_dir>/QCValue_sppband.txt：SPP 选段 + 频点统计汇总
     （中心频率、Δf、λW、MW 放大倍数、覆盖卫星数、有效观测样本数）。

文件命名
--------
  标签 = 各系统的「系统字母 + SPP 选段各码的频段键」组合，系统内双频用 '_'
  连接，各系统之间也用 '_' 连接，系统按字母序排列。例如
  BDS 选 C2I+C6I、Galileo 选 C1C+C5Q、GPS 选 C1C+C2W 时：

     QCValue_C2I_6I_E1C_5Q_G1C_2W.csv   全卫星逐历元时间序列
     QCValue_C19_C2I_6I.csv             单星（--sat C19）逐历元时间序列
     QCValue_G01_G1C_2W.csv             单星（--sat G01）逐历元时间序列
     QCValue_sppband.txt                SPP 频段统计汇总

单星表的统计口径
----------------
  单星表**按频段分列**，不做跨频段平均：每个 SPP 频段各出一组统计量。
  像 CSR（周跳）这类**与频段无关**的模块（周跳检测按卫星整体判定，不区分
  频段），只出一列不带频段后缀的列。

输出列
------
  时间/标识列：STATION（测站名，取观测文件名前 4 位）、Epoch（实际时间戳）、
                Idx（采样序号），不含相对秒数。

  每个模块只出**一列主统计量**，列名不带统计量后缀：
    CNR → 均值      MP / CMC / IOD / PN / HOD → RMS      CSR → 计数
  单星表在此基础上按频段加后缀：CNR_2I / CNR_6I / MP_2I / MP_6I …
  例外：COMP（完整率）与 SAT（卫星数）不按频段分列，见下。

  另附各系统可见卫星数（整站量，全星表与单星表一致）：
    SAT                    统计范围内各系统的可见卫星数之和

  完整率 COMP 为整段弧段的统计量（0~1 小数），不按频段分列（双频取值几乎相同）：
    全星表 = 跨全部卫星、跨其各频段取均值
    单星表 = 该星各频段完整率的均值

  SPP 逐历元解算结果置于最右（整站量，全星表与单星表一致）：
    SPP_E / SPP_N / SPP_U  相对参考原点的东/北/天残差 (m)
    SPP_3D                 sqrt(E²+N²+U²) 三维误差 (m)

时间轴
------
  - 行标签为各模块 CSV 中的**实际时间戳**；
  - STATION = 测站名，取观测文件名的**前 4 位**（大写），
    逐历元标记该行的数据来源（多文件跨站时据此区分）；
  - Idx   = 相对时间原点的整数索引，采样间隔由实际时间戳的中位差分推导
    （可用 --time_interval 覆盖），时间原点默认取**数据中最早的历元**
    （可用 --time_origin 指定其他基准）。
"""
import os
import csv
import math
import datetime
import statistics

from qc.qc_common import ensure_dir, parse_time
from com.readNav import readNav
from com.gnssParameter import getBandFreq, CLIGHT
from qc.qc_engine_exp import QCProcessor, _extract_date_key, _extract_station_from_filename

# 时间原点不设固定默认值：默认取数据中最早的历元（见 _derive_origin）

# 明细 CSV 后缀 → (指标名前缀, 统计类型)
#   level : 有量纲的观测量（CNR）→ 均值 / 最小值 / 样本数
#   resid : 零均值残差量（多路径/CMC/IOD/相位噪声/高阶差分）→ RMS / 峰值 / 样本数
#   flag  : 0/1 标记量（周跳、完整率）→ 计数 / 覆盖数
MODULE_SPECS = (
    ('cnr',           'CNR',  'level'),
    ('slip',          'CSR',  'flag'),
    ('multipath',     'MP',   'resid'),
    ('cmc',           'CMC',  'resid'),
    ('iod',           'IOD',  'resid'),
    ('phasenoise',    'PN',   'resid'),
    ('highOrderDiff', 'HOD',  'resid'),
    ('completeness',  'COMP', 'flag'),
)


# ========================= 累加器 =========================
class _Acc:
    """单历元、单模块的统计累加器（可跨文件合并）。"""

    __slots__ = ('n', 'nz', 's', 'ss', 'mn', 'mabs')

    def __init__(self):
        self.n = 0        # 有效（非空）样本数
        self.nz = 0       # 非零样本数
        self.s = 0.0      # 求和
        self.ss = 0.0     # 平方和
        self.mn = None    # 最小值
        self.mabs = 0.0   # 最大绝对值

    def add(self, v):
        self.n += 1
        if v != 0:
            self.nz += 1
        self.s += v
        self.ss += v * v
        if self.mn is None or v < self.mn:
            self.mn = v
        a = abs(v)
        if a > self.mabs:
            self.mabs = a

    def merge(self, o):
        self.n += o.n
        self.nz += o.nz
        self.s += o.s
        self.ss += o.ss
        if o.mn is not None and (self.mn is None or o.mn < self.mn):
            self.mn = o.mn
        if o.mabs > self.mabs:
            self.mabs = o.mabs

    # ---- 派生统计量 ----
    @property
    def mean(self):
        return self.s / self.n if self.n else None

    @property
    def rms(self):
        return math.sqrt(self.ss / self.n) if self.n else None


def _fmt(v, nd=4):
    if v is None:
        return ''
    try:
        if math.isnan(v) or math.isinf(v):
            return ''
    except TypeError:
        return str(v)
    return f'{v:.{nd}f}'


def _parse_col(col):
    """
    拆解明细 CSV 的列名。

    列名形如 'C19_C2I'（卫星 + 观测码）、'C05'（仅卫星，slip 模块）。
    返回 (prn, gsys, band_key)，band_key 为频点号+跟踪模式（如 '2I'），
    无频段信息时为 None。
    """
    if '_' in col:
        prn, code = col.split('_', 1)
    else:
        prn, code = col, ''
    gsys = prn[0] if prn else '?'
    key = code[1:] if len(code) >= 2 else None
    return prn, gsys, key


def _read_spp_txt(spp_file):
    """
    读取 SPP 结果文件 _spp.txt，提取逐历元的 ENU 残差。

    文件为定宽文本，列序：
        date mjd sod X[m] Y[m] Z[m] E[m] N[m] U[m] CLK[ns]
        HDOP VDOP PDOP GDOP

    返回 {epoch_str: (E, N, U)}；文件不存在或格式不符时返回空字典。
    """
    if not os.path.exists(spp_file):
        return {}
    out = {}
    try:
        with open(spp_file, 'r', encoding='utf-8') as f:
            lines = f.read().splitlines()
    except Exception:
        return {}

    for line in lines:
        t = line.split()
        # token 序：0-1 日期时间、2 mjd、3 sod、4-6 X/Y/Z、7-9 E/N/U、
        #           10 CLK[ns]、11-14 HDOP/VDOP/PDOP/GDOP
        if len(t) < 15:
            continue
        try:
            epoch = f'{t[0]} {t[1]}'
            e = float(t[7])       # E[m]
            n = float(t[8])       # N[m]
            u = float(t[9])       # U[m]
        except (ValueError, IndexError):
            continue
        out[epoch] = (e, n, u)
    return out


def _read_sat_csv(sat_file):
    """
    读取卫星数量明细 CSV（_sat.csv），提取逐历元各系统的可见卫星数。

    宽表格式：行=历元，列=各系统（'C'、'G' … 及 'ALL'）。
    返回 {epoch_str: {系统: 数量}}；文件不存在时返回空字典。
    """
    header, rows = _read_wide_csv(sat_file)
    if not header or not rows:
        return {}
    if 'epoch' not in header:
        return {}
    systems = [c for c in header if c and c != 'epoch']
    out = {}
    for row in rows:
        ep = (row.get('epoch') or '').strip()
        if not ep:
            continue
        counts = {}
        for sys in systems:
            raw = (row.get(sys) or '').strip()
            if raw == '':
                continue
            try:
                counts[sys] = int(float(raw))
            except ValueError:
                continue
        out[ep] = counts
    return out


def _read_wide_csv(csv_file):
    """
    读取宽表明细 CSV。

    Returns
    -------
    (header_cols, rows) : (list[str], list[dict])
        header_cols 含首列 'epoch'；rows 每行是 {col: str}
    """
    if not os.path.exists(csv_file):
        return None, None
    with open(csv_file, 'r', newline='', encoding='utf-8-sig') as f:
        reader = csv.reader(f)
        try:
            header = next(reader)
        except StopIteration:
            return None, None
        header = [h.strip() for h in header]
        rows = []
        for row in reader:
            if not row:
                continue
            rows.append(dict(zip(header, row)))
    return header, rows


def _band_freq_mhz(gsys, code):
    """伪距/载波观测码 → 中心频率 (MHz)，取不到返回 None。"""
    try:
        f = getBandFreq(gsys, code)
    except Exception:
        return None
    return f


# ========================= 引擎 =========================
class TimeSeriesQC:
    """time 模式引擎：多文件 → 逐历元时间序列 QCValue_<标签>.csv。"""

    def __init__(self, args):
        self.args = args
        # 时间原点：--time_origin 未指定时，取数据中最早的历元（自动推导）
        origin = getattr(args, 'time_origin', None)
        self.origin = parse_time(origin) if origin else None
        self.origin_auto = (self.origin is None)
        self.interval_override = getattr(args, 'time_interval', None)
        self.interval = None          # 采样间隔（秒），运行后确定

        # 选星参数 --sat：只对这些卫星单独出汇总表
        self.sat_filter = []
        for s in (getattr(args, 'sat', None) or []):
            s = str(s).strip().upper()
            if s and s not in self.sat_filter:
                self.sat_filter.append(s)

        self.acc = {}                 # (modkey, epoch_str) -> _Acc
        self.sat_acc = {}             # (prn, modname, band_key, epoch) -> _Acc
        self.sat_mod_bands = {}        # (prn, modname) -> 该模块涉及的频段键集合
        self.sat_epochs = {}          # prn -> set(epoch_str)
        # 完整率：(prn, 频段键) -> {'obs','arc','ratio'}
        # ratio = 有值历元数 / 可见弧段历元数，与 completeness 模块的
        # ratio_arc 口径一致（弧段完整率，整列恒定）
        self.comp_ratio = {}
        self.spp_enu = {}             # epoch_str -> (E, N, U) 逐历元 ENU 残差
        self.sat_count = {}           # epoch_str -> {系统: 可见卫星数}
        self.epoch_station = {}       # epoch_str -> 测站名（文件名前 4 位）
        self.epochs = set()           # 全部历元时间戳（字符串）
        self.spp_info = {}            # gsys -> {'codes': [...], 'keys': set()}
        self.spp_freq_counts = {}     # (prn, spp_code) -> freq 模块的观测次数
        self.file_summaries = []      # 每个文件的 SPP 选段记录
        self.band_usage = {}          # (gsys, key) -> {'prns': {prn: cnt}, 'nobs'}
        self.band_valid = {}          # (gsys, key) -> {prn: 有数据历元数}

    # ---------- 导航数据 ----------
    def _nav_for(self, obs_file):
        pattern = getattr(self.args, 'nav', None)
        if not pattern:
            return None
        import glob as _glob
        nav_files = _glob.glob(pattern)
        if not nav_files:
            print(f"   警告：未找到导航文件 '{pattern}'，跳过 SPP")
            return None
        if len(nav_files) == 1:
            key = nav_files[0]
        else:
            ok = _extract_date_key(obs_file)
            key = None
            if ok:
                for nf in nav_files:
                    if _extract_date_key(nf) == ok:
                        key = nf
                        break
            key = key or nav_files[0]
        try:
            return readNav(key)
        except Exception as e:
            print(f"   警告：读取导航文件失败 {key}: {e}")
            return None

    # ---------- SPP 选段记录 ----------
    @staticmethod
    def _spp_codes(band_choose):
        """{sys: {code: {}}} → {sys: [code, ...]}（保持选定顺序）"""
        out = {}
        for gsys, d in (band_choose or {}).items():
            out[gsys] = list(d.keys())
        return out

    def _write_spp_csv(self, sub_dir, base, codes_by_sys):
        """把本文件的 SPP 选段写入 <base>_sppband.csv（供断点续跑复用）。"""
        ensure_dir(sub_dir)
        path = os.path.join(sub_dir, f'{base}_sppband.csv')
        with open(path, 'w', newline='', encoding='utf-8') as f:
            w = csv.writer(f)
            w.writerow(['system', 'code1', 'code2', 'key1', 'key2',
                        'f1_MHz', 'f2_MHz', 'df_MHz', 'lamdaW_m', 'mw_amp'])
            for gsys in sorted(codes_by_sys):
                codes = codes_by_sys[gsys]
                if not codes:
                    continue
                c1 = codes[0]
                c2 = codes[1] if len(codes) > 1 else ''
                f1 = _band_freq_mhz(gsys, c1)
                f2 = _band_freq_mhz(gsys, c2) if c2 else None
                df = abs(f1 - f2) if (f1 and f2) else None
                lam = CLIGHT / (df * 1e6) if df else None
                amp = ((f1 + f2) / df) if df else None
                w.writerow([gsys, c1, c2,
                            c1[1:] if len(c1) >= 2 else '',
                            c2[1:] if len(c2) >= 2 else '',
                            _fmt(f1, 3), _fmt(f2, 3), _fmt(df, 3),
                            _fmt(lam, 4), _fmt(amp, 2)])
        return path

    def _read_spp_csv(self, path):
        """读取 <base>_sppband.csv → ({sys: [code1, code2]}, 明细行)"""
        if not os.path.exists(path):
            return {}, []
        codes = {}
        rows = []
        with open(path, 'r', newline='', encoding='utf-8-sig') as f:
            reader = csv.reader(f)
            try:
                next(reader)
            except StopIteration:
                return {}, []
            for row in reader:
                if not row or not row[0].strip():
                    continue
                rows.append(row)
                lst = codes.setdefault(row[0].strip(), [])
                # code1 / code2 两列都要取，缺一不可（否则频段过滤只保留一半）
                for idx in (1, 2):
                    if len(row) > idx:
                        c = row[idx].strip()
                        if c and c not in lst:
                            lst.append(c)
        return codes, rows

    def _write_spp_freq_csv(self, sub_dir, base, codes_by_sys, freq_data):
        """把 SPP 频段的逐 PRN 频段计数（来自 freq 模块）写入
        <base>_sppband_freq.csv。

        freq 模块统计的是**载波观测码**（如 'L2I'）的出现次数，
        SPP 选段用的是**伪距观测码**（如 'C2I'），二者共享同一个频段键
        （'2I'），据此对应。
        """
        ensure_dir(sub_dir)
        path = os.path.join(sub_dir, f'{base}_sppband_freq.csv')
        rows = []
        for prn in sorted(freq_data or {}):
            gsys = prn[0]
            codes = (codes_by_sys or {}).get(gsys) or []
            for code in codes:
                key = code[1:] if len(code) >= 2 else ''
                lcode = 'L' + key
                cnt = (freq_data.get(prn) or {}).get(lcode)
                if cnt is None:
                    continue
                rows.append((prn, gsys, code, lcode, cnt))
        with open(path, 'w', newline='', encoding='utf-8') as f:
            w = csv.writer(f)
            w.writerow(['prn', 'system', 'spp_code', 'phase_code', 'obs_count'])
            for r in rows:
                w.writerow(r)
        return path

    def _read_spp_freq_csv(self, path):
        """读取 <base>_sppband_freq.csv → {(prn, code): count}"""
        out = {}
        if not os.path.exists(path):
            return out
        with open(path, 'r', newline='', encoding='utf-8-sig') as f:
            reader = csv.reader(f)
            try:
                next(reader)
            except StopIteration:
                return out
            for row in reader:
                if len(row) < 5 or not row[0].strip():
                    continue
                try:
                    out[(row[0].strip(), row[2].strip())] = int(float(row[4]))
                except ValueError:
                    continue
        return out

    def _resolve_bands_from_head(self, obs_file, out_dir):
        """仅读观测头，重算 SPP 选段（用于已有 QCValue_station.csv 但缺
        _sppband.csv 的断点续跑场景，避免重跑全部 QC 模块）。"""
        from com.readObs import readObsHead

        try:
            head = readObsHead(obs_file, needSatList=False)
        except Exception as e:
            print(f'   警告：读取观测头失败 {obs_file}: {e}')
            return {}

        proc = QCProcessor(obs_file, out_dir, self.args)
        proc.obs_head = head
        try:
            band_choose = proc.resolve_spp_band_choose()
        except Exception as e:
            print(f'   警告：SPP 选段计算失败: {e}')
            return {}
        return self._spp_codes(band_choose)

    # ---------- 逐文件聚合 ----------
    def _aggregate_file(self, sub_dir, base, codes_by_sys, station=None):
        """
        读取本文件各模块的明细 CSV，按 SPP 频段过滤后逐历元累加。

        station 为该文件的测站名（文件名前 4 位），用于逐历元标记来源。
        """
        # 本文件允许参与统计的 (系统 → 频段键集合)
        allow = {}
        for gsys, codes in (codes_by_sys or {}).items():
            keys = {c[1:] for c in codes if len(c) >= 2}
            if keys:
                allow[gsys] = keys

        for suffix, mname, kind in MODULE_SPECS:
            path = os.path.join(sub_dir, f'{base}_{suffix}.csv')
            header, rows = _read_wide_csv(path)
            if not header or not rows:
                continue
            if 'epoch' not in header:
                continue
            cols = [c for c in header if c and c != 'epoch']

            # ---- 频段过滤：只保留 SPP 选中的频段 ----
            keep = []
            for c in cols:
                prn, gsys, key = _parse_col(c)
                if gsys not in allow:
                    continue
                if key is not None and key not in allow[gsys]:
                    continue
                keep.append((c, prn, gsys, key))

            if not keep:
                continue

            # ---- 频点覆盖统计（按 PRN 记观测次数）----
            for c, prn, gsys, key in keep:
                uk = (gsys, key) if key else (gsys, '*')
                u = self.band_usage.setdefault(uk, {'prns': {}, 'nobs': 0})
                u['prns'].setdefault(prn, 0)

            # 完整率模块：按 (prn, 频段) 记录「可见」与「有值」的历元，
            # 可见性以**该列非空**为准（空 = 该星该频段不可见），
            # 有值以标记为 1 为准（0 = 可见但该频段无值，属缺测）。
            comp_cols = ([(c, prn, key) for c, prn, gsys, key in keep]
                         if suffix == 'completeness' else [])

            for row in rows:
                ep = (row.get('epoch') or '').strip()
                if not ep:
                    continue
                self.epochs.add(ep)
                if station:
                    self.epoch_station.setdefault(ep, station)
                acc = self.acc.setdefault((mname, ep), _Acc())
                for c, prn, gsys, key in keep:
                    raw = (row.get(c) or '').strip()
                    if raw == '':
                        continue
                    try:
                        v = float(raw)
                    except ValueError:
                        continue
                    acc.add(v)
                    uk = (gsys, key) if key else (gsys, '*')
                    u = self.band_usage[uk]
                    u['nobs'] += 1
                    u['prns'][prn] = u['prns'].get(prn, 0) + 1
                    # 完整率模块的标记即「该频段该卫星有数据」，
                    # 单列出来作为有效数据历元数（不与其它模块混算）
                    if suffix == 'completeness':
                        vd = self.band_valid.setdefault(uk, {})
                        vd[prn] = vd.get(prn, 0) + 1
                    # ---- 选星参数 --sat：单独累加该星的统计量 ----
                    # 按频段分列累加：同一模块下不同频段各自独立统计，
                    # 不做跨频段平均；key 为 None 表示该模块与频段无关
                    # （如周跳检测按卫星整体判定），只统计一份。
                    if prn in self.sat_filter:
                        bk = key if key else ''
                        sacc = self.sat_acc.setdefault((prn, mname, bk, ep), _Acc())
                        sacc.add(v)
                        self.sat_epochs.setdefault(prn, set()).add(ep)
                        self.sat_mod_bands.setdefault((prn, mname), set()).add(bk)

            # ---- 完整率：按连续可见段统计 ----
            # 明细 CSV 的空值 = 该星该频段**不可见**（不在观测数据里），
            # 0 = 可见但该频段无值（真正的缺测），1 = 有值。
            # 分母只应计入「可见」的历元；若把不可见期也算进去，
            # 卫星中途长时间消失后再出现时完整率会被严重拉低
            # （实测 C19 可见 1165 历元，却因跨越 878 历元的中断只算出 57%）。
            if comp_cols:
                self._accum_comp_ratio(rows, comp_cols)

    def _accum_comp_ratio(self, rows, comp_cols):
        """按 (prn, 频段) 统计完整率：分母为连续可见段内的历元数。

        comp_cols: [(列名, prn, 频段键), ...]
        """
        # 每个 (prn, 频段) 收集：可见历元集合、有值历元集合
        visible = {}
        has_val = {}
        epochs_idx = {}
        for i, row in enumerate(rows):
            ep = (row.get('epoch') or '').strip()
            if ep:
                epochs_idx[ep] = i
        for row in rows:
            ep = (row.get('epoch') or '').strip()
            if not ep:
                continue
            for c, prn, key in comp_cols:
                raw = (row.get(c) or '').strip()
                if raw == '':
                    continue          # 不可见，不计入分母
                k = (prn, key or '')
                visible.setdefault(k, []).append(ep)
                if raw != '0':
                    has_val.setdefault(k, []).append(ep)

        for k, vis_eps in visible.items():
            # 按连续可见段切分：相邻可见历元在 epochs 中的索引间隔 > 1
            # 视为该星中途消失，期间的历元不计入分母
            segs = []
            cur = []
            prev_i = None
            for ep in vis_eps:
                i = epochs_idx.get(ep)
                if prev_i is not None and i is not None and i - prev_i > 1:
                    segs.append(cur)
                    cur = []
                cur.append(ep)
                prev_i = i
            if cur:
                segs.append(cur)
            arc_n = sum(len(s) for s in segs)
            obs_n = len(has_val.get(k, []))
            if arc_n <= 0:
                continue
            cr = self.comp_ratio.setdefault(k, {'obs': 0, 'arc': 0})
            cr['obs'] += obs_n
            cr['arc'] += arc_n

    # ---------- 标签（SPP 频段标识） ----------
    @staticmethod
    def _sys_tag(gsys, codes):
        """单系统的标签：系统字母 + SPP 选段各码的频段键，双频用 '_' 连接。

        系统字母只与**第一个**频段键直接相连（保持 C2I 这样的完整形式），
        第二个频段键以 '_' 追加：
          BDS 选 C2I+C6I    → 'C2I_6I'
          Galileo 选 C1C+C5Q → 'E1C_5Q'
          GPS 选 C1C+C2W    → 'G1C_2W'
        """
        if not codes:
            return gsys
        keys = []
        for code in codes:
            key = code[1:] if len(code) >= 2 else ''
            if key and key not in keys:
                keys.append(key)
        if not keys:
            return gsys
        return gsys + keys[0] + ''.join(f'_{k}' for k in keys[1:])

    def _band_tag(self):
        """全局标签：各系统标签用 '_' 连接，系统按字母序。"""
        if not self.spp_info:
            return ''
        parts = []
        for gsys in sorted(self.spp_info):
            tag = self._sys_tag(gsys, self.spp_info[gsys].get('codes'))
            if tag and tag not in parts:
                parts.append(tag)
        return '_'.join(parts)

    def _sat_tag(self, prn):
        """单星标签：取该星所属系统的 SPP 频段标识。"""
        gsys = prn[0] if prn else ''
        info = self.spp_info.get(gsys)
        if not info:
            return gsys
        return self._sys_tag(gsys, info.get('codes'))

    def _sat_band_keys(self, prn):
        """该星所属系统的 SPP 频段键列表（按选段顺序）。

        例：BDS 选 C2I+C6I → ['2I', '6I']。
        """
        gsys = prn[0] if prn else ''
        info = self.spp_info.get(gsys)
        if not info:
            return []
        keys = []
        for code in info.get('codes') or []:
            key = code[1:] if len(code) >= 2 else ''
            if key and key not in keys:
                keys.append(key)
        return keys

    # ---------- 采样间隔 / 时间原点 ----------
    def _derive_interval(self):
        """由实际时间戳推导采样间隔与时间原点。

        采样间隔 = 相邻历元时间差的中位数；
        时间原点 = --time_origin 指定值，未指定时取最早的历元。
        """
        eps = sorted(self.epochs)
        parsed = []
        for e in eps:
            try:
                parsed.append(datetime.datetime.strptime(e, '%Y-%m-%d %H:%M:%S'))
            except ValueError:
                continue

        if self.interval_override:
            try:
                self.interval = float(self.interval_override)
            except (TypeError, ValueError):
                self.interval = None
        else:
            diffs = []
            for a, b in zip(parsed, parsed[1:]):
                d = (b - a).total_seconds()
                if d > 0:
                    diffs.append(d)
            self.interval = statistics.median(diffs) if diffs else 30.0

        if self.origin is None:
            self.origin = parsed[0] if parsed else None
            self.origin_auto = True

    # ---------- 输出 ----------
    def _stat_columns(self, kind):
        """按统计类型给出该模块的输出列名。

        每个模块只出一列主统计量，列名**不带统计量后缀**：
          level (CNR)  → 均值
          resid (MP…) → RMS
          flag (CSR…) → 计数
        实际列名由 _stat_label() 给出。
        """
        if kind == 'level':
            return ['mean']
        if kind == 'resid':
            return ['rms']
        return ['cnt']                # flag：计数

    def _row_value(self, mname, kind, acc):
        cols = self._stat_columns(kind)
        if acc is None or acc.n == 0:
            return [''] * len(cols)
        if kind == 'level':
            return [_fmt(acc.mean, 3)]
        if kind == 'resid':
            return [_fmt(acc.rms, 4)]
        # flag：计数（非零样本数）
        return [str(acc.nz)]

    # ---------- 完整率（弧段口径，整列恒定） ----------
    def _comp_rate(self, prn=None, band_key=None):
        """
        完整率（0~1 小数）：有值历元数 / 可见弧段历元数。

        prn=None                → 跨全部卫星、跨其各频段取均值
                                   （先每星各频段均值，再跨星均值）
        prn 指定、band_key 指定 → 该星该频段的完整率
        prn 指定、band_key 为空 → 该星各频段完整率的均值

        该指标是**整段弧段**的统计量，不随历元变化，故同一列各行取值相同；
        也不按频段分列——双频完整率通常几乎相同，拆列无意义。
        """
        if not self.comp_ratio:
            return None
        if prn is not None:
            if band_key:
                d = self.comp_ratio.get((prn, band_key))
                if not d or not d['arc']:
                    return None
                return d['obs'] / d['arc']
            # 该星跨频段取均值
            vals = []
            for (p, _bk), d in self.comp_ratio.items():
                if p == prn and d['arc']:
                    vals.append(d['obs'] / d['arc'])
            if not vals:
                return None
            return sum(vals) / len(vals)

        # 跨卫星、跨频段取均值：先对每颗星的各频段取均值，再对卫星取均值
        per_prn = {}
        for (p, _bk), d in self.comp_ratio.items():
            if not d['arc']:
                continue
            per_prn.setdefault(p, []).append(d['obs'] / d['arc'])
        if not per_prn:
            return None
        prn_means = [sum(v) / len(v) for v in per_prn.values() if v]
        if not prn_means:
            return None
        return sum(prn_means) / len(prn_means)

    def _sat_module_bands(self, prn, mname):
        """单星表中某模块的频段列布局。

        返回 [(band_key, 列后缀), ...]，按 SPP 选段顺序排列。
        band_key 为 '' 表示该模块与频段无关（周跳等），只出一列。
        完整率（COMP）整段弧段统计、双频取值几乎相同，**不按频段分列**，
        统一出一列 COMP（值为该星各频段完整率的均值）。
        """
        if mname == 'COMP':
            return [('', '')] if self.sat_mod_bands.get((prn, mname)) else []
        bands = self.sat_mod_bands.get((prn, mname))
        if not bands:
            return []
        order = self._sat_band_keys(prn)
        out = []
        # 与频段无关的模块：band_key = ''
        if '' in bands:
            out.append(('', ''))
        for key in order:
            if key in bands:
                out.append((key, f'_{key}'))
        # SPP 选段之外的其它频段（保险起见也列出）
        for key in sorted(bands):
            if key and key not in order:
                out.append((key, f'_{key}'))
        return out

    def _write_qcvalue(self, out_file, prn=None):
        """写逐历元时间序列表。

        prn=None  → 全卫星汇总（QCValue_<标签>.csv），各频段合并统计
        prn='C19' → 只含该星（QCValue_C19_<标签>.csv），**按频段分列**：
                     双频各出一组统计量，不做跨频段平均；
                     与频段无关的模块（周跳等）只出一列。
        """
        ensure_dir(os.path.dirname(out_file) or '.')
        self._derive_interval()
        if self.origin is None:
            self.origin = datetime.datetime.now()

        # 逐星时只保留该星有数据的历元
        if prn is None:
            epochs = sorted(self.epochs)
        else:
            epochs = sorted(self.sat_epochs.get(prn) or ())

        # ---- 列布局 ----
        # layout: [(mname, kind, band_key, 频段后缀), ...]
        layout = []
        for _suffix, mname, kind in MODULE_SPECS:
            if prn is None:
                layout.append((mname, kind, None, ''))
            else:
                for band_key, sfx in self._sat_module_bands(prn, mname):
                    layout.append((mname, kind, band_key, sfx))

        # 列名 = 模块指标名 + 频段后缀（不带统计量后缀）
        # STATION = 测站名（观测文件名前 4 位），逐历元标记数据来源
        header_cols = ['STATION', 'Epoch', 'Idx']
        for mname, kind, band_key, sfx in layout:
            header_cols.append(f'{mname}{sfx}')
        # 卫星数量（整站量，与卫星无关，两表一致）
        # 只出一列 SAT：统计范围内各系统的卫星数之和
        # （--systems 限定单系统时，即为该系统的卫星数）
        has_sat = bool(self.sat_count)
        if has_sat:
            header_cols.append('SAT')
        # SPP 逐历元解算结果（整站量，与卫星无关，两表一致），置于最右
        has_spp = bool(self.spp_enu)
        if has_spp:
            header_cols += ['SPP_E', 'SPP_N', 'SPP_U', 'SPP_3D']

        with open(out_file, 'w', newline='', encoding='utf-8') as f:
            f.write('# Inf       : GNSS QC Time Series (per-epoch)\n')
            f.write(f'# Time      : {datetime.datetime.now():%Y-%m-%d %H:%M:%S}\n')
            f.write(f'# Origin    : {self.origin:%Y-%m-%d %H:%M:%S}'
                    '   (Idx 序号的时间原点；'
                    + ('由数据中最早历元自动推导' if self.origin_auto else '命令行指定')
                    + ')\n')
            f.write(f'# Interval  : {self.interval:g} s'
                    '   (由实际时间戳中位差分推导)\n')
            f.write('# Scope     : 仅统计「计算 SPP 时实际选用的频段」'
                    '（各文件 _sppband.csv）\n')
            if prn is None:
                f.write('# Rows      : 逐历元一行，时间戳取自各模块明细 CSV\n')
            else:
                f.write(f'# Satellite : {prn}（--sat 选星，仅该星数据）\n')
                f.write('# Rows      : 逐历元一行，仅含该星有数据的历元\n')
                f.write('# Band      : 指标按频段分列（后缀为频段键），'
                        '双频不合并；无后缀者为与频段无关的指标\n')
            for _suffix, mname, kind in MODULE_SPECS:
                if mname == 'COMP':
                    f.write('# COMP   : 弧段完整率(0~1)，整段弧段统计量，各行取值相同\n')
                    continue
                f.write(f'# {mname:<5} : '
                        f'{"均值" if kind == "level" else ("RMS" if kind == "resid" else "计数")}\n')
            if has_spp:
                f.write('# SPP    : 逐历元单点定位结果，E/N/U 为相对参考原点的'
                        '东/北/天残差(m)，3D 为 sqrt(E²+N²+U²)(m)\n')
            if has_sat:
                f.write('# SAT    : 逐历元统计范围内各系统的可见卫星数之和\n')
            f.write('#             END OF HEADER\n')

            w = csv.writer(f)
            w.writerow(header_cols)

            # 全卫星表的完整率：跨卫星、跨频段取均值（整列恒定）
            comp_rate_all = self._comp_rate() if prn is None else None

            for ep in epochs:
                try:
                    dt = datetime.datetime.strptime(ep, '%Y-%m-%d %H:%M:%S')
                except ValueError:
                    dt = None
                rel = (dt - self.origin).total_seconds() if dt else None
                idx = int(round(rel / self.interval)) if rel is not None and self.interval else ''
                row = [self.epoch_station.get(ep, ''), ep, idx]
                for mname, kind, band_key, sfx in layout:
                    if mname == 'COMP':
                        # 完整率是整段弧段的统计量，整列恒定
                        if prn is None:
                            rate = comp_rate_all
                        else:
                            rate = self._comp_rate(prn, band_key)
                        row.append(_fmt(rate, 4) if rate is not None else '')
                        continue
                    if prn is None:
                        acc = self.acc.get((mname, ep))
                    else:
                        acc = self.sat_acc.get((prn, mname, band_key, ep))
                    row.extend(self._row_value(mname, kind, acc))
                # 卫星数量（整站量）：统计范围内各系统之和
                # （--systems 限定单系统时，即为该系统的卫星数）
                if has_sat:
                    counts = self.sat_count.get(ep) or {}
                    row.append(str(sum(v for k, v in counts.items()
                                       if k != 'ALL')) if counts else '')
                # SPP 逐历元解算结果（整站量）
                if has_spp:
                    enu = self.spp_enu.get(ep)
                    if enu is None:
                        row += ['', '', '', '']
                    else:
                        e, n, u = enu
                        row += [_fmt(e, 3), _fmt(n, 3), _fmt(u, 3),
                                _fmt(math.sqrt(e * e + n * n + u * u), 3)]
                w.writerow(row)
        return out_file

    def _write_sppband_summary(self, out_file):
        """SPP 选段 + 频点统计汇总。"""
        with open(out_file, 'w', encoding='utf-8') as f:
            f.write('# Inf       : SPP Band Statistics Summary\n')
            f.write(f'# Time      : {datetime.datetime.now():%Y-%m-%d %H:%M:%S}\n')
            f.write('# SPP band  : 各文件实际用于 SPP 的伪距码对（合并去重）\n')
            f.write('#             END OF HEADER\n')
            f.write('+SYS\n')
            f.write(' gSys  Codes                f1(MHz)    f2(MHz)    DF(MHz)'
                    '    LamW(m)   MWamp   NKey  NPrn   NObs\n')
            f.write('-' * 96 + '\n')
            for gsys in sorted(self.spp_info):
                info = self.spp_info[gsys]
                codes = info['codes']
                c1 = codes[0] if codes else ''
                c2 = codes[1] if len(codes) > 1 else ''
                f1 = _band_freq_mhz(gsys, c1) if c1 else None
                f2 = _band_freq_mhz(gsys, c2) if c2 else None
                df = abs(f1 - f2) if (f1 and f2) else None
                lam = CLIGHT / (df * 1e6) if df else None
                amp = ((f1 + f2) / df) if df else None
                nkey = nprn = nobs = 0
                prn_set = set()
                for key in info['keys']:
                    u = self.band_usage.get((gsys, key))
                    if not u:
                        continue
                    nkey += 1
                    prn_set.update(u['prns'].keys())
                    nobs += u['nobs']
                nprn = len(prn_set)
                f.write(f'   {gsys:<4}{"+".join(codes):<20}'
                        f'{_fmt(f1, 3):>10}{_fmt(f2, 3):>11}{_fmt(df, 3):>11}'
                        f'{_fmt(lam, 4):>11}{_fmt(amp, 2):>9}'
                        f'{nkey:>6}{nprn:>6}{nobs:>7}\n')
            f.write('-' * 96 + '\n')
            f.write('-SYS\n')
            f.write('\n')

            # ---- 逐 PRN 的 SPP 频段观测统计 ----
            f.write('+SAT\n')
            f.write('  PRN    SPPBand  PhaseCode  Freq(MHz)   ObsNum'
                    '   Ratio(%)   ValidEpoch\n')
            f.write('-' * 78 + '\n')
            total_epochs = max(1, len(self.epochs))
            for (prn, code), cnt in sorted(self.spp_freq_counts.items()):
                gsys = prn[0]
                key = code[1:] if len(code) >= 2 else ''
                lcode = 'L' + key
                fq = _band_freq_mhz(gsys, lcode)
                vd = (self.band_valid.get((gsys, key)) or {}).get(prn, 0)
                f.write(f'  {prn:<6}{code:>9}{lcode:>11}{_fmt(fq, 3):>12}'
                        f'{cnt:>10}{cnt / total_epochs * 100:>11.2f}'
                        f'{vd:>12}\n')
            f.write('-' * 78 + '\n')
            f.write('-SAT\n')
        return out_file

    # ---------- 主流程 ----------
    def run(self, obs_files):
        total = len(obs_files)
        for idx, obs_file in enumerate(obs_files, 1):
            base = os.path.splitext(os.path.basename(obs_file))[0]
            sub_dir = os.path.join(self.args.output_dir, base)
            print(f'\n========== 处理文件 [{idx}/{total}]: {obs_file} ==========')

            spp_csv = os.path.join(sub_dir, f'{base}_sppband.csv')
            spp_freq_csv = os.path.join(sub_dir, f'{base}_sppband_freq.csv')
            station_csv = os.path.join(sub_dir, f'{base}_QCValue_station.csv')
            codes_by_sys = {}

            if os.path.exists(station_csv) and os.path.exists(spp_csv):
                # 断点续跑：复用已有结果与 SPP 选段
                print('   (检测到已有汇总结果，跳过重复计算)')
                codes_by_sys, _rows = self._read_spp_csv(spp_csv)
            elif os.path.exists(station_csv):
                # 已有 QC 结果但缺 SPP 选段记录：只读观测头重算选段
                print('   (检测到已有汇总结果，重算 SPP 选段)')
                codes_by_sys = self._resolve_bands_from_head(obs_file, sub_dir)
                if codes_by_sys:
                    self._write_spp_csv(sub_dir, base, codes_by_sys)
            else:
                nav_data = self._nav_for(obs_file)
                if nav_data is None and getattr(self.args, 'nav', None):
                    print('   (无导航数据，SPP 与完整率可见性口径将被跳过)')
                proc = QCProcessor(obs_file, sub_dir, self.args, nav_data=nav_data)
                try:
                    proc.process()
                except Exception as e:
                    print(f'文件 {obs_file} 处理失败，已跳过: {e}')
                    continue
                codes_by_sys = self._spp_codes(getattr(proc, 'spp_band_choose', None))
                self._write_spp_csv(sub_dir, base, codes_by_sys)
                self._write_spp_freq_csv(sub_dir, base, codes_by_sys,
                                         getattr(proc, 'freq_data', None))

            if not codes_by_sys:
                print('   (该文件无 SPP 选段记录，不参与频段过滤统计)')
                continue

            # 合并全局 SPP 信息
            for gsys, codes in codes_by_sys.items():
                info = self.spp_info.setdefault(gsys, {'codes': [], 'keys': set()})
                for c in codes:
                    if c not in info['codes']:
                        info['codes'].append(c)
                    if len(c) >= 2:
                        info['keys'].add(c[1:])
            self.file_summaries.append((base, codes_by_sys))

            # 合并 SPP 频段的 PRN 级频段计数（freq 模块口径）
            for (prn, code), cnt in self._read_spp_freq_csv(spp_freq_csv).items():
                self.spp_freq_counts[(prn, code)] = (
                    self.spp_freq_counts.get((prn, code), 0) + cnt)

            # 合并逐历元 SPP 解算结果（ENU 残差），SPP 是整站量，与卫星无关
            enu = _read_spp_txt(os.path.join(sub_dir, f'{base}_spp.txt'))
            if enu:
                self.spp_enu.update(enu)
            else:
                print('   (无 SPP 解算结果，跳过 SPP 列)')

            # 合并逐历元各系统可见卫星数（整站量，与卫星无关）
            # 各系统列已按 --systems 过滤，汇总表只出一列 SAT（范围内求和）
            sat = _read_sat_csv(os.path.join(sub_dir, f'{base}_sat.csv'))
            if sat:
                self.sat_count.update(sat)
            else:
                print('   (无卫星数量统计，跳过卫星数列)')

            # 逐历元汇总（station = 文件名前 4 位的测站名）
            station = _extract_station_from_filename(obs_file)
            self._aggregate_file(sub_dir, base, codes_by_sys, station=station)

        # ---- 输出 ----
        if not self.epochs:
            print('无有效历元数据，跳过 QCValue 汇总表输出')
            return

        # 标签 = 各系统 SPP 频段标识合并（BDS C2I + Galileo C1C + GPS C1C
        #      → C2I_E1C_G1C）
        tag = self._band_tag()
        suffix = f'_{tag}' if tag else ''

        out_qcvalue = os.path.join(self.args.output_dir, f'QCValue{suffix}.csv')
        self._write_qcvalue(out_qcvalue)
        print(f'\n-->时间序列质量表: {out_qcvalue}  ({len(self.epochs)} 个历元)')

        # ---- 选星参数 --sat：单星汇总表 ----
        for prn in self.sat_filter:
            eps = self.sat_epochs.get(prn)
            if not eps:
                print(f'   提示：{prn} 在 SPP 频段内无数据，跳过单星汇总表')
                continue
            st = self._sat_tag(prn)
            sfx = f'_{st}' if st else ''
            out_sat = os.path.join(self.args.output_dir,
                                   f'QCValue_{prn}{sfx}.csv')
            self._write_qcvalue(out_sat, prn=prn)
            print(f'-->单星时间序列: {out_sat}  ({len(eps)} 个历元)')

        out_band = os.path.join(self.args.output_dir, 'QCValue_sppband.txt')
        self._write_sppband_summary(out_band)
        print(f'-->SPP频段统计汇总: {out_band}')

        print(f'\n所有文件处理完成，结果保存在: {self.args.output_dir}')
