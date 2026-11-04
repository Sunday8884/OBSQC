import pandas as pd
import numpy as np
from pathlib import Path

# ========== 配置 ==========
INPUT_FILE = Path(r"F:\GNSS_data\O\2024\1\qc\C\QCValue_C2I_6I_7I.csv")
OUTPUT_FILE = Path(r"F:\OBSQC-eval\test\QCValue_C.csv")
# INPUT_FILE = Path(r"F:\GNSS_data\O\2024\1\qc\C\QCValue_C30_C2I_6I_7I.csv")
# OUTPUT_FILE = Path(r"F:\OBSQC-eval\test\QCValue_C30.csv")
MAX_GAP = 0         # 连续缺失 > MAX_GAP 的缺口不插值，保留空值
START_DATE =None    # 起始日期，形如 "2024-01-01"；None = 不限
END_DATE   =None    # 截止日期，形如 "2024-12-31"；None = 不限
NON_INTERP_COLS = ("STATION", "Epoch", "Idx",'SPP_E', 'SPP_N', 'SPP_U')  # 不参与插值


# ========== 绘图配置 ==========
PLOT_TIMESERIES = True       # 是否绘制时序曲线图（每列一个子图，共享 x 轴）
TS_COLS = ["CNR", "CSR", "MP", "CMC", "IOD", "PN", "HOD", "COMP", "SAT", "SPP_3D"] # 可自动展开为频段分列：CNR_2I
LOG_SCALE_COLS = ("")   # 需要对数纵轴的列
TS_TRIM_COLS = ("SPP_3D",)  # 需要截断极端离群值的列
TS_TRIM_Q = 0.99        # 截断分位数，None = 不截断
TS_OUT = Path(r"F:\OBSQC-eval\test\time.png")

PLOT_HIST = True             # 是否绘制 SPP_3D 频数统计图
HIST_COL = "SPP_3D"          # 频数统计目标列
HIST_BINS = 100              # 分箱数
HIST_LOGX = False            # SPP_3D 对数横轴
HIST_TRIM_Q = 0.99           # 线性横轴时按此分位数截断极端离群值
HIST_OUT = Path(r"F:\OBSQC-eval\test\spp3d.png")
# ==========================

def interpolate_gap_limited(series: pd.Series, max_gap: int) -> pd.Series:
    """对单列做线性插值，但以下两类缺口保持空值（不填充）：

    1) 连续缺失个数 > max_gap 的内部长缺口 —— 超过阈值不插值；
    2) 开头 / 结尾的边界缺口 —— 没有双侧已知点，不属于"插值"语义，一律留空。

    只保留"两端都有有效值且连续缺失 ≤ max_gap"的内部缺口做线性插值。
    做法：用 (~isna).cumsum() 给每段连续 NaN 编号并统计段长；用 cummax 标记
    首尾边界缺口；先插值整列，再把上述两类缺口 mask 回 NaN。

    注意：本函数只对**单站连续序列**有效。跨站调用会把一个站的缺口
    用相邻站的值填充，故必须配合 interpolate_by_station 使用。
    """
    is_na = series.isna()
    if not is_na.any():
        return series
    # 内部缺口每段长度
    grp = (~is_na).cumsum()                       # 每段连续 NaN 共享同一编号
    run_len = is_na.groupby(grp).transform("sum")   # 段内 NaN 个数
    # 边界缺口：开头直到首个有效值、结尾从末个有效值起
    leading = (~is_na).cummax()                     # 首个有效值之后才为 True
    trailing = (~is_na)[::-1].cummax()[::-1]        # 末个有效值之前才为 True
    edge = ~(leading & trailing)                    # 开头或结尾的 NaN 段
    protected = is_na & ((run_len > max_gap) | edge)
    filled = series.interpolate(method="linear", limit_direction="both")
    return filled.mask(protected)


def interpolate_by_station(df: pd.DataFrame, col: str, max_gap: int,
                            station_col: str = "STATION") -> pd.Series:
    """按测站分组做限长插值，杜绝跨站插值污染。

    多站数据里，各站的时间戳通常完全重合（如同一天并行观测）。
    若对整列直接插值，某站的缺口会被**相邻站**的数值填上，
    得到完全错误的数据。故必须逐站独立插值。

    无 STATION 列时退化为整列插值（单站场景）。
    """
    if station_col not in df.columns:
        return interpolate_gap_limited(df[col], max_gap)
    out = df[col].copy()
    for _st, idx in df.groupby(station_col).groups.items():
        out.loc[idx] = interpolate_gap_limited(df.loc[idx, col], max_gap)
    return out

def drop_rows_with_gaps(df: pd.DataFrame, cols):
    """剔除插值后仍含空值的行 —— **唯一的剔除规则**。

    cols 为**参与插值的列**（不含 SPP_* 等不插值列）：
    这些列插值完成后若仍有空值，说明该历元该指标确实没有数据
    （整列全空的指标会自然导致该站所有行被剔除，无需额外规则）。

    不参与插值的列（如 SPP_E/N/U/3D）不作为剔空依据，
    它们的空值原样保留。

    返回 (清理后 DataFrame, 移除行数)。
    """
    cols = [c for c in cols if c in df.columns]
    if not cols:
        return df, 0
    mask = df[cols].isna().any(axis=1)
    removed = int(mask.sum())
    return df.loc[~mask].reset_index(drop=True), removed


def read_csv_safe(path: Path) -> pd.DataFrame:
    """自动尝试常见编码读取 CSV，并跳过多余注释行（以 # 开头的行）。"""
    for enc in ("utf-8-sig", "utf-8", "gbk", "gb18030", "latin1"):
        try:
            return pd.read_csv(path, encoding=enc, comment="#")
        except UnicodeDecodeError:
            continue
    return pd.read_csv(path, encoding="utf-8", errors="replace", comment="#")

def filter_by_epoch(df: pd.DataFrame, start=None, end=None):
    """按 Epoch 的日期(YYYY-MM-DD)截断数据。

    start / end 形如 "2024-01-01" 或 None（不限制该端）。
    只比较"年-月-日"部分，忽略日内时分秒；区间两端均包含（闭区间）。
    返回 (过滤后 DataFrame, 被移除行数)；解析失败的 Epoch 视为不在区间内。
    """
    if start is None and end is None:
        return df, 0
    ep = pd.to_datetime(df["Epoch"], errors="coerce")
    day = ep.dt.floor("D")                       # 归一化到当日 00:00，仅保留日期
    mask = pd.Series(True, index=df.index)
    if start is not None:
        sd = pd.Timestamp(start).floor("D")
        mask &= day >= sd
    if end is not None:
        ed = pd.Timestamp(end).floor("D")
        mask &= day <= ed
    removed = int((~mask).sum())
    return df.loc[mask].reset_index(drop=True), removed

def setup_cjk_font():
    """尝试设置中文字体，失败则静默回退（避免缺字方块）。返回使用的字体名或 None。

    字体设为**列表**，但选取顺序有讲究：

    matplotlib 的 log 轴刻度标签是 mathtext 格式 `$\\mathdefault{10^{-2}}$`，
    其中指数里的负号会被 mathtext 渲染成 U+2212（数学减号），
    而 `\\mathdefault` 强制使用 rcParams 的默认字体。
    SimHei / SimSun / DengXian 等中文字体**缺 U+2212 字形**，
    于是报 "Font 'default' does not have a glyph for '\\u2212'" 并渲染成方块。

    注意 `font.sans-serif` 列表**不会**触发逐字形 fallback
    （把 DejaVu Sans 放前面会让中文变成方块），
    因此必须**选用本身就含 U+2212 的中文字体**：
    Microsoft YaHei / Noto Sans SC / STSong 均可，SimHei 不行。
    """
    from matplotlib import font_manager
    from matplotlib.ft2font import FT2Font
    # 按优先级排列；含 U+2212 的排在前面，缺字形的仅作后备
    candidates = ["Microsoft YaHei", "Noto Sans SC", "STSong",
                  "WenQuanYi Micro Hei", "SimHei"]
    avail = {f.name for f in font_manager.fontManager.ttflist}

    def _has_minus(name):
        """该字体是否含 U+2212（数学减号）字形。"""
        try:
            path = font_manager.findfont(
                font_manager.FontProperties(family=name))
            return FT2Font(path).get_char_index(0x2212) != 0
        except Exception:
            return False

    usable = [n for n in candidates if n in avail]
    # 含 U+2212 的排前面；否则 log 轴刻度会报缺字形警告
    usable.sort(key=lambda n: 0 if _has_minus(n) else 1)
    for name in usable:
        import matplotlib.pyplot as plt
        plt.rcParams["font.family"] = "sans-serif"
        plt.rcParams["font.sans-serif"] = [name] + [
            f for f in ("DejaVu Sans", "SimHei") if f != name]
        plt.rcParams["axes.unicode_minus"] = False
        return name
    return None


def resolve_ts_cols(df: pd.DataFrame, bases) -> list:
    """把基名（如 CNR）展开为实际列名：精确匹配，或前缀匹配 `base+'_'`
    （处理 CNR_2I/CNR_6I 等分频列）。去重保序，返回实际存在的列列表。
    """
    out, seen = [], set()
    for b in bases:
        if b in df.columns:
            matches = [b]
        else:
            matches = [c for c in df.columns if c.startswith(b + "_")]
        for c in matches:
            if c not in seen:
                seen.add(c)
                out.append(c)
    return out


def _trim_range(vals, trim_q):
    """按分位数求出用于显示的数值范围，截去极端离群值。

    vals   : 待统计的数值序列（可含 NaN）
    trim_q : 分位数（如 0.99）；None 或非 (0,1) 区间时不做截断

    **自动识别重尾方向**：比较上下两侧的尾部长度
    （上界侧 max-p99 vs 下界侧 p1-min），截掉更长的一侧。
    这样对正向重尾（SPP_3D 等距离类指标，主体在小值、偶发极大值）
    截上界；对全负或负向重尾的数据则截下界。
    单靠固定截上界在全负数据上会失效——那类数据的主体在负区间。

    返回 (lo, hi, n_cut)：lo/hi 为显示范围，n_cut 落在范围外的样本数。
    """
    v = pd.to_numeric(vals, errors="coerce").dropna()
    if v.empty:
        return 0.0, 1.0, 0
    lo, hi = float(v.min()), float(v.max())
    if trim_q is not None and 0 < trim_q < 1:
        p_hi = float(v.quantile(trim_q))
        p_lo = float(v.quantile(1 - trim_q))
        if (hi - p_hi) >= (p_lo - lo):
            hi = p_hi                      # 上侧尾部更长 → 截上界
        else:
            lo = p_lo                      # 下侧尾部更长 → 截下界
    n_cut = int(((v < lo) | (v > hi)).sum())
    return lo, hi, n_cut


def _trim_range_by_station(df, col, trim_q, station_col="STATION"):
    """**分站**计算分位数，取各站范围的并集作为全局显示范围。

    为什么不能直接用全局分位数：多站混合数据的分位数是按总体算的，
    水平偏低的测站会被误判成离群。实测 CNR 全体分位截断把 JFNG 站
    6.1%、ABPO 站 2.7% 的**正常**低值切掉了（这两站 CNR 整体偏低，
    并非异常）。

    分站计算：每站各按自己的 trim_q 分位取范围，
    全局 lo = 各站 lo 的最小值、hi = 各站 hi 的最大值（并集）。
    这样每站只截掉自己最极端的那部分，所有站的正常区间都完整保留。

    返回 (lo, hi, n_cut)。
    """
    v = pd.to_numeric(df[col], errors="coerce")
    if station_col not in df.columns or trim_q is None:
        return _trim_range(v, trim_q)
    los, his = [], []
    for _st, idx in df.groupby(station_col).groups.items():
        s = v.loc[idx]
        if s.notna().sum() == 0:
            continue
        lo, hi, _ = _trim_range(s, trim_q)
        los.append(lo)
        his.append(hi)
    if not los:
        return _trim_range(v, trim_q)
    lo, hi = min(los), max(his)
    n_cut = int(((v < lo) | (v > hi)).sum())
    return lo, hi, n_cut


def station_palette(stations, cmap_name="tab20"):
    """为各测站分配固定颜色：同一测站在所有子图中颜色一致，便于横向对比。

    测站数较少时用 tab20（色块区分度高）；超过 20 个改用 viridis 连续色带。
    """
    import matplotlib
    import matplotlib.pyplot as plt
    stations = list(stations)
    if len(stations) <= 20:
        cmap = matplotlib.colormaps.get_cmap(cmap_name)
        return {st: cmap(i % 20) for i, st in enumerate(stations)}
    cmap = matplotlib.colormaps.get_cmap("viridis")
    n = max(1, len(stations) - 1)
    return {st: cmap(i / n) for i, st in enumerate(stations)}


def plot_timeseries(df: pd.DataFrame, cols, out_path, title="时序曲线"):
    """绘制指定列的时序曲线图：横轴 Epoch，每列一个子图，共享 x 轴。

    各指标量纲差异大（如 CNR≈40、MP≈0、SPP_3D 可达数千），故每列独立面板、
    各自 y 轴，仅共享时间轴，避免互相压扁。cols 为基名，自动展开为频段分列。

    **按测站分线**：多站数据的同一时刻有多个不同值，若画成一条线会来回穿插
    严重混淆。此处每个测站一条线，**颜色在所有子图中固定**，便于横向对比；
    图例统一放在整图底部一行。
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.dates as mdates

    cols = resolve_ts_cols(df, cols)
    if not cols:
        print("[plot] 时序图：无可用列，跳过")
        return
    sub = df.copy()
    sub["_t"] = pd.to_datetime(sub["Epoch"], errors="coerce")
    sub = sub.dropna(subset=["_t"])
    if sub.empty:
        print("[plot] 时序图：时间列解析为空，跳过")
        return

    n = len(cols)
    ncols = 2
    nrows = (n + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols, figsize=(14, 3.0 * nrows), sharex=True)
    axes = np.atleast_1d(axes).ravel()

    # 按测站分组，每个测站一条线、颜色固定（各子图一致，便于横向对比）
    has_station = "STATION" in sub.columns
    if has_station:
        stations = sorted(sub["STATION"].astype(str).unique())
        colors = station_palette(stations)
    else:
        stations, colors = [], {}

    for ax, col in zip(axes, cols):
        if has_station:
            for st in stations:
                g = sub[sub["STATION"].astype(str) == st]
                if g.empty:
                    continue
                ax.plot(g["_t"], g[col], lw=0.7, color=colors[st],
                        label=st if ax is axes[0] else None)
        else:
            ax.plot(sub["_t"], sub[col], lw=0.7, color="#2a6f97")
        ax.set_ylabel(col)
        ax.grid(True, ls=":", alpha=0.4)
        ax.xaxis.set_major_locator(mdates.AutoDateLocator())
        ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y-%m-%d\n%H:%M"))
        # 纵轴处理是三条**独立**规则：
        #   对数轴：由 LOG_SCALE_COLS 指定（该列有正值时生效）
        #   截断  ：由 TS_TRIM_COLS 指定，按 TS_TRIM_Q **分站**截去极端离群值
        vals = pd.to_numeric(sub[col], errors="coerce")
        pos = vals[vals > 0]
        use_log = col in LOG_SCALE_COLS and not pos.empty
        if use_log:
            ax.set_yscale("log")
            ax.set_ylabel(f"{col} (log)")
        elif col in TS_TRIM_COLS:
            lo, hi, n_cut = _trim_range_by_station(sub, col, TS_TRIM_Q)
            # 常量列（lo == hi）不设 ylim，否则 matplotlib 会警告上下限相同
            if hi > lo:
                ax.set_ylim(lo, hi)
            if n_cut:
                ax.set_ylabel(f"{col} (截断 {n_cut} 个离群值)")
    for ax in axes[n:]:
        ax.set_visible(False)
    axes[0].set_title(title)

    # 图例统一放在整图底部一行（各子图不重复）
    handles, labels = axes[0].get_legend_handles_labels()
    fig.autofmt_xdate(rotation=30)
    if handles:
        # 图例按 ncol 分行，底部预留足够高度，避免最后一行被裁切
        ncol = min(8, len(handles))
        nrow = -(-len(handles) // ncol)          # 向上取整
        bottom = 0.025 + 0.035 * nrow
        fig.tight_layout(rect=(0, bottom, 1, 1))
        fig.legend(handles, labels, loc="lower center", ncol=ncol,
                   frameon=False, fontsize=8)
    else:
        fig.tight_layout()
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    print(f"[plot] 时序图已保存: {out_path}  (列: {cols}, 测站: {len(stations)})")


def plot_histogram(df: pd.DataFrame, col, out_path, bins=50, logx=True, title=None,
                   trim_q=None):
    """绘制指定列的频数统计图（直方图）。

    SPP_3D 呈重尾分布（常规 ~2 m，异常可达上万 m），logx=True 时用对数等距 bins，
    能同时看清常态与极端离群区间。并标注中位数/均值竖线。

    **线性横轴时的截断**：重尾数据的极端值会把坐标轴拉长，
    使 99% 的样本全挤在第一根柱子里，图完全不可读。
    故线性轴下按分位数 trim_q **截断上界**（默认 p99），
    只统计 [lo, hi] 区间内的样本，并在标题中注明被截断的数量。
    （只截上界：下界的小值是正常精度，不是离群值）

    trim_q=None 表示不截断。
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    if col not in df.columns:
        print(f"[plot] 频数图：列 {col} 不存在，跳过")
        return
    s = pd.to_numeric(df[col], errors="coerce").dropna()
    if s.empty:
        print(f"[plot] 频数图：列 {col} 无有效数值，跳过")
        return

    fig, ax = plt.subplots(figsize=(9, 5))
    trimmed = 0
    if logx:
        pos = s[s > 0]
        lo = pos.min() if len(pos) else s.min()
        hi = s.max()
        edges = np.logspace(np.log10(lo), np.log10(hi), bins + 1)
        ax.set_xscale("log")
        ax.hist(s, bins=edges, color="#2a6f97", edgecolor="white")
    else:
        # 线性横轴：双侧截断极端离群值，保证图可读
        lo, hi, trimmed = _trim_range(s, trim_q)
        kept = s[(s >= lo) & (s <= hi)]
        # 显式给边界，使分箱在 [lo, hi] 上等距，主体分布不被挤压
        edges = np.linspace(lo, hi, bins + 1)
        ax.hist(kept, bins=edges, color="#2a6f97", edgecolor="white")
        ax.set_xlim(lo, hi)
    ax.set_xlabel(col)
    ax.set_ylabel("频数 (count)")
    ax.grid(True, axis="y", ls=":", alpha=0.4)
    med, mean, p95 = s.median(), s.mean(), s.quantile(0.95)
    ax.axvline(med, color="#d95f02", ls="--", lw=1.2, label=f"median={med:.3f}")
    ax.axvline(mean, color="#1b9e77", ls="--", lw=1.2, label=f"mean={mean:.3f}")
    ax.legend()
    ttl = title or f"{col} 频数分布"
    if trimmed:
        ttl += f"（已截断 {trimmed} 个 >p{trim_q * 100:g} 的离群样本，"\
               f"最大 {s.max():.1f}）"
    ax.set_title(ttl)
    fig.tight_layout()
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    print(f"[plot] 频数图已保存: {out_path}  (n={len(s)}, median={med:.3f}, "
          f"mean={mean:.3f}, p95={p95:.3f})"
          + (f", 已截断 {trimmed} 个离群样本" if trimmed else ""))


def main():
    df = read_csv_safe(INPUT_FILE)
    total_before = len(df)
    print(f"读取 {INPUT_FILE.name}: {df.shape[0]} 行 × {df.shape[1]} 列")

    # 时间截断：按 Epoch 日期过滤（在开始/截止日期之间）
    if START_DATE is not None or END_DATE is not None:
        df, removed_time = filter_by_epoch(df, START_DATE, END_DATE)
        print(f"时间截断(Epoch∈[{START_DATE or '∞'}, {END_DATE or '∞'}]): "
              f"移除 {removed_time} 行 → 剩余 {len(df)} 行")

    # 数值列 = 数值类型且不在排除列表
    #   num_cols：参与插值的列（SPP_* 不插值），也是唯一的剔空判定依据
    num_cols = [c for c in df.columns
                if pd.api.types.is_numeric_dtype(df[c]) and c not in NON_INTERP_COLS]
    print(f"插值列({len(num_cols)}): {num_cols}")

    miss_before = int(df[num_cols].isna().sum().sum())
    # 按测站分组插值：多站数据各站时间戳重合，整列插值会跨站取值污染
    for c in num_cols:
        df[c] = interpolate_by_station(df, c, MAX_GAP)
    miss_after_interp = int(df[num_cols].isna().sum().sum())
    print(f"线性插值填充(按测站): 缺失 {miss_before} → 剩余 {miss_after_interp} "
          f"(已插值 {miss_before - miss_after_interp} 个；长缺口 {miss_after_interp} 个保留空)")

    # 剔空：唯一的剔除规则 —— 插值后仍有空值的行剔除
    df, removed = drop_rows_with_gaps(df, num_cols)
    print(f"剔空(插值列仍有空值的行): 原始 {total_before} → 剩余 {len(df)} (移除 {removed})")
    if "STATION" in df.columns:
        vc = df["STATION"].value_counts()
        print(f"   保留测站 {len(vc)} 个: " +
              ", ".join(f"{k}={v}" for k, v in vc.sort_index().items()))

    OUTPUT_FILE.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUTPUT_FILE, index=False, encoding="utf-8-sig")
    print(f"输出: {OUTPUT_FILE}  ({df.shape[0]} 行 × {df.shape[1]} 列)")

    # ---- 绘图 ----
    setup_cjk_font()
    if PLOT_TIMESERIES:
        plot_timeseries(df, TS_COLS, TS_OUT)
    if PLOT_HIST:
        plot_histogram(df, HIST_COL, HIST_OUT, bins=HIST_BINS, logx=HIST_LOGX,
                       trim_q=HIST_TRIM_Q)

if __name__ == "__main__":
    main()
