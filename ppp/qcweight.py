# -*- coding: utf-8 -*-
"""观测质量多准则定权（QC-informed weighting，实验性）

流水线（每历元、跨可见卫星），method 切换评分算法：
  gra    （默认）熵权 → 灰色关联系数（参考序列 = 全 1 理想星，分辨系数 ρ）
         → 加权灰色关联度 p_i ∈ (0,1]
  critic CRITIC 客观赋权（对比强度×冲突性，修复熵权忽略指标相关性）替换熵权，
         其余同 gra
  rank   秩变换合成：方向统一后每指标转跨星排名（对异常值免疫），加权平均
         名次归一化为 p_i（权仍用熵权，隔离"合成方式"单变量）
  topsis TOPSIS 贴近度：加权欧氏距离到理想星（全 1）与负理想星（全 0），
         c_i = D⁻/(D⁺+D⁻)（权用熵权；对指标极端值的排序比灰关联稳）
  fixed  AHP 固定权：专家经验固定权（FIXED_W）+ 灰色关联合成，无逐历元
         权重跳动，行为可复现
  game   博弈论组合赋权：熵权与固定权的几何均值组合（兼顾客观区分度与
         主观稳定性）+ 灰色关联合成
  abs    绝对阈值标定法（需调用方维护跨历元标定缓冲，见 abs_calibrate/
         qc_abs_weights）：指标对历史分位数 [q, 1−q] 打绝对分，p 跨历元
         可比——某历元全部星质量都差时相对评分（其余方法）会集体失真，
         绝对评分能正确反映"整体都差"

共用步骤：
  1. 逐星 QC 指标：|ΔMP| / |ΔCMC|（历元差分，消除模糊度常数）/
     SNR（两频段均值，dB-Hz）
  2. 列向 winsorize（鲁棒化，抗单星异常值压扁归一化域）
  3. min-max 正向化（越小越好的指标先反向平移）
  4. 输出 p_i → 观测权 w = p_i（σ² = σ0²/p_i 起点组合，见 ppp.py obsFactor）

说明：
  - PDOP 为历元级常量（同历元各星相同），跨卫星 min-max 归一化退化，
    不纳入指标集；
  - MP/CMC 原值含模糊度常数，逐星状态在周跳/中断后重置，差分自动消除；
  - 首历元（或全部缺失）差分不可用时以列中位数中性填充，全缺失时
    退化为 p=1（中性权，不改变 σ0）；
  - p 的历元间平滑在调用方（ppp.py）实现。
"""
import numpy as np

__all__ = ['qc_grey_weights', 'qc_abs_weights', 'abs_calibrate']

# 越小越好的指标名（min-max 前先反向）：
#   mpd=|ΔMP|  cmcd=|ΔCMC|  pn=相位噪声(GF三阶差分/√40)
#   hod=高阶差分噪声(L1四阶差分/√70)  iod=电离层变化率|GF/dt|
_INVERT = ('mpd', 'cmcd', 'pn', 'hod', 'iod')

# AHP 固定权（专家经验，缺省覆盖 6 指标；abs/fixed/game 共用）：
# SNR 最直接反映链路质量 > 多径类（mpd/cmcd）> 差分噪声类（pn/hod/iod）
FIXED_W = {'snr': 0.30, 'mpd': 0.25, 'cmcd': 0.10,
           'pn': 0.15, 'hod': 0.10, 'iod': 0.10}


def _minmax_positive(cols):
    """cols: 已同向（越大越好）的 1-D 数组列表 → X (n, m)，理想值 1"""
    n = len(cols[0])
    X = np.empty((n, len(cols)))
    for j, c in enumerate(cols):
        c = np.asarray(c, dtype=float)
        rng = c.max() - c.min()
        X[:, j] = 1.0 if rng < 1e-12 else (c - c.min()) / rng
    return X


def _entropy_weights(X):
    """熵权法：X 已正向化 (n, m) → w (m,)"""
    n, m = X.shape
    if n <= 1:
        return np.full(m, 1.0 / m)
    P = X / np.maximum(X.sum(axis=0, keepdims=True), 1e-12)
    P = np.clip(P, 1e-12, None)
    E = (P * np.log(P)).sum(axis=0) / np.log(n)     # e_j ∈ [0,1]
    d = 1.0 - E
    s = d.sum()
    return d / s if s > 1e-12 else np.full(m, 1.0 / m)


def _critic_weights(X):
    """CRITIC 客观赋权：X 已正向化 (n, m) → w (m,)

    w_j ∝ σ_j · Σ_k (1 − r_jk)：对比强度（列标准差）× 冲突性
    （与其他指标的相关系数之和）。熵权法假设指标独立，pn/hod/iod
    高度相关时权重被重复计入；CRITIC 相关性越强冲突性越小，自动去冗余。
    """
    n, m = X.shape
    if n <= 2 or m <= 1:            # 单指标无加权意义；n<3 相关不可估
        return np.full(m, 1.0 / m)
    sigma = X.std(axis=0)
    with np.errstate(invalid='ignore', divide='ignore'):
        R = np.corrcoef(X, rowvar=False)
    # 常数列/单列时 corrcoef 可能返回 0-d 或含 nan
    R = np.atleast_2d(np.asarray(R, dtype=float))
    R = np.where(np.isfinite(R), np.clip(R, -1.0, 1.0), 0.0)
    if R.shape[0] == m:
        np.fill_diagonal(R, 1.0)
    conflict = np.clip(1.0 - R.sum(axis=1), 1e-6, None)
    info = sigma * conflict
    s = info.sum()
    return info / s if s > 1e-12 else np.full(m, 1.0 / m)


def _avg_rank(c):
    """1-D 数组的平均秩（并列取平均），1 = 最小值。"""
    order = np.argsort(c, kind='stable')
    ranks = np.empty(len(c))
    ranks[order] = np.arange(1, len(c) + 1)
    for v in np.unique(c):
        msk = c == v
        if msk.sum() > 1:
            ranks[msk] = ranks[msk].mean()
    return ranks


def _rank_scores(X, w):
    """秩变换合成：X 已正向化 (n, m)（越大越好），w (m,) 指标权 → p (n,)

    每列转跨星平均秩并归一化到 [0,1]（1=最优），加权平均合成。
    秩对异常值单调不变 → 完全免疫离群值；代价是丢失量级信息。
    """
    n, m = X.shape
    if n <= 1:
        return np.ones(n)
    R = np.empty((n, m))
    for j in range(m):
        R[:, j] = (_avg_rank(X[:, j]) - 1.0) / (n - 1)
    return R @ w


def _grey_coeffs(X, rho):
    """灰色关联系数：参考序列 = 全 1（理想星）。返回 ξ (n, m)"""
    D = np.abs(1.0 - X)
    dmin, dmax = D.min(), D.max()
    if dmax < 1e-12:
        return np.ones_like(X)
    return (dmin + rho * dmax) / (D + rho * dmax)


def _fixed_weights(names):
    """AHP 固定权：names 中的指标取 FIXED_W 并归一化（缺失指标权重重分配）。"""
    w = np.array([FIXED_W.get(n, 0.10) for n in names], dtype=float)
    s = w.sum()
    return w / s if s > 1e-12 else np.full(len(names), 1.0 / len(names))


def _game_weights(X, names):
    """博弈论组合赋权：熵权 w1 与固定权 w2 的几何均值组合。

    w_j ∝ √(w1_j · w2_j)——最小化组合权与各权重方案偏差的常用简化
    （二人零和博弈 minimax 的 Nash 均衡解以乘性组合近似）。兼顾熵权的
    逐历元客观区分度与固定权的稳定性。
    """
    w1 = _entropy_weights(X)
    w2 = _fixed_weights(names)
    g = np.sqrt(np.clip(w1, 1e-12, None) * np.clip(w2, 1e-12, None))
    s = g.sum()
    return g / s if s > 1e-12 else np.full(len(names), 1.0 / len(names))


def _topsis_scores(X, w):
    """TOPSIS 贴近度：X 已正向化 (n, m)（越大越好，理想星=1，负理想=0）。

    加权欧氏距离 D⁺(到理想星)、D⁻(到负理想星)，
    c_i = D⁻/(D⁺+D⁻) ∈ [0,1]，越大质量越好。
    """
    Dp = np.sqrt((w[None, :] * (1.0 - X) ** 2).sum(axis=1))
    Dm = np.sqrt((w[None, :] * X ** 2).sum(axis=1))
    den = Dp + Dm
    return np.where(den > 1e-12, Dm / np.maximum(den, 1e-12), 1.0)


def abs_calibrate(hist, q=0.1, min_n=30):
    """跨历元指标历史 → 绝对标定表 calib {metric: (q_lo, q_hi)}。

    hist : {metric: [逐历元累积的指标值]}（调用方维护，见 ppp.py absHist）
    q    : 标定分位（q = 差档，1−q = 优档）
    min_n: 样本不足的指标不进入标定表（热身保护）
    """
    calib = {}
    for k, v in hist.items():
        if len(v) >= min_n:
            calib[k] = (float(np.quantile(v, q)), float(np.quantile(v, 1.0 - q)))
    return calib


def qc_abs_weights(satInd, calib):
    """绝对阈值标定法：指标对历史分位数打绝对分（p 跨历元可比）。

    satInd : {prn: {指标名: 值}}
    calib  : abs_calibrate 输出 {metric: (q_lo, q_hi)}
    返回   : {prn: p}，p ∈ (0,1]，越大越好。越小越好的指标
             s = (q_hi−x)/(q_hi−q_lo)，越大越好（snr）s = (x−q_lo)/(q_hi−q_lo)，
             截断到 [0,1] 后按 FIXED_W 加权合成（缺失指标权重重分配）。
    某历元全部星质量都差时，相对评分（gra 等）会集体接近 1 而失真；
    绝对评分对全历元给出真实低分。
    """
    prns = list(satInd.keys())
    if not prns or not calib:
        return {p: 1.0 for p in prns}
    out = {}
    for prn, ind in satInd.items():
        num = den = 0.0
        for k, x in ind.items():
            if k not in calib or x is None or not np.isfinite(x):
                continue
            lo, hi = calib[k]
            rng = hi - lo
            if rng < 1e-12:
                s = 1.0
            else:
                s = (hi - float(x)) / rng if k in _INVERT \
                    else (float(x) - lo) / rng
                s = min(max(s, 0.0), 1.0)
            wj = FIXED_W.get(k, 0.10)
            num += wj * s
            den += wj
        out[prn] = float(max(num / den, 1e-3)) if den > 1e-9 else 1.0
    return out


def _winsorize(M, q):
    """列向 winsorize：clip 到 [q, 1−q] 分位数。

    单星指标异常值（如 MP 突跳）不再压扁跨星 min-max 归一化域；
    q<=0 或 n 过小时不起作用。NaN 由调用方先填充。
    """
    if q <= 0.0 or M.shape[0] < 3:
        return M
    lo = np.nanquantile(M, q, axis=0)
    hi = np.nanquantile(M, 1.0 - q, axis=0)
    return np.clip(M, lo[None, :], hi[None, :])


def qc_grey_weights(satInd, rho=0.5, winsor=0.1, method='gra'):
    """计算每星的质量评分权（相对评分类方法；abs 走 qc_abs_weights）。

    satInd : {prn: {指标名: 值}}，指标缺失可直接缺键
    method : 'gra' 熵权+灰色关联 | 'critic' CRITIC权+灰色关联 |
             'rank' 秩变换合成 | 'topsis' TOPSIS 贴近度 |
             'fixed' AHP固定权+灰色关联 | 'game' 博弈论组合权+灰色关联
    返回   : {prn: p}，p ∈ (0,1]，越大表示该星观测质量相对越好。
    """
    prns = list(satInd.keys())
    if not prns:
        return {}
    names = sorted({k for d in satInd.values() for k in d})
    if not names:
        return {p: 1.0 for p in prns}
    M = np.array([[float(d[n]) if n in d and d[n] is not None
                   and np.isfinite(d[n]) else np.nan
                   for n in names] for d in (satInd[p] for p in prns)])
    # 缺失值以列中位数中性填充；整列缺失填 0（归一化后退化 → 全 1）
    colMed = np.nanmedian(M, axis=0)
    colMed = np.where(np.isfinite(colMed), colMed, 0.0)
    M = np.where(np.isfinite(M), M, colMed[None, :])
    # 鲁棒化：列向 winsorize（抗单星异常值）后再统一方向
    M = _winsorize(M, winsor)
    # 方向统一：越小越好的指标反向（max − x），其余（x − min）
    cols = [(M[:, j].max() - M[:, j]) if n in _INVERT else (M[:, j] - M[:, j].min())
            for j, n in enumerate(names)]
    X = _minmax_positive(cols)
    if method == 'fixed':
        w = _fixed_weights(names)
    elif method == 'game':
        w = _game_weights(X, names)
    else:
        w = _critic_weights(X) if method == 'critic' else _entropy_weights(X)
    if method == 'rank':
        p = _rank_scores(X, w)
    elif method == 'topsis':
        p = _topsis_scores(X, w)
    else:
        xi = _grey_coeffs(X, rho=rho)
        p = xi @ w
    return {prn: float(max(v, 1e-3)) for prn, v in zip(prns, p)}
