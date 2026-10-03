#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
madm/PCA.py
主成分分析：输入评价指标矩阵，输出主成分得分 / 综合得分。

用法：
    from madm.PCA import pca_reduce
    Z, info = pca_reduce(data, impacts=['max','min',...])          # 降维
    score, info = pca_reduce(data, score_mode='weighted')          # 综合打分
"""
import numpy as np


class PCAReducer:
    """PCA 降维器（z-score 标准化 + 特征分解 + 符号校正）"""

    def __init__(self, n_components=None, var_threshold=0.85,
                 standardize=True, sign_align=True):
        """
        :param n_components: 指定主成分个数；None 则按 var_threshold 自动确定
        :param var_threshold: 累计方差贡献率阈值（默认 0.85）
        :param standardize: 是否 z-score 标准化（消除量纲影响，默认 True）
        :param sign_align: 符号校正，使每个主成分的主要载荷为正
        """
        self.n_components = n_components
        self.var_threshold = var_threshold
        self.standardize = standardize
        self.sign_align = sign_align

        # 拟合后填充
        self.mean_ = None
        self.std_ = None
        self.components_ = None                # (n_features, k) 载荷
        self.explained_variance_ = None        # (k,) 各主成分方差
        self.explained_variance_ratio_ = None  # (k,) 各主成分贡献率
        self.cumulative_ratio_ = None          # (k,) 累计贡献率
        self.n_components_ = None              # 实际保留个数
        self.signs_ = None                     # (k,) 各主成分的符号

    # ---------- 拟合 ----------
    def fit(self, X):
        X = np.asarray(X, dtype=float)
        n_samples, n_features = X.shape
        if n_samples < 2:
            raise ValueError("样本数至少为 2")

        # 1. 中心化 + 可选标准化
        self.mean_ = X.mean(axis=0)
        if self.standardize:
            self.std_ = X.std(axis=0, ddof=1)
            self.std_[self.std_ < 1e-12] = 1.0
            Xs = (X - self.mean_) / self.std_
        else:
            self.std_ = np.ones(n_features)
            Xs = X - self.mean_

        # 2. 协方差矩阵特征分解（eigh 保证实对称矩阵稳定）
        cov = np.cov(Xs, rowvar=False)
        eigvals, eigvecs = np.linalg.eigh(cov)

        # 3. 降序排列
        idx = np.argsort(eigvals)[::-1]
        eigvals = np.clip(eigvals[idx], 0, None)
        eigvecs = eigvecs[:, idx]

        # 4. 符号校正（让每个主成分的绝对值最大的载荷为正）
        if self.sign_align:
            signs = np.ones(eigvecs.shape[1])
            for j in range(eigvecs.shape[1]):
                col = eigvecs[:, j]
                k = int(np.argmax(np.abs(col)))
                if col[k] < 0:
                    signs[j] = -1.0
            eigvecs = eigvecs * signs
        else:
            signs = np.ones(eigvecs.shape[1])
        self.signs_full_ = signs

        # 5. 排除接近零的特征值（n_samples < n_features 时必出现）
        thresh = 1e-10 * (eigvals[0] + 1e-12)
        nonzero_mask = eigvals > thresh
        if not np.any(nonzero_mask):
            nonzero_mask[0] = True
        k_max = int(np.sum(nonzero_mask))

        ratio = eigvals / (eigvals.sum() + 1e-12)
        cum = np.cumsum(ratio)

        # 6. 确定主成分个数
        if self.n_components is not None:
            k = min(int(self.n_components), k_max)
        else:
            k = int(np.searchsorted(cum, self.var_threshold) + 1)
            k = max(1, min(k, k_max))

        self.components_ = eigvecs[:, :k]
        self.explained_variance_ = eigvals[:k]
        self.explained_variance_ratio_ = ratio[:k]
        self.cumulative_ratio_ = cum[:k]
        self.n_components_ = k
        self.signs_ = signs[:k]
        return self

    # ---------- 变换 ----------
    def transform(self, X):
        X = np.asarray(X, dtype=float)
        Xs = (X - self.mean_) / self.std_
        return Xs @ self.components_          # (n_samples, k)

    def fit_transform(self, X):
        return self.fit(X).transform(X)

    # ---------- 简要信息 ----------
    def summary(self):
        lines = [f"主成分个数: {self.n_components_}"]
        for i in range(self.n_components_):
            lines.append(
                f"  PC{i+1}: 方差={self.explained_variance_[i]:.4f}  "
                f"贡献率={self.explained_variance_ratio_[i]*100:.2f}%  "
                f"累计={self.cumulative_ratio_[i]*100:.2f}%"
            )
        return "\n".join(lines)


# ============ 便捷函数 ============
def pca_reduce(data, impacts=None, n_components=None, var_threshold=0.95,
               standardize=True, sign_align=True, score_mode='none'):
    """
    对评价指标做主成分分析。

    :param data: 2D array (n_samples, n_features)，原始指标
    :param impacts: list of 'max'/'min'，指定每个指标方向；
                    'min' 指标会先正向化（max - x）
    :param n_components: 指定主成分个数，None 则按 var_threshold 自动
    :param var_threshold: 累计方差贡献率阈值（默认 0.85）
    :param standardize: 是否标准化
    :param sign_align: 符号校正
    :param score_mode:
        'none'     → 返回主成分得分矩阵 (n_samples, k)
        'weighted' → 返回以贡献率为权的主成分加权得分 (n_samples,)  ← 综合评分
        'first'    → 返回第一主成分得分 (n_samples,)

    :return: (result, info)
        info 包含：components, explained_variance, explained_variance_ratio,
                   cumulative_ratio, n_components, mean, std, signs
    """
    X = np.asarray(data, dtype=float).copy()

    # 1. 极小指标正向化
    if impacts is not None:
        if len(impacts) != X.shape[1]:
            raise ValueError("impacts 长度与指标个数不一致")
        for j, imp in enumerate(impacts):
            if imp == 'min':
                col = X[:, j]
                mx = np.nanmax(col)
                X[:, j] = mx - col

    # 2. 缺失值用列最小值填充
    for j in range(X.shape[1]):
        col = X[:, j]
        if np.any(np.isnan(col)):
            col[np.isnan(col)] = np.nanmin(col)
        X[:, j] = col

    # 3. 拟合 + 变换
    reducer = PCAReducer(n_components=n_components,
                         var_threshold=var_threshold,
                         standardize=standardize,
                         sign_align=sign_align)
    Z = reducer.fit_transform(X)

    info = {
        'components': reducer.components_,
        'explained_variance': reducer.explained_variance_,
        'explained_variance_ratio': reducer.explained_variance_ratio_,
        'cumulative_ratio': reducer.cumulative_ratio_,
        'n_components': reducer.n_components_,
        'mean': reducer.mean_,
        'std': reducer.std_,
        'signs': reducer.signs_,
    }

    # 4. 输出
    if score_mode == 'none':
        return Z, info
    elif score_mode == 'weighted':
        w = reducer.explained_variance_ratio_ / reducer.explained_variance_ratio_.sum()
        return Z @ w, info
    elif score_mode == 'first':
        return Z[:, 0], info
    else:
        raise ValueError("score_mode 仅支持 'none' / 'weighted' / 'first'")


# ============ 自测 ============
if __name__ == '__main__':
    np.set_printoptions(precision=4, suppress=True)

    # 模拟 6 个测站、8 个指标（8 个指标间有一定相关性）
    rng = np.random.default_rng(42)
    base = rng.normal(size=(6, 3))
    X = np.hstack([
        base[:, [0]] * 1.0 + 0.2 * rng.normal(size=(6, 1)),   # 指标1
        base[:, [0]] * 0.8 + 0.3 * rng.normal(size=(6, 1)),   # 指标2（与1相关）
        base[:, [1]] * 1.0 + 0.2 * rng.normal(size=(6, 1)),   # 指标3
        base[:, [1]] * 0.7 + 0.3 * rng.normal(size=(6, 1)),   # 指标4（与3相关）
        base[:, [2]] * 1.0 + 0.2 * rng.normal(size=(6, 1)),   # 指标5
        -base[:, [2]] + 0.3 * rng.normal(size=(6, 1)),        # 指标6（与5负相关）
        base[:, [0]] * 0.5 + base[:, [1]] * 0.5,
        base[:, [2]] * 0.5 + 0.2 * rng.normal(size=(6, 1)),
    ])
    impacts = ['max'] * 8

    print("========== 主成分得分（保留前几个主成分） ==========")
    Z, info = pca_reduce(X, impacts=impacts, var_threshold=0.85)
    print("主成分得分矩阵:\n", Z)
    print("\n各主成分贡献率:", np.round(info['explained_variance_ratio'], 4))
    print("累计贡献率:", np.round(info['cumulative_ratio'], 4))
    print("载荷矩阵 (行=原指标, 列=主成分):\n", info['components'])

    print("\n========== PCA 综合得分（加权） ==========")
    score, info = pca_reduce(X, impacts=impacts, score_mode='weighted')
    print("综合得分:", np.round(score, 4))
    print("排名 (1=最优):", np.argsort(-score) + 1)

    print("\n========== 只保留第一主成分 ==========")
    score1, _ = pca_reduce(X, impacts=impacts, score_mode='first')
    print("第一主成分得分:", np.round(score1, 4))
    print("排名:", np.argsort(-score1) + 1)