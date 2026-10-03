#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
GNSS 观测数据综合质量评估模块

供 prn 模式（按 PRN 评估）与 station 模式（按测站评估）共同复用。

对外接口:
    evaluate_matrix(labels, modules, data, ...)   # 直接对矩阵评估
    evaluate_qcvalue(csv_file, ...)               # 读取 _QCValue.csv 后评估
    DEFAULT_IMPACTS                               # 各 QC 指标方向
"""
import os
import csv
import numpy as np

from madm.GRA import (
    entropy_weight, ahp_weight, bwm_weight,
    grey_relational_analysis, topsis, waspas,
)
from madm.INCOME import INCOME
from madm.PCA import pca_reduce


# ========================= 评估常量 =========================
# 各 QC 指标方向：'max' 越大越好，'min' 越小越好
DEFAULT_IMPACTS = {
    'CSR':        'min',   # 周跳比
    'CNR':        'max',   # 载噪比
    'multipath':  'min',   # 多路径
    'phaseNoise': 'min',   # 相位噪声
    'CMC':        'min',   # 码减相位
    'IOD':        'min',   # 电离层变化率
    'highOrder':  'min',   # 高阶差分
    'freq':       'max',   # 频段利用率
    'spp3drms':   'min',   # SPP ENU 3D RMS
    'dop':        'min',   # SPP DOP 均值
}


# ========================= 辅助函数 =========================
def ensure_dir(path):
    if path and not os.path.exists(path):
        os.makedirs(path)


def _read_qcvalue_csv(csv_file):
    """读取 _QCValue.csv，返回 (labels, modules, data)。"""
    with open(csv_file, 'r', newline='') as f:
        reader = csv.reader(f)
        header = next(reader)
        modules = header[1:]
        labels, rows = [], []
        for row in reader:
            if not row or not row[0].strip():
                continue
            labels.append(row[0].strip())
            vals = []
            for v in row[1:]:
                v = v.strip()
                if v == '' or v.lower() in ('nan', 'none'):
                    vals.append(np.nan)
                elif v.lower() in ('inf', '+inf'):
                    vals.append(np.inf)
                elif v.lower() == '-inf':
                    vals.append(-np.inf)
                else:
                    try:
                        vals.append(float(v))
                    except ValueError:
                        vals.append(np.nan)
            rows.append(vals)
    return labels, modules, np.array(rows, dtype=float)


def _get_impacts(modules):
    """根据列名匹配方向，未匹配到的默认 max。"""
    impacts = []
    for m in modules:
        m_up = m.upper()
        imp = 'max'
        for k, v in DEFAULT_IMPACTS.items():
            if k.upper() in m_up:
                imp = v
                break
        impacts.append(imp)
    return impacts


def _preprocess_eval(data):
    """空值填充列最小值 + inf 替换 + 负值平移。"""
    data = data.copy()
    for j in range(data.shape[1]):
        col = data[:, j]
        finite_mask = np.isfinite(col)
        if not finite_mask.any():
            col[:] = 0.0
            data[:, j] = col
            continue
        finite_vals = col[finite_mask]
        col_min = np.min(finite_vals)
        col_max = np.max(finite_vals)

        if np.any(np.isposinf(col)):
            col[np.isposinf(col)] = col_max
        if np.any(np.isneginf(col)):
            col[np.isneginf(col)] = col_min

        col[np.isnan(col)] = col_min
        data[:, j] = col

    min_vals = data.min(axis=0)
    shift = np.where(min_vals < 0, -min_vals + 1e-6, 0.0)
    return data + shift


def _compute_weights(data_proc, modules, weight_method):
    n_criteria = data_proc.shape[1]
    if weight_method == 'entropy':
        return entropy_weight(data_proc)
    elif weight_method == 'ahp':
        n = n_criteria
        pairwise = np.ones((n, n))
        for i in range(n):
            for j in range(n):
                pairwise[i, j] = (n - i) / (n - j)
        return ahp_weight(pairwise)
    elif weight_method == 'bwm':
        n = n_criteria
        best_to_others = [1] + [3] * (n - 1)
        others_to_worst = [3] * (n - 1) + [1]
        return bwm_weight(modules, best_to_others, others_to_worst)
    else:
        return np.ones(n_criteria) / n_criteria


# ========================= 核心评估 =========================
def evaluate_matrix(labels, modules, data,
                    weight_method='entropy', eval_method='gra',
                    rho=0.5, lamb=0.5,
                    use_pca=False, pca_var_threshold=0.95,
                    pca_n_components=None, pca_print=True,
                    pca_report_file=None):
    """
    对 (labels, modules, data) 矩阵做综合评价。
    返回 (labels, modules, data_used, weights, scores, ranking) 或 None。
    """
    if len(labels) == 0 or data.size == 0:
        return None

    # ---- 1. 极小指标正向化 ----
    impacts = _get_impacts(modules)
    data_proc = data.copy()
    for j, imp in enumerate(impacts):
        if imp == 'min':
            col = data_proc[:, j]
            finite = col[np.isfinite(col)]
            if finite.size == 0:
                continue
            mx = np.nanmax(finite)
            data_proc[:, j] = mx - col

    # ---- 2. 空值/inf/负值预处理 ----
    data_proc = _preprocess_eval(data_proc)

    # ---- 3. PCA 降维 ----
    if use_pca:
        original_modules = list(modules)
        Z, pca_info = pca_reduce(
            data_proc,
            impacts=None,                    # 已统一正向化，无需再次处理
            n_components=pca_n_components,
            var_threshold=pca_var_threshold,
            standardize=True,
            sign_align=True,
            score_mode='none',
        )
        k = pca_info['n_components']
        comps = pca_info['components']            # (n_features, k)
        evr = np.asarray(pca_info['explained_variance_ratio']).ravel()
        cum = np.asarray(pca_info['cumulative_ratio']).ravel()

        # ---- 终端打印 ----
        if pca_print:
            print(f"    [PCA] 指标数 {len(original_modules)} → 主成分数 {k}，"
                  f"累计贡献率 {cum[-1] * 100:.2f}%")
            print(f"    [PCA] 各主成分贡献率: {np.round(evr, 4)}")
            print(f"    [PCA] 载荷矩阵（行=原指标，列=主成分）:")
            header = "        {:<16s}".format("指标") + \
                     "".join(f"{('PC' + str(i + 1)):>10s}" for i in range(k))
            print(header)
            for name, row in zip(original_modules, comps):
                print("        {:<16s}".format(name) +
                      "".join(f"{v:>10.4f}" for v in row))

        # ---- CSV 报告 ----
        if pca_report_file:
            try:
                ensure_dir(os.path.dirname(pca_report_file) or '.')
                with open(pca_report_file, 'w', newline='',
                          encoding='utf-8-sig') as f:
                    writer = csv.writer(f)

                    # 段 1：总体信息
                    writer.writerow(['# PCA Summary'])
                    writer.writerow(['指标数', len(original_modules)])
                    writer.writerow(['主成分数', k])
                    writer.writerow(['累计贡献率', f"{cum[-1]* 100:.2f}%"])
                    writer.writerow(['原始指标', '|'.join(original_modules)])
                    writer.writerow([])

                    # 段 2：各主成分贡献率
                    writer.writerow(['# Explained Variance Ratio'])
                    writer.writerow(['主成分', '贡献率', '累计贡献率'])
                    for i in range(k):
                        writer.writerow([f'PC{i + 1}',
                                         f"{evr[i] * 100:.2f}%",
                                         f"{cum[i] * 100:.2f}%"])
                    writer.writerow([])

                    # 段 3：载荷矩阵
                    writer.writerow(['# 载荷矩阵'])
                    writer.writerow(['指标'] +
                                    [f'PC{i + 1}' for i in range(k)])
                    for name, row in zip(original_modules, comps):
                        writer.writerow([name] +
                                        [f"{v:.4f}" for v in row])

                print(f"    [PCA]报告: {pca_report_file}")
            except Exception as e:
                print(f"    [PCA]报告写入失败: {e}")

        data_proc = Z
        modules = [f'PC{i + 1}' for i in range(k)]

    # ---- 4. 权重 ----
    weights = _compute_weights(data_proc, modules, weight_method)
    n_criteria = data_proc.shape[1]

    # ---- 5. 评价 ----
    if eval_method == 'gra':
        scores = grey_relational_analysis(data_proc, weights=weights, rho=rho)
    elif eval_method == 'topsis':
        scores = topsis(data_proc, weights=weights,
                        impacts=['max'] * n_criteria)
    elif eval_method == 'waspas':
        scores = waspas(data_proc, weights=weights, lamb=lamb)
    elif eval_method == 'income':
        income = INCOME(n_levels=3, k=3,
                        distance_weighted=True, binning='equal_width')
        income.fit_by_madm(data_proc, madm_method='topsis',
                           weights=weights, impacts=['max'] * n_criteria)
        scores = income.predict(data_proc)
    else:
        raise ValueError(f"未知评价方法: {eval_method}")

    ranking = np.argsort(np.argsort(-scores)) + 1
    return labels, modules, data_proc, weights, scores, ranking


def evaluate_qcvalue(csv_file, weight_method='entropy', eval_method='gra',
                     rho=0.5, lamb=0.5,
                     use_pca=False, pca_var_threshold=0.95,
                     pca_n_components=None,
                     pca_report_file=None):
    """读取 _QCValue.csv 并做综合评价。"""
    labels, modules, data = _read_qcvalue_csv(csv_file)
    if len(labels) == 0:
        print(f"文件 {csv_file} 无有效数据")
        return None
    return evaluate_matrix(
        labels, modules, data,
        weight_method=weight_method,
        eval_method=eval_method,
        rho=rho, lamb=lamb,
        use_pca=use_pca,
        pca_var_threshold=pca_var_threshold,
        pca_n_components=pca_n_components,
        pca_report_file=pca_report_file,
    )