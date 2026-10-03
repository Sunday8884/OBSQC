#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
k-means + KNN 融合的有监督质量分级模块（与 PCA.py 平级）

参考 Li et al. (2024, GPS Solutions) 的方法：
  - 无标签样本 → k-means 聚类；
  - 用每个簇的"物理精度指标"（本项目中为 SPP 3D RMS，代替论文的 PPP 精度）
    均值为簇打标签（簇均值越小质量越好，等级越高）；
  - 构造带标签训练集 → KNN 有监督分类（k 由交叉验证自动选取）。

关键约定：
  - label_module（默认 spp3drms）只作为监督标签，**不作为特征/指标**；
  - 标签缺失（NaN）的样本不参与聚类训练，最终等级由 KNN 预测补齐；
  - 等级范围 1 ~ n_clusters（1 最差，n_clusters 最好）。

依赖：
  - madm.EQCValue 的 _get_impacts / _preprocess_eval / ensure_dir；
  - madm.PCA 的 pca_reduce。

对外接口:
    evaluate_KNN(labels, modules, data, ...) -> dict | None
"""
import os
import csv
import numpy as np

from madm.EQCValue import _get_impacts, _preprocess_eval, ensure_dir
from madm.PCA import pca_reduce


def evaluate_KKNN(labels, modules, data,
                  label_module='spp3drms',
                  n_clusters=10, knn_k=None,
                  use_pca=False, pca_var_threshold=0.95,
                  pca_n_components=None,
                  pca_print=True, pca_report_file=None,
                  random_state=42):
    """
    k-means + KNN 融合的有监督质量分级。

    :param labels:        样本标签（station 名或日期）
    :param modules:       指标列名（须包含 label_module 列）
    :param data:          (n_samples, n_modules) 指标矩阵
    :param label_module:  监督标签列名（越小越好），不作为特征
    :param n_clusters:    k-means 聚类数（论文取 10）
    :param knn_k:         KNN 邻居数；None 时交叉验证自动选取
    :return: dict（见函数末尾）或 None（输入无效 / 缺少标签列 / sklearn 缺失）
    """
    try:
        from sklearn.cluster import KMeans
        from sklearn.neighbors import KNeighborsClassifier
        from sklearn.model_selection import cross_val_score
        from sklearn.preprocessing import StandardScaler
    except ImportError as e:
        print(f"[supervised] 缺少 scikit-learn：{e}")
        return None

    if len(labels) == 0 or data.size == 0:
        return None

    modules = list(modules)

    # ---- 1. 拆分标签列 / 特征列 ----
    label_idx = None
    for i, m in enumerate(modules):
        if m == label_module or m.lower() == label_module.lower():
            label_idx = i
            break
    if label_idx is None:
        print(f"[supervised] 未找到标签列 '{label_module}'，可用列：{modules}")
        return None
    if len(modules) < 2:
        print("[supervised] 特征数为 0，无法聚类")
        return None

    feature_modules = [m for i, m in enumerate(modules) if i != label_idx]
    X_raw   = np.delete(data, label_idx, axis=1).astype(float)
    y_label = data[:, label_idx].astype(float)

    # ---- 2. 特征正向化（min -> max，标签列不参与） ----
    impacts = _get_impacts(feature_modules)
    X = X_raw.copy()
    for j, imp in enumerate(impacts):
        if imp == 'min':
            col = X[:, j]
            finite = col[np.isfinite(col)]
            if finite.size == 0:
                continue
            mx = np.nanmax(finite)
            X[:, j] = mx - col

    # ---- 3. 空值 / inf / 负值预处理 ----
    X = _preprocess_eval(X)

    # ---- 4. 可选 PCA ----
    if use_pca:
        Z, pca_info = pca_reduce(
            X, impacts=None,
            n_components=pca_n_components,
            var_threshold=pca_var_threshold,
            standardize=True, sign_align=True, score_mode='none',
        )
        k_pc = pca_info['n_components']
        comps = pca_info['components']
        evr = np.asarray(pca_info['explained_variance_ratio']).ravel()
        cum = np.asarray(pca_info['cumulative_ratio']).ravel()

        if pca_print:
            print(f"    [PCA] 特征数 {len(feature_modules)} → 主成分数 {k_pc}，"
                  f"累计贡献率 {cum[-1] * 100:.2f}%")
            print(f"    [PCA] 各主成分贡献率: {np.round(evr, 4)}")
            print("    [PCA] 载荷矩阵（行=原特征，列=主成分）:")
            print("        {:<16s}".format("指标") +
                  "".join(f"{('PC' + str(i + 1)):>10s}" for i in range(k_pc)))
            for name, row in zip(feature_modules, comps):
                print("        {:<16s}".format(name) +
                      "".join(f"{v:>10.4f}" for v in row))

        if pca_report_file:
            try:
                ensure_dir(os.path.dirname(pca_report_file) or '.')
                with open(pca_report_file, 'w', newline='',
                          encoding='utf-8-sig') as f:
                    w = csv.writer(f)
                    w.writerow(['# PCA Summary (supervised)'])
                    w.writerow(['特征数', len(feature_modules)])
                    w.writerow(['主成分数', k_pc])
                    w.writerow(['累计贡献率', f"{cum[-1] * 100:.2f}%"])
                    w.writerow(['原始特征', '|'.join(feature_modules)])
                    w.writerow([])
                    w.writerow(['# Explained Variance Ratio'])
                    w.writerow(['主成分', '贡献率', '累计贡献率'])
                    for i in range(k_pc):
                        w.writerow([f'PC{i + 1}',
                                    f"{evr[i] * 100:.2f}%",
                                    f"{cum[i] * 100:.2f}%"])
                    w.writerow([])
                    w.writerow(['# 载荷矩阵'])
                    w.writerow(['指标'] + [f'PC{i + 1}' for i in range(k_pc)])
                    for name, row in zip(feature_modules, comps):
                        w.writerow([name] + [f"{v:.4f}" for v in row])
                print(f"    [PCA]报告: {pca_report_file}")
            except Exception as e:
                print(f"    [PCA]报告写入失败: {e}")

        X = Z
        feature_out = [f'PC{i + 1}' for i in range(k_pc)]
    else:
        feature_out = feature_modules

    # ---- 5. 有效标签样本过滤 ----
    valid_mask = np.isfinite(y_label)
    n_valid = int(valid_mask.sum())
    if n_valid < n_clusters:
        print(f"[supervised] 有效标签样本 {n_valid} < 聚类数 {n_clusters}，无法聚类")
        return None

    # ---- 6. 特征标准化（仅用有效样本拟合，避免缺失样本污染统计量） ----
    scaler = StandardScaler()
    scaler.fit(X[valid_mask])
    X_scaled_all   = scaler.transform(X)
    X_scaled_valid = X_scaled_all[valid_mask]
    y_valid        = y_label[valid_mask]

    # ---- 7. k-means 聚类 ----
    kmeans = KMeans(n_clusters=n_clusters, n_init=10,
                    random_state=random_state)
    cluster_ids_valid = kmeans.fit_predict(X_scaled_valid)

    # ---- 8. 簇内标签均值 → 等级（label 越小越好 → 均值越小等级越高） ----
    cluster_mean = np.full(n_clusters, np.nan, dtype=float)
    for c in range(n_clusters):
        m = (cluster_ids_valid == c)
        if m.sum() > 0:
            cluster_mean[c] = float(np.mean(y_valid[m]))

    sorted_clusters = sorted(
        range(n_clusters),
        key=lambda c: (-cluster_mean[c] if np.isfinite(cluster_mean[c])
                       else np.inf)   # NaN 簇排末尾，赋最低级
    )
    cluster_to_grade = {c: i + 1 for i, c in enumerate(sorted_clusters)}

    grades_valid = np.array([cluster_to_grade[c] for c in cluster_ids_valid],
                            dtype=int)

    # ---- 9. KNN 训练（k 交叉验证自动选择） ----
    auto_k = False
    if knn_k is None:
        auto_k = True
        best_k, best_acc = 5, -1.0
        cv_folds = min(10, n_valid)
        for kk in (3, 5, 7, 9, 11, 15):
            if kk >= n_valid or cv_folds < 2:
                continue
            try:
                knn = KNeighborsClassifier(n_neighbors=kk)
                sc  = cross_val_score(knn, X_scaled_valid, grades_valid,
                                      cv=cv_folds)
                acc = float(sc.mean())
                if acc > best_acc:
                    best_acc, best_k = acc, kk
            except Exception:
                continue
        if best_acc < 0:
            best_k, best_acc = min(5, n_valid - 1), float('nan')
        knn_k = best_k
    else:
        best_acc = float('nan')

    knn = KNeighborsClassifier(n_neighbors=knn_k)
    knn.fit(X_scaled_valid, grades_valid)

    # ---- 10. 全体样本等级（k-means 等级 + KNN 补齐缺失标签样本） ----
    predicted_all = knn.predict(X_scaled_all)
    final_grades = np.full(len(labels), -1, dtype=int)
    final_grades[valid_mask] = grades_valid
    final_grades[~valid_mask] = predicted_all[~valid_mask]

    cluster_ids_full = np.full(len(labels), -1, dtype=int)
    cluster_ids_full[valid_mask] = cluster_ids_valid

    print(f"    [supervised] k-means 聚类数 k={n_clusters}, "
          f"有效标签样本 {n_valid}/{len(labels)}")
    print(f"    [supervised] KNN n_neighbors={knn_k}"
          f"{' (auto)' if auto_k else ''}, "
          f"CV 平均准确率="
          f"{'N/A' if not np.isfinite(best_acc) else f'{best_acc * 100:.2f}%'}")
    print(f"    [supervised] 各簇平均 {label_module}: " +
          ", ".join(f"C{c}={cluster_mean[c]:.3f}"
                    for c in sorted_clusters if np.isfinite(cluster_mean[c])))

    return {
        'labels':             list(labels),
        'feature_modules':    feature_out,
        'data_used':          X,
        'label_module':       label_module,
        'label_values':       y_label,
        'cluster_ids':        cluster_ids_full,
        'valid_mask':         valid_mask,
        'grades':             final_grades,
        'knn_k':              knn_k,
        'knn_cv_acc':         best_acc,
        'n_clusters':         n_clusters,
        'cluster_mean_label': cluster_mean,
        'cluster_to_grade':   cluster_to_grade,
    }
