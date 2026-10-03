import numpy as np
from scipy.optimize import linprog


def entropy_weight(data):
    """
    熵权法计算客观权重
    :param data: 2D array, 行=方案, 列=指标 (需为正向化后的数据)
    :return: 权重数组
    """
    data = np.array(data, dtype = float)
    # 1. 归一化 (确保非负, 避免log(0))
    min_vals = np.min(data, axis = 0)
    max_vals = np.max(data, axis = 0)
    norm_data = (data - min_vals) / (max_vals - min_vals + 1e-10)

    # 2. 计算比重矩阵
    p_ij = norm_data / (np.sum(norm_data, axis = 0) + 1e-10)

    # 3. 计算熵值
    m = data.shape[0]
    k = 1 / np.log(m + 1e-10)
    e_j = -k * np.sum(p_ij * np.log(p_ij + 1e-10), axis = 0)

    # 4. 计算差异系数并归一化得到权重
    d_j = 1 - e_j
    weights = d_j / np.sum(d_j)
    return weights


def ahp_weight(pairwise_matrix):
    """
    AHP层次分析法
    :param pairwise_matrix: n x n 成对比较矩阵 (a_ij 表示指标i比指标j的重要程度)
    :return: 权重向量 (归一化特征向量)
    """
    matrix = np.array(pairwise_matrix, dtype = float)
    n = matrix.shape[0]

    # 1. 几何平均法求权重
    geom_mean = np.power(np.prod(matrix, axis = 1), 1 / n)
    weights = geom_mean / np.sum(geom_mean)

    # 2. 一致性检验 (CR < 0.1 通过)
    # 最大特征值
    eig_vals, _ = np.linalg.eig(matrix)
    lambda_max = np.max(eig_vals.real)
    CI = (lambda_max - n) / (n - 1)
    # 随机一致性指标 RI (n<=10)
    RI_dict = {1: 0, 2: 0, 3: 0.58, 4: 0.9, 5: 1.12, 6: 1.24, 7: 1.32, 8: 1.41, 9: 1.45, 10: 1.49}
    RI = RI_dict.get(n, 1.5)
    CR = CI / RI if RI != 0 else 0.0

    if CR >= 0.1:
        print(f"⚠️ 一致性比率 CR={CR:.3f} >= 0.1，请检查判断矩阵！")
    else:
        print(f"✅ AHP一致性检验通过 (CR={CR:.3f})")
    return weights


def bwm_weight(criteria_names, best_to_others, others_to_worst):
    """
    BWM 最优最劣法 (基于线性规划求解)
    :param criteria_names: 指标名称列表
    :param best_to_others: 最优指标相对于其他指标的偏好值 (如 [1, 3, 5] 表示最优比第j个指标重要3倍)
    :param others_to_worst: 其他指标相对于最劣指标的偏好值
    :return: 权重数组
    """
    n = len(criteria_names)
    # 构造线性规划: 最小化 ξ
    # 变量顺序: [w1, w2, ..., wn, ξ]
    # 目标函数: min ξ
    c = np.zeros(n + 1)
    c[-1] = 1

    # 约束: 对于最优指标 b 和每个指标 j
    # |w_b - a_bj * w_j| <= ξ  -> 线性化为两个不等式
    # 找到最优指标索引 b，最劣指标索引 w
    # 若最优指标与最劣指标未指定，通常取第一个为最优，最后一个为最劣
    # 假设 best_to_others 指向索引0， others_to_worst 指向索引n-1
    b_idx = 0  # 通常最优指标在第一位
    w_idx = n - 1  # 通常最劣指标在最后一位

    A_ub = []
    b_ub = []

    # 约束组1: w_b - a_bj * w_j <= ξ
    for j in range(n):
        if j == b_idx: continue
        row = np.zeros(n + 1)
        row[b_idx] = 1
        row[j] = -best_to_others[j]  # a_bj
        row[-1] = -1
        A_ub.append(row)
        b_ub.append(0)
        # 约束组2: -w_b + a_bj * w_j <= ξ
        row = np.zeros(n + 1)
        row[b_idx] = -1
        row[j] = best_to_others[j]
        row[-1] = -1
        A_ub.append(row)
        b_ub.append(0)

    # 约束组3: w_j - a_jw * w_w <= ξ
    for j in range(n):
        if j == w_idx: continue
        row = np.zeros(n + 1)
        row[j] = 1
        row[w_idx] = -others_to_worst[j]  # a_jw
        row[-1] = -1
        A_ub.append(row)
        b_ub.append(0)
        # 约束组4: -w_j + a_jw * w_w <= ξ
        row = np.zeros(n + 1)
        row[j] = -1
        row[w_idx] = others_to_worst[j]
        row[-1] = -1
        A_ub.append(row)
        b_ub.append(0)

    # 约束组5: sum(w) = 1
    A_eq = [np.ones(n + 1)]
    A_eq[-1][-1] = 0  # ξ 不参与求和
    b_eq = [1]

    # 约束组6: w >= 0, ξ >= 0
    bounds = [(0, None) for _ in range(n)] + [(0, None)]

    result = linprog(c, A_ub = A_ub, b_ub = b_ub, A_eq = A_eq, b_eq = b_eq, bounds = bounds, method = 'highs')
    if result.success:
        weights = result.x[:n]
        return weights / np.sum(weights)
    else:
        print("⚠️ BWM 未收敛，请检查输入向量")
        return np.ones(n) / n


def waspas(data, weights, lamb = 0.5):
    """
    WASPAS 方法 (集成了加权和 WSM 与加权积 WPM)
    :param data: 2D array (正向化后)
    :param weights: 权重
    :param lamb: 系数 (0~1), 一般取 0.5
    :return: 各方案综合得分
    """
    data = np.array(data, dtype = float)
    weights = np.array(weights)

    # 1. 线性归一化 (极大型)
    min_vals = np.min(data, axis = 0)
    max_vals = np.max(data, axis = 0)
    norm_data = (data - min_vals) / (max_vals - min_vals + 1e-10)

    # 2. WSM (加权和)
    wsm = np.sum(norm_data * weights, axis = 1)

    # 3. WPM (加权积)
    wpm = np.prod(norm_data ** weights, axis = 1)

    # 4. 综合得分
    Q = lamb * wsm + (1 - lamb) * wpm
    return Q


def topsis(data, weights=None, impacts=None):
    """
    TOPSIS 综合评价方法
    :param data: 2D numpy array, 行=方案, 列=指标
    :param weights: 1D array, 各指标权重（默认等权重）
    :param impacts: list of str, 各指标类型 ('max' 或 'min')，默认全为 'max'
    :return: 每个方案的贴近度 C 值 (排序依据)
    """
    data = np.array(data, dtype=float)
    n_obj, n_criteria = data.shape

    # 默认权重等权
    if weights is None:
        weights = np.ones(n_criteria) / n_criteria
    else:
        weights = np.array(weights) / np.sum(weights)  # 归一化

    # 默认所有指标为正向（越大越好）
    if impacts is None:
        impacts = ['max'] * n_criteria

    # ---- 1. 正向化处理 (极小型 → 极大型) ----
    for i, imp in enumerate(impacts):
        if imp == 'min':
            data[:, i] = np.max(data[:, i]) - data[:, i]   # 或 1/x，但此处用极差法
        # 若还有其他类型（如中间型、区间型）可自行扩展

    # ---- 2. 标准化 (向量归一化) ----
    norm_data = data / np.sqrt(np.sum(data**2, axis=0))

    # ---- 3. 加权 ----
    weighted_norm = norm_data * weights

    # ---- 4. 确定正理想解和负理想解 ----
    ideal_best = np.max(weighted_norm, axis=0)
    ideal_worst = np.min(weighted_norm, axis=0)

    # ---- 5. 计算到正负理想解的欧氏距离 ----
    dist_best = np.sqrt(np.sum((weighted_norm - ideal_best)**2, axis=1))
    dist_worst = np.sqrt(np.sum((weighted_norm - ideal_worst)**2, axis=1))

    # ---- 6. 计算相对贴近度 C ----
    C = dist_worst / (dist_best + dist_worst + 1e-10)  # 避免分母为0
    return C

def grey_relational_analysis(data, weights=None, reference=None, rho=0.5):
    """
    灰色关联度分析 (GRA)
    :param data: 2D array, 行=方案, 列=指标
    :param weights: 各指标权重 (默认等权)
    :param reference: 参考序列 (默认取各列最大值，即理想方案)
    :param rho: 分辨系数，通常取 0.5
    :return: 各方案的关联度 (排序依据)
    """
    data = np.array(data, dtype=float)
    n_obj, n_criteria = data.shape

    # 默认等权重
    if weights is None:
        weights = np.ones(n_criteria) / n_criteria
    else:
        weights = np.array(weights) / np.sum(weights)

    # ---- 1. 无量纲化 (均值化处理) ----
    # 也可以使用初值化 (除以第一行) 或极差标准化，这里选择均值化
    mean_cols = np.mean(data, axis=0)
    norm_data = data / mean_cols

    # ---- 2. 确定参考序列 (若未提供则取每列最优，因已均值化，最优为最大值) ----
    if reference is None:
        reference = np.max(norm_data, axis=0)   # 所有指标均为正向
    else:
        reference = np.array(reference) / mean_cols   # 同样对参考序列进行均值化

    # ---- 3. 计算关联系数 ----
    diff = np.abs(norm_data - reference)
    min_diff = np.min(diff)
    max_diff = np.max(diff)
    # 关联系数矩阵
    coef = (min_diff + rho * max_diff) / (diff + rho * max_diff + 1e-10)

    # ---- 4. 计算加权关联度 ----
    degree = np.sum(coef * weights, axis=1)
    return degree


if __name__ == "__main__":
    # ---- 原始正向化数据 ----
    data = np.array([
        [0.8, 0.9, 0.6],
        [0.7, 0.8, 0.9],
        [0.9, 0.5, 0.7],
        [0.6, 0.7, 0.8]
    ])
    criteria = ["指标A", "指标B", "指标C"]

    print("========== 1. 客观赋权：熵权法 ==========")
    w_entropy = entropy_weight(data)
    print(f"熵权法权重: {np.round(w_entropy, 3)}")

    print("\n--- 用熵权法权重评价 ---")
    score_topsis = topsis(data, weights = w_entropy)  # 沿用之前定义的函数
    score_gra = grey_relational_analysis(data, weights = w_entropy)
    score_waspas = waspas(data, weights = w_entropy)
    print(f"TOPSIS贴近度: {np.round(score_topsis, 3)} -> 排序: {np.argsort(-score_topsis) + 1}")
    print(f"灰色关联度: {np.round(score_gra, 3)} -> 排序: {np.argsort(-score_gra) + 1}")
    print(f"WASPAS得分: {np.round(score_waspas, 3)} -> 排序: {np.argsort(-score_waspas) + 1}")

    print("\n========== 2. 主观赋权：AHP ==========")
    # 构造 3x3 成对比较矩阵 (指标A最重要，指标B其次，指标C最次)
    pairwise = np.array([
        [1, 3, 5],
        [1 / 3, 1, 2],
        [1 / 5, 1 / 2, 1]
    ])
    w_ahp = ahp_weight(pairwise)
    print(f"AHP权重: {np.round(w_ahp, 3)}")

    # 与GRA结合 (示例)
    score_gra_ahp = grey_relational_analysis(data, weights = w_ahp)
    print(f"AHP+GRA关联度: {np.round(score_gra_ahp, 3)} -> 排序: {np.argsort(-score_gra_ahp) + 1}")

    print("\n========== 3. 主观赋权：BWM ==========")
    # 假设指标A是最优，指标C是最劣
    # best_to_others: [1, 3, 5] 表示 A比A重要1倍，A比B重要3倍，A比C重要5倍
    # others_to_worst: [5, 3, 1] 表示 A比C重要5倍，B比C重要3倍，C比C重要1倍
    best_to_others = [1, 3, 5]
    others_to_worst = [5, 3, 1]
    w_bwm = bwm_weight(criteria, best_to_others, others_to_worst)
    print(f"BWM权重: {np.round(w_bwm, 3)}")

    # 与WASPAS结合 (示例)
    score_waspas_bwm = waspas(data, weights = w_bwm)
    print(f"BWM+WASPAS得分: {np.round(score_waspas_bwm, 3)} -> 排序: {np.argsort(-score_waspas_bwm) + 1}")

    print("\n========== 4. 混合策略示例 ==========")
    # 例如: 熵权法权重 + TOPSIS，AHP权重 + WASPAS
    score_topsis_entropy = topsis(data, weights = w_entropy)
    score_waspas_ahp = waspas(data, weights = w_ahp)
    print(f"熵权TOPSIS排序: {np.argsort(-score_topsis_entropy) + 1}")
    print(f"AHP-WASPAS排序: {np.argsort(-score_waspas_ahp) + 1}")