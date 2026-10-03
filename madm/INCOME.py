import numpy as np
from itertools import product
from sklearn.preprocessing import MinMaxScaler
from collections import defaultdict


def topsis(data, weights = None, impacts = None):
    """
    TOPSIS 综合评价方法
    :param data: 2D numpy array, 行=方案, 列=指标
    :param weights: 1D array, 各指标权重（默认等权重）
    :param impacts: list of str, 各指标类型 ('max' 或 'min')，默认全为 'max'
    :return: 每个方案的贴近度 C 值 (排序依据)
    """
    data = np.array(data, dtype = float)
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
            data[:, i] = np.max(data[:, i]) - data[:, i]  # 或 1/x，但此处用极差法
        # 若还有其他类型（如中间型、区间型）可自行扩展

    # ---- 2. 标准化 (向量归一化) ----
    norm_data = data / np.sqrt(np.sum(data ** 2, axis = 0))

    # ---- 3. 加权 ----
    weighted_norm = norm_data * weights

    # ---- 4. 确定正理想解和负理想解 ----
    ideal_best = np.max(weighted_norm, axis = 0)
    ideal_worst = np.min(weighted_norm, axis = 0)

    # ---- 5. 计算到正负理想解的欧氏距离 ----
    dist_best = np.sqrt(np.sum((weighted_norm - ideal_best) ** 2, axis = 1))
    dist_worst = np.sqrt(np.sum((weighted_norm - ideal_worst) ** 2, axis = 1))

    # ---- 6. 计算相对贴近度 C ----
    C = dist_worst / (dist_best + dist_worst + 1e-10)  # 避免分母为0
    return C
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

class INCOME:
    """
    INCOME (INtelligent Characteristic Objects METhod)基于 kNN 与 COMET 的数据驱动多准则决策
    """
    def __init__(self, n_levels=3, k=3, distance_weighted=True, binning='equal_width'):
        """
        :param n_levels: 每个准则离散化的等级数（默认 3）
        :param k: kNN 中的邻居个数
        :param distance_weighted: 是否按距离倒数加权（True）或简单平均（False）
        :param binning: 离散化方式 'equal_width' 或 'equal_freq'
        """
        self.n_levels = n_levels
        self.k = k
        self.distance_weighted = distance_weighted
        self.binning = binning
        self.CO_values = None          # 特征对象对应的专家评价值 (dict: tuple -> score)
        self.bins = []                 # 每个准则的分箱边界
        self.scaler = MinMaxScaler()   # 用于归一化距离

    def _discretize(self, data):
        """
        根据训练数据确定每个准则的离散化边界
        :param data: 2D array (样本数 x 准则数)
        """
        n_criteria = data.shape[1]
        self.bins = []
        for j in range(n_criteria):
            col = data[:, j]
            if self.binning == 'equal_width':
                # 等宽分箱：按数值范围均匀切分
                bins = np.linspace(col.min(), col.max(), self.n_levels + 1)
            elif self.binning == 'equal_freq':
                # 等频分箱：按分位数切分
                percentiles = np.linspace(0, 100, self.n_levels + 1)
                bins = np.percentile(col, percentiles)
            else:
                raise ValueError("binning 仅支持 'equal_width' 或 'equal_freq'")
            self.bins.append(bins)

    def _map_to_level(self, value, col_idx):
        """
        将原始数值映射到离散等级（0 ~ n_levels-1）
        """
        bins = self.bins[col_idx]
        # 注意：bins[-1] 作为最大值，使用 np.digitize 并调整
        level = np.digitize(value, bins) - 1
        # 当 value 恰好等于 bins[-1] 时，level 会等于 n_levels，需要修正
        if level >= self.n_levels:
            level = self.n_levels - 1
        return level

    def _generate_COs(self):
        """
        生成所有可能的特征对象（笛卡尔积）
        返回：每个 CO 的等级组合（tuple of ints）
        """
        levels_per_criterion = [range(self.n_levels) for _ in range(len(self.bins))]
        all_COs = list(product(*levels_per_criterion))
        return all_COs

    def _get_CO_value(self, co_tuple):
        """
        获取某个特征对象的专家评价值 (需预先设置 self.CO_values)
        """
        return self.CO_values.get(co_tuple, None)

    def fit(self, data, expert_scores=None, expert_COs=None):
        """
        训练 INCOME 模型
        方式1：提供 expert_scores 与 data 一一对应（每个样本的专家评分）
        方式2：提供 expert_COs 及其评分 (expert_COs 为等级组合列表，expert_scores 对应评分)
        若两者均提供，优先使用方式2；若均未提供，则使用内置模拟函数生成。
        :param data: 训练数据（用于离散化）
        :param expert_scores: 与 data 行对应的专家评分列表
        :param expert_COs: 特征对象等级组合列表（tuple）
        """
        # 1. 离散化训练数据
        self._discretize(data)
        n_samples = data.shape[0]

        # 2. 生成所有特征对象
        all_COs = self._generate_COs()
        # 初始化 CO 评价值字典
        self.CO_values = {}

        # 3. 确定专家评价来源
        if expert_COs is not None and expert_scores is not None:
            # 直接使用给定的 CO 评价
            if len(expert_COs) != len(expert_scores):
                raise ValueError("expert_COs 与 expert_scores 长度不一致")
            for co, score in zip(expert_COs, expert_scores):
                self.CO_values[co] = score
            # 检查是否覆盖所有 CO？若未覆盖，后续 kNN 可能找不到邻居，但没关系，我们可用其他方式补充
            # 这里简单地用所有 CO 补全（若未评价，则使用函数模拟）
        else:
            # 若未提供专家评价，则根据 data 和 expert_scores（如果提供）生成
            if expert_scores is not None and len(expert_scores) == n_samples:
                # 将每个训练样本映射到等级组合，并将专家评分赋给相应 CO
                for i, row in enumerate(data):
                    co_tuple = tuple(self._map_to_level(row[j], j) for j in range(len(row)))
                    self.CO_values[co_tuple] = expert_scores[i]
            else:
                # 模拟专家评价：假设评分与训练数据的第一准则正相关（示例）
                print("⚠️ 未提供专家评价，使用模拟函数生成（评分 = 0.5 + 0.5 * 第一准则归一化值）")
                for i, row in enumerate(data):
                    co_tuple = tuple(self._map_to_level(row[j], j) for j in range(len(row)))
                    # 模拟评分：0~1，与第一个准则正相关
                    score = 0.5 + 0.5 * (row[0] - data[:,0].min()) / (data[:,0].max() - data[:,0].min() + 1e-10)
                    self.CO_values[co_tuple] = score

        # 若某些 CO 仍未被评价，补全为最邻近已评价 CO 的评分（或统一中值）
        # 此处使用简单补全：所有未赋值 CO 设为所有已评价 CO 评分的平均值
        if len(self.CO_values) < len(all_COs):
            mean_score = np.mean(list(self.CO_values.values()))
            for co in all_COs:
                if co not in self.CO_values:
                    self.CO_values[co] = mean_score
            print(f"⚠️ 部分 CO 未评价，已用平均值 {mean_score:.3f} 补全")

        # 将 CO 及其值存储为结构化数组以便 kNN 查询
        self.CO_list = list(self.CO_values.keys())
        self.CO_scores = np.array([self.CO_values[co] for co in self.CO_list])
        # 将等级组合转换为数值向量（用于距离计算）
        self.CO_vectors = np.array(self.CO_list, dtype=float)

        # 标准化特征对象向量（以便距离计算不受等级数量影响）
        self.scaler.fit(self.CO_vectors)

    # 依赖topsis() 函数
    def fit_by_madm(self, data, madm_method = 'topsis', weights = None, impacts = None):
        """
        不使用专家，而是使用传统MADM方法为所有CO生成评分
        """
        # 1. 离散化数据（确定边界）
        self._discretize(data)
        all_COs = self._generate_COs()

        # 2. 将离散等级还原为实际的数值（取每个等级区间的中位数）
        co_real_values = []
        for co_tuple in all_COs:
            row = []
            for j, level in enumerate(co_tuple):
                bins = self.bins[j]
                # 取该等级区间的中间值
                mid_val = (bins[level] + bins[level + 1]) / 2
                row.append(mid_val)
            co_real_values.append(row)
        co_real_values = np.array(co_real_values)

        # 3. 利用传统MADM计算得分（作为伪专家评分）
        if madm_method == 'topsis':
            scores = topsis(co_real_values, weights = weights, impacts = impacts)
        elif madm_method == 'waspas':
            scores = waspas(co_real_values, weights = weights)
        else:
            raise ValueError("请选择 topsis 或 waspas")

        # 4. 将得分存入CO字典
        self.CO_values = {}
        for co, score in zip(all_COs, scores):
            self.CO_values[co] = score

        # 5. 构建kNN索引（同之前的fit逻辑）
        self.CO_list = list(self.CO_values.keys())
        self.CO_scores = np.array([self.CO_values[co] for co in self.CO_list])
        self.CO_vectors = np.array(self.CO_list, dtype = float)
        self.scaler.fit(self.CO_vectors)

        print(f"✅ 已使用 {madm_method} 生成 {len(self.CO_values)} 个CO的伪专家评分")

    def predict(self, new_data):
        """
        对新样本进行评价
        :param new_data: 2D array (样本数 x 准则数)
        :return: 评分数组
        """
        new_data = np.array(new_data)
        scores = []
        for row in new_data:
            # 1. 将新样本映射到离散等级
            co_tuple = tuple(self._map_to_level(row[j], j) for j in range(len(row)))
            co_vec = np.array(co_tuple, dtype=float).reshape(1, -1)
            # 2. 标准化
            co_vec_norm = self.scaler.transform(co_vec)

            # 3. 计算与所有 CO 的距离 (欧氏距离)
            CO_vectors_norm = self.scaler.transform(self.CO_vectors)
            distances = np.linalg.norm(CO_vectors_norm - co_vec_norm, axis=1)
            # 4. 找出 k 个最近邻
            k_indices = np.argsort(distances)[:self.k]
            k_dists = distances[k_indices]
            k_scores = self.CO_scores[k_indices]

            # 5. 计算加权得分
            if self.distance_weighted:
                # 避免除零，添加小值
                weights = 1 / (k_dists + 1e-10)
                score = np.sum(weights * k_scores) / np.sum(weights)
            else:
                score = np.mean(k_scores)
            scores.append(score)
        return np.array(scores)


if __name__ == "__main__":
    # -------------------- 1. 准备数据 --------------------
    # 原始数据（未离散化）
    X_train = np.array([
        [0.2, 0.9],
        [0.5, 0.6],
        [0.8, 0.3],
        [0.1, 0.7],
        [0.9, 0.2],
        [0.3, 0.8]
    ])
    # 假设专家对每个训练样本给出了评分（0~1）
    y_train = np.array([0.3, 0.5, 0.7, 0.4, 0.8, 0.6])

    # 待评价
    X_test = np.array([
        [0.4, 0.5],
        [0.7, 0.4],
        [0.15, 0.85]
    ])

    # -------------------- 2. 创建 INCOME 模型 --------------------
    model = INCOME(n_levels=3, k=2, distance_weighted=True, binning='equal_width')
    # 训练：传入训练数据及其专家评分（也可以传入 CO 级评价，此处使用样本级评分）
    #model.fit(X_train, expert_scores=y_train)
    model.fit_by_madm(X_train)

    # -------------------- 3. 新方案 --------------------
    scores = model.predict(X_test)
    print("新方案得分:", np.round(scores, 3))
    print("排序 (从优到劣):", np.argsort(-scores) + 1)

    # -------------------- 4. 查看模型内部 --------------------
    print("离散化边界:", model.bins)
    print("特征对象数量:", len(model.CO_list))
    print("部分 CO 评分示例:", list(model.CO_values.items())[:5])