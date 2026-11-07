#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
pppfilter.py — 简明扩展卡尔曼滤波器（标量逐个量测更新）。

状态向量布局由 ppp.py 主引擎管理：
    [dx, dy, dz, clk_sys1, clk_sys2, ..., zwd, amb_track1, amb_track2, ...]
滤波器只负责 predict / update / 状态扩容，与具体模型解耦。
"""
import numpy as np


class KalmanFilter:
    def __init__(self, x0, P0):
        self.x = np.array(x0, dtype=np.float64).reshape(-1, 1)
        self.P = np.array(P0, dtype=np.float64).reshape(len(x0), len(x0))
        self.n = self.x.shape[0]

    # ---------------------------------------------------------- 状态扩容
    def add_state(self, xNew, varNew):
        """追加一个新状态（新弧段模糊度等），返回新状态的索引。"""
        self.x = np.vstack([self.x, [[float(xNew)]]])
        self.P = np.pad(self.P, ((0, 1), (0, 1)), mode='constant')
        self.P[self.n, self.n] = float(varNew)
        idx = self.n
        self.n += 1
        return idx

    # ---------------------------------------------------------- 时间更新
    def predict(self, Q):
        """
        时间更新（F = I，状态外推不变，方差累加过程噪声）。
        Q : n×n 过程噪声阵（增量）。
        """
        Q = np.array(Q, dtype=np.float64)
        if Q.shape != (self.n, self.n):
            raise ValueError('Q shape mismatch: %s vs %d' % (Q.shape, self.n))
        self.P = self.P + Q
        # 保持对称、正定
        self.P = 0.5 * (self.P + self.P.T)

    # ---------------------------------------------------------- 状态裁剪
    def remove_states(self, deadIdx):
        """
        删除指定索引的状态（已失效的模糊度弧段），返回 {旧索引: 新索引} 映射。
        调用方需按映射更新自己保存的状态索引。
        """
        dead = set(int(i) for i in deadIdx)
        if not dead:
            return {}
        if max(dead) >= self.n or min(dead) < 0:
            raise ValueError('dead index out of range')
        keep = [i for i in range(self.n) if i not in dead]
        mapping = {old: new for new, old in enumerate(keep)}
        self.x = self.x[keep, :].reshape(-1, 1)
        self.P = self.P[np.ix_(keep, keep)]
        self.n = len(keep)
        return mapping

    # ---------------------------------------------------------- 量测更新
    def update(self, H, z, R, hHat=None, guard=1e16, nsigma=None):
        """
        标量量测更新。
        H : 1×n 设计阵
        z : 量测值（obs）
        R : 量测方差（标量）
        hHat : 非线性量测模型在当前状态处的预测值 h(x)。
               量测模型非线性时（如距离 ρ=|rs−x|）必须显式传入，
               否则 innovation 按 z − H@x 计算（仅对线性模型成立）。
        nsigma : 粗差剔除阈值。若 |innov| > nsigma·sqrt(S)，
                 则不应用本次更新（先验拒绝），返回 (innov, S, True)。
        返回 (innov, innovVar, rejected)。
        滤波发散保护：若 S 或 Kalman 增益异常大则跳过本量测，返回 (None, None, False)。
        """
        H = np.array(H, dtype=np.float64).reshape(-1)
        if H.size < self.n:
            # 状态扩容（新弧段模糊度）后，旧量测对新状态的系数本就为 0，补零对齐
            H = np.pad(H, (0, self.n - H.size))
        elif H.size > self.n:
            raise ValueError('H size %d > state size %d' % (H.size, self.n))
        H = H.reshape(1, self.n)
        R = float(R)
        if hHat is None:
            hHat = float((H @ self.x).item())
        innov = float(np.asarray(z).reshape(-1)[0]) - float(hHat)
        HP = H @ self.P
        S = float((HP @ H.T).item()) + R
        if not np.isfinite(S) or S <= 0 or S > guard:
            return None, None, False
        if nsigma is not None and abs(innov) > float(nsigma) * np.sqrt(S):
            return innov, S, True
        K = (HP.T) / S
        self.x = self.x + K * innov
        IKH = np.eye(self.n) - K @ H
        # Joseph 形式，保证协方差对称正定
        self.P = IKH @ self.P @ IKH.T + K @ K.T * R
        self.P = 0.5 * (self.P + self.P.T)
        return innov, S, False
