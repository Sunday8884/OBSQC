#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
readSp3.py — 读取 IGS SP3 精密星历文件（sp3 / sp3-c / sp3-d 格式），并提供位置/钟差内插。

sp3Data 结构
------------
{
    'epochs'   : [datetime, ...]          # 全局历元列表（升序，所有卫星共用网格）
    'sampling' : float or None            # 采样间隔（秒），不规则时为 None
    'sat'      : {
        prn: {
            epoch(datetime): {
                'pos'    : np.array([x,y,z]) m,
                'clk'    : float s,
                'posStd' : float m   # 位置标准差（SP3-d P 行 / SP3-c 精度行；未知为 0）
                'clkStd' : float s   # 钟差标准差（未知为 0）
            }
        }
    },
    'bfact'    : (posBase, clkBase)       # std = base^exp（%f 行，缺省 2.0/2.0）
}

内插
----
- 位置：Lagrange 多项式内插（默认 10 阶，适用于 15min / 5min 采样）；
- 钟差：相邻两历元线性内插（SP3 钟差采样密，线性内插即可满足 dm 级 PPP）；
- 精度：取最近节点的 std（随历元缓变，无需内插）。
"""
import datetime
import numpy as np

# 缺省钟差标记（SP3 中以 999999.xxxxxxx 表示无钟差）
CLK_MISSING = 999999.0


def _parse_epoch_line(line):
    """解析 '*  YYYY MM DD HH MM SS.fffffff' 历元行，返回 datetime。"""
    try:
        year = int(line[3:7])
        month = int(line[8:10])
        day = int(line[11:13])
        hour = int(line[14:16])
        minute = int(line[17:19])
        sec = float(line[20:31])
        return datetime.datetime(year, month, day, hour, minute) + datetime.timedelta(seconds=sec)
    except Exception:
        return None


def readSp3(sp3File):
    """
    读取 SP3 文件，返回 sp3Data 字典（结构见模块 docstring）。
    支持 SP3-a/b/c/d；位置单位 km→m，钟差 μs→s。速度行（V*）忽略。
    精度（std）解析：
    - SP3-d：P 行 col 61-62/64-65/67-68 为 X/Y/Z std 指数（2 字符），col 70-72 为钟差 std 指数；
      std = base^exp，base 取头部 %f 行（缺省 2.0），位置 ×1e-3 (m)、钟差 ×1e-12 (s)
      （与 RTKLIB readsp3b 一致）；
    - SP3-c：位置记录后的精度行（每 4 字符一个整数指数，base=2），仅位置。
    """
    sp3Data = {'epochs': [], 'sampling': None, 'sat': {}, 'bfact': (2.0, 2.0)}
    nowEpoch = None
    lastWasP = False        # 上一非空行是否为 P 记录（识别 SP3-c 精度行）
    epochPsats = []         # 本历元 P 记录卫星顺序（SP3-c 精度行按此映射）

    with open(sp3File, 'r', encoding='utf-8', errors='ignore') as f:
        for line in f:
            if not line.strip():
                continue
            tag = line[0]

            # ---- 头部 ----
            if line.startswith('#') and not line.startswith('##'):
                # '#cV...' / '#dV...' 版本行
                continue
            if line.startswith('##'):
                # '##  week  sow  epoch_interval  ...' 采样间隔（秒）
                try:
                    tokens = line.split()
                    interval = float(tokens[3])
                    if interval > 0:
                        sp3Data['sampling'] = interval
                except Exception:
                    pass
                continue
            if line.startswith('%f'):
                # %f 行：std 底数（位置 base、钟差 base），与 RTKLIB readsp3h 一致
                try:
                    posBase = float(line[3:13])
                    clkBase = float(line[13:25])
                    if posBase > 0 and clkBase > 0:
                        sp3Data['bfact'] = (posBase, clkBase)
                except Exception:
                    pass
                continue

            # ---- 历元行 ----
            if tag == '*':
                nowEpoch = _parse_epoch_line(line)
                if nowEpoch is not None and nowEpoch not in sp3Data['epochs']:
                    sp3Data['epochs'].append(nowEpoch)
                epochPsats = []
                lastWasP = False
                continue

            # ---- 位置/钟差记录行 ----
            if tag == 'P' and nowEpoch is not None:
                prn = line[1:4]
                if len(prn) < 3 or not prn.strip():
                    continue
                try:
                    xkm = float(line[4:18])
                    ykm = float(line[18:32])
                    zkm = float(line[32:46])
                    clkus = float(line[46:60])
                except Exception:
                    continue
                # 缺省值判断
                if abs(xkm) < 1e-10 and abs(ykm) < 1e-10 and abs(zkm) < 1e-10:
                    continue
                clk = 0.0 if clkus >= CLK_MISSING else clkus * 1e-6

                # SP3-d：P 行尾部的 std 指数（62-63/65-66/68-69 位置，71-73 钟差；
                # 与 RTKLIB str2num(buff,61+j*3, j<3?2:3) 一致，0 基索引 61 起）
                posStd, clkStd = 0.0, 0.0
                posBase, clkBase = sp3Data['bfact']
                try:
                    exps = [float(line[61 + i * 3: 63 + i * 3]) for i in range(3)]
                    if posBase > 0 and min(exps) > 0 and max(exps) < 255:
                        posStd = (posBase ** float(np.mean(exps))) * 1e-3
                except Exception:
                    posStd = 0.0
                try:
                    expClk = float(line[70:73])
                    if clkBase > 0 and 0 < expClk < 255:
                        clkStd = (clkBase ** expClk) * 1e-12
                except Exception:
                    clkStd = 0.0

                sp3Data['sat'].setdefault(prn, {})[nowEpoch] = {
                    'pos': np.array([xkm * 1000.0, ykm * 1000.0, zkm * 1000.0]),
                    'clk': clk,
                    'posStd': posStd,
                    'clkStd': clkStd,
                }
                epochPsats.append(prn)
                lastWasP = True
                continue

            # ---- SP3-c 精度行：P 记录后紧跟的空格起始整数行 ----
            if lastWasP and tag == ' ' and nowEpoch is not None and epochPsats:
                vals = []
                ok = True
                for i in range(len(epochPsats)):
                    chunk = line[2 + i * 4: 6 + i * 4]
                    try:
                        v = int(chunk)
                    except Exception:
                        ok = False
                        break
                    vals.append(v)
                if ok and 0 < len(vals) <= 17:
                    for prn, v in zip(epochPsats, vals):
                        rec = sp3Data['sat'].get(prn, {}).get(nowEpoch)
                        if rec is not None and 0 < v < 255:
                            rec['posStd'] = float(2.0 ** v) * 1e-3
                    lastWasP = False
                continue

            lastWasP = False
            # 'V' 速度行、'EP'/'EV' 预报行、'/*' 注释、'++' 相关性等：忽略
            continue

    sp3Data['epochs'].sort()
    return sp3Data


def mergeSp3(sp3DataList):
    """合并多个 sp3Data（例如连续多天的精密星历），按卫星/历元去重（后读覆盖）。"""
    merged = {'epochs': [], 'sampling': None, 'sat': {}}
    seen_epochs = set()
    for sp3Data in sp3DataList:
        if sp3Data is None:
            continue
        if merged['sampling'] is None:
            merged['sampling'] = sp3Data.get('sampling')
        for ep in sp3Data['epochs']:
            if ep not in seen_epochs:
                seen_epochs.add(ep)
                merged['epochs'].append(ep)
        for prn, epDict in sp3Data['sat'].items():
            merged['sat'].setdefault(prn, {}).update(epDict)
    merged['epochs'].sort()
    return merged


def _lagrange(epochs, values, epoch, order):
    """
    Lagrange 多项式内插（通用标量版本，供外部/调试使用；引擎内使用向量化实现）。
    """
    n = len(epochs)
    if n < 2:
        return None
    order = min(order, n - 1)
    if order % 2 == 1:
        order -= 1
    if order < 2:
        order = 2 if n >= 3 else 1

    idx = 0
    best = abs((epochs[0] - epoch).total_seconds())
    for i in range(n):
        d = abs((epochs[i] - epoch).total_seconds())
        if d < best:
            best = d
            idx = i

    half = order // 2
    lo = max(0, idx - half)
    hi = min(n, lo + order + 1)
    lo = max(0, hi - order - 1)

    xs = np.array([(epochs[i] - epoch).total_seconds() for i in range(lo, hi)])
    ys = np.array(values[lo:hi], dtype=np.float64)

    result = 0.0
    m = len(xs)
    for j in range(m):
        w = 1.0
        for k in range(m):
            if k != j:
                w *= xs[k] / (xs[k] - xs[j])
        result += w * ys[j]
    return result


def getSp3PosClk(sp3Data, prn, epoch, posOrder=10, clkInterp='linear',
                 wantVel=False, wantStd=False):
    """
    内插得到 prn 在 epoch 时刻的位置与钟差。
    返回 (pos: np.array(3) m, clk: float s)；无数据/超出覆盖范围时返回 (None, None)。
    wantVel=True 时返回 (pos, clk, vel)——速度由同一批 Lagrange 节点求导得到，
    供相对论钟差改正（−2r·v/c²）与卫星姿态计算使用，几乎无额外开销。
    wantStd=True 时返回 (pos, clk, vel, posStd, clkStd)——std 取最近节点值
    （未知为 0），供量测方差 var_rs 加权使用（隐含 wantVel）。
    首次调用后按卫星缓存排序历元列表与时间戳（秒），后续用二分查找加速。
    """
    retNone = (None, None)
    if wantStd:
        retNone = (None, None, None, None, None)
    elif wantVel:
        retNone = (None, None, None)

    if prn not in sp3Data['sat']:
        return retNone
    satDict = sp3Data['sat'][prn]
    if not satDict:
        return retNone

    # ---- 缓存：排序列表 + 秒时间戳（相对首历元）----
    cache = sp3Data.setdefault('_cache', {})
    if prn not in cache:
        eps = sorted(satDict.keys())
        t0 = eps[0]
        cache[prn] = {
            'epochs': eps,
            'ts': np.array([(ep - t0).total_seconds() for ep in eps]),
            'pos': np.array([satDict[ep]['pos'] for ep in eps]),   # (n,3)
            'clk': np.array([satDict[ep]['clk'] for ep in eps]),   # (n,)
            'posStd': np.array([satDict[ep].get('posStd', 0.0) for ep in eps]),
            'clkStd': np.array([satDict[ep].get('clkStd', 0.0) for ep in eps]),
            't0': t0,
        }
    c = cache[prn]
    ts = c['ts']
    tTarget = (epoch - c['t0']).total_seconds()

    n = len(ts)
    if tTarget < ts[0] - 1.0 or tTarget > ts[-1] + 1.0:
        return retNone

    # 精确/邻近节点索引（二分）
    idx = int(np.searchsorted(ts, tTarget))
    if idx < n and ts[idx] == tTarget:
        pos, clk = c['pos'][idx].copy(), float(c['clk'][idx])
        vel = None
        if wantStd or wantVel:
            # 精确节点处：用相邻节点差分近似速度（仅用于改正量级 ~cm/s 项，足够）
            i0, i1 = max(idx - 1, 0), min(idx + 1, n - 1)
            vel = (c['pos'][i1] - c['pos'][i0]) / max(ts[i1] - ts[i0], 1e-9)
        if wantStd:
            return pos, clk, vel, float(c['posStd'][idx]), float(c['clkStd'][idx])
        if wantVel:
            return pos, clk, vel
        return pos, clk

    # ---- 位置：Lagrange（以索引切片参与计算）----
    order = min(posOrder, n - 1)
    if order % 2 == 1:
        order -= 1
    if order < 1:
        order = 1
    half = order // 2
    lo = max(0, idx - half)
    hi = min(n, lo + order + 1)
    lo = max(0, hi - order - 1)

    xs = ts[lo:hi] - tTarget
    ys = c['pos'][lo:hi]          # (m,3)
    m = len(xs)
    # 权重向量：l_j(0) = Π_{k≠j} x_k / (x_k - x_j)
    w = np.ones(m)
    for j in range(m):
        for k in range(m):
            if k != j:
                w[j] *= xs[k] / (xs[k] - xs[j])
    pos = w @ ys
    if not np.isfinite(pos).all():
        return retNone

    # 速度：Lagrange 基函数在 0 处的导数 l_j'(0) = −l_j(0)·Σ_{k≠j} 1/x_k
    vel = None
    if wantStd or wantVel:
        dl = np.array([-w[j] * np.sum(1.0 / np.delete(xs, j)) for j in range(m)])
        vel = dl @ ys
        if not np.isfinite(vel).all():
            return retNone

    # ---- 钟差：线性 ----
    clks = c['clk']
    if idx == 0:
        clk = float(clks[0])
    elif idx >= n:
        clk = float(clks[-1])
    else:
        t0e, t1e = ts[idx - 1], ts[idx]
        a = (tTarget - t0e) / (t1e - t0e)
        clk = float(clks[idx - 1] * (1.0 - a) + clks[idx] * a)

    if wantStd:
        # std 取最近节点（随历元缓变，无需内插）
        iNear = idx if idx < n else n - 1
        if idx > 0 and abs(ts[idx - 1] - tTarget) < abs(ts[iNear] - tTarget):
            iNear = idx - 1
        return pos, clk, vel, float(c['posStd'][iNear]), float(c['clkStd'][iNear])
    if wantVel:
        return pos, clk, vel
    return pos, clk


if __name__ == '__main__':
    import sys
    if len(sys.argv) > 1:
        sp3Data = readSp3(sys.argv[1])
        print('历元数:', len(sp3Data['epochs']), '采样间隔:', sp3Data['sampling'])
        print('卫星数:', len(sp3Data['sat']))
        if sp3Data['epochs']:
            ep = sp3Data['epochs'][len(sp3Data['epochs']) // 2]
            for prn in list(sp3Data['sat'])[:2]:
                pos, clk = getSp3PosClk(sp3Data, prn, ep)
                print(prn, ep, pos, clk)
