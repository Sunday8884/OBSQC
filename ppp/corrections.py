#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
corrections.py — PPP 误差改正模型：
    - gmst()            格林尼治恒星时
    - sunPosECEF()      太阳低精度星历（ECEF，m）
    - moonPosECEF()     月亮低精度星历（ECEF，m）
    - solidEarthTide()  固体潮位移改正（度-2主项，径向最大 ~0.3 m）
    - phaseWindUp()     相位缠绕改正（连续跟踪累计，单位：周）
    - tropMapping()     对流层湿延迟映射函数（与 spos.trop 同款）
"""
import math
import numpy as np

from com.gnssParameter import CLIGHT, coordSystem

OMEGA_E = coordSystem['WGS84'].omega


# ---------------------------------------------------------------- 恒星时
def gmst(epoch, dut1=0.0):
    """格林尼治平恒星时（弧度），输入 datetime（UTC）；dut1 = UT1−UTC（秒）。"""
    jd = epoch.toordinal() + 0.5 + \
        (epoch.hour * 3600 + epoch.minute * 60 + epoch.second + epoch.microsecond * 1e-6
         + dut1) / 86400.0
    tut = (jd - 2451545.0) / 36525.0
    # 度
    g = 280.46061837 + 360.98564736629 * (jd - 2451545.0) + 0.000387933 * tut * tut - tut ** 3 / 38710000.0
    return math.radians(g % 360.0)


# ---------------------------------------------------------------- 太阳
def sunPosECEF(epoch, dut1=0.0):
    """
    太阳低精度星历（精度 ~0.01°，满足固体潮/相位缠绕/卫星姿态计算）。
    返回 ECEF 坐标（m）。
    """
    jd = epoch.toordinal() + 0.5 + \
        (epoch.hour * 3600 + epoch.minute * 60 + epoch.second + epoch.microsecond * 1e-6) / 86400.0
    n = jd - 2451545.0
    L = math.radians(280.460 + 0.9856474 * n)          # 平黄经
    g = math.radians(357.528 + 0.9856003 * n)          # 平近点角
    lam = L + math.radians(1.915) * math.sin(g) + math.radians(0.020) * math.sin(2 * g)
    eps = math.radians(23.439 - 4.0e-7 * n)
    rAU = 1.00014 - 0.01671 * math.cos(g) - 0.00014 * math.cos(2 * g)
    AU = 1.495978707e11

    x = rAU * AU * math.cos(lam)
    y = rAU * AU * math.cos(eps) * math.sin(lam)
    z = rAU * AU * math.sin(eps) * math.sin(lam)

    # 惯性系 → ECEF
    theta = gmst(epoch, dut1)
    ct, st = math.cos(theta), math.sin(theta)
    return np.array([ct * x - st * y, st * x + ct * y, z])


# ---------------------------------------------------------------- 月亮
def moonPosECEF(epoch, dut1=0.0):
    """
    月亮低精度星历（精度 ~几角分，满足固体潮计算）。返回 ECEF 坐标（m）。
    （Meeus 主项级数，与 RTKLIB sun_moon 一致的主项子集）
    """
    mjd = epoch.toordinal() - 678576 + \
        (epoch.hour * 3600 + epoch.minute * 60 + epoch.second + epoch.microsecond * 1e-6) / 86400.0
    t = (mjd - 51544.5) / 36525.0

    f = [0.0] * 5
    f[0] = math.radians(218.3164477 + 481267.88123421 * t - 0.0015786 * t * t
                        + t ** 3 / 538841.0 - t ** 4 / 65194000.0)
    f[1] = math.radians(134.9633964 + 477198.8675055 * t + 0.0087414 * t * t
                        + t ** 3 / 69699.0 - t ** 4 / 14712000.0)
    f[2] = math.radians(93.2720950 + 483202.0175233 * t - 0.0036539 * t * t
                        - t ** 3 / 3526000.0 + t ** 4 / 863310000.0)
    f[3] = math.radians(297.8501954 + 445267.1114803 * t - 0.0019142 * t * t
                        + t ** 3 / 189475.0)
    f[4] = math.radians(125.0445479 - 1934.1362891 * t + 0.0020754 * t * t
                        + t ** 3 / 467441.0 - t ** 4 / 60616000.0)

    # 黄经（度）
    lon = (218.3164477 + 481267.88123421 * t
           + 6.288774 * math.sin(f[1])
           + 1.274027 * math.sin(2 * f[0] - f[1])
           + 0.658314 * math.sin(2 * f[0])
           + 0.213618 * math.sin(2 * f[1])
           - 0.185116 * math.sin(f[3])
           - 0.114332 * math.sin(2 * f[2])
           + 0.058793 * math.sin(2 * f[0] - 2 * f[1])
           + 0.057066 * math.sin(2 * f[0] - f[3] - f[1])
           + 0.053322 * math.sin(2 * f[0] + f[1])
           + 0.045758 * math.sin(2 * f[0] - f[2]))
    # 黄纬（度）
    lat = (5.128122 * math.sin(f[2])
           + 0.280602 * math.sin(f[0] + f[2])
           + 0.277693 * math.sin(f[1] + f[2])
           + 0.173237 * math.sin(2 * f[0] - f[2])
           + 0.055413 * math.sin(2 * f[1] - f[2])
           + 0.046271 * math.sin(2 * f[0] - 2 * f[2] + f[1]) * 0  # 保守取主项
           )
    # 地心距（m）
    r = (385000.56
         - 20905.355 * math.cos(f[1])
         - 3699.111 * math.cos(2 * f[0] - f[1])
         - 2955.968 * math.cos(2 * f[0])
         - 569.925 * math.cos(2 * f[1])) * 1000.0

    lon = math.radians(lon)
    lat = math.radians(lat)
    eps = math.radians(23.4392)

    xi = r * math.cos(lat) * math.cos(lon)
    yi = r * (math.cos(lat) * math.sin(lon) * math.cos(eps) - math.sin(lat) * math.sin(eps))
    zi = r * (math.cos(lat) * math.sin(lon) * math.sin(eps) + math.sin(lat) * math.cos(eps))

    theta = gmst(epoch, dut1)
    ct, st = math.cos(theta), math.sin(theta)
    return np.array([ct * xi - st * yi, st * xi + ct * yi, zi])


# ---------------------------------------------------------------- 固体潮
def solidEarthTide(epoch, siteXyz, dut1=0.0):
    """
    固体潮测站位移改正（度-2 主项，频率无关部分）。
    精度说明：忽略 K1/K2 频率相关项与度-3/4 项，最大遗漏 ~1-2 cm，
    对浮点模糊度 PPP 影响可忽略；如需 mm 级请外接 IERS 约令完整模型。

    siteXyz : 测站 ECEF 坐标（m）
    返回位移向量 [dx, dy, dz]（m，加到测站坐标上）。
    """
    GM_E = coordSystem['WGS84'].GM
    Re = coordSystem['WGS84'].a
    GM_M = GM_E * 0.0123000343   # 月亮引力常数（地月质量比 0.0123000343）
    GM_S = GM_E * 332946.0487    # 太阳引力常数
    H2 = 0.6078                  # Love 数（度-2，径向）
    L2 = 0.0847                  # Shida 数（度-2，切向）

    k = np.array(siteXyz, dtype=np.float64)
    rSite = np.linalg.norm(k)
    k /= rSite

    disp = np.zeros(3)
    for bodyPos, GMb in ((moonPosECEF(epoch, dut1), GM_M), (sunPosECEF(epoch, dut1), GM_S)):
        rb = np.linalg.norm(bodyPos)
        m = bodyPos / rb
        A = (GMb / GM_E) * (Re ** 4) / (rb ** 3)   # 特征幅度（月亮 ~0.358 m，太阳 ~0.165 m）
        km = float(np.dot(k, m))
        # 径向分量
        dr = A * H2 * 0.5 * (3.0 * km * km - 1.0)
        # 切向分量（指向天体）
        mPerp = m - km * k
        disp += dr * k + 3.0 * A * L2 * km * mPerp
    return disp


# ---------------------------------------------------------------- 相位缠绕
def _unit(v):
    return v / np.linalg.norm(v)


def phaseWindUp(prn, satPos, recPos, sunPosECEF, prevWindup):
    """
    计算卫星-接收机相位缠绕角（周）。

    prn         : 卫星号（保留接口）
    satPos/recPos : ECEF（m）
    sunPosECEF  : 太阳 ECEF（m）
    prevWindup  : 上一历元的累计缠绕角（周）；首历元传 0.0

    返回累计缠绕角（周，连续不跳变）。卫星姿态按标准偏航模型（EZ 指向地心）。
    """
    rs = np.array(satPos, dtype=np.float64)
    rr = np.array(recPos, dtype=np.float64)
    rSun = np.array(sunPosECEF, dtype=np.float64)

    # 卫星本体坐标基（标准姿态：EZ 指向地心）
    ez = _unit(-rs)
    es = _unit(rSun - rs)
    ey = _unit(np.cross(ez, es))
    ex = np.cross(ey, ez)
    if not (np.isfinite(ex).all() and np.isfinite(ey).all()):
        return prevWindup

    # 接机固定偶极子（ENU → ECEF）
    ekw = _unit(rr)
    exr = np.cross(np.array([0.0, 0.0, 1.0]), ekw)
    exr /= np.linalg.norm(exr)
    eyr = np.cross(ekw, exr)
    # GPS 右旋约定
    exg = np.array([-exr[0], -exr[1], exr[2]])
    eyg = np.array([-eyr[0], -eyr[1], eyr[2]])

    # 视线单位矢量（接收机 → 卫星）
    ekr = _unit(rs - rr)

    # 复数偶极子投影
    dr = complex(np.dot(exg, ekr), np.dot(eyg, ekr))
    ds = complex(np.dot(ex, ekr), np.dot(ey, ekr))
    if not (np.isfinite(dr.real) and np.isfinite(dr.imag)
            and np.isfinite(ds.real) and np.isfinite(ds.imag)):
        return prevWindup
    if abs(dr) < 1e-12 or abs(ds) < 1e-12:
        return prevWindup

    # 本历元缠绕角（弧度）
    ph = math.atan2(ds.imag, ds.real) - math.atan2(dr.imag, dr.real)
    if not np.isfinite(ph):
        return prevWindup
    # 相对上一历元解缠（选择离 prevWindup 最近的 2π 分支）
    prevRad = prevWindup * 2.0 * math.pi
    n = round((prevRad - ph) / (2.0 * math.pi))
    phNow = ph + 2.0 * math.pi * n
    return phNow / (2.0 * math.pi)


# ---------------------------------------------------------------- 对流层
def tropMapping(elevation):
    """湿延迟映射函数（与 spos/trop.py 同款简化 Niell 形式），elevation 为弧度。"""
    return 1.001 / np.sqrt(0.002001 + np.sin(elevation) ** 2)
