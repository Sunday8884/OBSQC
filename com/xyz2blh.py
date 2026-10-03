import math
from com.gnssParameter import CLIGHT, coordSystem
import numpy as np

def xyz2blh(X, Y, Z):
    """
    input :
            X           :X
            Y           :Y
            Z           :Z
    output:
            L           :经度(弧度)
            B           :纬度(弧度)
            H           :椭球高
    """
    a = 6378137.0
    e2 = 0.0066943799901413165
    b = 6356752.314245179
    L = np.arctan2(Y, X)
    N0 = a
    sqrtxy = np.sqrt(X * X + Y * Y)
    zps = Z / sqrtxy
    H0 = np.sqrt(X * X + Y * Y + Z * Z) - np.sqrt(a * b)
    B0 = np.arctan(zps / (1.0 - e2 * N0 / (N0 + H0)))
    while True:
        N1 = a / np.sqrt(1.0 - e2 * (np.sin(B0) ** 2.0))
        H = sqrtxy / np.cos(B0) - N1
        B = np.arctan(zps / (1.0 - e2 * N1 / (N1 + H)))
        if np.abs(H - H0) > 1e-3 and np.abs(B - B0) > 1.0e-4:
            N0 = N1
            H0 = H
            B0 = B
        else:
            break
    return B, L, H

# ECEF坐标系到东北天（NEU）局部坐标系的旋转矩阵
def bl2enuRot(lat, lon):
    sp = np.sin(lat)
    cp = np.cos(lat)
    sl = np.sin(lon)
    cl = np.cos(lon)
    enuRot = np.array([[-sl, cl, 0],
                  [-sp*cl, -sp*sl, cp],
                  [cp*cl, cp*sl, sp]])
    return enuRot

# 计算卫星到接收机的距离及方向余弦，
def geodist(rs, rr):
    """ calculate geometric distance """
    e = rr-rs
    r = np.linalg.norm(e)
    e = e/r
    r += 7.2921151467E-5*(rs[0]*rr[1]-rs[1]*rr[0])/CLIGHT
    return r, e