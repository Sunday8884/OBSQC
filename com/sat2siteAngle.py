import numpy as np
"""计算卫星相对于地面测站的几何角度，包括天顶角（Zenith）、方位角（Azimuth）和高度角（Elevation），输出单位均为弧度。"""
# 方位角
def sat2siteAngle(satX, satY, satZ, siteX, siteY, siteZ):
    """
    input :
            satX        :卫星X
            satY        :卫星Y
            satZ        :卫星Z
            siteX       :站点X
            siteX       :站点Y
            siteX       :站点Z
    output:
            Zenith      :天顶角(弧度)
            Azimuth     :方位角(弧度)
            Elevation   :高度角(弧度)
    """
    import math
    from com.xyz2neu import xyz2neu
    from com.xyz2blh import xyz2blh
    from com.xy2azi import xy2azi
    System = 'WGS84'
    north, east, up = xyz2neu(siteX, siteY, siteZ, satX, satY, satZ, System)
    siteB, siteL, H = xyz2blh(siteX, siteY, siteZ, System)
    satB, satL, H = xyz2blh(satX, satY, satZ, System)
    Azimuth = xy2azi(siteL, siteB, satL, satB)
    Elevation = math.atan(up / math.sqrt(north * north + east * east))
    Zenith = math.pi / 2 - Elevation
    return Zenith, Azimuth, Elevation

# 高度角
def getEle(rs, rr, enuRot):
    """
    input :
            satX        :卫星X
            satY        :卫星Y
            satZ        :卫星Z
            siteX       :站点X
            siteX       :站点Y
            siteX       :站点Z
    output:
            Elevation   :高度角(弧度)
    """
    deltaPos = rs - rr
    east, north, up = enuRot@deltaPos
    Azimuth = np.arctan2(east,north)
    Elevation = np.arctan(up / np.sqrt(north * north + east * east))
    Zenith = np.pi / 2 - Elevation
    return Zenith, Azimuth, Elevation