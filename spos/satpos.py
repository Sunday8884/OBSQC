from itertools import islice
import numpy as np
import datetime

from com.gnssParameter import CLIGHT, bds3MEOList
from com.nav2posclk import getClkInNav, getPosInNav

"""
给定一个历元 epoch，遍历该历元下所有卫星，利用双频伪距估算信号传播时间，迭代修正发射时刻，
然后从导航文件中获取该时刻的卫星位置、速度和钟差。
"""
def satPos(bandChoose, obsData, navData, epoch):
    """
    bandChoose：双频配置字典，结构如 {'G': {'C1C': {...}, 'C2W': {...}}}。

    obsData：观测数据字典，{epoch: {prn: {band: value}}}。

    navData：导航数据字典（由 readNav 读取）。

    epoch：当前接收时刻（datetime 对象）。
    """
    satPosInEpoch = {}
    satClkInEpoch = {}
    for prn in obsData[epoch]: # 遍历当前历元的所有卫星
        gnssSys = prn[0]
        if gnssSys not in bandChoose:
            continue
        if gnssSys == 'C' and (prn not in bds3MEOList): # 仅处理北斗三号 MEO
            continue
        if prn not in navData:
            continue
        BAND1, BAND2 = list(islice(bandChoose[gnssSys].keys(), 2))
        C1, C2 = obsData[epoch][prn][BAND1], obsData[epoch][prn][BAND2] # 获取双频伪距
        if C1 is None or C2 is None:
            continue
        bandDif = abs(C1 - C2)
        if bandDif > 100:  # 剔除粗差
            continue
        aveP = (C1 + C2) / 2
        # 发射时间 = 接收时间 - 传播时间(距离/光速)
        epochFind = epoch - datetime.timedelta(seconds=aveP/CLIGHT)
        satClk = getClkInNav(navData, prn, epochFind)

        if satClk is None:
            continue

        # 最终发射时刻 = 发射时刻 - 卫星钟差
        epochFind = epochFind - datetime.timedelta(seconds=satClk)
        satX, satY, satZ, satVX, satVY, satVZ = getPosInNav(navData, prn, epochFind, epochFind)
        if satX is not None and satClk is not None:
            satPosInEpoch[prn] = np.array([satX, satY, satZ, satVX, satVY, satVZ])
            satClkInEpoch[prn] = satClk
    return satPosInEpoch, satClkInEpoch