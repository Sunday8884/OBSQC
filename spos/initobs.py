from itertools import islice
from com.gnssTime import datetime2doys, datetime2mjds

"""
SPP初始化：
初始化每个历元观测数据的函数。它从原始观测数据中筛选出有效卫星（具有双频伪距观测值且与星历匹配），提取伪距观测值并组织成便于定位解算的结构。
"""
def initObs(bandChoose, obsData, epoch, satPosInEpoch, satClkInEpoch):
    # time
    doys = datetime2doys(epoch)
    mjds = datetime2mjds(epoch)

    # init data
    obsDataInTime = {}
    obsDataInTime['gSysInf'] = {}
    obsDataInTime['obs'] = {}
    obsDataInTime['prnIndex'] = {}
    obsDataInTime['satUse'] = []
    obsDataInTime['doys'] = doys
    obsDataInTime['mjds'] = mjds

    # check data
    prni = 0
    for prn in obsData[epoch]:
        gnssSys = prn[0]
        if gnssSys not in bandChoose:
            continue
        BAND1, BAND2 = list(islice(bandChoose[gnssSys].keys(), 2))
        P1, P2 = obsData[epoch][prn][BAND1], obsData[epoch][prn][BAND2] #读取双频伪距观测值

        # 有效性检查：有星历位置和钟差，P不空，剔除粗差 > 30。
        # satPosInEpoch 和 satClkInEpoch 存储每个历元的卫星位置和钟差数据
        if prn not in satPosInEpoch or prn not in satClkInEpoch:
            continue
        if None in [P1, P2]:
            continue
        if abs(P1 - P2) > 30:
            continue

        # 存储可用卫星
        obsDataInTime['prnIndex'][prn] = prni
        prni += 1
        obsDataInTime['satUse'].append(prn)
        obsDataInTime['obs'][prn] = {}
        obsDataInTime['obs'][prn]['P1'] = P1
        obsDataInTime['obs'][prn]['P2'] = P2
        if gnssSys not in obsDataInTime['gSysInf']:
            obsDataInTime['gSysInf'][gnssSys] = {}
            obsDataInTime['gSysInf'][gnssSys]['lambda1'] = bandChoose[gnssSys][BAND1]['lambda']
            obsDataInTime['gSysInf'][gnssSys]['lambda2'] = bandChoose[gnssSys][BAND2]['lambda']
            obsDataInTime['gSysInf'][gnssSys]['f1'] = bandChoose[gnssSys][BAND1]['freq']
            obsDataInTime['gSysInf'][gnssSys]['f2'] = bandChoose[gnssSys][BAND2]['freq']
    isLess = False

    if len(obsDataInTime['satUse']) < 4: # 四颗卫星
        isLess = True
    return obsDataInTime, isLess