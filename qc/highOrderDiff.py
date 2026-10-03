from com.gnssParameter import getBandFreq, CLIGHT, bandComb, obsFreq
import numpy as np
import datetime
"""
实现了 高阶差分（High‑Order Difference） 的计算，用于评估 GNSS 载波相位观测值的噪声水平。
通过对连续弧段内载波相位观测值（或几何距离组合）进行多次差分，可以消除低阶项（如常数、线性趋势），
保留高阶噪声，从而衡量观测值的平滑度或高频抖动，RMS 值越小，表明载波相位越平滑，数据质量越高。
"""
#计算一个序列的 N 阶差分
def culHighOrderDiff(epochOrderData):
    diffCount = len(epochOrderData)
    backDiff = epochOrderData
    while diffCount > 1:
        updateDiff = []
        for diffIndex in range(len(backDiff) - 1, 0, -1):
            updateDiff.append(backDiff[diffIndex] - backDiff[diffIndex-1])
        backDiff = updateDiff[::-1]
        diffCount += -1
    return backDiff[0]

"""
对观测文件中的每个卫星‑波段组合，提取连续无周跳的载波相位弧段，计算其 高阶差分值，并输出为时间序列。
highData = {
    'G': {
        'G01': {
            'L1C': {
                'epoch': [datetime, ...],
                'highOrder': [float, ...]   # 高阶差分值 (单位: m)
            },
            ...
        },
        ...
    }
}
"""
def getObsHighOrderDiff(obsHead, obsData, diffNum=4, start_datetime=None, end_datetime=None,
                        sat_list=None, sys_list=None, band_choose=None):
    """
    Calculate high-order differences for carrier phase observations.

    Parameters:
    ----------
    obsHead : dict
        Observation file header.
    obsData : dict
        Observation data {epoch: {prn: {band: value}}}.
    diffNum : int
        Order of difference (default 4).
    start_datetime : datetime, optional
        Start time for processing (inclusive).
    end_datetime : datetime, optional
        End time for processing (inclusive).
    sat_list : list of str, optional
        List of PRNs to process. If None, use all from obsHead['prn'].
    sys_list : list of str, optional
        List of GNSS system characters (e.g., ['G','C','E']). If None, use all systems.
    band_choose : dict, optional
        Dictionary specifying which bands to use per system, e.g., {'G': ['1C','2W']}.

    Returns:
    -------
    highData : dict
        {sys: {prn: {band: {'epoch': list, 'highOrder': list}}}}
    """
    obsType = obsHead['OBS TYPES']
    bandChoose = {}
    for gSys in obsType:
        if gSys not in bandChoose:
            bandChoose[gSys] = {}
        for band in obsType[gSys]:
            if 'C' == band[0] or 'P' == band[0]:
                if band[1] not in bandChoose[gSys]:
                    bandChoose[gSys][band[1]] = []
                bandChoose[gSys][band[1]].append(band)
    highBandChoose = {}

    # Determine satellite list with system filter
    if sat_list is None:
        sat_list = obsHead['prn']
    if sys_list is not None:
        sat_list = [prn for prn in sat_list if prn[0] in sys_list]
    # Band filter
    band_choose_filter = band_choose if band_choose is not None else {}

    for gSys in obsType:
        # Skip system if not in sys_list
        if sys_list is not None and gSys not in sys_list:
            continue
        if len(bandChoose[gSys]) < 2 or gSys not in obsFreq:
            continue
        if band_choose_filter and gSys not in band_choose_filter:
            continue
        if gSys not in highBandChoose:
            highBandChoose[gSys] = {}
        for band in obsType[gSys]:
            if 'C' == band[0] or 'P' == band[0]:
                if band_choose_filter and gSys in band_choose_filter:
                    if band[1:] not in band_choose_filter[gSys]:
                        continue
                if band in bandChoose[gSys][list(bandChoose[gSys])[0]]:
                    band2 = bandChoose[gSys][list(bandChoose[gSys])[1]][0]
                else:
                    band2 = bandChoose[gSys][list(bandChoose[gSys])[0]][0]

                bandL = 'L' + band[1:]
                band2L = 'L' + band2[1:]
                if bandL not in obsType[gSys] or band2 not in obsType[gSys] or band2L not in obsType[gSys]:
                    continue

                band1freq = getBandFreq(gSys, band)
                band2freq = getBandFreq(gSys, band2)

                if band1freq is None or band2freq is None:
                    continue

                band1freq = band1freq * 1.e6
                band1lamda = CLIGHT / band1freq

                band2freq = band2freq * 1.e6
                band2lamda = CLIGHT / band2freq

                highBandChoose[gSys][band] = {}
                highBandChoose[gSys][band]['band1C'] = band
                highBandChoose[gSys][band]['band2C'] = band2
                highBandChoose[gSys][band]['band1L'] = bandL
                highBandChoose[gSys][band]['band2L'] = band2L
                highBandChoose[gSys][band]['band1freq'] = band1freq
                highBandChoose[gSys][band]['band1lambda'] = band1lamda
                highBandChoose[gSys][band]['band2freq'] = band2freq
                highBandChoose[gSys][band]['band2lambda'] = band2lamda

                highBandChoose[gSys][band]['COEFL1'] = -(band1freq * band1freq + band2freq*band2freq)/(band1freq*band1freq - band2freq*band2freq)
                highBandChoose[gSys][band]['COEFL2'] = 2*band2freq*band2freq/(band1freq*band1freq - band2freq*band2freq)

                lamdaW = CLIGHT / (band1freq - band2freq)
                highBandChoose[gSys][band]['lamdaW'] = lamdaW

    # Get sampling interval from first two epochs
    sorted_epochs = sorted(obsData.keys())
    if len(sorted_epochs) < 2:
        interval = 30.0
    else:
        interval = (sorted_epochs[1] - sorted_epochs[0]).total_seconds()

    highData = {}
    for prn in sat_list:
        gSys = prn[0]
        if gSys not in obsFreq or gSys not in highBandChoose:
            continue

        for band in highBandChoose[gSys]:
            bandEpochList = []
            highList = []
            highGoodList = []
            epochGoodList = []
            band1L = highBandChoose[gSys][band]['band1L']
            band2L = highBandChoose[gSys][band]['band2L']
            band1C = highBandChoose[gSys][band]['band1C']
            band2C = highBandChoose[gSys][band]['band2C']
            freq1 = highBandChoose[gSys][band]['band1freq']
            freq2 = highBandChoose[gSys][band]['band2freq']
            lamda1 = highBandChoose[gSys][band]['band1lambda']
            lamda2 = highBandChoose[gSys][band]['band2lambda']
            lamdaW = highBandChoose[gSys][band]['lamdaW']
            COEFL1 = highBandChoose[gSys][band]['COEFL1']
            COEFL2 = highBandChoose[gSys][band]['COEFL2']

            for epoch in sorted_epochs:
                if prn not in obsData[epoch].keys():
                    continue
                if start_datetime is not None and epoch < start_datetime:
                    continue
                if end_datetime is not None and epoch > end_datetime:
                    continue
                L1 = obsData[epoch][prn].get(band1L)
                L2 = obsData[epoch][prn].get(band2L)
                P1 = obsData[epoch][prn].get(band1C)
                P2 = obsData[epoch][prn].get(band2C)
                if L1 is None or L2 is None or P1 is None or P2 is None:
                    continue
                if abs(P2 - P1) > 100:
                    continue

                N_MW = L1 - L2 - (freq1 * P1 + freq2 * P2) / (freq1 + freq2) / lamdaW
                L_GF = lamda1 * L1 - lamda2 * L2
                P_GF = P2 - P1
                nowL = L1 * lamda1

                MW_CLIP = False
                GF_CLIP = False
                if bandEpochList == []:
                    bandEpochList = [epoch]
                    N_MW_list = [[N_MW]]
                    MEAN_MW_list = [[N_MW]]
                    SIGMA_list = [[0.3]]
                    L_GF_list = [[L_GF]]
                    P_GF_list = [[P_GF]]
                    slipList = [0]
                    highList = [[nowL]]
                    continue
                if (epoch - bandEpochList[-1]).total_seconds() > interval:
                    bandEpochList.append(epoch)
                    N_MW_list.append([N_MW])
                    MEAN_MW_list.append([N_MW])
                    SIGMA_list = [[0.3]]
                    L_GF_list.append([L_GF])
                    P_GF_list.append([P_GF])
                    slipList.append(0)
                    highList.append([nowL])
                    continue

                if len(N_MW_list[-1]) < 2:
                    SIGMA_I_BACK = 0.15
                else:
                    SIGMA_I_BACK = SIGMA_list[-1][-1]

                N_MW_MEAN_BACK = MEAN_MW_list[-1][-1]
                N_MW_MEAN = N_MW_MEAN_BACK + (N_MW - N_MW_MEAN_BACK)/(len(MEAN_MW_list[-1]) + 1)
                SIGMA2 = SIGMA_I_BACK * SIGMA_I_BACK + ((N_MW - N_MW_MEAN_BACK) ** 2 - SIGMA_I_BACK*SIGMA_I_BACK) / (len(MEAN_MW_list[-1]) + 1)
                SIGMA = np.sqrt(SIGMA2)
                if abs(N_MW - N_MW_MEAN) >= 4 * SIGMA_I_BACK:
                    MW_CLIP = True

                # GF-CHECK
                if not MW_CLIP and len(L_GF_list[-1]) > 1:
                    if abs(L_GF - L_GF_list[-1][-1]) > 0.15:
                        GF_CLIP = True

                if MW_CLIP or GF_CLIP:
                    bandEpochList.append(epoch)
                    N_MW_list.append([N_MW])
                    MEAN_MW_list.append([N_MW])
                    SIGMA_list.append([0.15])
                    L_GF_list.append([L_GF])
                    P_GF_list.append([P_GF])
                    slipList.append(1)
                    highList.append([nowL])
                else:
                    bandEpochList.append(epoch)
                    N_MW_list[-1].append(N_MW)
                    MEAN_MW_list[-1].append(N_MW_MEAN)
                    L_GF_list[-1].append(L_GF)
                    SIGMA_list[-1].append(SIGMA)
                    P_GF_list[-1].append(P_GF)
                    slipList.append(0)
                    highList[-1].append(nowL)
                    if len(highList[-1]) > diffNum * 2:
                        difIndex = (diffNum + 1)*(-1)
                        nowDiff = culHighOrderDiff(highList[-1][difIndex:])
                        if diffNum == 4:
                            nowDiff = nowDiff / np.sqrt(70)
                        else:
                            nowDiff = nowDiff / np.sqrt(20)
                        if abs(nowDiff) > 100:
                            continue
                        highGoodList.append(nowDiff)
                        epochGoodList.append(epoch)
            if len(highGoodList) == 0:
                continue
            if gSys not in highData:
                highData[gSys] = {}
            if prn not in highData[gSys]:
                highData[gSys][prn] = {}
            highData[gSys][prn][band1L] = {}
            highData[gSys][prn][band1L]['epoch'] = epochGoodList
            highData[gSys][prn][band1L]['highOrder'] = highGoodList

    return highData

# 简化版（无 MW/GF 周跳检测），同样增加 sys_list
def getObsHighOrderDiffcc(obsHead, obsData, diffNum=4, start_datetime=None, end_datetime=None,
                          sat_list=None, sys_list=None, band_choose=None):
    """
    Simplified version of getObsHighOrderDiff (no MW/GF cycle slip detection).
    """
    obsType = obsHead['OBS TYPES']
    bandChoose = {}
    for gSys in obsType:
        if gSys not in bandChoose:
            bandChoose[gSys] = {}
        for band in obsType[gSys]:
            if 'C' == band[0] or 'P' == band[0]:
                if band[1] not in bandChoose[gSys]:
                    bandChoose[gSys][band[1]] = []
                bandChoose[gSys][band[1]].append(band)
    highBandChoose = {}

    if sat_list is None:
        sat_list = obsHead['prn']
    if sys_list is not None:
        sat_list = [prn for prn in sat_list if prn[0] in sys_list]
    band_choose_filter = band_choose if band_choose is not None else {}

    for gSys in obsType:
        if sys_list is not None and gSys not in sys_list:
            continue
        if len(bandChoose[gSys]) < 2 or gSys not in obsFreq:
            continue
        if band_choose_filter and gSys not in band_choose_filter:
            continue
        if gSys not in highBandChoose:
            highBandChoose[gSys] = {}
        for band in obsType[gSys]:
            if 'C' == band[0] or 'P' == band[0]:
                if band_choose_filter and gSys in band_choose_filter:
                    if band[1:] not in band_choose_filter[gSys]:
                        continue
                if band in bandChoose[gSys][list(bandChoose[gSys])[0]]:
                    band2 = bandChoose[gSys][list(bandChoose[gSys])[1]][0]
                else:
                    band2 = bandChoose[gSys][list(bandChoose[gSys])[0]][0]

                bandL = 'L' + band[1:]
                band2L = 'L' + band2[1:]
                if bandL not in obsType[gSys] or band2 not in obsType[gSys] or band2L not in obsType[gSys]:
                    continue

                band1freq = getBandFreq(gSys, band)
                band2freq = getBandFreq(gSys, band2)

                if band1freq is None or band2freq is None:
                    continue

                band1freq = band1freq * 1.e6
                band1lamda = CLIGHT / band1freq

                band2freq = band2freq * 1.e6
                band2lamda = CLIGHT / band2freq

                highBandChoose[gSys][band] = {}
                highBandChoose[gSys][band]['band1C'] = band
                highBandChoose[gSys][band]['band2C'] = band2
                highBandChoose[gSys][band]['band1L'] = bandL
                highBandChoose[gSys][band]['band2L'] = band2L
                highBandChoose[gSys][band]['band1freq'] = band1freq
                highBandChoose[gSys][band]['band1lambda'] = band1lamda
                highBandChoose[gSys][band]['band2freq'] = band2freq
                highBandChoose[gSys][band]['band2lambda'] = band2lamda

                highBandChoose[gSys][band]['COEFL1'] = -(band1freq * band1freq + band2freq*band2freq)/(band1freq*band1freq - band2freq*band2freq)
                highBandChoose[gSys][band]['COEFL2'] = 2*band2freq*band2freq/(band1freq*band1freq - band2freq*band2freq)

                lamdaW = CLIGHT / (band1freq - band2freq)
                highBandChoose[gSys][band]['lamdaW'] = lamdaW

    sorted_epochs = sorted(obsData.keys())
    if len(sorted_epochs) < 2:
        interval = 30.0
    else:
        interval = (sorted_epochs[1] - sorted_epochs[0]).total_seconds()

    highData = {}
    for prn in sat_list:
        gSys = prn[0]
        if gSys not in obsFreq or gSys not in highBandChoose:
            continue

        for band in highBandChoose[gSys]:
            bandEpochList = []
            highList = []
            highGoodList = []
            epochGoodList = []
            band1L = highBandChoose[gSys][band]['band1L']
            band2L = highBandChoose[gSys][band]['band2L']
            band1C = highBandChoose[gSys][band]['band1C']
            band2C = highBandChoose[gSys][band]['band2C']
            freq1 = highBandChoose[gSys][band]['band1freq']
            freq2 = highBandChoose[gSys][band]['band2freq']
            lamda1 = highBandChoose[gSys][band]['band1lambda']
            lamda2 = highBandChoose[gSys][band]['band2lambda']
            lamdaW = highBandChoose[gSys][band]['lamdaW']
            COEFL1 = highBandChoose[gSys][band]['COEFL1']
            COEFL2 = highBandChoose[gSys][band]['COEFL2']

            for epoch in sorted_epochs:
                if prn not in obsData[epoch].keys():
                    continue
                if start_datetime is not None and epoch < start_datetime:
                    continue
                if end_datetime is not None and epoch > end_datetime:
                    continue
                L1 = obsData[epoch][prn].get(band1L)
                L2 = obsData[epoch][prn].get(band2L)
                P1 = obsData[epoch][prn].get(band1C)
                P2 = obsData[epoch][prn].get(band2C)
                if L1 is None or L2 is None or P1 is None or P2 is None:
                    continue
                if abs(P2 - P1) > 100:
                    continue
                nowL = L1 * lamda1

                if bandEpochList == []:
                    bandEpochList.append(epoch)
                    highList = [[nowL]]
                    continue
                if (epoch - bandEpochList[-1]).total_seconds() > interval:
                    highList.append([nowL])
                    bandEpochList.append(epoch)
                    continue
                highList[-1].append(nowL)
                bandEpochList.append(epoch)
                if len(highList[-1]) > diffNum * 2:
                    difIndex = (diffNum + 1)*(-1)
                    nowDiff = culHighOrderDiff(highList[-1][difIndex:])
                    if diffNum == 4:
                        nowDiff = nowDiff / np.sqrt(70)
                    else:
                        nowDiff = nowDiff / np.sqrt(20)
                    if abs(nowDiff) > 100:
                        continue
                    highGoodList.append(nowDiff)
                    epochGoodList.append(epoch)
            if len(highGoodList) == 0:
                continue
            if gSys not in highData:
                highData[gSys] = {}
            if prn not in highData[gSys]:
                highData[gSys][prn] = {}
            highData[gSys][prn][band1L] = {}
            highData[gSys][prn][band1L]['epoch'] = epochGoodList
            highData[gSys][prn][band1L]['highOrder'] = highGoodList

    return highData

#写入汇总和详细CSV（已包含）
def writeHod(highData, highFile):
    import math
    import datetime
    summary_dict = {}
    highFileWrite = open(highFile, 'w+')
    highFileWrite.write('# Inf       : Fourth-order diff\n')
    highFileWrite.write('# Time      : ' + datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S') + '\n')
    highFileWrite.write('#             END OF HEADER\n')
    highFileWrite.write('\n')
    highLines = ['\n']
    highFileWrite.write('+SAT\n')
    highFileWrite.write(' prnFreq    RMS[m]    MIN[m]    MAX[m]\n')
    highFileWrite.write(38*'-' + '\n')
    for gnssSys in highData:
        for prn in highData[gnssSys]:
            for band in highData[gnssSys][prn]:
                if len(highData[gnssSys][prn][band]['highOrder']) == 0:
                    continue
                prnFreqLine = ' ' + prn + '_' + band
                square_sum = sum(x ** 2 for x in highData[gnssSys][prn][band]['highOrder'])
                rms = math.sqrt(square_sum / len(highData[gnssSys][prn][band]['highOrder']))
                prnFreqLine += '%10.6f' % rms + '%10.6f' % min(highData[gnssSys][prn][band]['highOrder']) + '%10.6f' % max(highData[gnssSys][prn][band]['highOrder']) + '\n'
                highFileWrite.write(prnFreqLine)
                highLines.append('+' + prn + '_' + band + '\n')
                highLines.append('               epoch   diff(m)\n')
                for epoch, highOrder in zip(highData[gnssSys][prn][band]['epoch'], highData[gnssSys][prn][band]['highOrder']):
                    highLines.append(' ' + epoch.strftime('%Y-%m-%d %H:%M:%S') + '%10.6f' % highOrder + '\n')
                highLines.append('-' + prn + '_' + band + '\n')
                highLines.append('\n')
                summary_dict.setdefault('highOrder', {}).setdefault(prn, {})[band] = {'val': rms}
    highFileWrite.write(38*'-' + '\n')
    highFileWrite.write('-SAT\n\n')
    highFileWrite.writelines(highLines)
    return summary_dict

def writeHodDetailed(highData, csvFile):
    """
    将高阶差分详细序列输出为宽表 CSV 格式：行=历元，列=prn_band，单元格=highOrder 值（缺失留空）
    """
    import csv
    all_prn_bands = set()
    all_epochs = set()
    for gnssSys in highData:
        for prn in highData[gnssSys]:
            for band, data in highData[gnssSys][prn].items():
                if not data['highOrder']:
                    continue
                prn_band = f"{prn}_{band}"
                all_prn_bands.add(prn_band)
                for ep in data['epoch']:
                    all_epochs.add(ep)

    if not all_prn_bands or not all_epochs:
        print("警告：没有有效的高阶差分数据，未生成 CSV 文件。")
        return

    sorted_epochs = sorted(all_epochs)
    sorted_prn_bands = sorted(all_prn_bands)

    epoch_data = {ep: {} for ep in sorted_epochs}
    for gnssSys in highData:
        for prn in highData[gnssSys]:
            for band, data in highData[gnssSys][prn].items():
                if not data['highOrder']:
                    continue
                prn_band = f"{prn}_{band}"
                for ep, val in zip(data['epoch'], data['highOrder']):
                    epoch_data[ep][prn_band] = val

    with open(csvFile, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['epoch'] + sorted_prn_bands)
        for ep in sorted_epochs:
            row = [ep.strftime('%Y-%m-%d %H:%M:%S')]
            for prn_band in sorted_prn_bands:
                val = epoch_data[ep].get(prn_band, '')
                row.append(f"{val:.6f}" if val != '' else '')
            writer.writerow(row)