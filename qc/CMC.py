# -*- coding: utf-8 -*-
from com.gnssParameter import getBandFreq, CLIGHT, obsFreq
import numpy as np
import datetime
"""
历元间差分码减相位（Time‑differenced Code Minus Phase, CMC）
CMC是评估GNSS伪距观测噪声和多路径效应的重要指标，直接用伪距减去载波相位，最终项中剩下的是两倍的电离层误差、整周模糊度、伪距多路径、载波相位多路径、伪距噪声、载波相位噪声。
计算历元间差分的CMC序列，可以消除模糊度等常数项，保留伪距噪声与多路径的变化，反映了伪距观测的时间变化噪声。
cmcData = 
{
    'G': {
        'G01': {
            'C1C': {
                'epoch': [datetime, ...],
                'cmc':   [float, ...]   # 历元间差分的CMC (单位: m)
            },
            ...
        },
        ...
    },
    ...
}

SumData = 
{
    'CMC': {
        'prn': {
            'band': {
                'val': [datetime, ...],
            },
            ...
        },
        ...
    },
    ...
}
"""
def CMC(obsHead, obsData, start_datetime=None, end_datetime=None,
        sat_list=None, sys_list=None, band_choose=None):
    """
    Calculate time-differenced Code Minus Phase (CMC) for quality assessment.

    Parameters:
    ----------
    obsHead : dict
        Observation file header (contains 'OBS TYPES', 'prn', etc.)
    obsData : dict
        Observation data {epoch: {prn: {band: value}}}
    start_datetime : datetime, optional
        Start time for processing (inclusive)
    end_datetime : datetime, optional
        End time for processing (inclusive)
    sat_list : list of str, optional
        List of PRNs to process (e.g., ['G01','G02']). If None, use all from obsHead['prn'].
    sys_list : list of str, optional
        List of GNSS system characters (e.g., ['G','C','E']). If None, use all systems.
    band_choose : dict, optional
        Dictionary specifying which bands to use per system, e.g., {'G': ['1C','2W'], 'C': ['2I','6I']}.
        If None, use all available band combinations.

    Returns:
    -------
    cmcData : dict
        Structure: {sys: {prn: {band: {'epoch': list, 'cmc': list}}}}
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

    # Determine satellite list
    if sat_list is None:
        sat_list = obsHead['prn']
    # Apply system filter
    if sys_list is not None:
        sat_list = [prn for prn in sat_list if prn[0] in sys_list]

    # Build band selection for each system (if provided)
    band_choose_filter = band_choose if band_choose is not None else {}

    cmcBandChoose = {}
    for gSys in obsType:
        # Skip if system not in sys_list
        if sys_list is not None and gSys not in sys_list:
            continue
        if len(bandChoose[gSys]) < 2 or gSys not in obsFreq:
            continue
        if band_choose_filter and gSys not in band_choose_filter:
            continue
        if gSys not in cmcBandChoose:
            cmcBandChoose[gSys] = {}
        for band in obsType[gSys]:
            if 'C' == band[0] or 'P' == band[0]:
                # Apply band filtering if provided
                if band_choose_filter and gSys in band_choose_filter:
                    if band[1:] not in band_choose_filter[gSys]:
                        continue
                # Determine the companion band (the other frequency)
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

                cmcBandChoose[gSys][band] = {}
                cmcBandChoose[gSys][band]['band1freq'] = band1freq
                cmcBandChoose[gSys][band]['band1lambda'] = band1lamda
                cmcBandChoose[gSys][band]['band2freq'] = band2freq
                cmcBandChoose[gSys][band]['band2lambda'] = band2lamda
                cmcBandChoose[gSys][band]['band1C'] = band
                cmcBandChoose[gSys][band]['band2C'] = band2
                cmcBandChoose[gSys][band]['band1L'] = bandL
                cmcBandChoose[gSys][band]['band2L'] = band2L

                cmcBandChoose[gSys][band]['COEFL1'] = -(band1freq * band1freq + band2freq * band2freq) / (
                            band1freq * band1freq - band2freq * band2freq)
                cmcBandChoose[gSys][band]['COEFL2'] = 2 * band2freq * band2freq / (
                            band1freq * band1freq - band2freq * band2freq)

                lamdaW = CLIGHT / (band1freq - band2freq)
                cmcBandChoose[gSys][band]['lamdaW'] = lamdaW

    # Sort epochs for consistent interval calculation
    sorted_epochs = sorted(obsData.keys())
    if len(sorted_epochs) < 2:
        interval = 30.0  # default if not enough epochs
    else:
        interval = (sorted_epochs[1] - sorted_epochs[0]).total_seconds()

    cmcData = {}
    for prn in sat_list:
        gSys = prn[0]
        if gSys not in obsFreq or gSys not in cmcBandChoose:
            continue

        for band in cmcBandChoose[gSys]:
            bandEpochList = []
            cmcList = []
            epochGoodList = []
            cmcGoodList = []

            band1L = cmcBandChoose[gSys][band]['band1L']
            band2L = cmcBandChoose[gSys][band]['band2L']
            band1C = cmcBandChoose[gSys][band]['band1C']
            band2C = cmcBandChoose[gSys][band]['band2C']
            freq1 = cmcBandChoose[gSys][band]['band1freq']
            freq2 = cmcBandChoose[gSys][band]['band2freq']
            lamda1 = cmcBandChoose[gSys][band]['band1lambda']
            lamda2 = cmcBandChoose[gSys][band]['band2lambda']
            lamdaW = cmcBandChoose[gSys][band]['lamdaW']

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
                if abs(P2 - P1) > 30:
                    continue

                N_MW = L1 - L2 - (freq1 * P1 + freq2 * P2) / (freq1 + freq2) / lamdaW
                L_GF = lamda1 * L1 - lamda2 * L2
                P_GF = P2 - P1
                nowCmc = P1 - lamda1 * L1  # Code minus phase, time-differenced later

                MW_CLIP = False
                GF_CLIP = False

                if not bandEpochList:
                    # Start a new segment
                    bandEpochList.append(epoch)
                    N_MW_list = [[N_MW]]
                    MEAN_MW_list = [[N_MW]]
                    SIGMA_list = [[0.3]]
                    L_GF_list = [[L_GF]]
                    P_GF_list = [[P_GF]]
                    slipList = [0]
                    cmcList = [[nowCmc]]
                    continue

                # Check for data gap
                if (epoch - bandEpochList[-1]).total_seconds() > interval:
                    bandEpochList.append(epoch)
                    N_MW_list.append([N_MW])
                    MEAN_MW_list.append([N_MW])
                    SIGMA_list.append([0.3])
                    L_GF_list.append([L_GF])
                    P_GF_list.append([P_GF])
                    slipList.append(0)
                    cmcList.append([nowCmc])
                    continue

                # Update MW statistics
                if len(N_MW_list[-1]) < 2:
                    SIGMA_I_BACK = 0.15
                else:
                    SIGMA_I_BACK = SIGMA_list[-1][-1]

                N_MW_MEAN_BACK = MEAN_MW_list[-1][-1]
                N_MW_MEAN = N_MW_MEAN_BACK + (N_MW - N_MW_MEAN_BACK) / (len(MEAN_MW_list[-1]) + 1)
                SIGMA2 = SIGMA_I_BACK * SIGMA_I_BACK + ((N_MW - N_MW_MEAN_BACK) ** 2 - SIGMA_I_BACK * SIGMA_I_BACK) / (
                            len(MEAN_MW_list[-1]) + 1)
                SIGMA = np.sqrt(SIGMA2)
                if abs(N_MW - N_MW_MEAN) >= 4 * SIGMA_I_BACK:
                    MW_CLIP = True

                # GF check
                if not MW_CLIP and len(L_GF_list[-1]) > 1:
                    if abs(L_GF - L_GF_list[-1][-1]) > 0.15:
                        GF_CLIP = True

                # Reset flags for CMC decision (note: original code sets them False again)
                MW_CLIP = False
                GF_CLIP = False
                deltCMC = (nowCmc - cmcList[-1][-1]) / np.sqrt(2)

                if MW_CLIP or GF_CLIP or abs(deltCMC) > 1000:
                    # Cycle slip or large jump, start new segment
                    bandEpochList.append(epoch)
                    N_MW_list.append([N_MW])
                    MEAN_MW_list.append([N_MW])
                    SIGMA_list.append([0.15])
                    L_GF_list.append([L_GF])
                    P_GF_list.append([P_GF])
                    slipList.append(1)
                    cmcList.append([nowCmc])
                else:
                    # Good epoch, keep in current segment
                    epochGoodList.append(epoch)
                    cmcGoodList.append(deltCMC)
                    bandEpochList.append(epoch)
                    N_MW_list[-1].append(N_MW)
                    MEAN_MW_list[-1].append(N_MW_MEAN)
                    L_GF_list[-1].append(L_GF)
                    SIGMA_list[-1].append(SIGMA)
                    P_GF_list[-1].append(P_GF)
                    cmcList[-1].append(nowCmc)
                    slipList.append(0)

            # Store results for this satellite and band
            if gSys not in cmcData:
                cmcData[gSys] = {}
            if prn not in cmcData[gSys]:
                cmcData[gSys][prn] = {}
            cmcData[gSys][prn][band] = {
                'epoch': epochGoodList,
                'cmc': cmcGoodList
            }

    return cmcData

def writeCmc(cmcData, cmcFile):
    """仅输出汇总表：每个卫星‑波段的 RMS、最小值、最大值（单位：米），精度 4 位小数。同时返回统计字典。"""
    import math
    import datetime
    summary_dict = {}
    with open(cmcFile, 'w+') as f:
        f.write('# Inf       : Time-differenced code minus phase\n')
        f.write('# Time      : ' + datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S') + '\n')
        f.write('#             END OF HEADER\n\n')
        f.write('+SAT\n')
        f.write(' prnFreq    RMS[m]    MIN[m]    MAX[m]\n')
        f.write(38 * '-' + '\n')
        for gnssSys in cmcData:
            for prn in cmcData[gnssSys]:
                for band, data in cmcData[gnssSys][prn].items():
                    vals = data['cmc']
                    if not vals:
                        continue
                    square_sum = sum(x ** 2 for x in vals)
                    rms = math.sqrt(square_sum / len(vals))
                    min_val = min(vals)
                    max_val = max(vals)
                    # 写入文件
                    f.write(f" {prn}_{band} {rms:10.4f} {min_val:10.4f} {max_val:10.4f}\n")
                    summary_dict.setdefault('CMC', {}).setdefault(prn, {})[band] = {'val': rms}
        f.write(38 * '-' + '\n')
        f.write('-SAT\n')
    return summary_dict   # 返回统计列表

def writeCmcDetailed(cmcData, csvFile):
    """
    将 CMC 详细序列输出为宽表 CSV 格式：行=历元，列=prn_band，单元格=cmc 值（缺失留空）
    """
    import csv
    # 收集所有 prn_band 和所有历元
    all_prn_bands = set()
    all_epochs = set()
    for gnssSys in cmcData:
        for prn in cmcData[gnssSys]:
            for band, data in cmcData[gnssSys][prn].items():
                if not data['cmc']:
                    continue
                prn_band = f"{prn}_{band}"
                all_prn_bands.add(prn_band)
                for ep in data['epoch']:
                    all_epochs.add(ep)

    if not all_prn_bands or not all_epochs:
        print("警告：没有有效的 CMC 数据，未生成 CSV 文件。")
        return

    # 排序
    sorted_epochs = sorted(all_epochs)
    sorted_prn_bands = sorted(all_prn_bands)

    # 构建每个历元对应的 {prn_band: cmc} 字典
    epoch_data = {ep: {} for ep in sorted_epochs}
    for gnssSys in cmcData:
        for prn in cmcData[gnssSys]:
            for band, data in cmcData[gnssSys][prn].items():
                if not data['cmc']:
                    continue
                prn_band = f"{prn}_{band}"
                for ep, cmc_val in zip(data['epoch'], data['cmc']):
                    epoch_data[ep][prn_band] = cmc_val

    # 写入 CSV
    with open(csvFile, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['epoch'] + sorted_prn_bands)
        for ep in sorted_epochs:
            row = [ep.strftime('%Y-%m-%d %H:%M:%S')]
            for prn_band in sorted_prn_bands:
                val = epoch_data[ep].get(prn_band, '')
                if val != '':
                    row.append(f"{val:.4f}")   # 与汇总表精度一致，保留4位小数
                else:
                    row.append('')
            writer.writerow(row)
