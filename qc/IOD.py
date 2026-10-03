from com.gnssParameter import getBandFreq, CLIGHT, obsFreq
import numpy as np
import datetime
"""
电离层延迟变化率（Ionospheric Delay Rate, IOD）
通过双频载波相位的几何无关（GF）组合得到电离层延迟，再对其做时间差分，得到电离层延迟的变化速率（单位：m/s）。
该指标可用于评估电离层活跃程度、探测电离层扰动，也可作为数据质量检验的辅助参数。
iodData = {
    'G': {
        'G01': {
            'C1C': {               # 伪距波段代码
                'epoch': [datetime, ...],   # 历元列表
                'iod':   [float, ...]       # 电离层延迟变化率 (m/s)
            },
            ...
        },
        ...
    },
    ...
}
"""
def IOD(obsHead, obsData, start_datetime=None, end_datetime=None,
        sat_list=None, sys_list=None, band_choose=None):
    """
    Computes the Ionospheric Delay Rate (IOD) for GNSS observations.

    Args:
        obsHead (dict): Observation header.
        obsData (dict): Observation data.
        start_datetime, end_datetime (datetime, optional): Time range.
        sat_list (list, optional): List of PRNs.
        sys_list (list, optional): List of GNSS system characters (e.g., ['G','C','E']).
        band_choose (dict, optional): Band selection.

    Returns:
        dict: IOD data {sys: {prn: {band: {'epoch': list, 'iod': list}}}}
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

    # Determine satellite list with system filter
    if sat_list is None:
        sat_list = obsHead['prn']
    if sys_list is not None:
        sat_list = [prn for prn in sat_list if prn[0] in sys_list]

    # Band filter
    band_choose_filter = band_choose if band_choose is not None else {}

    iodBandChoose = {}
    for gSys in obsType:
        # Skip if system not in sys_list
        if sys_list is not None and gSys not in sys_list:
            continue
        if len(bandChoose[gSys]) < 2 or gSys not in obsFreq:
            continue
        if band_choose_filter and gSys not in band_choose_filter:
            continue
        if gSys not in iodBandChoose:
            iodBandChoose[gSys] = {}
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

                iodBandChoose[gSys][band] = {}
                iodBandChoose[gSys][band]['band1freq'] = band1freq
                iodBandChoose[gSys][band]['band1lambda'] = band1lamda
                iodBandChoose[gSys][band]['band2freq'] = band2freq
                iodBandChoose[gSys][band]['band2lambda'] = band2lamda
                iodBandChoose[gSys][band]['band1C'] = band
                iodBandChoose[gSys][band]['band2C'] = band2
                iodBandChoose[gSys][band]['band1L'] = bandL
                iodBandChoose[gSys][band]['band2L'] = band2L

                iodBandChoose[gSys][band]['COEFL1'] = -(band1freq*band1freq + band2freq*band2freq) / (band1freq*band1freq - band2freq*band2freq)
                iodBandChoose[gSys][band]['COEFL2'] = 2*band2freq*band2freq / (band1freq*band1freq - band2freq*band2freq)

                lamdaW = CLIGHT / (band1freq - band2freq)
                iodBandChoose[gSys][band]['lamdaW'] = lamdaW

    # Sort epochs and compute sampling interval
    sorted_epochs = sorted(obsData.keys())
    if len(sorted_epochs) < 2:
        interval = 30.0
    else:
        interval = (sorted_epochs[1] - sorted_epochs[0]).total_seconds()

    iodData = {}
    for prn in sat_list:
        gSys = prn[0]
        if gSys not in obsFreq or gSys not in iodBandChoose:
            continue

        for band in iodBandChoose[gSys]:
            bandEpochList = []
            iodList = []
            epochGoodList = []
            iodGoodList = []
            band1L = iodBandChoose[gSys][band]['band1L']
            band2L = iodBandChoose[gSys][band]['band2L']
            band1C = iodBandChoose[gSys][band]['band1C']
            band2C = iodBandChoose[gSys][band]['band2C']
            freq1 = iodBandChoose[gSys][band]['band1freq']
            freq2 = iodBandChoose[gSys][band]['band2freq']
            lamda1 = iodBandChoose[gSys][band]['band1lambda']
            lamda2 = iodBandChoose[gSys][band]['band2lambda']
            lamdaW = iodBandChoose[gSys][band]['lamdaW']
            COEFL2 = iodBandChoose[gSys][band]['COEFL2']

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
                nowIOD = COEFL2 / 2 * L_GF

                MW_CLIP = False
                GF_CLIP = False
                if not bandEpochList:
                    bandEpochList = [epoch]
                    N_MW_list = [[N_MW]]
                    MEAN_MW_list = [[N_MW]]
                    SIGMA_list = [[0.3]]
                    L_GF_list = [[L_GF]]
                    P_GF_list = [[P2 - P1]]
                    iodList = [[nowIOD]]
                    continue
                if (epoch - bandEpochList[-1]).total_seconds() > interval:
                    bandEpochList.append(epoch)
                    N_MW_list.append([N_MW])
                    MEAN_MW_list.append([N_MW])
                    SIGMA_list = [[0.3]]
                    L_GF_list.append([L_GF])
                    P_GF_list.append([P2 - P1])
                    iodList.append([nowIOD])
                    continue

                if len(N_MW_list[-1]) < 2:
                    SIGMA_I_BACK = 0.15
                else:
                    SIGMA_I_BACK = SIGMA_list[-1][-1]

                N_MW_MEAN_BACK = MEAN_MW_list[-1][-1]
                N_MW_MEAN = N_MW_MEAN_BACK + (N_MW - N_MW_MEAN_BACK) / (len(MEAN_MW_list[-1]) + 1)
                SIGMA2 = SIGMA_I_BACK * SIGMA_I_BACK + ((N_MW - N_MW_MEAN_BACK) ** 2 - SIGMA_I_BACK * SIGMA_I_BACK) / (len(MEAN_MW_list[-1]) + 1)
                SIGMA = np.sqrt(SIGMA2)
                if abs(N_MW - N_MW_MEAN) >= 4 * SIGMA_I_BACK:
                    MW_CLIP = True

                if not MW_CLIP and len(L_GF_list[-1]) > 1:
                    if abs(L_GF - L_GF_list[-1][-1]) > 0.15:
                        GF_CLIP = True

                if MW_CLIP or GF_CLIP:
                    bandEpochList.append(epoch)
                    N_MW_list.append([N_MW])
                    MEAN_MW_list.append([N_MW])
                    SIGMA_list.append([0.15])
                    L_GF_list.append([L_GF])
                    P_GF_list.append([P2 - P1])
                    iodList.append([nowIOD])
                else:
                    epochGoodList.append(epoch)
                    iodForFreq1 = (nowIOD - iodList[-1][-1]) / interval
                    iodGoodList.append(iodForFreq1)
                    bandEpochList.append(epoch)
                    N_MW_list[-1].append(N_MW)
                    MEAN_MW_list[-1].append(N_MW_MEAN)
                    L_GF_list[-1].append(L_GF)
                    SIGMA_list[-1].append(SIGMA)
                    P_GF_list[-1].append(P2 - P1)
                    iodList[-1].append(nowIOD)

            if gSys not in iodData:
                iodData[gSys] = {}
            if prn not in iodData[gSys]:
                iodData[gSys][prn] = {}
            iodData[gSys][prn][band] = {
                'epoch': epochGoodList,
                'iod': iodGoodList
            }

    return iodData

def writeIon(ionData, ionFile):
    """输出汇总表：每个卫星‑波段的 RMS、最小值、最大值（单位：m/s），精度 6 位小数。"""
    import math
    summary_dict = {}
    with open(ionFile, 'w+') as f:
        f.write('# Inf       : Derivative of the ionospheric delay\n')
        f.write('# Time      : ' + datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S') + '\n')
        f.write('#             END OF HEADER\n\n')
        f.write('+SAT\n')
        f.write(' prnFreq    RMS[m]    MIN[m]    MAX[m]\n')
        f.write(38 * '-' + '\n')
        for gnssSys in ionData:
            for prn in ionData[gnssSys]:
                for band, data in ionData[gnssSys][prn].items():
                    if not data['iod']:
                        continue
                    square_sum = sum(x ** 2 for x in data['iod'])
                    rms = math.sqrt(square_sum / len(data['iod']))
                    line = f" {prn}_{band} {rms:10.6f} {min(data['iod']):10.6f} {max(data['iod']):10.6f}\n"
                    f.write(line)
                    summary_dict.setdefault('IOD', {}).setdefault(prn, {})[band] = {'val': rms}
        f.write(38 * '-' + '\n')
        f.write('-SAT\n')
    return summary_dict

def writeIonDetailed(ionData, csvFile):
    """
    将 IOD 详细序列输出为宽表 CSV 格式：行=历元，列=prn_band，单元格=iod 值（m/s，缺失留空）
    """
    import csv
    # 收集所有 prn_band 和所有历元
    all_prn_bands = set()
    all_epochs = set()
    for gnssSys in ionData:
        for prn in ionData[gnssSys]:
            for band, data in ionData[gnssSys][prn].items():
                if not data['iod']:
                    continue
                prn_band = f"{prn}_{band}"
                all_prn_bands.add(prn_band)
                for ep in data['epoch']:
                    all_epochs.add(ep)

    if not all_prn_bands or not all_epochs:
        print("警告：没有有效的 IOD 数据，未生成 CSV 文件。")
        return

    # 排序
    sorted_epochs = sorted(all_epochs)
    sorted_prn_bands = sorted(all_prn_bands)

    # 构建每个历元对应的 {prn_band: iod} 字典
    epoch_data = {ep: {} for ep in sorted_epochs}
    for gnssSys in ionData:
        for prn in ionData[gnssSys]:
            for band, data in ionData[gnssSys][prn].items():
                if not data['iod']:
                    continue
                prn_band = f"{prn}_{band}"
                for ep, iod_val in zip(data['epoch'], data['iod']):
                    epoch_data[ep][prn_band] = iod_val

    # 写入 CSV
    with open(csvFile, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['epoch'] + sorted_prn_bands)
        for ep in sorted_epochs:
            row = [ep.strftime('%Y-%m-%d %H:%M:%S')]
            for prn_band in sorted_prn_bands:
                val = epoch_data[ep].get(prn_band, '')
                if val != '':
                    row.append(f"{val:.6f}")   # 与汇总表精度一致，保留6位小数
                else:
                    row.append('')
            writer.writerow(row)
