from com.gnssParameter import getBandFreq, CLIGHT, obsFreq
import numpy as np
import datetime
"""
多径效应（Multipath）通过特定的伪距‑载波相位组合来提取。
采用经典的双频多路径组合（MP组合）公式（伪距P1、P2和载波相位L1、L2进行线性组合，消除几何距离、钟差、电离层一阶项，保留伪距多径和噪声）。
同时集成了 Turboedit 周跳检测方法，用于划分连续弧段，并对每个弧段内的多径序列减去其均值，以消除模糊度等常数项。
mpData = 
{
    'G': {
        'G01': {
            'C1C': {                        # 伪距波段代码
                'epoch': [datetime, ...],   # 所有历元
                'mp':    [float, ...]       # 多路径值序列（已去除非整段均值）
            },
            ...
        },
        ...
    }
}
"""
def multipath(obsHead, obsData, start_datetime=None, end_datetime=None,
              sat_list=None, sys_list=None, band_choose=None):
    """
    This subroutine calculates the multipath for each satellite from GNSS
    observation data, using Turboedit method to detect cycle slips.

    Parameters:
    ----------
    obsHead :
        Observation File Header in Python Dictionary Format
    obsData :
        Observation File Data in Python Dictionary Format
    start_datetime : datetime, optional
        Start time for processing (inclusive)
    end_datetime : datetime, optional
        End time for processing (inclusive)
    sat_list : list, optional
        List of PRNs to process (e.g., ['G01','G02']). If None, use all from obsHead['prn'].
    sys_list : list, optional
        List of GNSS system characters (e.g., ['G','C','E']). If None, use all systems.
    band_choose : dict, optional
        Dictionary specifying which bands to use per system, e.g., {'G': ['1C','2W']}

    Returns:
    ----------
    mpData :
        mp DATA in Python Dictionary Format {sys: {prn: {band: {'epoch': list, 'mp': list}}}}

    Notes
    ----------
        Modified for Python by Chuntao Chang
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
    mpBandChoose = {}

    # Determine satellite list
    if sat_list is None:
        sat_list = obsHead['prn']
    # Apply system filter to satellite list
    if sys_list is not None:
        sat_list = [prn for prn in sat_list if prn[0] in sys_list]
    # Band filter
    band_choose_filter = band_choose if band_choose is not None else {}

    for gSys in obsType:
        # Skip if system not in sys_list (if sys_list is provided)
        if sys_list is not None and gSys not in sys_list:
            continue
        if len(bandChoose[gSys]) < 2 or gSys not in obsFreq:
            continue
        if band_choose_filter and gSys not in band_choose_filter:
            continue
        if gSys not in mpBandChoose:
            mpBandChoose[gSys] = {}
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
                if band1freq is None:
                    continue
                band1freq = band1freq * 1.e6
                band2freq = getBandFreq(gSys, band2)
                if band2freq is None:
                    continue
                band2freq = band2freq * 1.e6
                mpBandChoose[gSys][band] = {}
                mpBandChoose[gSys][band]['band1C'] = band
                mpBandChoose[gSys][band]['band2C'] = band2
                mpBandChoose[gSys][band]['band1L'] = bandL
                mpBandChoose[gSys][band]['band2L'] = band2L
                # print(gSys, band)  # optional debug print
                band1lamda = CLIGHT / band1freq
                mpBandChoose[gSys][band]['band1freq'] = band1freq
                mpBandChoose[gSys][band]['band1lambda'] = band1lamda

                band2lamda = CLIGHT / band2freq
                mpBandChoose[gSys][band]['band2freq'] = band2freq
                mpBandChoose[gSys][band]['band2lambda'] = band2lamda

                mpBandChoose[gSys][band]['COEFL1'] = -(band1freq * band1freq + band2freq*band2freq)/(band1freq*band1freq - band2freq*band2freq)
                mpBandChoose[gSys][band]['COEFL2'] = 2*band2freq*band2freq/(band1freq*band1freq - band2freq*band2freq)

                lamdaW = CLIGHT / (band1freq - band2freq)
                mpBandChoose[gSys][band]['lamdaW'] = lamdaW

    # Compute sampling interval from sorted epochs
    sorted_epochs = sorted(obsData.keys())
    if len(sorted_epochs) < 2:
        interval = 30.0  # fallback
    else:
        interval = (sorted_epochs[1] - sorted_epochs[0]).total_seconds()

    mpData = {}
    for prn in sat_list:
        gSys = prn[0]
        if gSys not in obsFreq or gSys not in mpBandChoose:
            continue

        for band in mpBandChoose[gSys]:
            bandEpochList = []
            mpList = []
            band1L = mpBandChoose[gSys][band]['band1L']
            band2L = mpBandChoose[gSys][band]['band2L']
            band1C = mpBandChoose[gSys][band]['band1C']
            band2C = mpBandChoose[gSys][band]['band2C']
            freq1 = mpBandChoose[gSys][band]['band1freq']
            freq2 = mpBandChoose[gSys][band]['band2freq']
            lamda1 = mpBandChoose[gSys][band]['band1lambda']
            lamda2 = mpBandChoose[gSys][band]['band2lambda']
            lamdaW = mpBandChoose[gSys][band]['lamdaW']
            COEFL1 = mpBandChoose[gSys][band]['COEFL1']
            COEFL2 = mpBandChoose[gSys][band]['COEFL2']

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
                # if abs(P2 - P1) > 100:   # optional, currently commented in original
                #     continue

                N_MW = L1 - L2 - (freq1 * P1 + freq2 * P2) / (freq1 + freq2) / lamdaW
                L_GF = lamda1 * L1 - lamda2 * L2
                P_GF = P2 - P1
                nowMP = P1 + COEFL1 * CLIGHT / freq1 * L1 + COEFL2 * CLIGHT/freq2 * L2

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
                    mpList = [[nowMP]]
                    continue
                if (epoch - bandEpochList[-1]).total_seconds() > interval:
                    bandEpochList.append(epoch)
                    N_MW_list.append([N_MW])
                    MEAN_MW_list.append([N_MW])
                    SIGMA_list = [[0.3]]
                    L_GF_list.append([L_GF])
                    P_GF_list.append([P_GF])
                    slipList.append(0)
                    meanMP_N = sum(mpList[-1]) / len(mpList[-1])
                    mpList[-1] = [element - meanMP_N for element in mpList[-1]]
                    mpList.append([nowMP])
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
                    meanMP_N = sum(mpList[-1]) / len(mpList[-1])
                    mpList[-1] = [element - meanMP_N for element in mpList[-1]]
                    mpList.append([nowMP])
                else:
                    bandEpochList.append(epoch)
                    N_MW_list[-1].append(N_MW)
                    MEAN_MW_list[-1].append(N_MW_MEAN)
                    L_GF_list[-1].append(L_GF)
                    SIGMA_list[-1].append(SIGMA)
                    P_GF_list[-1].append(P_GF)
                    slipList.append(0)
                    mpList[-1].append(nowMP)

            # After loop, subtract mean from last segment if any
            if mpList:
                meanMP_N = sum(mpList[-1]) / len(mpList[-1])
                mpList[-1] = [element - meanMP_N for element in mpList[-1]]

            if gSys not in mpData:
                mpData[gSys] = {}
            if prn not in mpData[gSys]:
                mpData[gSys][prn] = {}
            mpData[gSys][prn][band] = {}
            mpData[gSys][prn][band]['epoch'] = bandEpochList
            mpData[gSys][prn][band]['mp'] = [item for sublist in mpList for item in sublist]

    return mpData

"""
汇总表：每个卫星‑波段的 RMS、最小值、最大值（单位：米），精度 4 位小数。
"""
def writeMp(mpData, mpFile):
    import math
    summary_dict = {}
    with open(mpFile, 'w+') as f:
        f.write('# Inf       : Multipath\n')
        f.write('# Time      : ' + datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S') + '\n')
        f.write('#             END OF HEADER\n\n')
        f.write('+SAT\n')
        f.write(' prnFreq    RMS[m]    MIN[m]    MAX[m]\n')
        f.write(38*'-' + '\n')
        for gnssSys in mpData:
            for prn in mpData[gnssSys]:
                for band, data in mpData[gnssSys][prn].items():
                    if not data['mp']:
                        continue
                    square_sum = sum(x**2 for x in data['mp'])
                    rms = math.sqrt(square_sum / len(data['mp']))
                    line = f" {prn}_{band} {rms:10.4f} {min(data['mp']):10.4f} {max(data['mp']):10.4f}\n"
                    f.write(line)
                    summary_dict.setdefault('multipath', {}).setdefault(prn, {})[band] = {'val':rms}
        f.write(38*'-' + '\n')
        f.write('-SAT\n')
    return summary_dict


def writeMpDetailed(mpData, csvFile):
    """
    将多路径详细序列输出为宽表 CSV 格式：行=历元，列=prn_band，单元格=mp 值（缺失留空）
    """
    import csv
    # 收集所有 prn_band 和所有历元
    all_prn_bands = set()
    all_epochs = set()
    for gnssSys in mpData:
        for prn in mpData[gnssSys]:
            for band, data in mpData[gnssSys][prn].items():
                if not data['mp']:
                    continue
                prn_band = f"{prn}_{band}"
                all_prn_bands.add(prn_band)
                for ep in data['epoch']:
                    all_epochs.add(ep)

    if not all_prn_bands or not all_epochs:
        print("警告：没有有效的多路径数据，未生成 CSV 文件。")
        return

    # 排序
    sorted_epochs = sorted(all_epochs)
    sorted_prn_bands = sorted(all_prn_bands)

    # 构建每个历元对应的 {prn_band: mp} 字典
    epoch_data = {ep: {} for ep in sorted_epochs}
    for gnssSys in mpData:
        for prn in mpData[gnssSys]:
            for band, data in mpData[gnssSys][prn].items():
                if not data['mp']:
                    continue
                prn_band = f"{prn}_{band}"
                for ep, mp_val in zip(data['epoch'], data['mp']):
                    epoch_data[ep][prn_band] = mp_val

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