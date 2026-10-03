from com.gnssParameter import getBandFreq, CLIGHT, obsFreq
import numpy as np
import datetime
"""
几何无关（GF）组合(双频载波相位差分)的三阶时间差分来估计 GNSS 载波相位噪声。
GF组合消除了几何距离和钟差误差的影响，保留了电离层延迟和噪声，通过取第三个差分消除了电离层趋势项，提取了载波相位噪声。
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
对每个卫星‑波段组合，在连续无周跳的弧段内，对 GF 组合值 L_GF = λ1·L1 - λ2·L2 滑动计算三阶差分，
并除以缩放因子 sqrt(40)，得到相位噪声估计值（单位：米），用于评估载波相位的平滑度。可视为载波相位观测值的等效噪声水平
phaseNoiseData = 
{
    'G': {
        'G01': {
            'C1C': {                       # 伪距波段代码
                'epoch': [datetime, ...],    # 三阶差分对应的历元（窗口最后一个历元）
                'phaseNoise': [float, ...]   # 三阶差分 / sqrt(40)，单位米
            },
            ...
        },
        ...
    }
}
"""
def phaseNoise(obsHead, obsData, start_datetime=None, end_datetime=None,
               sat_list=None, sys_list=None, band_choose=None):
    """
    Analyzes phase noise using the third-order time difference of the Geometry-Free (GF) phase combination.

    This function processes GNSS observation data to calculate phase noise by analyzing the third-order
    difference of the GF phase combination. The GF combination helps isolate ionospheric effects, and
    the third-order difference provides insights into high-frequency noise components. The function
    supports multiple GNSS systems (e.g., GPS, BDS).

    Args:
        obsHead (dict): Observation header containing metadata (e.g., OBS TYPES, prn).
        obsData (dict): Observation data dictionary with epochs as keys.
        start_datetime (datetime, optional): Start time for processing (inclusive). Defaults to None.
        end_datetime (datetime, optional): End time for processing (inclusive). Defaults to None.
        sat_list (list, optional): List of PRNs to process. If None, use all from obsHead['prn'].
        sys_list (list, optional): List of GNSS system characters (e.g., ['G','C','E']). If None, use all systems.
        band_choose (dict, optional): Dictionary specifying which bands to use per system,
                                      e.g., {'G': ['1C','2W']}. Defaults to None.

    Returns:
        dict: Phase noise data dictionary containing the third-order differences for each satellite and band.
              Structure: {sys: {prn: {band: {'epoch': list, 'phaseNoise': list}}}}
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
    iodBandChoose = {}

    # Determine satellite list
    if sat_list is None:
        sat_list = obsHead['prn']
    # Apply system filter
    if sys_list is not None:
        sat_list = [prn for prn in sat_list if prn[0] in sys_list]
    # Band filter
    band_choose_filter = band_choose if band_choose is not None else {}

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
                # print('band1freq:', band1freq, 'band2freq:', band2freq)  # optional
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

                iodBandChoose[gSys][band]['COEFL1'] = -(band1freq * band1freq + band2freq*band2freq)/(band1freq*band1freq - band2freq*band2freq)
                iodBandChoose[gSys][band]['COEFL2'] = 2*band2freq*band2freq/(band1freq*band1freq - band2freq*band2freq)

                lamdaW = CLIGHT / (band1freq - band2freq)
                iodBandChoose[gSys][band]['lamdaW'] = lamdaW

    # Sort epochs and compute sampling interval
    sorted_epochs = sorted(obsData.keys())
    if len(sorted_epochs) < 2:
        interval = 30.0
    else:
        interval = (sorted_epochs[1] - sorted_epochs[0]).total_seconds()

    phaseNoiseData = {}
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
                nowGF = L_GF

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
                    iodList = [[nowGF]]
                    continue
                if (epoch - bandEpochList[-1]).total_seconds() > interval:
                    bandEpochList.append(epoch)
                    N_MW_list.append([N_MW])
                    MEAN_MW_list.append([N_MW])
                    SIGMA_list = [[0.3]]
                    L_GF_list.append([L_GF])
                    P_GF_list.append([P_GF])
                    slipList.append(0)
                    iodList.append([nowGF])
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
                    iodList.append([nowGF])
                else:
                    bandEpochList.append(epoch)
                    N_MW_list[-1].append(N_MW)
                    MEAN_MW_list[-1].append(N_MW_MEAN)
                    L_GF_list[-1].append(L_GF)
                    SIGMA_list[-1].append(SIGMA)
                    P_GF_list[-1].append(P_GF)
                    iodList[-1].append(nowGF)
                    slipList.append(0)
                    if len(iodList[-1]) > 3:
                        epochGoodList.append(epoch)
                        nowDiff = culHighOrderDiff(iodList[-1][-4:]) / np.sqrt(40)
                        iodGoodList.append(nowDiff)

            if gSys not in phaseNoiseData:
                phaseNoiseData[gSys] = {}
            if prn not in phaseNoiseData[gSys]:
                phaseNoiseData[gSys][prn] = {}
            phaseNoiseData[gSys][prn][band] = {}
            phaseNoiseData[gSys][prn][band]['epoch'] = epochGoodList
            phaseNoiseData[gSys][prn][band]['phaseNoise'] = iodGoodList

    return phaseNoiseData


"""
汇总表：每个卫星‑波段的 RMS、最小值、最大值（单位：米），精度 6 位小数。
详细序列：每个卫星‑波段的所有好历元及其相位噪声值（精度 6 位小数）。
"""
def writePhaseNoise(phaseNoiseData, phaseNoiseFile):
    """
    仅输出汇总表：每个卫星‑波段的 RMS、最小值、最大值（单位：米），精度 6 位小数。
    """
    import math
    summary_dict = {}
    with open(phaseNoiseFile, 'w+') as f:
        f.write('# Inf       : PhaseNoise\n')
        f.write('# Time      : ' + datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S') + '\n')
        f.write('#             END OF HEADER\n\n')
        f.write('+SAT\n')
        f.write(' prnFreq    RMS[m]    MIN[m]    MAX[m]\n')
        f.write(38*'-' + '\n')
        for gnssSys in phaseNoiseData:
            for prn in phaseNoiseData[gnssSys]:
                for band, data in phaseNoiseData[gnssSys][prn].items():
                    if not data['phaseNoise']:
                        continue
                    square_sum = sum(x**2 for x in data['phaseNoise'])
                    rms = math.sqrt(square_sum / len(data['phaseNoise']))
                    line = f" {prn}_{band} {rms:10.6f} {min(data['phaseNoise']):10.6f} {max(data['phaseNoise']):10.6f}\n"
                    f.write(line)
                    summary_dict.setdefault('phaseNoise', {}).setdefault(prn, {})[band] = {'val': rms}
        f.write(38*'-' + '\n')
        f.write('-SAT\n')
    return summary_dict

def writePhaseNoiseDetailed(phaseNoiseData, csvFile):
    """
    将详细序列输出为宽表 CSV 格式：行=历元，列=prn_band，单元格=phaseNoise 值（缺失留空）
    """
    import csv
    # 1. 收集所有出现过的 prn_band 和所有历元
    all_prn_bands = set()
    all_epochs = set()
    for gnssSys in phaseNoiseData:
        for prn in phaseNoiseData[gnssSys]:
            for band, data in phaseNoiseData[gnssSys][prn].items():
                if not data['phaseNoise']:
                    continue
                prn_band = f"{prn}_{band}"
                all_prn_bands.add(prn_band)
                for ep in data['epoch']:
                    all_epochs.add(ep)

    if not all_prn_bands or not all_epochs:
        print("警告：没有有效的相位噪声数据，未生成 CSV 文件。")
        return

    # 2. 排序
    sorted_epochs = sorted(all_epochs)
    sorted_prn_bands = sorted(all_prn_bands)

    # 3. 构建每个历元对应的 {prn_band: noise} 字典
    epoch_data = {ep: {} for ep in sorted_epochs}
    for gnssSys in phaseNoiseData:
        for prn in phaseNoiseData[gnssSys]:
            for band, data in phaseNoiseData[gnssSys][prn].items():
                if not data['phaseNoise']:
                    continue
                prn_band = f"{prn}_{band}"
                for ep, noise in zip(data['epoch'], data['phaseNoise']):
                    epoch_data[ep][prn_band] = noise

    # 4. 写入 CSV
    with open(csvFile, 'w', newline='') as f:
        writer = csv.writer(f)
        # 表头：epoch, 然后是所有 prn_band
        writer.writerow(['epoch'] + sorted_prn_bands)
        # 数据行
        for ep in sorted_epochs:
            row = [ep.strftime('%Y-%m-%d %H:%M:%S')]
            for prn_band in sorted_prn_bands:
                val = epoch_data[ep].get(prn_band, '')
                # 保留 6 位小数，但确保空值留空
                if val != '':
                    row.append(f"{val:.6f}")
                else:
                    row.append('')
            writer.writerow(row)

