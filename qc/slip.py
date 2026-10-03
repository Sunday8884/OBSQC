"""
输入：obsHead, obsData
    │
    ▼
确定卫星列表 & 时间范围
    │
    ▼
选择最佳双频组合（优先 L1/L2，或 L1/L5 等，由 bandComb 定义）
    │
    ▼
对每颗卫星：
    │
    ▼
逐历元处理（按时间顺序）
    │
    ├─ 若该历元无此卫星或观测缺失 → 跳过
    │
    ├─ 若时间间隔 > 采样间隔 → 开始新段
    │
    ├─ 计算 N_MW、L_GF、P_GF
    │
    ├─ MW 检测：
    │   └─ 若 |N_MW - mean| ≥ 4*sigma → 周跳
    │
    ├─ GF 检测（若 MW 未检测到且段长度>1）：
    │   └─ 若 |ΔL_GF| > 0.15 m → 周跳
    │
    ├─ 若检测到周跳 → 结束当前段，开始新段，标记 1
    │
    └─ 否则 → 更新段内统计量，标记 0
    │
    ▼
返回 slipData
slipData = {
    'G': {
        'G01': {
            'slip': [0, 0, 1, ...],   # 周跳标记列表
            'epochs': [datetime, ...]  # 对应的历元列表
        },
        ...
    },
    ...
}
start_datetime: 开始时间（datetime对象），只处理此时间及之后的历元。
end_datetime: 结束时间（datetime对象），只处理此时间及之前的历元。
sat_list: 卫星列表（字符串列表，如['G01','G02']），只处理列表中的卫星。
sys_list: 系统列表（字符列表，如['G','C','E']），只处理这些GNSS系统的卫星。
"""
from com.gnssParameter import getBandFreq, CLIGHT, bandComb
import numpy as np
import time
import datetime

def turboedit(obsHead, obsData, start_datetime = None, end_datetime = None, sat_list = None, sys_list = None):
    """
    This subroutine calculates the Cycle Slips Ratio (CSR) for each satellite from GNSS
    observation data, using Turboedit method to detect cycle slips.

    Parameters:
    ----------
    obsHead :
        Observation File Header in Python Dictionary Format

    obsData :
        Observation File Data in Python Dictionary Format

    start_datetime : datetime, optional
        Start time for processing (inclusive). If None, process from the beginning.

    end_datetime : datetime, optional
        End time for processing (inclusive). If None, process to the end.

    sat_list : list, optional
        List of PRN strings to process (e.g., ['G01','G02']). If None, use all satellites in obsHead.

    sys_list : list, optional
        List of GNSS system characters (e.g., ['G','C','E']). If None, use all systems.

    Returns:
    ----------
    slipData :
        slip DATA in Python Dictionary Format {sys: {prn: [0/1, ...]}}

    Notes
    ----------
        Modified for Python by Chuntao Chang
    """

    # Determine satellite list
    if sat_list is not None:
        satList = sat_list
    else:
        satList = obsHead['prn']

    # Apply system filter if provided
    if sys_list is not None:
        satList = [prn for prn in satList if prn[0] in sys_list]

    obsType = obsHead['OBS TYPES']
    bandChoose = {}
    for gSys in obsType:
        if gSys not in bandChoose:
            bandChoose[gSys] = {}
        for band in obsType[gSys]:
            if 'L' == band[0]:
                if band[1] not in bandChoose[gSys]:
                    bandChoose[gSys][band[1]] = []
                bandC = 'C' + band[1:]
                if bandC in obsType[gSys]:
                    bandChoose[gSys][band[1]].append(band)
                elif 'P' + band[1:] in obsType[gSys]:
                    bandChoose[gSys][band[1]].append(band)
    cycleSlipBandChoose = {}
    for gSys in bandComb:
        if gSys in bandChoose:
            for comb in bandComb[gSys]:
                band1 = str(comb[0])
                band2 = str(comb[1])
                if band1 in bandChoose[gSys] and band2 in bandChoose[gSys]:
                    if len(bandChoose[gSys][band1]) == 0 or len(bandChoose[gSys][band2]) == 0:
                        continue
                    if bandChoose[gSys][band1][0] is None or bandChoose[gSys][band2][0] is None:
                        continue
                    freq1 = getBandFreq(gSys, bandChoose[gSys][band1][0])
                    freq2 = getBandFreq(gSys, bandChoose[gSys][band2][0])
                    if freq1 is None or freq2 is None:
                        continue
                    cycleSlipBandChoose[gSys] = {}
                    freq1 = freq1 * 1.e6
                    freq2 = freq2 * 1.e6
                    cycleSlipBandChoose[gSys][bandChoose[gSys][band1][0]] = {}
                    cycleSlipBandChoose[gSys][bandChoose[gSys][band2][0]] = {}
                    cycleSlipBandChoose[gSys][bandChoose[gSys][band1][0]]['freq'] = freq1
                    cycleSlipBandChoose[gSys][bandChoose[gSys][band2][0]]['freq'] = freq2
                    cycleSlipBandChoose[gSys][bandChoose[gSys][band1][0]]['lambda'] = CLIGHT / freq1
                    cycleSlipBandChoose[gSys][bandChoose[gSys][band2][0]]['lambda'] = CLIGHT / freq2
                    cycleSlipBandChoose[gSys]['lambdaW'] = CLIGHT / (freq1 - freq2)
                    break

    slipData = {}
    # Get sampling interval from first two epochs
    epoch_list = sorted(obsData.keys())
    if len(epoch_list) < 2:
        interval = 30.0
    else:
        interval = (epoch_list[1] - epoch_list[0]).total_seconds()

    for prn in satList:
        slipList = []
        epochs_this_prn = []
        prnEpochList = []  # 用于分段
        # 分段内统计变量
        prn_N_MW_list = []
        prn_M_MW_list = []
        SIGMA_list = []
        L_GF_list = []
        P_GF_list = []

        for epoch in epoch_list:
            if prn not in obsData[epoch]:
                continue
            gSys = prn[0]
            if gSys not in cycleSlipBandChoose:
                continue
            if start_datetime and epoch < start_datetime:
                continue
            if end_datetime and epoch > end_datetime:
                continue

            band1L = list(cycleSlipBandChoose[gSys])[0]
            band2L = list(cycleSlipBandChoose[gSys])[1]
            band1C = 'C' + band1L[1:]
            if band1C not in obsData[epoch][prn]:
                band1C = 'P' + band1L[1:]
            band2C = 'C' + band2L[1:]
            if band2C not in obsData[epoch][prn]:
                band2C = 'P' + band2L[1:]
            L1 = obsData[epoch][prn].get(band1L)
            L2 = obsData[epoch][prn].get(band2L)
            P1 = obsData[epoch][prn].get(band1C)
            P2 = obsData[epoch][prn].get(band2C)
            if L1 is None or L2 is None or P1 is None or P2 is None:
                continue

            # 记录有效历元
            epochs_this_prn.append(epoch)

            freq1 = cycleSlipBandChoose[gSys][band1L]['freq']
            freq2 = cycleSlipBandChoose[gSys][band2L]['freq']
            lambdaW = cycleSlipBandChoose[gSys]['lambdaW']
            lambda1 = cycleSlipBandChoose[gSys][band1L]['lambda']
            lambda2 = cycleSlipBandChoose[gSys][band2L]['lambda']
            N_MW = L1 - L2 - (freq1 * P1 + freq2 * P2) / (freq1 + freq2) / lambdaW
            L_GF = lambda1 * L1 - lambda2 * L2

            MW_CLIP = False
            GF_CLIP = False
            if prnEpochList == []:
                # 开始新段
                prnEpochList.append(epoch)
                prn_N_MW_list.append([N_MW])
                prn_M_MW_list.append([N_MW])
                SIGMA_list.append([0.3])
                L_GF_list.append([L_GF])
                P_GF_list.append([P2 - P1])
                slipList.append(0)
                continue
            if (epoch - prnEpochList[-1]).total_seconds() > interval:
                prnEpochList.append(epoch)
                prn_N_MW_list.append([N_MW])
                prn_M_MW_list.append([N_MW])
                SIGMA_list.append([0.3])
                L_GF_list.append([L_GF])
                P_GF_list.append([P2 - P1])
                slipList.append(0)
                continue

            if len(prn_N_MW_list[-1]) < 2:
                SIGMA_I_BACK = 0.15
            else:
                SIGMA_I_BACK = SIGMA_list[-1][-1]

            N_MW_MEAN_BACK = prn_M_MW_list[-1][-1]
            N_MW_MEAN = N_MW_MEAN_BACK + (N_MW - N_MW_MEAN_BACK) / (len(prn_M_MW_list[-1]) + 1)
            SIGMA2 = SIGMA_I_BACK * SIGMA_I_BACK + ((N_MW - N_MW_MEAN_BACK) ** 2 - SIGMA_I_BACK * SIGMA_I_BACK) / (
                        len(prn_M_MW_list[-1]) + 1)
            SIGMA = np.sqrt(SIGMA2)
            if abs(N_MW - N_MW_MEAN) >= 4 * SIGMA_I_BACK:
                MW_CLIP = True

            if not MW_CLIP and len(L_GF_list[-1]) > 1:
                if abs(L_GF - L_GF_list[-1][-1]) > 0.15:
                    GF_CLIP = True

            if MW_CLIP or GF_CLIP:
                prnEpochList.append(epoch)
                prn_N_MW_list.append([N_MW])
                prn_M_MW_list.append([N_MW])
                SIGMA_list.append([0.15])
                L_GF_list.append([L_GF])
                P_GF_list.append([P2 - P1])
                slipList.append(1)
            else:
                prnEpochList.append(epoch)
                prn_N_MW_list[-1].append(N_MW)
                prn_M_MW_list[-1].append(N_MW_MEAN)
                L_GF_list[-1].append(L_GF)
                SIGMA_list[-1].append(SIGMA)
                P_GF_list[-1].append(P2 - P1)
                slipList.append(0)

        if slipList:
            gnssSys = prn[0]
            if gnssSys not in slipData:
                slipData[gnssSys] = {}
            slipData[gnssSys][prn] = {
                'slip': slipList,
                'epochs': epochs_this_prn
            }

    return slipData

def writeSlip(slipData, slipFile):
    summary_dict = {}
    slipFileWrite = open(slipFile, 'w+')
    slipFileWrite.write('# Inf       : Slip\n')
    slipFileWrite.write('# Time      : ' + datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S') + '\n')
    slipFileWrite.write('#             END OF HEADER\n')
    slipFileWrite.write('\n')
    satLines = ['+SAT\n']
    satLines.append('  PRN  ObsNum   Slip   CSR\n')   # CSR = 1000 / (O / Slip)
    satLines.append('  ---------------------------\n')
    for gnssSys in slipData:
        for prn, data in slipData[gnssSys].items():
            slip_list = data['slip']
            if sum(slip_list) == 0:
                CSR = '%9.2f' % 0
                # CSR = '   NoSlip'
                summary_val = 0       # <-- NoSlip 时在汇总中记为 0
            else:
                CSR = '%9.2f' %  ( 1000 / ( len(slip_list) / sum(slip_list) ) ) #  CSR = 1000/(总历元数 / 周跳次数)
                summary_val = CSR if CSR != float('inf') else float('inf')
            line = str(prn).rjust(5) + '%8d' % len(slip_list) + '%7d' % sum(slip_list) + CSR + '\n'
            satLines.append(line)
            summary_dict.setdefault('CSR', {}).setdefault(prn, {})['CSR'] = {'val': summary_val}
    satLines.append('  ---------------------------\n')
    satLines.append('-SAT\n')
    slipFileWrite.writelines(satLines)
    slipFileWrite.close()
    return summary_dict

def writeSlipDetailed(slipData, slipFile, all_epochs=None, start_datetime=None, end_datetime=None):
    import csv
    # 收集所有卫星 PRN
    all_prns = []
    for gnssSys in slipData:
        for prn in slipData[gnssSys]:
            all_prns.append(prn)
    all_prns.sort()

    # 确定所有历元（从所有卫星的 epochs 中收集）
    if all_epochs is None:
        all_epochs_set = set()
        for gnssSys in slipData:
            for prn, data in slipData[gnssSys].items():
                all_epochs_set.update(data['epochs'])
        all_epochs = sorted(all_epochs_set)
    else:
        all_epochs = sorted(all_epochs)
        if start_datetime:
            all_epochs = [ep for ep in all_epochs if ep >= start_datetime]
        if end_datetime:
            all_epochs = [ep for ep in all_epochs if ep <= end_datetime]

    # 构建映射：每颗卫星的 {历元: 周跳标记}
    sat_epoch_flag = {}
    for gnssSys in slipData:
        for prn, data in slipData[gnssSys].items():
            epochs_list = data['epochs']
            slip_list = data['slip']
            if len(epochs_list) != len(slip_list):
                print(f"警告：卫星 {prn} 的历元数与周跳标记数不一致，将按较短者截断。")
                min_len = min(len(epochs_list), len(slip_list))
                epochs_list = epochs_list[:min_len]
                slip_list = slip_list[:min_len]
            sat_epoch_flag[prn] = {ep: flag for ep, flag in zip(epochs_list, slip_list)}

    # 写入 CSV
    with open(slipFile, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['epoch'] + all_prns)
        for ep in all_epochs:
            row = [ep.strftime('%Y-%m-%d %H:%M:%S')]
            for prn in all_prns:
                flag = sat_epoch_flag.get(prn, {}).get(ep, '')
                row.append(str(flag) if flag != '' else '')
            writer.writerow(row)
