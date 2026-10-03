# -*- coding: utf-8 -*-
import datetime
"""
实现将载噪比（Carrier‑to‑Noise Ratio, CNR, 单位 dBHz）数据写入文本文件的功能。
writeCnr 函数会统计每个系统、每个卫星‑波段组合的平均值、最小值和最大值，并以固定格式输出。
cnrData =
{
    'G': {                     # GNSS 系统字符（G,C,E,...）
        'G01': {               # 卫星 PRN
            'L1C': {           # 波段代码（载波相位观测类型）
                datetime1: 45.2,   # 历元 → CNR 值（dBHz）
                datetime2: 46.1,
                ...
            },
            'L2W': {...}
        },
        ...
    },
    'C': {...},
    ...
}
"""
# Cnr统计.txt
def writeCnr(cnrData, cnrFile):
    import math
    import datetime
    summary_dict = {}
    cnrFileWrite = open(cnrFile, 'w+')
    cnrFileWrite.write('# Inf       : Carrier-to-Noise Ratio [dBHz]\n')
    cnrFileWrite.write('# Time      : ' + datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S') + '\n')
    cnrFileWrite.write('#             END OF HEADER\n')
    cnrFileWrite.write('\n')
    cnrSysData = {}
    cnrSatData = {}
    for gnssSys in cnrData:
        cnrSysData[gnssSys] = {}
        cnrSysData[gnssSys]['AVE'] = 0.0
        cnrSysData[gnssSys]['MIN'] = 999.0
        cnrSysData[gnssSys]['MAX'] = 0.0
        bandSysSum = 0.0
        bandSysNum = 0.0

        for prn in cnrData[gnssSys]:
            if prn not in cnrSatData:
                cnrSatData[prn] = {}
            for band in cnrData[gnssSys][prn]:
                if band not in cnrSatData[prn]:
                    cnrSatData[prn][band] = {}
                cnrSatData[prn][band]['AVE'] = 0.0
                cnrSatData[prn][band]['MIN'] = 999.0
                cnrSatData[prn][band]['MAX'] = 0.0
                bandSum = 0.0
                bandNum = 0.0
                for epoch, nowcnr in cnrData[gnssSys][prn][band].items():
                    if nowcnr > cnrSatData[prn][band]['MAX']:
                        cnrSatData[prn][band]['MAX'] = nowcnr
                    if nowcnr < cnrSatData[prn][band]['MIN']:
                        cnrSatData[prn][band]['MIN'] = nowcnr
                    if nowcnr > cnrSysData[gnssSys]['MAX']:
                        cnrSysData[gnssSys]['MAX'] = nowcnr
                    if nowcnr < cnrSysData[gnssSys]['MIN']:
                        cnrSysData[gnssSys]['MIN'] = nowcnr
                    bandSum += nowcnr
                    bandNum += 1
                    bandSysSum += nowcnr
                    bandSysNum += 1
                if bandNum == 0:
                    continue
                cnrSatData[prn][band]['AVE'] = bandSum / bandNum
        if bandSysNum != 0:
            cnrSysData[gnssSys]['AVE'] = bandSysSum / bandSysNum

    cnrFileWrite.write('+SYS\n')
    cnrFileWrite.write('    gSys   AVE   MIN   MAX\n')
    cnrFileWrite.write('--------------------------\n')
    for gnssSys in cnrSysData:
        line = str(gnssSys).rjust(8) + '%6.1f' % cnrSysData[gnssSys]['AVE'] + '%6.1f' % cnrSysData[gnssSys][
            'MIN'] + '%6.1f' % cnrSysData[gnssSys]['MAX'] + '\n'
        cnrFileWrite.write(line)
    cnrFileWrite.write('--------------------------\n')
    cnrFileWrite.write('+SYS\n')
    cnrFileWrite.write('\n')

    cnrFileWrite.write('+SAT\n')
    cnrFileWrite.write(' prnFreq   AVE   MIN   MAX\n')
    cnrFileWrite.write('--------------------------\n')
    for prn in cnrSatData:
        for band in cnrSatData[prn]:
            prn_band = prn + '_' + band
            line = prn_band.rjust(8) + '%6.1f' % cnrSatData[prn][band]['AVE'] + '%6.1f' % cnrSatData[prn][band][
                'MIN'] + '%6.1f' % cnrSatData[prn][band]['MAX'] + '\n'
            cnrFileWrite.write(line)
            summary_dict.setdefault('CNR', {}).setdefault(prn, {})[band] = {'val': cnrSatData[prn][band]['AVE']}
    cnrFileWrite.write('--------------------------\n')
    cnrFileWrite.write('-SAT\n')
    cnrFileWrite.close()
    return summary_dict
# Cnr原始数据.csv
def writeCnrDetailed(cnrData, detailFile):
    """
    将原始 CNR 数据输出为宽表格式：
    - 行：每个历元
    - 列：prn_band（卫星+频段组合，如 G01_L1C）
    - 单元格：CNR 值 (dBHz)，缺失值留空（空字符串）
    """
    # 1. 收集所有出现过的 prn_band 和所有历元
    all_prn_bands = set()
    all_epochs = set()

    for gnssSys in cnrData:
        for prn in cnrData[gnssSys]:
            for band in cnrData[gnssSys][prn]:
                prn_band = f"{prn}_{band}"
                all_prn_bands.add(prn_band)
                for epoch in cnrData[gnssSys][prn][band]:
                    all_epochs.add(epoch)

    # 2. 排序（历元升序，prn_band 按自然顺序排序）
    sorted_epochs = sorted(all_epochs)
    sorted_prn_bands = sorted(all_prn_bands)

    # 3. 为每个历元构建一个字典 {prn_band: cnr}
    epoch_data = {epoch: {} for epoch in sorted_epochs}
    for gnssSys in cnrData:
        for prn in cnrData[gnssSys]:
            for band in cnrData[gnssSys][prn]:
                prn_band = f"{prn}_{band}"
                for epoch, cnr_val in cnrData[gnssSys][prn][band].items():
                    epoch_data[epoch][prn_band] = cnr_val

    # 4. 写入 CSV
    with open(detailFile, 'w+', newline='') as f:
        # 写入表头
        header = ['epoch'] + sorted_prn_bands
        f.write(','.join(header) + '\n')
        # 写入数据行
        for epoch in sorted_epochs:
            epoch_str = epoch.strftime('%Y-%m-%d %H:%M:%S') if hasattr(epoch, 'strftime') else str(epoch)
            row = [epoch_str]
            for prn_band in sorted_prn_bands:
                val = epoch_data[epoch].get(prn_band, '')
                row.append(str(val) if val != '' else '')
            f.write(','.join(row) + '\n')