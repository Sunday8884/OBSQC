# -*- coding: utf-8 -*-
import datetime

"""
将 GNSS 观测数据中每个卫星各频率的出现次数统计写入文本文件，同时按系统汇总各频率的总出现次数。
freqData = {
    'G01': {'L1C': 120, 'L2W': 115, 'L5X': 80},
    'G02': {'L1C': 118, 'L2W': 112, 'L5X': 75},
    'C06': {'B1I': 100, 'B3I': 98},
    ...
}
"""

def writeFreqData(freqData, freqFile, total_epochs = 2880):
    """
    写入频段统计文件，同时返回每个卫星各频段的出现频率（计数 / total_epochs）。

    Parameters:
        freqData (dict): 频段计数数据 {prn: {freq: count}}
        freqFile (str): 输出文件路径
        total_epochs (int): 总历元数，默认为 2880（一天 30 秒采样）

    Returns:
        dict: 统计字典 {'freq': {prn: {freq: {'val': frequency}}}}
    """
    freqFileWrite = open(freqFile, 'w+')
    freqFileWrite.write('# Inf       : Frequency Count Statistics\n')
    freqFileWrite.write('# Time      : ' + datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S') + '\n')
    freqFileWrite.write('#             END OF HEADER\n')
    freqSysData = {}
    freqFileWrite.write('\n')
    satLines = []
    satLines.append('+SAT\n')

    summary_dict = {}

    for prn in freqData:
        gSys = prn[0]
        if gSys not in freqSysData:
            freqSysData[gSys] = {}
        line = ' ' + prn + ' -'
        for freq in freqData[prn]:
            count = freqData[prn][freq]
            if freq not in freqSysData[gSys]:
                freqSysData[gSys][freq] = 0
            freqSysData[gSys][freq] += count
            line += ' ' + freq + ': ' + '%7d' % count
            # 出现频率 = count / total_epochs
            freq_val = count / total_epochs
            summary_dict.setdefault('freq', {}).setdefault(prn, {})[freq] = {'val': freq_val}
        line += '\n'
        satLines.append(line)
    satLines.append('-SAT\n')

    freqFileWrite.write('+SYS\n')
    for gSys in freqSysData:
        line = ' ' + gSys + '   -'
        for freq in freqSysData[gSys]:
            line += ' ' + freq + ': ' + '%7d' % freqSysData[gSys][freq]
        line += '\n'
        freqFileWrite.write(line)
    freqFileWrite.write('-SYS\n')
    freqFileWrite.write('\n')
    freqFileWrite.writelines(satLines)
    freqFileWrite.close()

    return summary_dict