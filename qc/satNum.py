# -*- coding: utf-8 -*-
import datetime
"""
该模块用于将 GNSS 观测数据中每个历元的卫星数量统计写入文本文件，并提供按 GNSS 系统的统计摘要（平均值、最小值、最大值）。
satNumData = {
    'stat': {
        'G': {'AVE': 9.2, 'MIN': 4, 'MAX': 12},
        'C': {'AVE': 8.5, 'MIN': 3, 'MAX': 11},
        ...
    },
    'data': {
        datetime(2023,1,1,0,0,0): {'G': 9, 'C': 8, ...},
        datetime(2023,1,1,0,30,0): {'G': 10, 'C': 7, ...},
        ...
    }
}
"""
def writeSatNum(satNumData, satNumFile):
    # SatNum[prn][satNum]
    satNumFileWrite = open(satNumFile, 'w+')
    satNumFileWrite.write('# Inf       : satNum Count Statistics\n')
    satNumFileWrite.write('# Time      : ' + datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S') + '\n')
    satNumFileWrite.write('#             END OF HEADER\n')
    satNumFileWrite.write('\n')

    satLines = ['+SAT\n']

    firstLine = 'epoch'.rjust(20)
    firstLine += ' '
    for gSys in satNumData['stat']:
        firstLine += str(gSys).rjust(4)

    satLines.append(firstLine + '\n')

    for epoch in satNumData['data']:
        line = ' ' + epoch.strftime('%Y-%m-%d %H:%M:%S') + ' '
        for gSys in satNumData['stat']:
            line += '%4d' % satNumData['data'][epoch][gSys]
        line += '\n'
        satLines.append(line)
    satLines.append('-SAT\n')

    satNumFileWrite.write('+SYS\n')
    satNumFileWrite.write(' gSys  AVE  MIN  MAX\n')
    satNumFileWrite.write('--------------------\n')
    for gSys in satNumData['stat']:
        line = str(gSys).rjust(5) + '%5.1f' % satNumData['stat'][gSys]['AVE'] + '%5d' % satNumData['stat'][gSys][
            'MIN'] + '%5d' % satNumData['stat'][gSys]['MAX'] + '\n'
        satNumFileWrite.write(line)

    satNumFileWrite.write('--------------------\n')
    satNumFileWrite.write('-SYS\n')

    satNumFileWrite.write('\n')
    satNumFileWrite.writelines(satLines)
