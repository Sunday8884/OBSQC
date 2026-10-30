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
def countSatBySystem(obsData, systems=None, start_datetime=None, end_datetime=None):
    """
    逐历元统计各 GNSS 系统的可见卫星数。

    systems 为允许统计的系统列表（来自 --systems / ALL_GNSS_SYSTEMS）：
      - 列表中**具体系统**（G、C、E …）决定统计范围，范围外的系统不计入；
      - 'ALL' 只是合计标记，不参与范围判定；当列表里除 'ALL' 外没有
        具体系统时（如 --systems ALL），各系统列才统计所有出现的系统；
      - 'ALL' 列**始终统计所有 PRN**（不过滤系统），与 _sat.txt 中
        plotSatNum 的 ALL 口径一致。

    返回 {epoch(datetime): {系统: 数量}}。
    """
    allow = None
    if systems is not None:
        concrete = {s for s in systems if s and s != 'ALL'}
        # 只有 'ALL'（或空）时不限系统；否则按具体系统过滤
        allow = concrete if concrete else None

    out = {}
    for ep in sorted(obsData.keys()):
        if start_datetime and ep < start_datetime:
            continue
        if end_datetime and ep > end_datetime:
            continue
        counts = {}
        total = 0
        for prn in obsData[ep]:
            total += 1
            sys_char = prn[0]
            if allow is not None and sys_char not in allow:
                continue
            counts[sys_char] = counts.get(sys_char, 0) + 1
        if total:
            counts['ALL'] = total
            out[ep] = counts
    return out


def writeSatNumDetailed(satCount, detailFile):
    """
    将逐历元各系统卫星数写成宽表 CSV：
        行 = 历元，列 = 各系统（如 C / G / E ...）及 ALL（所有 PRN 合计）
    """
    if not satCount:
        return None

    systems = []
    for ep in satCount:
        for sys in satCount[ep]:
            if sys not in systems:
                systems.append(sys)
    # ALL 放在最后
    systems = [s for s in systems if s != 'ALL']
    if 'ALL' in satCount[next(iter(satCount))]:
        systems.append('ALL')

    with open(detailFile, 'w+', newline='') as f:
        f.write(','.join(['epoch'] + systems) + '\n')
        for ep in sorted(satCount.keys()):
            epoch_str = ep.strftime('%Y-%m-%d %H:%M:%S') if hasattr(ep, 'strftime') else str(ep)
            row = [epoch_str]
            for sys in systems:
                v = satCount[ep].get(sys, '')
                row.append(str(v) if v != '' else '')
            f.write(','.join(row) + '\n')
    return detailFile


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
