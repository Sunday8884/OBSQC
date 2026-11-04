# -*- coding: utf-8 -*-
"""
卫星数据完整率可视化
=====================

两类图：
1. 完整率柱状图（按系统分面）——  直观看出哪些卫星/频段存在数据缺口
2. 缺测历元时间线图（按系统分面）—— 定位缺口发生在什么时段
"""
import matplotlib.pyplot as plt
import matplotlib.dates as mdate
import numpy as np
import datetime


def plotCompleteness(compData, pngFile=None, gapData=None):
    """
    绘制卫星数据完整率柱状图。

    Parameters
    ----------
    compData : dict
        {gSys: {prn: {band: {'ratio_arc': float, 'ratio_span': float,
                              'obs': int, 'arc': int, 'gap_count': int,
                              'max_gap': int, ...}}}}
    pngFile : str, optional
        图片保存路径；None 则直接显示
    gapData : dict, optional
        缺口明细 {gSys: {prn: {band: [{'start','end','n','sec'}, ...]}}}，
        保留参数以兼容调用方；柱状图不再在柱顶标注缺口数

    Returns
    -------
    fig : matplotlib Figure
    """
    if not compData:
        print("警告：无完整率数据，无法绘图。")
        return None

    sys_list = sorted(compData.keys())
    n_sys = len(sys_list)

    fig, axes = plt.subplots(n_sys, 1, figsize=(13, max(3.5, 2.6 * n_sys)),
                             squeeze=False)
    axes = axes.ravel()

    bar_colors = {'G': '#4878CF', 'C': '#EE6666', 'E': '#6ACC65',
                  'R': '#D5BB67', 'J': '#B07AA1', 'I': '#58A8C8',
                  'S': '#C0703A'}

    for si, gSys in enumerate(sys_list):
        ax = axes[si]
        sats = sorted(compData[gSys].keys())
        if not sats:
            continue

        # 收集该系统所有 (prn, band) 标签与数值
        labels = []
        ratios = []
        for prn in sats:
            for band in sorted(compData[gSys][prn]):
                r = compData[gSys][prn][band]
                labels.append(f'{prn}_{band}')
                ratios.append(r['ratio_arc'] * 100)

        x = np.arange(len(labels))
        color = bar_colors.get(gSys, '#999999')
        # 低于 90% 的条标红，便于快速定位问题卫星
        colors = [color if v >= 90 else '#D93A3A' for v in ratios]

        ax.bar(x, ratios, color=colors, width=0.72, zorder=3)

        # 90% 参考线
        ax.axhline(90, color='#555555', linestyle='--', linewidth=0.9, zorder=4)

        ax.set_ylim(0, 105)
        ax.set_ylabel(f'{gSys} Completeness [%]', fontsize='medium')
        ax.grid(axis='y', zorder=0, alpha=0.3)
        ax.tick_params(axis='x', labelsize='x-small', rotation=90)
        ax.tick_params(axis='y', labelsize='small')

        # 频率过高时抽稀标签
        step = max(1, len(labels) // 60)
        ax.set_xticks(x[::step])
        ax.set_xticklabels(labels[::step])

        n_bad = sum(1 for v in ratios if v < 90)
        ax.set_title(f'{gSys} system: {len(labels)} sat-band pairs, '
                     f'{n_bad} below 90%   (dashed line = 90%)',
                     fontsize='small', loc='left')

    fig.subplots_adjust(left=0.08, right=0.99, bottom=0.06, top=0.97, hspace=0.45)
    if pngFile:
        plt.savefig(pngFile, dpi=200, bbox_inches='tight')
    else:
        plt.show()
    plt.close(fig)
    return fig


def plotCompletenessTimeline(gapData, pngFile=None):
    """
    绘制各卫星数据缺口的时间线图。

    每个缺口段画一条横线，长度为该缺口持续时长，红色。
    纵轴为卫星（PRN_BAND）。

    Parameters
    ----------
    gapData : dict
        {gSys: {prn: {band: [{'start','end','n','sec'}, ...]}}}
    pngFile : str, optional
        图片保存路径
    """
    if not gapData:
        print("警告：无缺口数据，无法绘图。")
        return None

    sys_list = sorted(gapData.keys())
    n_sys = len(sys_list)

    fig, axes = plt.subplots(n_sys, 1, figsize=(13, max(3.5, 2.2 * n_sys)),
                             squeeze=False)
    axes = axes.ravel()

    for si, gSys in enumerate(sys_list):
        ax = axes[si]
        entries = []   # (label, start, end, n)
        for prn in sorted(gapData[gSys]):
            for band in sorted(gapData[gSys][prn]):
                for seg in gapData[gSys][prn][band]:
                    entries.append((f'{prn}_{band}', seg['start'], seg['end'],
                                    seg['n']))
        if not entries:
            ax.set_visible(False)
            continue

        entries.sort(key=lambda e: e[0])
        labels = [e[0] for e in entries]
        y = np.arange(len(entries))

        for i, (lab, s, e, n) in enumerate(entries):
            ax.plot([s, e], [i, i], color='#D93A3A', linewidth=2.2,
                    solid_capstyle='butt', zorder=3)
            ax.text(e, i, f'  {n}', va='center', ha='left',
                    fontsize='x-small', color='#D93A3A', zorder=4)

        ax.set_yticks(y)
        ax.set_yticklabels(labels, fontsize='x-small')
        ax.set_ylim(-0.8, len(entries) - 0.2)
        ax.grid(axis='x', zorder=0, alpha=0.3)
        ax.tick_params(axis='y', labelsize='x-small')
        xfmt = mdate.DateFormatter('%dD-%H:%M')
        ax.xaxis.set_major_formatter(xfmt)
        ax.tick_params(axis='x', labelsize='small')
        ax.set_title(f'{gSys} system: {len(entries)} data gaps',
                     fontsize='small', loc='left')

    fig.subplots_adjust(left=0.10, right=0.97, bottom=0.06, top=0.97, hspace=0.45)
    if pngFile:
        plt.savefig(pngFile, dpi=200, bbox_inches='tight')
    else:
        plt.show()
    plt.close(fig)
    return fig
