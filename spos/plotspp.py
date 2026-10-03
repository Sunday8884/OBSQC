#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
SPP 绘图模块：将 spp.py 中的绘图逻辑独立出来，便于复用与维护。
"""
import numpy as np
import matplotlib.dates as mdate
import matplotlib.pyplot as plt

from com.xyz2neu import xyz2neu
from com.pub import rms


def plot_spp(posData, XMEAN, YMEAN, ZMEAN, plot_show=False, plot_save=None):
    """
    绘制 SPP 结果（E/N/U 误差 + DOP）。

    参数
    ----
    posData : dict
        {epoch(datetime): {'x','y','z','clk','HDOP','VDOP','PDOP','GDOP'}}
    XMEAN, YMEAN, ZMEAN : float
        参考坐标（用于将 XYZ 转成 E/N/U）。
    plot_show : bool
        是否显示绘图窗口。
    plot_save : str or None
        保存路径（png）。为 None 时不保存。
    """
    epochs_list = list(posData.keys())
    if not epochs_list:
        print("无有效数据，无法绘图。")
        return

    Elist, NList, UList = [], [], []
    HDOPList, VDOPList, PDOPList, GDOPList = [], [], [], []
    for epoch in epochs_list:
        N, E, U = xyz2neu(XMEAN, YMEAN, ZMEAN,
                          posData[epoch]['x'], posData[epoch]['y'], posData[epoch]['z'],
                          'WGS84')
        NList.append(N)
        Elist.append(E)
        UList.append(U)
        HDOPList.append(posData[epoch]['HDOP'])
        VDOPList.append(posData[epoch]['VDOP'])
        PDOPList.append(posData[epoch]['PDOP'])
        GDOPList.append(posData[epoch]['GDOP'])

    fig = plt.figure(figsize=(12, 10))
    ax_e = fig.add_subplot(4, 1, 1)
    ax_n = fig.add_subplot(4, 1, 2)
    ax_u = fig.add_subplot(4, 1, 3)
    ax_dop = fig.add_subplot(4, 1, 4)

    for ax in [ax_e, ax_n, ax_u, ax_dop]:
        ax.grid(True, zorder=0)

    ax_e.scatter(epochs_list, Elist, marker='+', s=5, zorder=10, color='blue')
    ax_e.plot(epochs_list, Elist, color='gray', alpha=0.7, zorder=10)
    ax_e.set_ylabel('E [m]', fontsize=16)
    ax_e.text(0.02, 0.92, f'RMS: {rms(Elist):.3f} m',
              transform=ax_e.transAxes, verticalalignment='top',
              horizontalalignment='left', color='black',
              bbox=dict(facecolor='white', edgecolor='gray', boxstyle='round,pad=0.4'),
              fontsize=15, zorder=30)

    ax_n.scatter(epochs_list, NList, marker='+', s=5, zorder=10, color='green')
    ax_n.plot(epochs_list, NList, color='gray', alpha=0.7, zorder=10)
    ax_n.set_ylabel('N [m]', fontsize=16)
    ax_n.text(0.02, 0.92, f'RMS: {rms(NList):.3f} m',
              transform=ax_n.transAxes, verticalalignment='top',
              horizontalalignment='left', color='black',
              bbox=dict(facecolor='white', edgecolor='gray', boxstyle='round,pad=0.4'),
              fontsize=15, zorder=30)

    ax_u.scatter(epochs_list, UList, marker='+', s=5, zorder=10, color='red')
    ax_u.plot(epochs_list, UList, color='gray', alpha=0.7, zorder=10)
    ax_u.set_ylabel('U [m]', fontsize=16)
    ax_u.text(0.02, 0.92, f'RMS: {rms(UList):.3f} m',
              transform=ax_u.transAxes, verticalalignment='top',
              horizontalalignment='left', color='black',
              bbox=dict(facecolor='white', edgecolor='gray', boxstyle='round,pad=0.4'),
              fontsize=15, zorder=30)

    ax_dop.plot(epochs_list, HDOPList, label='HDOP', color='blue', linewidth=1.5)
    ax_dop.plot(epochs_list, VDOPList, label='VDOP', color='green', linewidth=1.5)
    ax_dop.plot(epochs_list, PDOPList, label='PDOP', color='red', linewidth=1.5)
    ax_dop.plot(epochs_list, GDOPList, label='GDOP', color='purple', linewidth=1.5)
    ax_dop.legend(loc='upper right', fontsize=12, ncol=2)
    ax_dop.set_ylabel('DOP', fontsize=16)
    ax_dop.xaxis.set_major_formatter(mdate.DateFormatter('%dD-%HH'))

    avg_hdop = np.mean(HDOPList)
    avg_vdop = np.mean(VDOPList)
    avg_pdop = np.mean(PDOPList)
    avg_gdop = np.mean(GDOPList)
    ax_dop.text(0.02, 0.92,
                f'Avg HDOP: {avg_hdop:.2f}  Avg VDOP: {avg_vdop:.2f}\n'
                f'Avg PDOP: {avg_pdop:.2f}  Avg GDOP: {avg_gdop:.2f}',
                transform=ax_dop.transAxes, verticalalignment='top',
                horizontalalignment='left', color='black',
                bbox=dict(facecolor='white', edgecolor='gray', boxstyle='round,pad=0.4'),
                fontsize=12, zorder=30)

    ax_e.set_xticklabels([])
    ax_n.set_xticklabels([])
    ax_u.set_xticklabels([])

    fig.suptitle(f'POS: [{XMEAN:.2f}, {YMEAN:.2f}, {ZMEAN:.2f}]', fontsize=20)
    fig.subplots_adjust(left=0.09, right=0.99, bottom=0.06, top=0.90, hspace=0.3)

    if plot_save:
        plt.savefig(plot_save, dpi=300, bbox_inches='tight')
    if plot_show:
        plt.show()
    else:
        plt.close(fig)