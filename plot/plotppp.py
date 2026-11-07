#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
plotppp.py — PPP 结果绘图。

与 spos.plotspp.plot_spp 的区别：
    - E/N/U 三条误差曲线合并到同一坐标系（原为三幅子图）；
    - 自动统计 PPP 收敛时间：自首个解算历元起，到 3D 误差
      sqrt(E²+N²+U²) 首次“永久”低于阈值（默认 0.3 m）的时刻，
      图中以竖直虚线标出，未收敛段以浅色阴影显示；
    - 文本框给出收敛时间、全天与收敛后 ENU RMS、平均 DOP。

收敛判据（可在 ppp.toml 用 ppp_conv_thresh 配置，单位 m）：
    convIdx = 第一个 i，使得 err3d[i:] 全部 < 阈值。
    若不存在这样的 i，则判定未收敛。
"""
import numpy as np
import matplotlib.dates as mdate
import matplotlib.pyplot as plt

from com.xyz2neu import xyz2neu
from com.pub import rms

# 中文字体（Windows 优先微软雅黑/黑体），并修复负号显示
plt.rcParams['font.sans-serif'] = ['Microsoft YaHei', 'SimHei', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False


def convergence_time(epochs_list, err3d, thresh=0.3, hold_min=5.0):
    """
    计算 PPP 收敛时间（滑动窗判据）。

    判据：首个历元 i，使得其后 hold_min 分钟内的所有历元 3D 误差
    sqrt(E²+N²+U²) 均低于 thresh——即“降入阈值并持续保持”的时刻。

    epochs_list : datetime 列表（升序）
    err3d       : 与 epochs_list 对应的 3D 误差序列
    thresh      : 收敛阈值（m）
    hold_min    : 需要持续保持低于阈值的时长（分钟）

    返回 (convTime_min, convIdx)：未收敛返回 (None, None)。
    """
    n = len(err3d)
    if n < 2:
        return None, None
    dts = np.array([(epochs_list[i + 1] - epochs_list[i]).total_seconds()
                    for i in range(n - 1)])
    dt = float(np.median(dts))
    if dt <= 0:
        dt = 30.0
    k = max(1, int(round(hold_min * 60.0 / dt)))

    for i in range(n - k + 1):
        if err3d[i:i + k].max() < thresh:
            convTime = (epochs_list[i] - epochs_list[0]).total_seconds() / 60.0
            return convTime, i
    return None, None


def plot_ppp(posData, refX, refY, refZ, plot_show=False, plot_save=None,
             conv_thresh=0.3, conv_hold_min=5.0):
    """
    绘制 PPP 结果：上幅为合并的 E/N/U 误差（单坐标系，含收敛时间标注），
    下幅为 DOP 序列。

    参数
    ----
    posData : dict
        {epoch(datetime): {'x','y','z','clk','HDOP','VDOP','PDOP','GDOP', ...}}
    refX, refY, refZ : float
        参考坐标（用于将 XYZ 转成 N/E/U）。
    plot_show : bool
        是否显示绘图窗口。
    plot_save : str or None
        保存路径（png）。为 None 时不保存。
    conv_thresh : float
        收敛阈值（3D 误差，米）。
    conv_hold_min : float
        需持续保持低于阈值的时长（分钟）。
    """
    epochs_list = list(posData.keys())
    if not epochs_list:
        print("无有效数据，无法绘图。")
        return

    Elist, NList, UList = [], [], []
    HDOPList, VDOPList, PDOPList, GDOPList = [], [], [], []
    for epoch in epochs_list:
        N, E, U = xyz2neu(refX, refY, refZ,
                          posData[epoch]['x'], posData[epoch]['y'], posData[epoch]['z'],
                          'WGS84')
        NList.append(N)
        Elist.append(E)
        UList.append(U)
        HDOPList.append(posData[epoch]['HDOP'])
        VDOPList.append(posData[epoch]['VDOP'])
        PDOPList.append(posData[epoch]['PDOP'])
        GDOPList.append(posData[epoch]['GDOP'])

    E = np.asarray(Elist)
    N = np.asarray(NList)
    U = np.asarray(UList)
    err3d = np.sqrt(E ** 2 + N ** 2 + U ** 2)
    convTime, convIdx = convergence_time(epochs_list, err3d, conv_thresh,
                                         conv_hold_min)

    # ---------------- 统计文本框 ----------------
    lines = []
    if convTime is not None:
        lines.append('收敛时间: %.1f min (3D < %.2f m, 持续 %.0f min)'
                     % (convTime, conv_thresh, conv_hold_min))
    else:
        lines.append('未收敛 (阈值 %.2f m, 持续 %.0f min)'
                     % (conv_thresh, conv_hold_min))
    if convIdx is not None and convIdx > 0:
        lines.append('收敛后 RMS  E %.3f  N %.3f  U %.3f m'
                     % (rms(list(E[convIdx:])), rms(list(N[convIdx:])),
                        rms(list(U[convIdx:]))))
    lines.append('全天 RMS     E %.3f  N %.3f  U %.3f m'
                 % (rms(Elist), rms(NList), rms(UList)))
    fig = plt.figure(figsize=(13, 9))
    ax_enu = fig.add_subplot(2, 1, 1)
    ax_dop = fig.add_subplot(2, 1, 2)
    for ax in (ax_enu, ax_dop):
        ax.grid(True, zorder=0)

    # ---------------- 收敛段阴影 ----------------
    if convIdx is not None and convIdx > 0:
        ax_enu.axvspan(epochs_list[0], epochs_list[convIdx],
                       color='0.85', alpha=0.5, zorder=1)

    # ---------------- 合并的 ENU 曲线 ----------------
    for vals, color, label in ((Elist, 'blue', 'E'),
                               (NList, 'green', 'N'),
                               (UList, 'red', 'U')):
        ax_enu.plot(epochs_list, vals, color='gray', alpha=0.35,
                    linewidth=0.8, zorder=5)
        ax_enu.scatter(epochs_list, vals, marker='+', s=5, zorder=10,
                       color=color, label=label)
    ax_enu.axhline(0, color='black', linewidth=0.8, zorder=6)
    if convIdx is not None:
        ax_enu.axvline(epochs_list[convIdx], color='darkorange',
                       linestyle='--', linewidth=1.5, zorder=20)

    ax_enu.set_ylabel('ENU error [m]', fontsize=15)
    ax_enu.set_xlabel('Time', fontsize=15)
    ax_enu.legend(loc='lower right', fontsize=13, ncol=3, framealpha=0.9)
    ax_enu.xaxis.set_major_formatter(mdate.DateFormatter('%dD-%HH:%MM'))

    # ---------------- 统计文本框 ----------------
    ax_enu.text(0.02, 0.96, '\n'.join(lines),
                transform=ax_enu.transAxes, verticalalignment='top',
                horizontalalignment='left', color='black', fontsize=13,
                bbox=dict(facecolor='white', edgecolor='gray',
                          boxstyle='round,pad=0.4'), zorder=30)

    # ---------------- DOP ----------------
    ax_dop.plot(epochs_list, HDOPList, label='HDOP', color='blue', linewidth=1.5)
    ax_dop.plot(epochs_list, VDOPList, label='VDOP', color='green', linewidth=1.5)
    ax_dop.plot(epochs_list, PDOPList, label='PDOP', color='red', linewidth=1.5)
    ax_dop.plot(epochs_list, GDOPList, label='GDOP', color='purple', linewidth=1.5)
    ax_dop.legend(loc='upper right', fontsize=12, ncol=2)
    ax_dop.set_ylabel('DOP', fontsize=15)
    ax_dop.xaxis.set_major_formatter(mdate.DateFormatter('%dD-%HH'))

    fig.suptitle('PPP POS: [%s, %s, %s]' % tuple('%.2f' % v for v in (refX, refY, refZ)),
                 fontsize=18)
    fig.subplots_adjust(left=0.07, right=0.99, bottom=0.07, top=0.91, hspace=0.25)

    if plot_save:
        plt.savefig(plot_save, dpi=300, bbox_inches='tight')
    if plot_show:
        plt.show()
    else:
        plt.close(fig)
