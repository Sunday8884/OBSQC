# from com.pub import rms
# import matplotlib.pyplot as plt
# import matplotlib.dates as mdate
#
# def plotCMC(cmcData, pngFile=None, fig=None, figsize=(12, 6), ax=None, start_datetime=None, end_datetime=None, sat_list=None,
#             band_choose=None):
#     """
#     Plot time-differenced Code Minus Phase (CMC) data.
#
#     Parameters
#     ----------
#     cmcData : dict
#         CMC data structure as returned by CMC() function.
#     pngFile : str, optional
#         If provided, save the figure to this path. Otherwise display interactively.
#     fig : matplotlib.figure.Figure, optional
#         Existing figure to plot on. If None, a new figure is created.
#     ax : matplotlib.axes.Axes, optional
#         Existing axes to plot on (if fig is also provided). If None, a new axes is created.
#     start_datetime, end_datetime : datetime, optional
#         Time range to display (not used for filtering data, only for x‑axis limits).
#     sat_list, band_choose : not used, kept for API consistency.
#     """
#     # Create figure/axes if not given
#     if figsize is None:
#         figsize = (12, 6)
#     if fig is None:
#         fig = plt.figure(figsize=figsize)
#     # if fig is None:
#     #     fig = plt.figure()
#     # else:
#     #     fig = plt.figure(figsize=figsize)
#     if ax is None:
#         ax = fig.add_subplot(111)
#
#     gnssSysNum = len(cmcData)
#     nowAxNum = 1
#     for gnssSys in cmcData:
#         if gnssSysNum > 1:
#             ax = fig.add_subplot(gnssSysNum, 1, nowAxNum)
#         else:
#             ax = ax  # reuse single axes
#         ax.set_ylabel(gnssSys + '_P [m]', fontsize='medium')
#         ax.grid(zorder=0)
#         ax.tick_params(axis='y', labelsize='medium')
#         cmcall = []
#         bandData = {}
#         for prn in cmcData[gnssSys]:
#             for band in cmcData[gnssSys][prn]:
#                 epochList = cmcData[gnssSys][prn][band]['epoch']
#                 cmcList = cmcData[gnssSys][prn][band]['cmc']
#                 cmcall += cmcList
#                 ax.scatter(epochList, cmcList, marker='+', s=1, zorder=10, alpha=0.8)
#                 if band not in bandData:
#                     bandData[band] = []
#                 bandData[band] += cmcList
#         rmsLine = 'RMS:'
#         bandIndex = 0
#         for band in bandData:
#             if len(bandData[band]) == 0:
#                 continue
#             if bandIndex == 3:
#                 rmsLine += '\n         '
#             rmsLine += band + ' (' + '%.2f' % rms(bandData[band]) + ') '
#             bandIndex += 1
#         ax.text(0.02, 0.92, rmsLine,
#                 transform=ax.transAxes, verticalalignment='top',
#                 horizontalalignment='left', color='black',
#                 bbox=dict(facecolor='white', edgecolor='gray', boxstyle='round,pad=0.4'),
#                 zorder=30, fontsize='medium')
#         sortcmc = sorted(cmcall)
#         if len(sortcmc) > 10:
#             y_min = sortcmc[-5]
#             y_max = sortcmc[5]
#         else:
#             y_min = min(cmcall) if cmcall else -5
#             y_max = max(cmcall) if cmcall else 5
#         y_min_rounded = (y_min // 5) * 5
#         y_max_rounded = ((y_max + 4) // 5) * 5
#         y_max_abs = max(abs(y_max_rounded), abs(y_min_rounded))
#         if y_max_abs < 5:
#             y_max_abs = 5
#         ax.set_ylim(-y_max_abs, y_max_abs)
#         if nowAxNum != gnssSysNum:
#             ax.set_xticklabels([])
#         else:
#             xfmt = mdate.DateFormatter('%dD-%HH')
#             ax.xaxis.set_major_formatter(xfmt)
#             ax.tick_params(axis='x', labelsize='medium')
#         nowAxNum += 1
#
#     fig.subplots_adjust(left=0.12, right=0.99, bottom=0.04, top=0.99)
#
#     if pngFile is not None:
#         plt.savefig(pngFile)
#     else:
#         plt.show()
#
#
#
from com.pub import rms
import matplotlib.pyplot as plt
import matplotlib.dates as mdate
import numpy as np

def plotCMC(cmcData, pngFile=None, fig=None, figsize=(12, 6), ax=None,
            start_datetime=None, end_datetime=None, sat_list=None, band_choose=None):
    """
    绘制历元间差分码减相位（CMC）数据。

    Parameters
    ----------
    cmcData : dict
        CMC 数据（由 CMC() 返回）。
    pngFile : str, optional
        保存图片路径。
    fig : matplotlib.figure.Figure, optional
        已有 Figure 对象，默认新建。
    figsize : tuple, optional
        图片尺寸（英寸），默认 (12, 6)。
    ax : matplotlib.axes.Axes, optional
        已有 Axes（仅在单系统时使用），默认自动创建。
    start_datetime, end_datetime : datetime, optional
        时间范围（仅用于设置 x 轴限制，数据已在 CMC 中过滤）。
    sat_list, band_choose : 保留接口，未使用。
    """
    if fig is None:
        fig = plt.figure(figsize=figsize)

    gnssSysNum = len(cmcData)
    if gnssSysNum == 0:
        print("警告：无 CMC 数据可绘图。")
        return

    # 如果是多系统，动态调整 figure 高度
    if gnssSysNum > 1:
        fig.set_figheight(figsize[1] * (0.5 + 0.5 * gnssSysNum))  # 随子图数量增加高度

    # 为每个系统创建子图
    axes = []
    for i, gnssSys in enumerate(cmcData.keys(), 1):
        if gnssSysNum > 1:
            ax = fig.add_subplot(gnssSysNum, 1, i)
        else:
            if ax is None:
                ax = fig.add_subplot(111)
            # 否则复用传入的 ax
        axes.append((gnssSys, ax))

    # 逐个系统绘图
    for idx, (gnssSys, ax) in enumerate(axes):
        ax.set_ylabel(f'{gnssSys}_P [m]', fontsize='medium')
        ax.grid(zorder=0, linestyle='--', alpha=0.6)
        ax.tick_params(axis='y', labelsize='medium')

        all_cmc = []
        band_vals = {}

        # 遍历该系统的所有卫星和波段
        for prn in cmcData[gnssSys]:
            for band, data in cmcData[gnssSys][prn].items():
                epochs = data['epoch']
                cmc_vals = data['cmc']
                if not epochs:
                    continue
                all_cmc.extend(cmc_vals)
                ax.scatter(epochs, cmc_vals, marker='+', s=1, zorder=10, alpha=0.8)
                band_vals.setdefault(band, []).extend(cmc_vals)

        # 添加 RMS 统计文本框
        rms_line = 'RMS:'
        b_idx = 0
        for band, vals in band_vals.items():
            if not vals:
                continue
            if b_idx == 3:
                rms_line += '\n         '
            rms_line += f'{band} ({rms(vals):.2f}) '
            b_idx += 1
        ax.text(0.02, 0.92, rms_line,
                transform=ax.transAxes, verticalalignment='top',
                horizontalalignment='left', color='black',
                bbox=dict(facecolor='white', edgecolor='gray', boxstyle='round,pad=0.4'),
                zorder=30, fontsize='medium')

        # 动态设置 y 轴范围（使用 5%~95% 百分位数，避免离群点拉宽范围）
        if all_cmc:
            sorted_vals = np.sort(all_cmc)
            lower = np.percentile(sorted_vals, 5)
            upper = np.percentile(sorted_vals, 95)
            # 以 5 为步长向上取整
            y_min = np.floor(lower / 5) * 5
            y_max = np.ceil(upper / 5) * 5
            y_abs = max(abs(y_min), abs(y_max))
            y_abs = max(y_abs, 5.0)  # 至少 ±5m
            ax.set_ylim(-y_abs, y_abs)
        else:
            ax.set_ylim(-5, 5)

        # x 轴格式：仅对最后一个子图显示 x 轴标签
        if idx == len(axes) - 1:
            xfmt = mdate.DateFormatter('%H:%M')
            ax.xaxis.set_major_formatter(xfmt)
            ax.tick_params(axis='x', labelsize='medium', rotation=30)
            # 自动调整刻度数量避免重叠
            ax.xaxis.set_major_locator(mdate.AutoDateLocator(maxticks=12))
        else:
            ax.set_xticklabels([])

    fig.subplots_adjust(left=0.12, right=0.99, bottom=0.08, top=0.95, hspace=0.3)

    if pngFile is not None:
        plt.savefig(pngFile, dpi=300, bbox_inches='tight')
    else:
        plt.show()
    plt.close(fig)