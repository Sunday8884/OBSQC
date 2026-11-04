import numpy as np
import matplotlib.dates as mdate
import matplotlib.pyplot as plt

def plotFreq(obsHead, obsData, pngFile=None, start_datetime=None, end_datetime=None,
             sat_list=None, sys_list=None, band_choose=None):
    """
    绘制每个卫星的频段可用性散点图，并返回频段统计。

    参数
    ----------
    obsHead, obsData : 观测文件头和数据
    pngFile : 保存图片路径，可选
    start_datetime, end_datetime : 时间范围
    sat_list : 卫星列表，默认使用 obsHead['prn']，再根据 sys_list 过滤
    sys_list : 系统列表，例如 ['G','C']
    band_choose : 波段选择字典，例如 {'G': ['L1C','L2W']}，若为 None 或空，则自动使用所有 L 波段
    """
    # 确定卫星列表
    if sat_list is None:
        sat_list = obsHead['prn']
    if sys_list is not None:
        sat_list = [prn for prn in sat_list if prn[0] in sys_list]
    if not sat_list:
        print("警告：没有符合条件的卫星，无法绘制频段图。")
        return {}

    sat_list.sort(reverse=True)
    satBZ = np.arange(len(sat_list))

    obsType = obsHead['OBS TYPES']
    plotList = {}
    plotListIndex = {}
    prnIndex = 0
    for prn in sat_list:
        plotList[prn] = {}
        plotListIndex[prn] = {}
        plotListIndex[prn]['index'] = prnIndex
        if prn[0] in obsType:
            freqIndex = 0
            for freq in obsType[prn[0]]:
                if freq[0] == 'L':
                    plotList[prn][freq] = np.empty((0, 2))
                    plotListIndex[prn][freq] = freqIndex
                    freqIndex += 1
        prnIndex += 1

    # 确定要绘制的频段
    freqChoose = {}
    if band_choose is None or band_choose == {}:
        # 自动使用所有 L 波段
        for gSys in obsType:
            freqChoose[gSys] = []
            for freq in obsType[gSys]:
                if freq[0] == 'L':
                    freqChoose[gSys].append(freq)
    else:
        for gSys, bands in band_choose.items():
            if gSys in obsType:
                freqChoose[gSys] = [f for f in obsType[gSys] if f in bands]

    # 如果仍然没有可用的频段，返回空
    if not any(freqChoose.values()):
        print("警告：没有找到可绘制的 L 波段。")
        return {}

    # 创建图形和轴
    fig, ax = plt.subplots(figsize=(15, 20))
    ax.set_yticks(satBZ)
    ax.set_yticklabels(list(plotList))
    ax.set_ylim(-0.5, satBZ[-1] + 0.5)
    if start_datetime and end_datetime:
        ax.set_xlim(start_datetime, end_datetime)
    else:
        epochs = sorted(obsData.keys())
        if epochs:
            ax.set_xlim(epochs[0], epochs[-1])
    ax.tick_params(axis='y', labelsize='small')
    ax.grid(zorder=0, alpha=0.1)

    colorList = ['orange', 'royalblue', 'tomato', 'g', 'lightpink', 'darkcyan', 'indigo', 'teal']
    markerList = ['o', 'o', '.', '.', '.', '.', '.', '.', '.', '.']
    markerSizeList = [5, 2, 2, 2, 1, 1, 1, 1, 1, 1]
    layerList = [0.0, 0.0, -0.1, 0.1, -0.1, 0.1, -0.1, 0.1]
    freqListTmp = set()
    freqData = {}

    for idx, prn in enumerate(sat_list):
        gSys = prn[0]
        if gSys not in freqChoose:
            continue
        freqIndex = 0
        freqData[prn] = {}
        for freq in freqChoose[gSys]:
            freqLabel = gSys + '_' + freq
            epochList = []
            yvals = []
            count = 0
            for epoch in obsData:
                if start_datetime and epoch < start_datetime:
                    continue
                if end_datetime and epoch > end_datetime:
                    continue
                if prn in obsData[epoch] and obsData[epoch][prn].get(freq) is not None:
                    epochList.append(epoch)
                    yvals.append(idx + layerList[freqIndex % len(layerList)])
                    count += 1
            freqData[prn][freq] = count
            if epochList:
                if freqLabel not in freqListTmp:
                    freqListTmp.add(freqLabel)
                    ax.scatter(epochList, yvals,
                               color=colorList[freqIndex % len(colorList)],
                               marker=markerList[freqIndex % len(markerList)],
                               s=markerSizeList[freqIndex % len(markerSizeList)],
                               zorder=10, label=freqLabel, alpha=0.5)
                else:
                    ax.scatter(epochList, yvals,
                               color=colorList[freqIndex % len(colorList)],
                               marker=markerList[freqIndex % len(markerList)],
                               s=markerSizeList[freqIndex % len(markerSizeList)],
                               zorder=10, alpha=0.5)
            freqIndex += 1

    ax.tick_params(axis='x', labelsize='small')
    if freqListTmp:
        ax.legend(loc='upper right', bbox_to_anchor=(1.16, 1), markerscale=3)
    xfmt = mdate.DateFormatter('%dD-%H:%M')
    ax.xaxis.set_major_formatter(xfmt)
    fig.subplots_adjust(left=0.08, right=0.87, bottom=0.03, top=0.99)
    plt.tight_layout()

    if pngFile is not None:
        plt.savefig(pngFile)
    else:
        plt.show()

    plt.close(fig)

    return freqData
