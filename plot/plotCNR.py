import matplotlib.pyplot as plt
import numpy as np
import datetime

def plotCnr(obsHead, obsData, pngFile=None, start_datetime=None, end_datetime=None,
            sat_list=None, sys_list=None, band_choose=None):
    """
    Plot CNR (Carrier‑to‑Noise Ratio) as boxplots per frequency band.

    Parameters
    ----------
    obsHead, obsData : as from readObs
    pngFile : str, optional
    start_datetime, end_datetime : datetime, optional
        Time range to include.
    sat_list : list of str, optional
        PRN list to plot.
    sys_list : list of str, optional
        GNSS systems to include.
    band_choose : dict, optional
        e.g. {'G': ['S1C', 'S2W']}  (note 'S' prefix for CNR)
    """
    # Build cnrData from obsData
    cnrData = {}
    if sat_list is None:
        sat_list = obsHead['prn']
    # Determine epochs to use
    sorted_epochs = sorted(obsData.keys())
    for epoch in sorted_epochs:
        if start_datetime and epoch < start_datetime:
            continue
        if end_datetime and epoch > end_datetime:
            continue
        for prn in obsData[epoch]:
            gnssSys = prn[0]
            if sys_list and gnssSys not in sys_list:
                continue
            if sat_list and prn not in sat_list:
                continue
            if gnssSys not in cnrData:
                cnrData[gnssSys] = {}
            if prn not in cnrData[gnssSys]:
                cnrData[gnssSys][prn] = {}
            for band in obsData[epoch][prn]:
                if band[0] != 'S':   # CNR values are stored with 'S' prefix
                    continue
                if obsData[epoch][prn][band] is None:
                    continue
                if band_choose and gnssSys in band_choose and band not in band_choose[gnssSys]:
                    continue
                if band not in cnrData[gnssSys][prn]:
                    cnrData[gnssSys][prn][band] = {}
                cnrData[gnssSys][prn][band][epoch] = obsData[epoch][prn][band]

    gnssSysNum = len(cnrData)
    fig = plt.figure()
    nowAxNum = 1
    for gnssSys in cnrData:
        ax = fig.add_subplot(gnssSysNum, 1, nowAxNum)
        ax.set_ylabel(f'{gnssSys} [dBHz]', fontsize='medium')
        ax.grid(zorder=0)
        ax.tick_params(axis='y', labelsize='medium')
        allCnrDict = {}
        for prn in cnrData[gnssSys]:
            for band in cnrData[gnssSys][prn]:
                cnrList = list(cnrData[gnssSys][prn][band].values())
                if band not in allCnrDict:
                    allCnrDict[band] = []
                allCnrDict[band] += cnrList
        if not allCnrDict:
            nowAxNum += 1
            continue
        values = list(allCnrDict.values())
        labels = list(allCnrDict.keys())
        bp = ax.boxplot(values, positions=range(1, len(values) + 1))
        for flier in bp['fliers']:
            flier.set(marker='o', markeredgecolor='red', alpha=0.3, markersize=4)
        ax.set_xticks(range(1, len(labels) + 1))
        ax.set_xticklabels(labels, fontsize='small')
        y_min, y_max = ax.get_ylim()
        lower_bound = np.floor(y_min / 10) * 10
        upper_bound = np.ceil(y_max / 10) * 10
        ax.set_ylim(lower_bound, upper_bound)
        ax.set_ylim(0, 60)
        ax.yaxis.set_major_locator(plt.MultipleLocator(20))
        nowAxNum += 1

    fig.subplots_adjust(left=0.08, right=0.99, bottom=0.03, top=0.99, hspace=0.3)
    if pngFile:
        plt.savefig(pngFile, dpi=300, bbox_inches='tight')
    else:
        plt.show()
    plt.close(fig)
    return cnrData