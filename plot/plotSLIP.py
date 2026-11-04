import matplotlib.pyplot as plt

def plotCycleSlip(slipData, pngFile=None):
    """
    Plot Cycle Slip Ratio (CSR) as bar chart per satellite.
    Satellites are sorted naturally (e.g., C01, C02...).

    Parameters
    ----------
    slipData : dict
        Structure: {sys: {prn: {'slip': [0/1, ...], 'epochs': [...]}}}
    pngFile : str, optional
    """
    import matplotlib.pyplot as plt

    fig = plt.figure()
    plt.subplots_adjust(wspace=1, hspace=1)
    gnssSysNum = len(slipData)
    nowAxNum = 1
    for gnssSys in slipData:
        ax = fig.add_subplot(gnssSysNum, 1, nowAxNum)
        ax.set_ylabel(gnssSys + '_CSR')
        ax.grid(zorder=0)
        ax.tick_params(axis='y', labelsize='medium')

        # 获取该系统下所有卫星并按 PRN 排序
        prn_list = sorted(slipData[gnssSys].keys())
        X = []
        X_PRN = []
        CSR_LIST = []
        prnNum = 1
        for prn in prn_list:
            slip_list = slipData[gnssSys][prn]['slip']
            if sum(slip_list) == 0:                     # 无周跳
                CSR = 0
            else:
                CSR =  1000 / ( len(slip_list) / sum(slip_list) )  #  1000/(总历元数 / 周跳次数)
            CSR_LIST.append(CSR)
            X.append(prnNum)
            X_PRN.append(prn)
            prnNum += 1

        ax.bar(X, CSR_LIST, width=0.5, color='tomato', zorder=10)
        ax.set_xticks(X, X_PRN)
        ax.tick_params(axis='x', labelsize='small', labelrotation=90)
        nowAxNum += 1

    fig.subplots_adjust(hspace=0.35, wspace=0.5, left=0.12, right=0.99, bottom=0.07, top=0.99)
    plt.tight_layout()
    if pngFile:
        plt.savefig(pngFile)
    else:
        plt.show()
    plt.close(fig)
