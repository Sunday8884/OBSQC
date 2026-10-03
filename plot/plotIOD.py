import matplotlib.pyplot as plt
import matplotlib.dates as mdate

def plotIOD(iodData, pngFile=None):
    """
    Plot Ionospheric Delay Rate (IOD) data.

    Parameters
    ----------
    iodData : dict
        As returned by IOD().
    pngFile : str, optional
    """
    fig = plt.figure()
    gnssSysNum = len(iodData)
    nowAxNum = 1
    for gnssSys in iodData:
        ax = fig.add_subplot(gnssSysNum, 1, nowAxNum)
        ax.set_ylabel(gnssSys + '_L [m/s]', fontsize='medium')
        ax.grid(zorder=0)
        ax.tick_params(axis='y', labelsize='medium')
        for prn in iodData[gnssSys]:
            for band in iodData[gnssSys][prn]:
                epochList = iodData[gnssSys][prn][band]['epoch']
                iodList = iodData[gnssSys][prn][band]['iod']
                ax.scatter(epochList, iodList, marker='+', s=1, zorder=10, alpha=0.8)
        if nowAxNum != gnssSysNum:
            ax.set_xticklabels([])
        else:
            xfmt = mdate.DateFormatter('%dD-%HH')
            ax.xaxis.set_major_formatter(xfmt)
            ax.tick_params(axis='x', labelsize='medium')
        nowAxNum += 1
    fig.subplots_adjust(left=0.12, right=0.99, bottom=0.08, top=0.96)
    if pngFile:
        plt.savefig(pngFile)
    else:
        plt.show()
    plt.close(fig)