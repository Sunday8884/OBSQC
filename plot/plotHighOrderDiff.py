import matplotlib.pyplot as plt
import matplotlib.dates as mdate

def plotHighOrderDiff(highOrderDiffData, pngFile=None):
    """
    Plot high‑order difference data.

    Parameters
    ----------
    highOrderDiffData : dict
        As returned by getObsHighOrderDiff().
    pngFile : str, optional
    """
    fig = plt.figure()
    gnssSysNum = len(highOrderDiffData)
    nowAxNum = 1
    for gnssSys in highOrderDiffData:
        ax = fig.add_subplot(gnssSysNum, 1, nowAxNum)
        ax.set_ylabel(gnssSys + '_L [m]')
        ax.grid(zorder=0)
        ax.tick_params(axis='y', labelsize=10)
        for prn in highOrderDiffData[gnssSys]:
            for band in highOrderDiffData[gnssSys][prn]:
                epochList = highOrderDiffData[gnssSys][prn][band]['epoch']
                highList = highOrderDiffData[gnssSys][prn][band]['highOrder']
                ax.scatter(epochList, highList, marker='+', s=1, zorder=10, alpha=0.8)
        if nowAxNum != gnssSysNum:
            ax.set_xticklabels([])
        else:
            xfmt = mdate.DateFormatter('%dD-%H:%M')
            ax.xaxis.set_major_formatter(xfmt)
            ax.tick_params(axis='x', labelsize=10)
        nowAxNum += 1
    fig.subplots_adjust(left=0.12, right=0.99, bottom=0.1, top=0.99)
    if pngFile:
        plt.savefig(pngFile)
    else:
        plt.show()

    plt.close(fig)
