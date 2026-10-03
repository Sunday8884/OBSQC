import matplotlib.pyplot as plt
import matplotlib.dates as mdate
import matplotlib
from com.pub import rms

def plotMultipath(mpData, pngFile=None, figsize=(84/25.4*2.0, 84/25.4), dpi=500):
    """
    Plot pseudorange multipath.

    Parameters
    ----------
    mpData : dict
        As returned by multipath().
    pngFile : str, optional
    figsize : tuple, optional
        Figure size in inches.
    dpi : int, optional
        Resolution for saved figure.
    """
    matplotlib.rcParams['font.family'] = 'Arial'
    matplotlib.rcParams['font.size'] = 10
    fig = plt.figure(figsize=figsize, dpi=dpi)
    gnssSysNum = len(mpData)
    nowAxNum = 1
    for gnssSys in mpData:
        ax = fig.add_subplot(gnssSysNum, 1, nowAxNum)
        ax.set_ylabel(gnssSys + ' [m]', fontsize='medium')
        ax.yaxis.labelpad = 0
        ax.grid(zorder=0)
        ax.tick_params(axis='y', labelsize='medium')
        ax.yaxis.set_major_locator(matplotlib.ticker.MultipleLocator(5))
        bandData = {}
        for prn in mpData[gnssSys]:
            for band in mpData[gnssSys][prn]:
                epochList = mpData[gnssSys][prn][band]['epoch']
                mpList = mpData[gnssSys][prn][band]['mp']
                ax.scatter(epochList, mpList, marker='+', s=1, zorder=10)
                if band not in bandData:
                    bandData[band] = []
                bandData[band] += mpList
        rmsLine = 'RMS:'
        bandIndex = 0
        for band in bandData:
            if not bandData[band]:
                continue
            if bandIndex == 3:
                rmsLine += '\n         '
            rmsLine += f'{band} ({rms(bandData[band]):.2f}) '
            bandIndex += 1
        ax.text(0.02, 0.92, rmsLine,
                transform=ax.transAxes, verticalalignment='top',
                horizontalalignment='left', color='black',
                bbox=dict(facecolor='white', edgecolor='gray', boxstyle='round,pad=0.4'),
                zorder=30, fontsize='medium')
        y_min, y_max = ax.get_ylim()
        y_min_rounded = (y_min // 5) * 5
        y_max_rounded = ((y_max + 4) // 5) * 5
        y_lim = max(abs(y_max_rounded), abs(y_min_rounded))
        ax.set_ylim(-y_lim, y_lim)
        if nowAxNum != gnssSysNum:
            ax.set_xticklabels([])
        else:
            xfmt = mdate.DateFormatter('%dD-%HH')
            ax.xaxis.set_major_formatter(xfmt)
            ax.tick_params(axis='x', labelsize='medium')
        nowAxNum += 1
    fig.subplots_adjust(left=0.15, right=0.99, bottom=0.10, top=0.95, hspace=0.1)
    #fig.tight_layout(pad = 2.0, h_pad = 2.0, w_pad = 2.0)
    if pngFile:
        plt.savefig(pngFile)
    else:
        plt.show()

    plt.close(fig)