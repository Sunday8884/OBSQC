import matplotlib.pyplot as plt
import matplotlib.dates as mdate
from com.pub import rms

def plotPhaseNoise(phaseNoiseData, pngFile=None):
    """
    Plot phase noise estimated from third‑order difference of GF combination.

    Parameters
    ----------
    phaseNoiseData : dict
        As returned by phaseNoise().
    pngFile : str, optional
    """
    fig = plt.figure()
    gnssSysNum = len(phaseNoiseData)
    nowAxNum = 1
    for gnssSys in phaseNoiseData:
        ax = fig.add_subplot(gnssSysNum, 1, nowAxNum)
        ax.set_ylabel(gnssSys + ' [cm]', fontsize='medium')
        ax.grid(zorder=0)
        ax.tick_params(axis='y', labelsize='medium')
        bandData = {}
        for prn in phaseNoiseData[gnssSys]:
            for band in phaseNoiseData[gnssSys][prn]:
                epochList = phaseNoiseData[gnssSys][prn][band]['epoch']
                noiseList = [x * 100.0 for x in phaseNoiseData[gnssSys][prn][band]['phaseNoise']]
                ax.scatter(epochList, noiseList, marker='+', s=1, zorder=10, alpha=0.8)
                if band not in bandData:
                    bandData[band] = []
                bandData[band] += noiseList
        rmsLine = 'RMS:'
        for band in bandData:
            if bandData[band]:
                rmsLine += f'L{band[1:]} ({rms(bandData[band]):.2f}) '
        ax.text(0.02, 0.92, rmsLine,
                transform=ax.transAxes, verticalalignment='top',
                horizontalalignment='left', color='black',
                bbox=dict(facecolor='white', edgecolor='gray', boxstyle='round,pad=0.4'),
                zorder=30, fontsize='medium')
        if nowAxNum != gnssSysNum:
            ax.set_xticklabels([])
        else:
            xfmt = mdate.DateFormatter('%dD-%HH')
            ax.xaxis.set_major_formatter(xfmt)
            ax.tick_params(axis='x', labelsize='medium')
        y_min, y_max = ax.get_ylim()
        y_min_rounded = (y_min // 5) * 5
        y_max_rounded = ((y_max + 4) // 5) * 5
        y_lim = max(abs(y_max_rounded), abs(y_min_rounded))
        if y_lim < 5:
            y_lim = 5
        ax.set_ylim(-y_lim, y_lim)
        nowAxNum += 1
    fig.subplots_adjust(left=0.12, right=0.99, bottom=0.06, top=0.93)
    if pngFile:
        plt.savefig(pngFile)
    else:
        plt.show()

    plt.close(fig)