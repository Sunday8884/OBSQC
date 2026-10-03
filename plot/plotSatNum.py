import matplotlib.pyplot as plt
import matplotlib.dates as mdate


def plotSatNum(obsData, pngFile = None, start_datetime = None, end_datetime = None,
               gnssSystem = ['G', 'C', 'R', 'E', 'L', 'S', 'J', '0', 'W', 'ALL']):
    """
    绘制每个系统（以及总卫星数 ALL）的卫星数量随时间变化图，并返回统计数据。

    参数
    ----------
    obsData : dict
        从 readObs() 获得的观测数据字典。
    pngFile : str, 可选
        保存图片的路径。若为 None 则交互显示。
    start_datetime, end_datetime : datetime, 可选
        时间范围过滤。
    gnssSystem : list of str
        需要绘制的 GNSS 系统列表（例如 ['G','C']）。'ALL' 始终会被强制添加。

    返回
    -------
    satNumData : dict
        结构：{
            'stat': {系统: {'AVE': float, 'MIN': int, 'MAX': int}},
            'data': {历元: {系统: 计数}}
        }
    """
    if pngFile is not None:
        plt.switch_backend('agg')

    # 强制包含 'ALL' 曲线
    if 'ALL' not in gnssSystem:
        gnssSystem = gnssSystem + ['ALL']

    # 初始化数据容器
    epochList = []
    sys_data = {sys: [] for sys in gnssSystem}  # 每个系统的绘图数据列表
    sys_has_data = {sys: 0 for sys in gnssSystem}  # 标记系统是否有数据

    satNumData = {'stat': {}, 'data': {}}
    for sys in gnssSystem:
        satNumData['stat'][sys] = {'AVE': 0, 'MIN': 0, 'MAX': 0}

    # 获取排序后的历元和采样间隔
    sorted_epochs = sorted(obsData.keys())
    if len(sorted_epochs) < 2:
        interval = 30.0
    else:
        interval = (sorted_epochs[1] - sorted_epochs[0]).total_seconds()

    backEp = sorted_epochs[0]

    for ep in sorted_epochs:
        if start_datetime and ep < start_datetime:
            continue
        if end_datetime and ep > end_datetime:
            continue

        # 统计当前历元各系统的卫星数量
        counts = {sys: 0 for sys in gnssSystem}
        for satPrn in obsData[ep]:
            # 总是统计 ALL
            counts['ALL'] += 1
            sys_has_data['ALL'] += 1
            sys_char = satPrn[0]
            if sys_char in counts:
                counts[sys_char] += 1
                sys_has_data[sys_char] += 1

        # 存储数据用于报告
        satNumData['data'][ep] = counts.copy()

        # 处理数据中断（间隔 > 2倍采样间隔时插入零点断开连线）
        if (ep - backEp).total_seconds() <= interval * 2:
            epochList.append(ep)
            for sys in gnssSystem:
                sys_data[sys].append(counts[sys])
        else:
            # 插入断点：上一个历元位置置0，当前历元位置先置0再置真实值
            epochList.append(backEp)
            for sys in gnssSystem:
                sys_data[sys].append(0)
            epochList.append(ep)
            for sys in gnssSystem:
                sys_data[sys].append(0)
            epochList.append(ep)
            for sys in gnssSystem:
                sys_data[sys].append(counts[sys])
        backEp = ep

    # 创建图形
    fig, ax = plt.subplots(figsize = (12, 5))

    # 绘制有数据的每个系统（图例使用英文）
    for sys in gnssSystem:
        if sys_has_data[sys] == 0:
            continue
        data_list = sys_data[sys]
        if not data_list:
            continue
        ave = sum(data_list) / len(data_list)
        mn = min(data_list)
        mx = max(data_list)
        label = f"{sys} AVE: {ave:.1f} MIN: {mn} MAX: {mx}"
        ax.plot(epochList, data_list, linewidth = 1.5, label = label, zorder = 20)
        satNumData['stat'][sys] = {'AVE': ave, 'MIN': mn, 'MAX': mx}

    ax.set_ylabel('Number Of Sat', fontsize = 'medium')
    ax.legend(loc = "best")
    xfmt = mdate.DateFormatter('%dD-%H:%M')
    ax.xaxis.set_major_formatter(xfmt)
    ax.tick_params(axis = 'x', labelsize = 'small')
    ax.tick_params(axis = 'y', labelsize = 'medium')
    ax.grid(zorder = 10)
    fig.subplots_adjust(left = 0.08, right = 0.99, bottom = 0.05, top = 0.99)

    if pngFile is not None:
        plt.savefig(pngFile)
    else:
        plt.show()

    plt.close(fig)

    return satNumData
