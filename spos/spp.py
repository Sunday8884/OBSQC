from com.sat2siteAngle import getEle, sat2siteAngle
from com.gnssParameter import CLIGHT, bandDictGetFreq, bandDictGetFreqPos, bds3MEOList
from com.xyz2neu import xyz2neu
from com.gnssTime import datetime2allgnssTime
from spos.initobs import initObs
from spos.satpos import satPos
from spos.corr import earthRotation, relativistic, satPCO
from com.xyz2blh import bl2enuRot, geodist, xyz2blh
from com.pub import printPanda, rms
from spos.trop import tropColins
from com.readNav import readNav
from com.nav2posclk import getTgdInNav
from spos.plotspp import plot_spp
import datetime
import numpy as np
import time

def read_snx_coord(snx_file, station = None):
    """
    从 SINEX 文件的 +SOLUTION/ESTIMATE 块中读取指定测站的 STAX/STAY/STAZ。
    返回 numpy 数组 [X, Y, Z]，单位 m。
    """
    true_xyz = np.full(3, np.nan)

    with open(snx_file, 'r', encoding='utf-8') as f:
        for line in f:
            if '+SOLUTION/ESTIMATE' in line:
                break
        else:
            raise RuntimeError('未在 .snx 文件中找到 +SOLUTION/ESTIMATE 标记')

        for line in f:
            if line.startswith('-SOLUTION/ESTIMATE'):
                break
            if line.startswith('*') or not line.strip():
                continue

            tokens = line.split()
            if len(tokens) < 9:
                continue

            typ  = tokens[1]
            code = tokens[2]

            if code == station and typ in ('STAX', 'STAY', 'STAZ'):
                val = float(tokens[8])
                if typ == 'STAX':
                    true_xyz[0] = val
                elif typ == 'STAY':
                    true_xyz[1] = val
                elif typ == 'STAZ':
                    true_xyz[2] = val

    if np.any(np.isnan(true_xyz)):
        raise RuntimeError(f'未找到测站 {station} 的完整 STAX/STAY/STAZ 坐标')

    return true_xyz

def spp(obsHead, obsData, navData, bandChoose, startDatetime, endDatetime,
        posFile=None, plot_show=False, plot_save=None, snx_file=None, station=None):

    start_time = time.time()
    # 内部策略
    posCoord = 'ALL'
    cutoff = 5
    fisrtRun = True

    # 初始坐标
    if obsHead['Approx Position'][0] is not None:
        refCoord = obsHead['Approx Position']

    # cutOff
    cutoff = np.radians(cutoff)

    # 获取频率
    bandChoose = bandDictGetFreqPos(bandChoose)
    gSysNumInEpoch = len(list(bandChoose))

    # 初值
    initDict = [1000.0, 1000.0, 1000.0]
    for gSys in bandChoose:
        initDict.append(0)

    X0 = np.asmatrix(initDict).T

    gSysIndex = {}
    gnssCt = 0
    for gSys in bandChoose:
        gSysIndex[gSys] = gnssCt
        gnssCt += 1

    # 结果
    posData = {}
    X = []
    Y = []
    Z = []

    # 历元循环 最小二乘求解 无电离层组合
    for epoch in obsData:
        if startDatetime is not None and startDatetime > epoch:  # Datetime为空
            continue
        if endDatetime is not None and endDatetime < epoch:
            continue

        satPosInEpoch, satClkInEpoch = satPos(bandChoose, obsData, navData, epoch)
        obsDataInTime, isLess = initObs(bandChoose, obsData, epoch, satPosInEpoch, satClkInEpoch)
        if isLess:
            continue

        less4 = False
        # 历元初始化
        posData[epoch] = {}
        # 时间转换为GNSS
        nowGnssTime = datetime2allgnssTime(epoch)

        # init
        satUse = obsDataInTime['satUse']
        mjds = obsDataInTime['mjds']
        doys = obsDataInTime['doys']
        prnIndex = obsDataInTime['prnIndex']
        satNumInEpoch = len(satUse)

        # init matrix
        gInit = np.zeros((satNumInEpoch, 3 + gSysNumInEpoch))
        wInit = np.zeros((satNumInEpoch, satNumInEpoch))
        yInit = np.zeros((satNumInEpoch, 1))

        iter = 0
        diverged = False
        while iter < 12:
            iter += 1
            if less4:
                break
            satNum = 0
            G = np.copy(gInit)
            W = np.copy(wInit)
            y = np.copy(yInit)

            # init spos clk
            sitX_init = float(X0[0, 0])
            sitY_init = float(X0[1, 0])
            sitZ_init = float(X0[2, 0])
            sitClkAll_init = []
            for gSys in gSysIndex:
                sitClkAll_init.append(float(X0[3 + gSysIndex[gSys], 0]))

            # site spos
            rr = np.array([sitX_init, sitY_init, sitZ_init])
            lat, lon, h = xyz2blh(sitX_init, sitY_init, sitZ_init)
            enuRot = bl2enuRot(lat, lon)
            blh = np.array([lat, lon, h])

            # 卫星循环
            for prn in satUse:
                # 判断系统是否被选择
                gnssSys = prn[0]
                prnIndexInC = prnIndex[prn]

                # init corr
                relativisticCorr = 0.0
                tropCorr = 0.0
                weight = 1
                TGD = 0.0
                shapiroCorr = 0.0

                # sat poss
                rs = np.copy(satPosInEpoch[prn])

                # 获取频率
                P1, P2 = obsDataInTime['obs'][prn]['P1'], obsDataInTime['obs'][prn]['P2']
                f1, f2 = obsDataInTime['gSysInf'][gnssSys]['f1'], obsDataInTime['gSysInf'][gnssSys]['f2']

                # 无电离层组合观测值
                P_IF = (f1 ** 2) / (f1 ** 2 - f2 ** 2) * P1 - (f2 ** 2) / (f1 ** 2 - f2 ** 2) * P2

                # dis satDir
                approxDis, satDir = geodist(rs[:3], rr)

                # sat clock
                satClk = satClkInEpoch[prn]

                # 跳过星历/观测/钟差异常的卫星（NaN 不参与最小二乘，否则污染整个历元）
                if not (np.isfinite(rs[:3]).all()
                        and np.isfinite(P_IF)
                        and np.isfinite(satClk)):
                    continue

                # sit clock
                sitClk_init = np.copy(sitClkAll_init[gSysIndex[gnssSys]])

                # TGD
                TGD = getTgdInNav(navData, f1, f2, prn, epoch)

                # trop
                if not fisrtRun:
                    Zenith, Azimuth, Elevation = getEle(rs[:3], rr, enuRot)
                    if Elevation < cutoff:
                        continue
                    if h < 200000:
                        tropCorr = tropColins(Elevation, lat, lon, h, doys)

                    # weight
                    weight = 1.0

                # power
                W[prnIndexInC, prnIndexInC] = weight

                # direction cosines
                G[prnIndexInC, :3] = satDir
                G[prnIndexInC, gSysIndex[gnssSys] + 3] = 1.0

                # OMC
                approxDis += relativisticCorr - satClk * CLIGHT + sitClk_init + TGD + tropCorr - shapiroCorr
                ls = P_IF - (approxDis)
                y[prnIndexInC, 0] = ls

                # sat num count
                satNum += 1

            less4 = False
            if satNum < 4:
                printPanda(
                    'epoch - ' + epoch.strftime('%Y-%m-%d %H:%M:%S') + ' -> The number of satellites is less than ' + str(
                        satNum))
                less4 = True
                continue

            # lsq
            try:
                X_MAT = sppLsq(G, W, y)
                # 计算DOP值
                Q = np.linalg.pinv(G.T @ G)
                HDOP = np.sqrt(Q[0, 0] + Q[1, 1])
                VDOP = np.sqrt(Q[2, 2])
                PDOP = np.sqrt(Q[0, 0] + Q[1, 1] + Q[2, 2])
                GDOP = np.sqrt(Q[0, 0] + Q[1, 1] + Q[2, 2] + Q[3, 3])
                posData[epoch]['HDOP'] = HDOP
                posData[epoch]['VDOP'] = VDOP
                posData[epoch]['PDOP'] = PDOP
                posData[epoch]['GDOP'] = GDOP
            except:
                continue

            # 估计值
            Xi = X0 + X_MAT
            if not np.isfinite(Xi).all():
                # 解算发散（NaN/inf）：丢弃本历元，且不更新初值，避免污染后续历元
                diverged = True
                break
            fisrtRun = False
            # next epoch init
            X0 = np.copy(Xi)
            # del bad sat
            if iter > 3:
                delSat = {}
                to_remove = None
                for lIndex, l_value in enumerate(y):
                    if abs(l_value[0]) > 10000:
                        delSat[satUse[lIndex]] = abs(l_value[0])
                        X_MAT[0] = 1
                if delSat:
                    max_delSat_key = max(delSat, key = delSat.get)
                    to_remove = satUse.index(max_delSat_key)
                if to_remove is not None:
                    del satUse[to_remove]

            if abs(X_MAT[0]) <= 0.0001 and abs(X_MAT[1]) <= 0.0001 and abs(X_MAT[2]) <= 0.0001:
                break

        if less4 or diverged:
            continue

        posData[epoch]['x'] = Xi[0][0].item()
        posData[epoch]['y'] = Xi[1][0].item()
        posData[epoch]['z'] = Xi[2][0].item()
        posData[epoch]['clk'] = Xi[3][0].item()
        X.append(Xi[0][0].item())
        Y.append(Xi[1][0].item())
        Z.append(Xi[2][0].item())

    if len(X) == 0:
        return None, None, None, None
    XMEAN = (sum(X) / len(X))
    YMEAN = (sum(Y) / len(Y))
    ZMEAN = (sum(Z) / len(Z))

    print('SPP POS -> ' + '%15.3f' % XMEAN + '%15.3f' % YMEAN + '%15.3f' % ZMEAN)
    end_time = time.time()
    execution_time = end_time - start_time

    # === 选择 ENU 参考坐标：优先使用 SNX 周解，否则回退到坐标均值 ===
    refX, refY, refZ = XMEAN, YMEAN, ZMEAN
    if snx_file is not None and station is not None:
        try:
            refX, refY, refZ = read_snx_coord(snx_file, station)
            print('SNX REF -> ' + '%15.3f' % refX + '%15.3f' % refY + '%15.3f' % refZ)
        except Exception as e:
            print(f"读取 SNX 失败（{e}），改用坐标均值作为 ENU 参考。")
            refX, refY, refZ = XMEAN, YMEAN, ZMEAN

    #  pos 文件（ENU 相对 refX/refY/refZ）
    if posFile is not None:
        writePosData(posData, posFile, refX, refY, refZ)

    # 绘图（ENU与标题均使用 refX/refY/refZ）
    if plot_show or plot_save:
        plot_spp(posData, refX, refY, refZ,
                 plot_show=plot_show, plot_save=plot_save)

    # print("spos Time : ", execution_time, "s")
    return posData, XMEAN, YMEAN, ZMEAN

def writePosData(posData, posFile, XMEAN=None, YMEAN=None, ZMEAN=None):
    """
    写出 SPP 结果。
    输出列：
        date  mjd  sod   X[m]  Y[m]  Z[m]  CLK[ns]  E[m]  N[m]  U[m]  HDOP  VDOP  PDOP  GDOP
    其中 E/N/U 为相对所有有效历元坐标均值（XMEAN,YMEAN,ZMEAN）的东、北、天误差。
    """
    if XMEAN is None or YMEAN is None or ZMEAN is None:
        if len(posData) == 0:
            return
        XMEAN = float(np.mean([posData[ep]['x'] for ep in posData]))
        YMEAN = float(np.mean([posData[ep]['y'] for ep in posData]))
        ZMEAN = float(np.mean([posData[ep]['z'] for ep in posData]))

    posOpen = open(posFile, 'w+')
    posOpen.write('\n')
    header = ('%25s' % 'date' +
              '%8s' % 'mjd' +
              '%11s' % 'sod' +
              '%20s' % 'X[m]' + '%20s' % 'Y[m]' + '%20s' % 'Z[m]' +
              '%13s' % 'E[m]' + '%13s' % 'N[m]' + '%13s' % 'U[m]' +
              '%20s' % 'CLK[ns]' +
              '%10s' % 'HDOP' + '%10s' % 'VDOP' + '%10s' % 'PDOP' + '%10s' % 'GDOP' +
              '\n')
    posOpen.write(header)

    for epoch in posData:
        N, E, U = xyz2neu(XMEAN, YMEAN, ZMEAN,
                          posData[epoch]['x'], posData[epoch]['y'], posData[epoch]['z'],
                          'WGS84')
        nowGnssTime = datetime2allgnssTime(epoch)
        posOpen.write(
            ' ' + epoch.strftime('%Y-%m-%d %H:%M:%S') +
            '%6d' % nowGnssTime.mjd +
            '%10.2f' % nowGnssTime.sod +
            '%20.3f' % posData[epoch]['x'] +
            '%20.3f' % posData[epoch]['y'] +
            '%20.3f' % posData[epoch]['z'] +
            '%13.3f' % E +
            '%13.3f' % N +
            '%13.3f' % U +
            '%20.5f' % (posData[epoch]['clk'] / CLIGHT * 1.e9) +
            '%10.3f' % posData[epoch]['HDOP'] +
            '%10.3f' % posData[epoch]['VDOP'] +
            '%10.3f' % posData[epoch]['PDOP'] +
            '%10.3f' % posData[epoch]['GDOP'] +
            '\n')
    posOpen.close()


def sppLsq(G_MAT, W_MAT, y_MAT):
    G_MAT = G_MAT.astype(np.float64)
    W_MAT = W_MAT.astype(np.float64)
    y_MAT = y_MAT.astype(np.float64)
    P_MAT_inv = G_MAT.T @ W_MAT @ G_MAT
    P_MAT = np.linalg.pinv(P_MAT_inv)
    X_MAT = P_MAT @ G_MAT.T @ W_MAT @ y_MAT
    return X_MAT


if __name__ == '__main__':
    obsFile = r'F:\OBSQC\test\wuh20050.24o'
    brdFile = r'F:\OBSQC\test\brdc0050.24p'

    bandChoose = {
                  'C': {'C2I': {'freq': 1561098000.0, 'lambda': 0.19203948631027648},
                        'C6I': {'freq': 1268520000.0, 'lambda': 0.2363324646044209}}
                  }

    from com.readObs import readObs, readObsHead
    import datetime

    startDatetime = datetime.datetime(2024, 1, 5)
    endDatetime = datetime.datetime(2024, 1, 6)
    obsData = readObs(obsFile)
    obsHead = readObsHead(obsFile)
    brdData = readNav(brdFile)
    posFile = r'F:\OBSQC\test\wuh20050.txt'
    posData, XMEAN, YMEAN, ZMEAN = spp(obsHead, obsData, brdData, bandChoose,
                                       startDatetime, endDatetime,
                                       posFile=posFile,
                                       plot_show=False, plot_save='sppresult.png')