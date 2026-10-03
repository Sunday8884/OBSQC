"""
用户调用 readObs()
    ↓
判断 RINEX 版本（2.x / 3.x）
    ↓
调用 readObs2() / readObs3()
    ↓
若未提供 obsHead → 调用 readObs2Head() / readObs3Head()
    ↓
扫描文件体，逐历元解析：
    - 解析时间行 → 生成 datetime 键
    - 解析卫星 PRN → 初始化该历元下各卫星的空观测值字典
    - 解析观测值行 → 按文件头观测类型顺序填充数值
    (3.x 额外处理 LLI)
    ↓
返回 obsData字典
obsData = {                  # 嵌套字典
    datetime(2025,1,1,0,0,0): {
        'C01': {
            'C1X': 23012345.6,   # 伪距观测值 (米)
            'L1X': 123456789.0,  # 载波相位观测值 (周)
            'L1XLLI': 0,         # 失锁标志 (仅 RINEX 3.x)
            'C2I': 23012346.2,
            'L2I': 123456790.5,
            'L2ILLI': 0,
            ...
        },
        'C02': {...},
        ...
    },
    datetime(2025,1,1,0,30,0): {...},
    ...
}

"""

# 读取rinex文件时间转为datetime
def readRinex3LineTime(line):
    '''
    读取rinex3.x文件时间转为datetime
    '''
    import datetime
    format_string = "> %Y %m %d %H %M %S.%f"
    lineDatetime = datetime.datetime.strptime(line[:28], format_string)
    return lineDatetime
def readRinex2LineTime(line):
    '''
    读取rinex2.x文件时间转为datetime
    '''
    import datetime
    line = str(line).split()
    year = int(line[0])
    month = int(line[1])
    day = int(line[2])
    hour = int(line[3])
    minute = int(line[4])
    second = round(float(line[5]))
    if year < 80:
        year += 2000
    else:
        year += 1900
    lineDatetime = datetime.datetime(year, month, day, hour, minute)
    lineDatetime += datetime.timedelta(seconds = second)
    return lineDatetime
# 读文件头
def readObs2Head(obsFile, needSatList = True):
    '''
    读取rinex2.x文件的文件头
    '''

    obsFileOpen = open(obsFile, 'r+')
    obsFileLines = obsFileOpen.readlines()
    obsFileOpen.close()
    gnssSysList = ['C', 'R', 'G', 'E', 'I', 'J', 'S', 'W', '0']
    obsHead = {'version': None,
               'PGM': 'gnssbox',
               'RUN': 'WHU',
               'MARKER NAME': None,
               'MARKER NUMBER': None,
               'Receiver': None,
               'Receiver Type': None,
               'Antenna': None,
               'Antenna Type': None,
               'Agency': None,
               'Approx Position': [None, None, None],
               'Antenna Delta': [None, None, None],
               'interval': None,
               'Leap Seconds': None,
               'epoch': [],
               'prn': [],
               'OBS TYPES': None}
    satList = []
    endHeadLineNum = 0
    for line in obsFileLines:
        endHeadLineNum += 1
        if 'END OF HEADER' in line:
            endHeadLineNum += 1
            break
    satNum = 0
    if needSatList:
        for lineIndex in range(len(obsFileLines[endHeadLineNum:])):
            nowLineIndex = endHeadLineNum - 1 + lineIndex
            nowLine = obsFileLines[nowLineIndex]
            satListInEpoch = []
            try:
                readRinex2LineTime(nowLine)
                satNum = int(nowLine[30:32])
                while len(satListInEpoch) < satNum:
                    for prn in [nowLine[32:-1][i:i + 3] for i in range(0, len(nowLine[32:-1]), 3)]:
                        if prn == '   ':
                            break
                        if prn[0] == ' ':
                            prn = 'G' + prn[1:]
                        satListInEpoch.append(prn)
                    nowLineIndex += 1
                    nowLine = obsFileLines[nowLineIndex]
                for prn in satListInEpoch:
                    if prn not in satList:
                        satList.append(prn)
            except:
                pass
    obsHead['prn'] = satList
    obs = {}
    for line in obsFileLines:
        if 'RINEX VERSION / TYPE' in line:
            obsHead['version'] = float(line.split()[0])
        elif 'PGM / RUN BY / DATE' in line:
            obsHead['PGM'] = line[:20].strip()
            obsHead['RUN'] = line[20:40].strip()
        elif 'MARKER NAME' in line:
            obsHead['MARKER NAME'] = line[:20].strip()
        elif 'MARKER NUMBER' in line:
            obsHead['MARKER NUMBER'] = line[:20].strip()
        elif 'OBSERVER / AGENCY' in line:
            obsHead['Agency'] = line[20:40].strip()
        elif 'REC # / TYPE / VERS' in line:
            obsHead['Receiver'] = line[:20].strip()
            obsHead['Receiver Type'] = line[20:40].strip()
        elif 'ANT # / TYPE' in line:
            obsHead['Antenna'] = line[:20].strip()
            obsHead['Antenna Type'] = line[20:40].strip()
        elif 'APPROX POSITION XYZ' in line:
            x = float(line.split()[0])
            y = float(line.split()[1])
            z = float(line.split()[2])
            obsHead['Approx Position'] = [x, y, z]
        elif 'ANTENNA: DELTA H/E/N' in line:
            h = float(line.split()[0])
            e = float(line.split()[1])
            n = float(line.split()[2])
            obsHead['Antenna Delta'] = [h, e, n]
        elif 'INTERVAL' in line:
            obsHead['interval'] = int(float(line.split()[0]))
        elif 'LEAP SECONDS' in line:
            obsHead['Leap Seconds'] = int(line.split()[0])
        elif '# / TYPES OF OBSERV' in line:
            gnssSys = line[0]
            if gnssSys == ' ':
                gnssSys = 'G'
            if gnssSys in gnssSysList:
                if gnssSys in obs:
                    obs[gnssSys] += line[10:60].split()
                else:
                    obs[gnssSys] = line[10:60].split()
            else:
                obs[gnssSys] += line[10:60].split()

    for prn in satList:
        if prn[0] not in obs:
            obs[prn[0]] = obs['G']

    for gSys in obs:
        newBandList = []
        for bandIndex in range(len(obs[gSys])):
            newBandList.append(obsType2to3(gSys, obs[gSys][bandIndex]))
        obs[gSys] = newBandList

    obsHead['OBS TYPES'] = obs
    return obsHead
def readObs3Head(obsFile, needSatList = False):
    '''
    读取rinex3.x文件的文件头
    '''

    obsHead = {'version': None,
               'PGM': 'gnssbox',
               'RUN': 'WHU',
               'MARKER NAME': None,
               'MARKER NUMBER': None,
               'Receiver': None,
               'Receiver Type': None,
               'Antenna': None,
               'Antenna Type': None,
               'Agency': None,
               'Approx Position': [None, None, None],
               'Antenna Delta': [None, None, None],
               'interval': None,
               'Leap Seconds': None,
               'epoch': [],
               'prn': [],
               'OBS TYPES': None}
    gnssSysList = ['C', 'R', 'G', 'E', 'I', 'J', 'S', 'W', '0', 'X']
    obs = {}
    obsFileOpen = open(obsFile, 'r+')
    line = obsFileOpen.readline()
    readBegin = False
    satList = []
    while line != '':
        if needSatList and readBegin:
            if '>' not in line:
                sat = line[:3]
                gnssSys = line[0]
                if gnssSys in ['W', '0', 'X', 'L', '5', '1']:
                    gnssSys = 'W'
                    sat = gnssSys + sat[1:]
                if sat not in satList:
                    satList.append(sat)
        if 'RINEX VERSION / TYPE' in line:
            obsHead['version'] = float(line.split()[0])
        elif 'PGM / RUN BY / DATE' in line:
            obsHead['PGM'] = line[:20].strip()
            obsHead['RUN'] = line[20:40].strip()
        elif 'MARKER NAME' in line:
            obsHead['MARKER NAME'] = line[:20].strip()
        elif 'MARKER NUMBER' in line:
            obsHead['MARKER NUMBER'] = line[:20].strip()
        elif 'OBSERVER / AGENCY' in line:
            obsHead['Agency'] = line[20:40].strip()
        elif 'REC # / TYPE / VERS' in line:
            obsHead['Receiver'] = line[:20].strip()
            obsHead['Receiver Type'] = line[20:40].strip()
        elif 'ANT # / TYPE' in line:
            obsHead['Antenna'] = line[:20].strip()
            obsHead['Antenna Type'] = line[20:40].strip()
        elif 'APPROX POSITION XYZ' in line:
            x = float(line.split()[0])
            y = float(line.split()[1])
            z = float(line.split()[2])
            obsHead['Approx Position'] = [x, y, z]
        elif 'ANTENNA: DELTA H/E/N' in line:
            h = float(line.split()[0])
            e = float(line.split()[1])
            n = float(line.split()[2])
            obsHead['Antenna Delta'] = [h, e, n]
        elif 'INTERVAL' in line:
            obsHead['interval'] = int(float(line.split()[0]))
        elif 'LEAP SECONDS' in line:
            obsHead['Leap Seconds'] = int(line.split()[0])
        elif 'SYS / # / OBS TYPES' in line:
            if line[0] in gnssSysList:
                gnssSys = line[0]
                if gnssSys in ['W', '0', 'X', 'L', '5', '1']:
                    gnssSys = 'W'
                if gnssSys in obs:
                    obs[gnssSys] = obs[gnssSys] + line[7:60].split()
                else:
                    obs[gnssSys] = line[7:60].split()
            else:
                obs[gnssSys] += line[7:60].split()
        elif 'END OF HEADER' in line:
            if needSatList:
                readBegin = True
            else:
                break
        line = obsFileOpen.readline()
    obsHead['prn'] = satList
    obsFileOpen.close()

    for gSys in obs:
        newBandList = []
        for bandIndex in range(len(obs[gSys])):
            band = obs[gSys][bandIndex]
            if len(band) == 2:
                newBandList.append(obsType2to3(gSys, band))
            else:
                newBandList.append(band)
        obs[gSys] = newBandList

    obsHead['OBS TYPES'] = obs
    return obsHead
# 读文件体
def readObs2(obsFile, obsHead = None):
    '''
    读取rinex2.x文件的文件体
    存储方式 -> obsData[epoch][prn][band] -> band value
    '''
    import time
    start_time = time.time()

    if obsHead is None:
        obsHead = readObs2Head(obsFile, needSatList = True)
    obsFileOpen = open(obsFile, 'r+')
    obsFileLines = obsFileOpen.readlines()
    obsFileOpen.close()

    obsData = {}
    endHeadLineNum = 0
    for line in obsFileLines:
        endHeadLineNum += 1
        if 'END OF HEADER' in line:
            break

    isHeadLine = True
    isHead2Line = False
    satNum = 0
    bandIndex = 0
    nowSatNum = 0
    for lineIndex in range(len(obsFileLines[endHeadLineNum:])):
        line = obsFileLines[endHeadLineNum + lineIndex]

        try:
            if 'COMMENT' in line:
                satIndex = 0
                bandIndex = 0
                continue
            if isHeadLine:
                obsDatetime = readRinex2LineTime(line)
                obsData[obsDatetime] = {}
                satNum = int(line[30:32])
            if isHeadLine or isHead2Line:
                nowSatNum = 0
                for prn in [line[32:-1][i:i + 3] for i in range(0, len(line[32:-1]), 3)]:
                    if prn == '   ':
                        break
                    if prn[0] == ' ':
                        prn = 'G' + prn[1:]
                    obsData[obsDatetime][prn] = {}
                    nowSys = prn[0]
                    for band in obsHead['OBS TYPES'][nowSys]:
                        obsData[obsDatetime][prn][band] = None
                if len(obsData[obsDatetime]) < satNum:
                    isHead2Line = True
                    isHeadLine = False
                else:
                    isHead2Line = False
                    isHeadLine = False
            elif not isHeadLine and not isHead2Line:
                nowPrn = list(obsData[obsDatetime])[nowSatNum]
                nowSys = nowPrn[0]
                # obsHead['OBS TYPES'][nowSys]
                line = "{:<80}".format(line[:-1])
                for bandStr in [line[i:i + 16] for i in range(0, len(line), 16)]:
                    bandStr = bandStr[:14].strip()
                    try:
                        bandValue = float(bandStr)
                    except:
                        bandValue = None
                    nowBandName = obsHead['OBS TYPES'][nowSys][bandIndex]
                    if bandIndex == len(obsHead['OBS TYPES'][nowSys]) - 1:
                        obsData[obsDatetime][nowPrn][nowBandName] = bandValue
                        bandIndex = 0
                        nowSatNum += 1
                        break
                    obsData[obsDatetime][nowPrn][nowBandName] = bandValue
                    bandIndex += 1
                if nowPrn == list(obsData[obsDatetime])[-1]:
                    isHeadLine = True
        except:
            isHeadLine = True
            isHead2Line = False
            bandIndex = 0
            satNum = 0
            continue

    end_time = time.time()
    execution_time = end_time - start_time
    print("Read OBS Time : ", execution_time, "s")
    return obsData
def readObs3(obsFile, obsHead = None):
    '''
    读取rinex3.x文件的文件体, 存储方式 -> obsDat[epoch][prn][band] -> band value
    '''
    import time
    start_time = time.time()

    if obsHead is None:
        obsHead = readObs3Head(obsFile)
    obsFileOpen = open(obsFile, 'r+')
    obsFileLines = obsFileOpen.readlines()
    obsFileOpen.close()

    obsData = dict()
    endHeadLineNum = 0
    for line in obsFileLines:
        endHeadLineNum += 1
        if 'END OF HEADER' in line:
            endHeadLineNum += 1
            break
    obsDatetime = None
    lineIndex = endHeadLineNum - 1
    for line in obsFileLines[endHeadLineNum - 1:]:
        lineIndex += 1
        if line[0] == '>':
            if obsDatetime is not None:
                obsData[obsDatetime].update(prnData)
            obsDatetime = readRinex3LineTime(line)
            if obsDatetime not in obsData.keys():
                obsData[obsDatetime] = dict()
                prnData = dict()
        else:
            prn = line[:3]
            gnssSys = line[0]
            if gnssSys in ['W', '0', 'X', 'L', '5', '1']:
                gnssSys = 'W'
                prn = gnssSys + prn[1:]
            nowLineIndex = 3
            if gnssSys not in obsHead['OBS TYPES']:
                continue
            if prn[1] == ' ':
                prn = prn[0] + '0' + prn[2]
            prnData[prn] = dict()
            for band in obsHead['OBS TYPES'][gnssSys]:
                bandStr = line[nowLineIndex:nowLineIndex + 14].strip()
                nowLineIndex += 16
                bandValue = None
                if bandStr:
                    try:
                        bandValue = float(bandStr)
                        if bandValue < -999999999.0 or bandValue == 0.0:
                            bandValue = None
                    except:
                        pass

                    if band[0] == 'L':
                        lli_idx = nowLineIndex - 2
                        if lli_idx < len(line) and line[lli_idx].strip():
                            prnData[prn][band + 'LLI'] = int(line[lli_idx])
                        else:
                            prnData[prn][band + 'LLI'] = 0
                else:
                    if band[0] == 'L': prnData[prn][band + 'LLI'] = 0
                prnData[prn][band] = bandValue
    obsData[obsDatetime].update(prnData)

    end_time = time.time()
    execution_time = end_time - start_time
    print("Read OBS Time : ", execution_time, "s")
    return obsData

# 读观测文件
def readObs(obsFile, obsHead = None):
    '''
    obsData = {                  # 嵌套字典
    datetime(2025,1,1,0,0,0): {
        'C01': {
            'C1X': 23012345.6,   # 伪距观测值 (米)
            'L1X': 123456789.0,  # 载波相位观测值 (周)
            'L1XLLI': 0,         # 失锁标志 (仅 RINEX 3.x)
            'C2I': 23012346.2,
            'L2I': 123456790.5,
            'L2ILLI': 0,
            ...
        },
        'C02': {...},
        ...
    },
    datetime(2025,1,1,0,30,0): {...},
    ...
}
    '''
    obsFileOpen = open(obsFile, 'r+')
    obsFileLine = obsFileOpen.readline()
    obsFileOpen.close()
    obsVersion = 0.0
    if 'RINEX VERSION / TYPE' in obsFileLine:
        obsVersion = float(obsFileLine.split()[0])
    else:
        return None

    if 2.0 <= obsVersion < 3.0:
        obsData = readObs2(obsFile, obsHead)
    elif 4.0 > obsVersion >= 3.0:
        obsData = readObs3(obsFile, obsHead)
    else:
        print('Higher versions of Rinex are not currently supported!')
        return None
    return obsData

def readObsHead(obsFile, needSatList = False):
    '''
    obsHead = {                  #字典
    'version': 3.05,
    'MARKER NAME': 'WUH200CHN',
    'prn': ['C01', 'C02', ...],  #文件中出现过的所有卫星 PRN
    'OBS TYPES': {               #观测类型列表
        'C': ['C1X', 'C2I', 'C5X', 'L1X', 'L2I', 'L5X', ...],
        'G': ['C1C', 'C2W', 'L1C', 'L2W', ...],
        ...
    },
    ...
}
    '''
    obsFileOpen = open(obsFile, 'r+')
    obsFileLine = obsFileOpen.readline()
    obsFileOpen.close()
    obsVersion = 0.0
    # for line in obsFileLines:
    if 'RINEX VERSION / TYPE' in obsFileLine:
        obsVersion = float(obsFileLine.split()[0])

    if 2.0 <= obsVersion < 3.0:
        obsHead = readObs2Head(obsFile, needSatList)
    elif 4.0 > obsVersion >= 3.0:
        obsHead = readObs3Head(obsFile, needSatList)
    else:
        print('Higher versions of Rinex are not currently supported!')
        return None
    return obsHead

#===========辅助函数===============
# 检查是否为obs
def is_obs(obs_file):
    '''
    Check whether the file is obs
    '''
    try:
        obs_file_open = open(obs_file, 'r+')
        obs_file_line = obs_file_open.readline()
    except:
        return False
    obs_file_open.close()
    if 'RINEX VERSION / TYPE' in obs_file_line and 'NAV' not in obs_file_line and obs_file_line[20] != 'C':
        return True
    else:
        return False
# 获取给定路径内所有的obs文件
def get_obs_in_path(rinex_path):
    '''
    获取给定路径内所有的obs文件
    '''
    import os
    output_filelist = []
    for root, dirs, files in os.walk(rinex_path):
        for filename in files:
            abs_filename = os.path.join(root, filename)
            if is_obs(abs_filename):
                output_filelist.append(abs_filename)
    return output_filelist
# 获取obs文件内所有的观测卫星
def getSatInObs3(obsFile):
    '''
    获取obs文件内所有的观测卫星
    '''
    obsOpen = open(obsFile, 'r+')
    obsLines = obsOpen.readlines()
    obsOpen.close()
    readBegin = False
    satList = []
    for line in obsLines:
        if readBegin:
            if '>' not in line:
                line = line.split()
                sat = line[0]
                if sat not in satList:
                    satList.append(sat)
        if 'END OF HEADER' in line:
            readBegin = True
    return satList

# 读取观测文件中卫星对应的时间
def getRinex3SatTime(rinexFile):
    '''
    读取观测文件中卫星对应的时间
    rinexFile         : 观测数据文件
    rinexSatTime      : 存储形式
    {prn1 : [datetime1, datetim2, ..., datetimen],
     prn2 : [datetime1, datetim2, ..., datetimen],
     ...,
     prnn : [datetime1, datetim2, ..., datetimen]}
    '''
    rinexFileLineOpen = open(rinexFile, 'r+')
    rinexFileLine = rinexFileLineOpen.readlines()
    rinexFileLineOpen.close()
    satList = []
    rinexSatTime = {}
    endHeadLineNum = 0
    for line in rinexFileLine:
        endHeadLineNum += 1
        if 'END OF HEADER' in line:
            endHeadLineNum += 1
            break
    for line in rinexFileLine[endHeadLineNum:]:
        satPrn = line[:3]
        if satPrn not in satList and line[0] != '>':
            satList.append(satPrn)
            rinexSatTime[satPrn] = []
    timeLine = []
    lineNum = 0
    for line in rinexFileLine:
        if line[0] == '>':
            timeLine.append(lineNum)
        lineNum += 1
    for lineNumIndex in range(0, len(timeLine) - 1):
        nowDateTime = readRinex3LineTime(rinexFileLine[timeLine[lineNumIndex]])
        for line in rinexFileLine[timeLine[lineNumIndex] + 1:timeLine[lineNumIndex + 1]]:
            satPrn = line[:3]
            rinexSatTime[satPrn].append(nowDateTime)
    return rinexSatTime
# 返回排序后的所有历元时间列表
def getRinex3TimeSeries(rinexFile):
    rinexFileLineOpen = open(rinexFile, 'r+')
    rinexFileLine = rinexFileLineOpen.readlines()
    rinexFileLineOpen.close()
    readBegin = False
    timeList = []
    for line in rinexFileLine:
        if readBegin:
            if '>' in line:
                nowTime = readRinex3LineTime(line)
                timeList.append(nowTime)
        if 'END OF HEADER' in line:
            readBegin = True

    def get_list(date):
        return date.timestamp()

    timeList = sorted(timeList, key = lambda date: get_list(date))
    return timeList

# 获取开始与结束时间
def get_obsfile_ST_ET_quick(rinexFile):
    '''
    get clkfile start_datetime & end_datetime quick
    '''
    obs_file_open = open(rinexFile, 'r+')
    obs_file_line = obs_file_open.readlines()
    obs_file_open.close()
    for line in obs_file_line:
        if '>' in line[0]:
            start_datetime = readRinex3LineTime(line)
            break
    for line in obs_file_line[-200:]:
        if '>' in line[0]:
            end_line = line
    end_datetime = readRinex3LineTime(end_line)
    return start_datetime, end_datetime
#计算文件中的采样间隔
def getIntervalInRinex3(rinexFile):
    rinexFileLineOpen = open(rinexFile, 'r+')
    rinexFileLine = rinexFileLineOpen.readlines()
    rinexFileLineOpen.close()
    readBegin = False
    timeList = []
    for line in rinexFileLine:
        if readBegin:
            if '>' in line:
                nowTime = readRinex3LineTime(line)
                timeList.append(nowTime)
        if 'END OF HEADER' in line:
            readBegin = True

    def get_list(date):
        return date.timestamp()

    timeList = sorted(timeList, key = lambda date: get_list(date))
    interval = 999999999999
    for nowTimeIndex in range(len(timeList[1:])):
        intervalTemp = (timeList[nowTimeIndex + 1] - timeList[nowTimeIndex]).total_seconds()
        if intervalTemp < interval:
            interval = intervalTemp
    return interval
#将 RINEX 2.x 的观测类型转换为 3.x 命名规范
def obsType2to3(gSys, band):
    if band[0] == 'C' or band[0] == 'L' or band[0] == 'S' or band[0] == 'D':
        if band[1] == '1' or band[1] == '2':
            if gSys == 'G':
                if band == 'L2':
                    band = band[:2] + 'P'
                else:
                    band = band[:2] + 'C'
            else:
                band = band[:2] + 'P'
        else:
            if gSys == 'G':
                band = band[:2] + 'I'
            else:
                band = band[:2] + 'P'
    elif band[0] == 'P':
        band = 'C' + band[1] + 'P'
    else:
        return band
    return band