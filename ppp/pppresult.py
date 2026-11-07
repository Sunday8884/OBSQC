#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
pppresult.py — PPP 结果写出（与 spos.writePosData 同款列格式，末尾追加 PPP 扩展列）。

输出列：
    date  mjd  sod  X[m] Y[m] Z[m]  E[m] N[m] U[m]  CLK[ns]
    HDOP VDOP PDOP GDOP  NSAT  ZWD[m]  SIG3D[m]  DX[m] DY[m] DZ[m]

说明：前 14 列与 *_spp.txt 完全一致，保证 QC 下游模块（如 qc_time._read_spp_txt）
可直接读取 PPP 结果；NSAT/ZWD/SIG3D 为 PPP 扩展列；
DX/DY/DZ 为坐标改正数（滤波状态相对初始坐标的 ECEF 差值，m）。
"""
import numpy as np

from com.gnssParameter import CLIGHT
from com.xyz2neu import xyz2neu
from com.gnssTime import datetime2allgnssTime


def writePosData(posData, posFile, refX=None, refY=None, refZ=None):
    """
    写出 PPP 逐历元结果。
    posData : {epoch: {'x','y','z','clk','HDOP',...,'NSAT','ZWD','SIG3D'}}
    refX/Y/Z: ENU 参考坐标（None 时用有效历元坐标均值）。
    """
    valid = [ep for ep in posData if 'x' in posData[ep]]
    if not valid:
        print('警告：无有效 PPP 历元，未写出 ' + str(posFile))
        return
    if refX is None or refY is None or refZ is None:
        refX = float(np.mean([posData[ep]['x'] for ep in valid]))
        refY = float(np.mean([posData[ep]['y'] for ep in valid]))
        refZ = float(np.mean([posData[ep]['z'] for ep in valid]))

    Elist, Nlist, Ulist = [], [], []
    with open(posFile, 'w') as f:
        f.write('\n')
        header = ('%25s' % 'date' +
                  '%8s' % 'mjd' +
                  '%11s' % 'sod' +
                  '%20s' % 'X[m]' + '%20s' % 'Y[m]' + '%20s' % 'Z[m]' +
                  '%13s' % 'E[m]' + '%13s' % 'N[m]' + '%13s' % 'U[m]' +
                  '%20s' % 'CLK[ns]' +
                  '%10s' % 'HDOP' + '%10s' % 'VDOP' + '%10s' % 'PDOP' + '%10s' % 'GDOP' +
                  '%8s' % 'NSAT' + '%12s' % 'ZWD[m]' + '%12s' % 'SIG3D[m]' +
                  '%12s' % 'DX[m]' + '%12s' % 'DY[m]' + '%12s' % 'DZ[m]' +
                  '\n')
        f.write(header)

        for epoch in valid:
            rec = posData[epoch]
            N, E, U = xyz2neu(refX, refY, refZ, rec['x'], rec['y'], rec['z'], 'WGS84')
            Elist.append(E)
            Nlist.append(N)
            Ulist.append(U)
            nowGnssTime = datetime2allgnssTime(epoch)
            f.write(
                ' ' + epoch.strftime('%Y-%m-%d %H:%M:%S') +
                '%6d' % nowGnssTime.mjd +
                '%10.2f' % nowGnssTime.sod +
                '%20.4f' % rec['x'] +
                '%20.4f' % rec['y'] +
                '%20.4f' % rec['z'] +
                '%13.4f' % E +
                '%13.4f' % N +
                '%13.4f' % U +
                '%20.5f' % (rec['clk'] / CLIGHT * 1.e9) +
                '%10.3f' % rec['HDOP'] +
                '%10.3f' % rec['VDOP'] +
                '%10.3f' % rec['PDOP'] +
                '%10.3f' % rec['GDOP'] +
                '%8d' % rec.get('NSAT', 0) +
                '%12.4f' % rec.get('ZWD', 0.0) +
                '%12.4f' % rec.get('SIG3D', 0.0) +
                '%12.4f' % rec.get('dx', 0.0) +
                '%12.4f' % rec.get('dy', 0.0) +
                '%12.4f' % rec.get('dz', 0.0) +
                '\n')

    def _rms(v):
        v = np.asarray(v, dtype=np.float64)
        return float(np.sqrt(np.mean(v ** 2))) if len(v) else 0.0

    print('  ENU RMS (ref=%s): E %.4f  N %.4f  U %.4f  [m]' %
          ('%.2f,%.2f,%.2f' % (refX, refY, refZ), _rms(Elist), _rms(Nlist), _rms(Ulist)))
