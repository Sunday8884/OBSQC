#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
satposSp3.py — 基于精密星历（SP3 + 可选精密钟差）计算卫星位置与钟差。
接口与 spos/satpos.py 保持一致：给定接收历元，返回该历元所有可用卫星的位置、速度与钟差。

改正说明
--------
- 周期相对论钟差：IGS 精密钟差按约定**不含**周期相对论项，
  此处按 dtrel = −2·(r⃗·v⃗)/c² 补回（等价于 −2√(GM·a)·e·sinE/c²，最大 ~±2e-7 s ≈ ±7 m）；
  广播星历路径（getClkInNav）已含该项，勿重复施加。
- 速度：由 SP3 内插节点直接求 Lagrange 导数。
"""
import datetime
import numpy as np

from com.gnssParameter import CLIGHT
from com.readSp3 import getSp3PosClk
from com.readClk import getClkFromRinex


def satPosSp3(bandChoose, obsData, sp3Data, clkData, epoch, config=None):
    """
    计算接收历元 epoch 各卫星的发射时刻位置/速度/钟差。

    bandChoose : 双频配置 {'G': {'C1C': {...}, 'C2W': {...}}}（需已填入 freq/lambda）
    obsData    : {epoch: {prn: {band: value}}}
    sp3Data    : readSp3.readSp3 的返回
    clkData    : readClk.readRinexClk 的返回（可为 None，此时用 SP3 内嵌钟差）
    epoch      : 接收时刻（datetime）
    config     : dict（用到 sp3_pos_order）

    返回 (satPosInEpoch, satClkInEpoch, satVarInEpoch)：
        satPosInEpoch = {prn: np.array([x,y,z,vx,vy,vz])}
        satClkInEpoch = {prn: clk(s)}            # 已含周期相对论项
        satVarInEpoch = {prn: var(m^2)}          # 卫星轨道+钟差方差（SP3 精度信息，
                                                 #  posStd²+(c·clkStd)²，未知为 0），
                                                 #  供量测 R 加权（RTKLIB var_rs 同源）
    """
    if config is None:
        config = {}
    posOrder = int(config.get('sp3_pos_order', 10))

    satPosInEpoch = {}
    satClkInEpoch = {}
    satVarInEpoch = {}

    for prn in obsData[epoch]:
        gnssSys = prn[0]
        if gnssSys not in bandChoose:
            continue
        bands = list(bandChoose[gnssSys].keys())[:2]
        if len(bands) < 2:
            continue
        C1 = obsData[epoch][prn].get(bands[0])
        C2 = obsData[epoch][prn].get(bands[1])
        if C1 is None or C2 is None:
            continue
        if abs(C1 - C2) > 100:  # 粗差剔除
            continue

        # 发射时间：t_tx = tag − pr/c − dts_sat（pr 含 c·dtr − c·dts，接收机钟差抵消；
        # 先用粗时刻查钟差，再迭代一次，RTKLIB satposs 同源思想）
        aveP = (C1 + C2) / 2.0
        epochFind = epoch - datetime.timedelta(seconds=aveP / CLIGHT)

        pos, clk, vel, posStd, clkStd = getSp3PosClk(sp3Data, prn, epochFind,
                                                     posOrder=posOrder,
                                                     clkInterp=config.get('sp3_clk_interp', 'linear'),
                                                     wantStd=True)
        if pos is None or clk is None:
            continue
        # 一次迭代：扣除卫星钟差对应的发射时刻偏差（|dts|可达 ±ms → 位置误差 ~m 级）
        epochFind = epoch - datetime.timedelta(seconds=aveP / CLIGHT + clk)
        pos, clk, vel, posStd, clkStd = getSp3PosClk(sp3Data, prn, epochFind,
                                                     posOrder=posOrder,
                                                     clkInterp=config.get('sp3_clk_interp', 'linear'),
                                                     wantStd=True)
        if pos is None or clk is None:
            continue

        # 周期相对论钟差改正（IGS 精密钟差约定不含该项；若钟差产品已含，配置 sp3_apply_rel=false）
        # dtrel = −2·(r⃗·v⃗)/c² ≈ −2√(GM·a)·e·sinE/c²
        dtrel = 0.0
        if config.get('sp3_apply_rel', True):
            dtrel = -2.0 * float(np.dot(pos, vel)) / (CLIGHT ** 2)
        clk += dtrel

        # 精密钟差文件优先（同样补相对论项）
        if clkData is not None:
            clk2 = getClkFromRinex(clkData, prn, epochFind)
            if clk2 is not None:
                clk = clk2 + dtrel

        satPosInEpoch[prn] = np.array([pos[0], pos[1], pos[2], vel[0], vel[1], vel[2]])
        satClkInEpoch[prn] = clk
        # 卫星轨道/钟差方差（量测 R 加权）
        satVarInEpoch[prn] = posStd ** 2 + (CLIGHT * clkStd) ** 2

    return satPosInEpoch, satClkInEpoch, satVarInEpoch
