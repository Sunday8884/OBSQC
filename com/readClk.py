#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
readClk.py — 读取 IGS 精密钟差文件（RINEX Clock 格式，*.clk / *.clk_30s），并提供线性内插。

clkData 结构
------------
{
    'epochs' : [datetime, ...]     # 全局历元列表（升序）
    'sat'    : {prn: {epoch: clk}} # 卫星钟差（秒，已含相对论项，按 IGS 约定）
}

说明
----
- 只解析 'AS' 记录（卫星钟差），'AR'（接收机钟差）忽略；
- 文件中每个历元每颗卫星一条记录，可直接线性内插到观测历元。
"""
import datetime
import numpy as np


def _parse_clk_epoch(tokens):
    """tokens: [YYYY MM DD HH MM SS.ssssss ...]（datetime 部分），返回 datetime。"""
    try:
        year = int(tokens[0])
        if year < 100:
            year += 2000
        month = int(tokens[1])
        day = int(tokens[2])
        hour = int(tokens[3])
        minute = int(tokens[4])
        sec = float(tokens[5])
        return datetime.datetime(year, month, day, hour, minute) + datetime.timedelta(seconds=sec)
    except Exception:
        return None


def readRinexClk(clkFile):
    """
    读取 RINEX Clock 文件，返回 clkData 字典（结构见模块 docstring）。
    """
    clkData = {'epochs': [], 'sat': {}}

    with open(clkFile, 'r', encoding='utf-8', errors='ignore') as f:
        for line in f:
            if not line.strip() or line[0] == '*':
                continue
            recType = line[:2]
            if recType not in ('AS',):
                continue
            # 'AS  G01  2024 01 05 00 00 00.000000  5  3  0.123456789012 ...'
            prn = line[3:7].strip()
            if not prn:
                continue
            tokens = line[8:].split()
            if len(tokens) < 8:
                continue
            epoch = _parse_clk_epoch(tokens)
            if epoch is None:
                continue
            # 钟差为历元后第一个科学计数/小数字段
            # （历元与钟差之间可能存在 1~2 个整数字段：数据个数/标志，随文件来源而异）
            clk = None
            for tok in tokens[6:]:
                try:
                    val = float(tok)
                except Exception:
                    continue
                if abs(val) < 1e-2:
                    clk = val
                    break
            if clk is None:
                continue
            clkData['sat'].setdefault(prn, {})[epoch] = clk
            if epoch not in clkData['epochs']:
                clkData['epochs'].append(epoch)

    clkData['epochs'].sort()
    return clkData


def getClkFromRinex(clkData, prn, epoch, maxGap=660.0):
    """
    线性内插 prn 在 epoch 时刻的钟差（秒）。
    maxGap : 允许的最大相邻历元间隔（秒），超出视为无数据。
             默认 660 s，兼容 5min（IGS final）与 15min 采样钟差产品的线性内插。
    无数据时返回 None。
    """
    if clkData is None or prn not in clkData['sat']:
        return None
    satDict = clkData['sat'][prn]
    if not satDict:
        return None

    if epoch in satDict:
        return satDict[epoch]

    epochs = sorted(satDict.keys())
    for i in range(len(epochs) - 1):
        if epochs[i] <= epoch <= epochs[i + 1]:
            gap = (epochs[i + 1] - epochs[i]).total_seconds()
            if gap > maxGap:
                return None
            a = (epoch - epochs[i]).total_seconds() / gap
            return satDict[epochs[i]] * (1.0 - a) + satDict[epochs[i + 1]] * a
    return None


if __name__ == '__main__':
    import sys
    if len(sys.argv) > 1:
        clkData = readRinexClk(sys.argv[1])
        print('历元数:', len(clkData['epochs']), '卫星数:', len(clkData['sat']))
        if clkData['epochs']:
            ep = clkData['epochs'][len(clkData['epochs']) // 2]
            for prn in list(clkData['sat'])[:3]:
                print(prn, ep, getClkFromRinex(clkData, prn, ep))
