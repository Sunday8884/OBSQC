#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
PPP — 精密单点定位（Precise Point Positioning）模块包。

公开接口：
    from ppp import ppp
    posData, refX, refY, refZ = ppp(obsHead, obsData, sp3Data, clkData,
                                    bandChoose, startDatetime, endDatetime,
                                    config=cfg_or_json_path, posFile=..., plot_save=...)
"""
from ppp.ppp import ppp
