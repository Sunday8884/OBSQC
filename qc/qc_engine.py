#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
GNSS 观测数据质量检查引擎 —— PRN 模式
==================================================

单文件 QC → PRN 级汇总（_QCValue.csv）。

公共骨架（读取观测、执行 9 个 QC 模块、PRN 级聚合）见 qc/qc_common.py，
本模块只负责 prn 模式特有的输出组织。

综合评价（MADM）不在本仓库范围内，由独立项目读取 _QCValue.csv /
QCValue_comparison.csv 完成。
"""
from qc.qc_common import (
    QCProcessorBase,
    parse_time,
    parse_band_string,
    ensure_dir,
)


class QCProcessor(QCProcessorBase):
    """单个观测文件的 QC 处理器，封装了读取、各 QC 模块和结果汇总。"""

    def __init__(self, obs_file, out_dir, args, nav_data=None):
        super().__init__(obs_file, out_dir, args, nav_data=nav_data)


# 兼容旧引用：qc_engine.parse_time / parse_band_string / ensure_dir
__all__ = ['QCProcessor', 'parse_time', 'parse_band_string', 'ensure_dir']
