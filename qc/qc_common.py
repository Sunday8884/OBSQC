#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
GNSS 观测数据质量检查 —— 公共基础设施
==================================================

本模块抽取 prn 模式（qc_engine.py）与 station/Sstation 模式（qc_engine_exp.py）
两个引擎中完全重复的部分，供二者共同继承：

1. 工具函数：parse_time / parse_band_string / ensure_dir
2. QCProcessorBase：单文件 QC 的公共骨架
   - 观测文件读取（load_obs）
   - 9 个 QC 指标模块的执行（run_satnum ... run_highOrderDiff）
   - PRN 级聚合与 _QCValue.csv 输出
3. 波段选择归一化（normalize_band_choose / band_choose_variants）

两个引擎的差异仅在于：
- prn 模式       ：只需 PRN 级 _QCValue.csv
- station 模式   ：额外需要 SPP（ENU 3D RMS、DOP 均值）与测站级 _QCValue_station.csv

【波段选择语义（重要）】
命令行 --bands 使用「频段键」格式，即去掉观测类型前缀的形式：
    --bands "G:1C,2W;C:2I,6I"
其中 '1C' 表示 L1 频点 + C 跟踪模式，与 RINEX 观测码 'L1C' / 'C1C' / 'S1C' 一一对应。

但不同 QC 模块对波段的用法不同，直接透传会导致 --bands 在部分模块失效：
  * freq / CNR   ：按完整观测码匹配（'L1C' / 'S1C'）
  * 其余 5 个模块：按频段键匹配（band[1:]，即 '1C'）
因此这里统一在基类中完成归一化，再按模块需要分发出不同形式，
保证 --bands 对所有模块一致生效。
"""
import os
import math
import datetime
import csv

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from com.readObs import readObsHead, readObs
from qc.satNum import writeSatNum
from qc.freq import writeFreqData
from qc.CNR import writeCnr, writeCnrDetailed
from qc.slip import turboedit, writeSlip, writeSlipDetailed
from qc.phasenoise import phaseNoise, writePhaseNoise, writePhaseNoiseDetailed
from qc.multipath import multipath, writeMp, writeMpDetailed
from qc.CMC import CMC, writeCmc, writeCmcDetailed
from qc.IOD import IOD, writeIon, writeIonDetailed
from qc.highOrderDiff import getObsHighOrderDiff, writeHod, writeHodDetailed
from qc.completeness import (computeCompleteness, writeCompleteness,
                              writeCompletenessDetailed)
from plot.plotSatNum import plotSatNum
from plot.plotCnr import plotCnr
from plot.plotCycleSlip import plotCycleSlip
from plot.plotnoise import plotPhaseNoise
from plot.plotMultipath import plotMultipath
from plot.plotCMC import plotCMC
from plot.plotIOD import plotIOD
from plot.plotFreq import plotFreq
from plot.plotHighOrderDiff import plotHighOrderDiff
from plot.plotCompleteness import (plotCompleteness, plotCompletenessTimeline)

# RINEX 观测类型前缀：用于把频段键还原为完整观测码
OBS_TYPE_PREFIXES = ('L', 'C', 'P', 'D', 'S')

# 各系统缺省启用的 GNSS 系统（--systems 未指定时）
ALL_GNSS_SYSTEMS = ['G', 'C', 'R', 'E', 'L', 'S', 'J', 'W', 'ALL']


# ========================= 工具函数 =========================
def parse_time(timestr):
    """将时间字符串（YYYY-MM-DD HH:MM:SS）解析为 datetime 对象，空字符串返回 None"""
    if not timestr or timestr.strip() == "":
        return None
    return datetime.datetime.strptime(timestr, '%Y-%m-%d %H:%M:%S')


def parse_band_string(band_str):
    """
    解析波段字符串，如 "G:1C,2W;C:2I,6I" 转换为字典：
    {'G': ['1C','2W'], 'C': ['2I','6I']}

    同时兼容用户直接书写完整观测码的形式（如 'L1C'），会自动去掉前缀，
    统一为频段键形式，保证下游各模块语义一致。
    """
    if not band_str:
        return {}
    band_dict = {}
    for part in band_str.split(';'):
        if ':' not in part:
            continue
        sys, bands = part.split(':', 1)
        keys = []
        for b in bands.split(','):
            b = b.strip()
            if not b:
                continue
            # 'L1C' -> '1C'；已是频段键则原样保留
            if len(b) >= 3 and b[0] in OBS_TYPE_PREFIXES:
                b = b[1:]
            if b not in keys:
                keys.append(b)
        if keys:
            band_dict[sys.strip()] = keys
    return band_dict


def ensure_dir(path):
    """确保目录存在，不存在则创建"""
    if not os.path.exists(path):
        os.makedirs(path)


def band_choose_variants(band_choose):
    """
    将归一化的频段键字典分发为各模块需要的形式。

    Parameters
    ----------
    band_choose : dict
        频段键字典，如 {'G': ['1C','2W']}

    Returns
    -------
    (keys, phase_codes, snr_codes) : tuple of dict
        keys       : {'G': ['1C','2W']}      —— 频段键，供 multipath/CMC/IOD/
                                            phasenoise/highOrderDiff 使用
        phase_codes: {'G': ['L1C','L2W']}    —— 完整载波观测码，供 freq 使用
        snr_codes  : {'G': ['S1C','S2W']}    —— 完整载噪比观测码，供 CNR 使用
    """
    keys = band_choose or {}
    phase_codes = {s: ['L' + b for b in bands] for s, bands in keys.items()}
    snr_codes = {s: ['S' + b for b in bands] for s, bands in keys.items()}
    return keys, phase_codes, snr_codes


# ========================= QC 处理器基类 =========================
class QCProcessorBase:
    """
    单个观测文件的 QC 处理器骨架。

    子类负责实现 write_summary()，决定输出哪些汇总表：
      - qc_engine.QCProcessor       → PRN 级 _QCValue.csv
      - qc_engine_exp.QCProcessor   → PRN 级 + 测站级 _QCValue_station.csv
    """

    # 模块名 -> 方法名（顺序即执行顺序）
    MODULES = (
        'satnum', 'freq', 'completeness', 'cnr', 'slip', 'phasenoise',
        'multipath', 'cmc', 'iod', 'highOrderDiff',
    )

    def __init__(self, obs_file, out_dir, args, nav_data=None):
        self.obs_file = obs_file
        self.out_dir = out_dir
        self.args = args
        self.nav_data = nav_data
        self.base_name = os.path.splitext(os.path.basename(obs_file))[0]

        # 解析过滤条件
        self.start_dt = parse_time(args.start) if args.start else None
        self.end_dt = parse_time(args.end) if args.end else None
        self.sat_list = args.satellites if args.satellites else None
        self.sys_list = args.systems if args.systems else None

        # 波段归一化：三种形式供不同模块使用
        self.band_choose = parse_band_string(args.bands) if args.bands else {}
        (self.band_choose_keys,
         self.band_choose_phase,
         self.band_choose_snr) = band_choose_variants(self.band_choose)

        # 观测数据与汇总结果
        self.obs_head = None
        self.obs_data = None
        self.total_epochs = 0
        self.qc_value = {}

    # ---------- 内部工具 ----------
    def _outpath(self, suffix):
        return os.path.join(self.out_dir, f'{self.base_name}{suffix}')

    def _systems(self):
        """--systems 未指定时返回全部系统"""
        return self.sys_list if self.sys_list else ALL_GNSS_SYSTEMS

    # ---------- 读取观测文件 ----------
    def load_obs(self):
        ensure_dir(self.out_dir)
        print(f"-->读取观测文件: {self.obs_file}")
        self.obs_head = readObsHead(self.obs_file, needSatList=True)
        self.obs_data = readObs(self.obs_file, obsHead=self.obs_head)
        self.total_epochs = len(self.obs_data)
        print(f"历元数量: {self.total_epochs}")

    # ---------- 各 QC 模块 ----------
    def run_satnum(self):
        print("-->统计卫星数量(sat)")
        satNumData = plotSatNum(
            self.obs_data,
            pngFile=self._outpath('_sat.png'),
            start_datetime=self.start_dt, end_datetime=self.end_dt,
            gnssSystem=self._systems(),
        )
        if satNumData:
            writeSatNum(satNumData, self._outpath('_sat.txt'))

    def run_freq(self):
        print("-->统计频段使用情况(freq)")
        freqData = plotFreq(
            self.obs_head, self.obs_data,
            pngFile=self._outpath('_freq.png'),
            start_datetime=self.start_dt, end_datetime=self.end_dt,
            sat_list=self.sat_list, sys_list=self.sys_list,
            band_choose=self.band_choose_phase,
        )
        if freqData:
            d = writeFreqData(freqData, self._outpath('_freq.txt'),
                              total_epochs=self.total_epochs)
            self.qc_value.update(d)

    def run_completeness(self):
        """统计卫星数据完整率（弧段口径），并识别连续缺测缺口。

        完整率以「卫星可见弧段」为分母，而非文件总历元数，
        因此不会被卫星本身不可见的时间稀释（详见 qc/completeness.py）。

        若提供了广播星历（station 模式），额外计算「理论可见」口径的完整率。
        """
        print("-->统计数据完整率(completeness)")
        compData = computeCompleteness(
            self.obs_head, self.obs_data,
            start_datetime=self.start_dt, end_datetime=self.end_dt,
            sat_list=self.sat_list, sys_list=self.sys_list,
            band_choose=self.band_choose_keys,
            nav_data=self.nav_data,
            cutoff=getattr(self.args, 'cutoff', 10.0),
        )
        if not compData:
            print("    无可统计的卫星数据")
            return
        d = writeCompleteness(compData, self._outpath('_completeness.txt'),
                              total_epochs=self.total_epochs)
        self.qc_value.update(d)
        gaps = writeCompletenessDetailed(
            compData, self.obs_data, self.obs_head,
            self._outpath('_completeness.csv'),
            band_choose=self.band_choose_keys,
            start_datetime=self.start_dt, end_datetime=self.end_dt,
        )
        plotCompleteness(compData, pngFile=self._outpath('_completeness.png'),
                         gapData=gaps)
        if gaps:
            plotCompletenessTimeline(gaps,
                                     pngFile=self._outpath('_gaps.png'))

    def run_cnr(self):
        print("-->统计载噪比(CNR)")
        cnrData = plotCnr(
            self.obs_head, self.obs_data,
            pngFile=self._outpath('_cnr.png'),
            start_datetime=self.start_dt, end_datetime=self.end_dt,
            sat_list=self.sat_list, sys_list=self.sys_list,
            band_choose=self.band_choose_snr,
        )
        if cnrData:
            d = writeCnr(cnrData, self._outpath('_cnr.txt'))
            self.qc_value.update(d)
            writeCnrDetailed(cnrData, self._outpath('_cnr.csv'))

    def run_slip(self):
        print("-->检测周跳(slip)")
        slipData = turboedit(
            self.obs_head, self.obs_data,
            start_datetime=self.start_dt, end_datetime=self.end_dt,
            sat_list=self.sat_list, sys_list=self.sys_list,
        )
        if slipData:
            d = writeSlip(slipData, self._outpath('_slip.txt'))
            self.qc_value.update(d)
            writeSlipDetailed(slipData, self._outpath('_slip.csv'),
                              self.obs_data, self.start_dt, self.end_dt)
            plotCycleSlip(slipData, pngFile=self._outpath('_slip.png'))

    def run_phasenoise(self):
        print("-->计算相位噪声(phasenoise)")
        phaseNoiseData = phaseNoise(
            self.obs_head, self.obs_data,
            start_datetime=self.start_dt, end_datetime=self.end_dt,
            sat_list=self.sat_list, sys_list=self.sys_list,
            band_choose=self.band_choose_keys,
        )
        if phaseNoiseData:
            d = writePhaseNoise(phaseNoiseData, self._outpath('_phasenoise.txt'))
            self.qc_value.update(d)
            writePhaseNoiseDetailed(phaseNoiseData, self._outpath('_phasenoise.csv'))
            plotPhaseNoise(phaseNoiseData, pngFile=self._outpath('_phasenoise.png'))

    def run_multipath(self):
        print("-->计算多路径效应(multipath)")
        mpData = multipath(
            self.obs_head, self.obs_data,
            start_datetime=self.start_dt, end_datetime=self.end_dt,
            sat_list=self.sat_list, sys_list=self.sys_list,
            band_choose=self.band_choose_keys,
        )
        if mpData:
            d = writeMp(mpData, self._outpath('_multipath.txt'))
            self.qc_value.update(d)
            writeMpDetailed(mpData, self._outpath('_multipath.csv'))
            plotMultipath(mpData, pngFile=self._outpath('_multipath.png'))

    def run_cmc(self):
        print("-->计算历元间差分码减相位(CMC)...")
        cmcData = CMC(
            self.obs_head, self.obs_data,
            start_datetime=self.start_dt, end_datetime=self.end_dt,
            sat_list=self.sat_list, sys_list=self.sys_list,
            band_choose=self.band_choose_keys,
        )
        if cmcData:
            d = writeCmc(cmcData, self._outpath('_cmc.txt'))
            self.qc_value.update(d)
            writeCmcDetailed(cmcData, self._outpath('_cmc.csv'))
            plotCMC(cmcData, pngFile=self._outpath('_cmc.png'))

    def run_iod(self):
        print("-->计算电离层延迟变化率(IOD)")
        iodData = IOD(
            self.obs_head, self.obs_data,
            start_datetime=self.start_dt, end_datetime=self.end_dt,
            sat_list=self.sat_list, sys_list=self.sys_list,
            band_choose=self.band_choose_keys,
        )
        if iodData:
            d = writeIon(iodData, self._outpath('_iod.txt'))
            self.qc_value.update(d)
            writeIonDetailed(iodData, self._outpath('_iod.csv'))
            plotIOD(iodData, pngFile=self._outpath('_iod.png'))

    def run_highOrderDiff(self):
        print("-->计算载波相位高阶差分(highOrderDiff)")
        highData = getObsHighOrderDiff(
            self.obs_head, self.obs_data,
            diffNum=4,
            start_datetime=self.start_dt, end_datetime=self.end_dt,
            sat_list=self.sat_list, sys_list=self.sys_list,
            band_choose=self.band_choose_keys,
        )
        if highData:
            d = writeHod(highData, self._outpath('_highOrderDiff.txt'))
            self.qc_value.update(d)
            writeHodDetailed(highData, self._outpath('_highOrderDiff.csv'))
            plotHighOrderDiff(highData, pngFile=self._outpath('_highOrderDiff.png'))

    # ---------- 聚合：PRN 级 ----------
    def _aggregate_prn_module_val(self):
        """
        把各模块返回的 {模块: {PRN: {波段: {val}}}} 聚合为 {PRN: {模块: 值}}。
        同一 PRN 跨波段取均值；无法转为数值的保留原字符串。
        """
        prn_module_val = {}
        modules_all = []
        for module, prn_dict in self.qc_value.items():
            if module not in modules_all:
                modules_all.append(module)
            for prn, band_dict in prn_dict.items():
                nums = []
                str_val = None
                for _band, val_dict in band_dict.items():
                    val = val_dict['val']
                    if isinstance(val, str):
                        try:
                            nums.append(float(val))
                        except ValueError:
                            str_val = val
                    else:
                        nums.append(float(val))

                if nums:
                    avg = sum(nums) / len(nums)
                    if math.isinf(avg):
                        val_str = 'inf' if avg > 0 else '-inf'
                    else:
                        val_str = f"{avg:.4f}"
                elif str_val is not None:
                    val_str = str_val
                else:
                    val_str = ''

                prn_module_val.setdefault(prn, {})[module] = val_str
        return prn_module_val, modules_all

    def _write_prn_summary(self, prn_module_val, modules_all):
        """写出 PRN 级明细 _QCValue.csv"""
        detail_file = self._outpath('_QCValue.csv')
        prns_all = sorted(prn_module_val.keys())
        with open(detail_file, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(['prn'] + modules_all)
            for prn in prns_all:
                writer.writerow([prn] +
                                [prn_module_val[prn].get(m, '') for m in modules_all])
        return detail_file

    # ---------- 主流程 ----------
    def write_summary(self):
        """写出 PRN 级汇总表（子类可扩展）"""
        if not self.qc_value:
            return None
        prn_module_val, modules_all = self._aggregate_prn_module_val()
        if not prn_module_val:
            return None
        detail_file = self._write_prn_summary(prn_module_val, modules_all)
        print(f"-->PRN级统计值: {detail_file}")
        return prn_module_val, modules_all

    def run_modules(self):
        """按开关顺序执行各 QC 模块，每模块结束释放绘图内存"""
        for name in self.MODULES:
            if getattr(self.args, name, False):
                getattr(self, f'run_{name}')()
                plt.close('all')

    def process(self):
        """执行完整 QC 流程：读数据 → 逐模块计算 → 汇总"""
        self.load_obs()
        self.run_modules()
        return self.write_summary()
