#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
GNSS 观测数据质量检查引擎
"""
import os
import datetime
import math
import csv
import numpy as np
import matplotlib
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
from plot.plotSatNum import plotSatNum
from plot.plotCnr import plotCnr
from plot.plotCycleSlip import plotCycleSlip
from plot.plotnoise import plotPhaseNoise
from plot.plotMultipath import plotMultipath
from plot.plotCMC import plotCMC
from plot.plotIOD import plotIOD
from plot.plotFreq import plotFreq
from plot.plotHighOrderDiff import plotHighOrderDiff
from madm.EQCValue import evaluate_qcvalue

matplotlib.use('Agg')


# ========================= 辅助函数 =========================
def parse_time(timestr):
    """将时间字符串（YYYY-MM-DD HH:MM:SS）解析为 datetime 对象，空字符串返回 None"""
    if not timestr or timestr.strip() == "":
        return None
    return datetime.datetime.strptime(timestr, '%Y-%m-%d %H:%M:%S')
def parse_band_string(band_str):
    """
    解析波段字符串，如 "G:1C,2W;C:2I,6I" 转换为字典：
    {'G': ['1C','2W'], 'C': ['2I','6I']}
    """
    if not band_str:
        return {}
    band_dict = {}
    for part in band_str.split(';'):
        if ':' not in part:
            continue
        sys, bands = part.split(':', 1)
        band_dict[sys] = [b.strip() for b in bands.split(',') if b.strip()]
    return band_dict
def ensure_dir(path):
    """确保目录存在，不存在则创建"""
    if not os.path.exists(path):
        os.makedirs(path)
# ========================= QC 处理器 =========================
class QCProcessor:
    """单个观测文件的 QC 处理器，封装了读取、各 QC 模块和结果汇总。"""
    # 模块名 -> 方法名（顺序即执行顺序）
    MODULES = (
        'satnum', 'freq', 'cnr', 'slip', 'phasenoise',
        'multipath', 'cmc', 'iod', 'highOrderDiff',
    )
    def __init__(self, obs_file, out_dir, args):
        self.obs_file = obs_file
        self.out_dir = out_dir
        self.args = args
        self.base_name = os.path.splitext(os.path.basename(obs_file))[0]

        # 解析过滤条件
        self.start_dt = parse_time(args.start) if args.start else None
        self.end_dt = parse_time(args.end) if args.end else None
        self.sat_list = args.satellites if args.satellites else None
        self.sys_list = args.systems if args.systems else None
        self.band_choose = parse_band_string(args.bands) if args.bands else {}

        # 观测数据与汇总结果
        self.obs_head = None
        self.obs_data = None
        self.total_epochs = 0
        self.qc_value = {}

    # ---------- 内部工具 ----------
    def _outpath(self, suffix):
        return os.path.join(self.out_dir, f'{self.base_name}{suffix}')

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
            gnssSystem=self.sys_list if self.sys_list else
                       ['G', 'C', 'R', 'E', 'L', 'S', 'J', 'W', 'ALL'],
        )
        if satNumData:
            writeSatNum(satNumData, self._outpath('_sat.txt'))

    def run_freq(self):
        print("-->统计频段使用情况(freq)")
        freqData = plotFreq(
            self.obs_head, self.obs_data,
            pngFile=self._outpath('_freq.png'),
            start_datetime=self.start_dt, end_datetime=self.end_dt,
            sat_list=self.sat_list, sys_list=self.sys_list, band_choose=self.band_choose,
        )
        if freqData:
            d = writeFreqData(freqData, self._outpath('_freq.txt'),
                              total_epochs=self.total_epochs)
            self.qc_value.update(d)

    def run_cnr(self):
        print("-->统计载噪比(CNR)")
        cnrData = plotCnr(
            self.obs_head, self.obs_data,
            pngFile=self._outpath('_cnr.png'),
            start_datetime=self.start_dt, end_datetime=self.end_dt,
            sat_list=self.sat_list, sys_list=self.sys_list, band_choose=self.band_choose,
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
            sat_list=self.sat_list, sys_list=self.sys_list, band_choose=self.band_choose,
        )
        if phaseNoiseData:
            d = writePhaseNoise(phaseNoiseData, self._outpath('_phasenoise.txt'))
            self.qc_value.update(d)
            writePhaseNoiseDetailed(phaseNoiseData, self._outpath('_phasenoise.csv'))
            plotPhaseNoise(phaseNoiseData, pngFile=self._outpath('_phasenoise.png'))

    def run_multipath(self):
        print("-->计算多径效应(multipath)")
        mpData = multipath(
            self.obs_head, self.obs_data,
            start_datetime=self.start_dt, end_datetime=self.end_dt,
            sat_list=self.sat_list, sys_list=self.sys_list, band_choose=self.band_choose,
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
            sat_list=self.sat_list, sys_list=self.sys_list, band_choose=self.band_choose,
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
            sat_list=self.sat_list, sys_list=self.sys_list, band_choose=self.band_choose,
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
            sat_list=self.sat_list, sys_list=self.sys_list, band_choose=self.band_choose,
        )
        if highData:
            d = writeHod(highData, self._outpath('_highOrderDiff.txt'))
            self.qc_value.update(d)
            writeHodDetailed(highData, self._outpath('_highOrderDiff.csv'))
            plotHighOrderDiff(highData, pngFile=self._outpath('_highOrderDiff.png'))

    # ---------- 汇总 ----------
    def write_summary(self):
        if not self.qc_value:
            return
        qc_value_file = self._outpath('_QCValue.csv')

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

        prns_all = sorted(prn_module_val.keys())
        with open(qc_value_file, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(['prn'] + modules_all)
            for prn in prns_all:
                row = [prn] + [prn_module_val[prn].get(m, '') for m in modules_all]
                writer.writerow(row)

        print(f"-->汇总统计已保存至: {qc_value_file}")

    # ---------- 综合评价 ----------
    def run_evaluate(self):
        """读取 _QCValue.csv，调用评估模块做 GRA/TOPSIS/WASPAS/INCOME 评估。
        PCA 参数从命令行 args 读取，由 evaluate_qcvalue 内部执行降维。
        """
        csv_file = self._outpath('_QCValue.csv')
        if not os.path.exists(csv_file):
            print("未找到 _QCValue.csv，跳过评估")
            return

        print("-->综合评价(evaluate)")
        weight_method = getattr(self.args, 'eval_weight', 'entropy')
        eval_method = getattr(self.args, 'eval_method', 'gra')
        rho = getattr(self.args, 'eval_rho', 0.5)
        lamb = getattr(self.args, 'eval_lambda', 0.5)

        # ----  PCA 参数 ----
        use_pca = getattr(self.args, 'pca', False)
        pca_var = getattr(self.args, 'pca_var', 0.85)
        pca_n = getattr(self.args, 'pca_n', None)
        pca_report_file = None
        if use_pca:
            pca_report_file = self._outpath(
                f'_eval_{weight_method}_{eval_method}_pca_report.csv'
            )
        result = evaluate_qcvalue(
            csv_file,
            weight_method=weight_method,
            eval_method=eval_method,
            rho=rho, lamb=lamb,
            use_pca=use_pca,
            pca_var_threshold=pca_var,
            pca_n_components=pca_n,
            pca_report_file=pca_report_file,
        )

        if result is None:
            return

        prns, modules, data_used, weights, scores, ranking = result
        order = np.argsort(-scores)

        # ---- 保存评估结果 CSV（PCA 时列名自动变为 PC1、PC2 ...） ----
        out_file = self._outpath(
            f'_eval_{weight_method}_{eval_method}'
            f'{"_pca" if use_pca else ""}.csv'
        )
        with open(out_file, 'w', newline = '') as f:
            writer = csv.writer(f)
            writer.writerow(['PRN'] + list(modules) + ['Score', 'Rank'])
            for i in order:
                raw_row = []
                for v in data_used[i]:
                    if np.isnan(v):
                        raw_row.append('')
                    elif np.isposinf(v):
                        raw_row.append('inf')
                    elif np.isneginf(v):
                        raw_row.append('-inf')
                    else:
                        raw_row.append(f"{v:.4f}")
                writer.writerow([prns[i]] + raw_row +
                                [f"{scores[i]:.4f}", int(ranking[i])])
            writer.writerow([])
            writer.writerow(['Weight'] + [f"{w:.4f}" for w in weights] + ['', ''])

        # ---- 终端打印 ----
        print(f"\n===== 评估结果: {self.base_name} "
              f"[{weight_method} + {eval_method}"
              f"{' + PCA' if use_pca else ''}] =====")
        print("指标权重:")
        for m, w in zip(modules, weights):
            print(f"    {m:<14s} {w:.4f}")
        print(f"    {'PRN':<10s}{'Score':>12s}{'Rank':>8s}")
        print("    " + "-" * 30)
        for i in order:
            print(f"    {prns[i]:<10s}{scores[i]:>12.4f}{int(ranking[i]):>8d}")
        print(f"-->评估结果已保存至: {out_file}")

    # ---------- 主流程 ----------
    def process(self):
        """执行完整 QC 流程：读数据 → 逐模块计算 → 汇总 → 综合评估"""
        self.load_obs()
        for name in self.MODULES:
            if getattr(self.args, name, False):
                getattr(self, f'run_{name}')()
                plt.close('all')
        self.write_summary()

        # 汇总完成后运行综合评价
        if getattr(self.args, 'evaluate', False):
            self.run_evaluate()