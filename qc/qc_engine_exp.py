#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
GNSS 观测数据站间质量检查引擎
- 每个文件保留 PRN 级明细
- 每个文件聚合成station QCValue
- station 模式下额外计算 SPP：ENU3DRMS、DOP 平均值
- 所有文件汇总为一张多站对比表 QCValue_comparison.csv
- 可选跨站综合评价（GRA/TOPSIS/WASPAS/INCOME）
"""
import os
import re
import datetime
import math
import csv
import glob
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

from com.readObs import readObsHead, readObs
from com.readNav import readNav
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
from madm.EQCValue import  evaluate_matrix
from madm.KKNN import evaluate_KKNN


matplotlib.use('Agg')

# SPP 默认波段
DEFAULT_SPP_BANDS = {
    'G': ['C1C', 'C2W'],
    'C': ['C2I', 'C6I'],
    'E': ['C1X', 'C5X'],
}

# ========================= 辅助函数 =========================
def parse_time(timestr):
    if not timestr or timestr.strip() == "":
        return None
    return datetime.datetime.strptime(timestr, '%Y-%m-%d %H:%M:%S')
def parse_band_string(band_str):
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
    if not os.path.exists(path):
        os.makedirs(path)
def _extract_date_key(filename):
    """从文件名提取 '年份_日期编码' 组合键，用于匹配观测与导航文件。"""
    basename = os.path.basename(filename)
    match = re.search(r'(\d{4})\.(\d{2})[opn]', basename)
    if match:
        date_code = match.group(1)
        year_code = match.group(2)
        return f"{year_code}_{date_code}"
    return None
def _extract_date_label(filename):
    """从 RINEX 文件名提取可排序的日期标签，如 'wuh20050.24o' -> '2024-005'。
    用于 Sstation 模式构建单站时间序列。无法解析时返回 None。"""
    basename = os.path.basename(filename)
    match = re.search(r'(\d{3})\d\.(\d{2})[opn]', basename)
    if match:
        doy = match.group(1)
        year = int(match.group(2))
        year += 2000 if year < 80 else 1900
        return f"{year}-{doy}"
    return None
def _parse_snx_date(filename):
    """从 IGS 长命名日解文件提取 (year, doy)。
    例: 'IGS0OPSSNX_20240070000_01D_01D_SOL.SNX' -> (2024, 7)"""
    basename = os.path.basename(filename)
    m = re.search(r'_(\d{4})(\d{3})', basename)
    if not m:
        return None
    year, doy = int(m.group(1)), int(m.group(2))
    return (year, doy) if 1 <= doy <= 366 else None

def _resolve_snx_file(snx_pattern, obs_file=None):
    """对 --snx 参数做通配符展开并返回匹配文件。
    匹配到多个文件时，按年积日与观测文件对应（IGS 日解命名 '_YYYYDDD'，
    如 IGS0OPSSNX_20240070000_01D_01D_SOL.SNX）；匹配不到或无法解析日期时
    返回 None，由上层回退到 SPP 坐标均值。
    """
    if not snx_pattern:
        return None
    matches = glob.glob(snx_pattern)
    if not matches:
        # 大小写不敏感回退（Windows/MINGW 常见）
        import fnmatch
        d = os.path.dirname(snx_pattern) or '.'
        pat = os.path.basename(snx_pattern)
        if os.path.isdir(d):
            for f in os.listdir(d):
                if fnmatch.fnmatch(f.lower(), pat.lower()):
                    matches.append(os.path.join(d, f))
    if not matches:
        print(f"警告：--snx 未匹配到任何文件（pattern={snx_pattern}）")
        return None

    if len(matches) == 1:
        return matches[0]

    # 多文件：必须按年积日精确对应，匹配不到则放弃 SNX 参考
    if obs_file:
        label = _extract_date_label(obs_file)   # '2024-310'
        if label:
            year, doy = (int(x) for x in label.split('-'))
            for f in matches:
                if _parse_snx_date(f) == (year, doy):
                    return f
            print(f"警告：{os.path.basename(obs_file)} ({label}) "
                  f"未找到对应 SNX 日解文件，跳过 SNX 参考（回退坐标均值）")
        else:
            print(f"警告：{os.path.basename(obs_file)} 无法解析日期，"
                  f"跳过 SNX 参考（回退坐标均值）")
    else:
        print(f"警告：--snx 匹配到多个文件但未提供观测文件，跳过 SNX 参考")
    return None
def _extract_station_from_filename(obs_file):
    """从观测文件名前 4 位提取测站名（大写）。"""
    basename = os.path.basename(obs_file)
    if len(basename) < 4:
        return None
    station = basename[:4].upper()
    if not station.isalnum():
        return None
    return station
# ========================= QC 处理器 =========================
class QCProcessor:
    """单个观测文件的 QC 处理器。输出 station 级一行的 QCValue。
    station 模式下额外计算 SPP 的 ENU 3D RMS 与 DOP 平均值。
    """
    MODULES = (
        'satnum', 'freq', 'cnr', 'slip', 'phasenoise',
        'multipath', 'cmc', 'iod', 'highOrderDiff',
    )

    def __init__(self, obs_file, out_dir, args, nav_data=None):
        self.obs_file = obs_file
        self.out_dir = out_dir
        self.args = args
        self.nav_data = nav_data
        self.base_name = os.path.splitext(os.path.basename(obs_file))[0]

        self.start_dt = parse_time(args.start) if args.start else None
        self.end_dt = parse_time(args.end) if args.end else None
        self.sat_list = args.satellites if args.satellites else None
        self.sys_list = args.systems if args.systems else None
        self.band_choose = parse_band_string(args.bands) if args.bands else {}

        self.obs_head = None
        self.obs_data = None
        self.total_epochs = 0
        self.qc_value = {}

    # ---------- 内部工具 ----------
    def _outpath(self, suffix):
        return os.path.join(self.out_dir, f'{self.base_name}{suffix}')

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
            sat_list=self.sat_list, sys_list=self.sys_list,
            band_choose=self.band_choose,
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
            sat_list=self.sat_list, sys_list=self.sys_list,
            band_choose=self.band_choose,
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
            band_choose=self.band_choose,
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
            sat_list=self.sat_list, sys_list=self.sys_list,
            band_choose=self.band_choose,
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
            band_choose=self.band_choose,
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
            band_choose=self.band_choose,
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
            band_choose=self.band_choose,
        )
        if highData:
            d = writeHod(highData, self._outpath('_highOrderDiff.txt'))
            self.qc_value.update(d)
            writeHodDetailed(highData, self._outpath('_highOrderDiff.csv'))
            plotHighOrderDiff(highData, pngFile=self._outpath('_highOrderDiff.png'))

    # ---------- SPP 模块 ----------
    def _build_spp_band_choose(self):
        """构建 spp() 所需的 bandChoose 字典 {sys: {band: {}}}。"""
        if self.args.bands:
            band_dict = {}
            for part in self.args.bands.split(';'):
                if ':' not in part:
                    continue
                sys, bands = part.split(':', 1)
                band_dict[sys] = {b.strip(): {} for b in bands.split(',') if b.strip()}
            if band_dict:
                return band_dict
        # 使用默认波段
        systems = self.sys_list if self.sys_list else list(DEFAULT_SPP_BANDS.keys())
        return {s: {c: {} for c in DEFAULT_SPP_BANDS[s]}
                for s in systems if s in DEFAULT_SPP_BANDS}

    @staticmethod
    def _extract_spp_metrics(posData, XMEAN, YMEAN, ZMEAN):
        """
        从 spp() 返回的 posData（dict: {epoch: {'x','y','z','clk','HDOP','VDOP','PDOP','GDOP'}}）
        中提取：
          - ENU 3D RMS（相对 XMEAN/YMEAN/ZMEAN 的 ENU 残差）
          - DOP 平均值（优先 PDOP，回退 GDOP）
        """
        if not posData or XMEAN is None or YMEAN is None or ZMEAN is None:
            return None, None

        from com.xyz2neu import xyz2neu

        e_vals, n_vals, u_vals, dop_vals = [], [], [], []
        dop_name_used = None

        for epoch, rec in posData.items():
            try:
                x = rec.get('x')
                y = rec.get('y')
                z = rec.get('z')
                if x is None or y is None or z is None:
                    continue

                # xyz2neu 返回 (N, E, U)
                N, E, U = xyz2neu(XMEAN, YMEAN, ZMEAN,
                                  float(x), float(y), float(z), 'WGS84')
                if np.isfinite([E, N, U]).all():
                    e_vals.append(float(E))
                    n_vals.append(float(N))
                    u_vals.append(float(U))

                # DOP：优先 PDOP，回退 GDOP，再回退 HDOP/VDOP 组合
                dop = rec.get('PDOP')
                if dop is None:
                    dop = rec.get('GDOP')
                if dop is None:
                    hd = rec.get('HDOP')
                    vd = rec.get('VDOP')
                    if hd is not None and vd is not None:
                        dop = float(hd) + float(vd)
                if dop is not None and np.isfinite(dop):
                    dop_vals.append(float(dop))
                    if dop_name_used is None:
                        dop_name_used = 'PDOP' if 'PDOP' in rec else (
                            'GDOP' if 'GDOP' in rec else 'H+V')
            except Exception:
                continue

        rms_3d = None
        if e_vals:
            arr_e = np.asarray(e_vals, dtype=float)
            arr_n = np.asarray(n_vals, dtype=float)
            arr_u = np.asarray(u_vals, dtype=float)
            rms_e = float(np.sqrt(np.mean(arr_e ** 2)))
            rms_n = float(np.sqrt(np.mean(arr_n ** 2)))
            rms_u = float(np.sqrt(np.mean(arr_u ** 2)))
            rms_3d = float(np.sqrt(rms_e ** 2 + rms_n ** 2 + rms_u ** 2))

        dop_mean = float(np.mean(dop_vals)) if dop_vals else None
        return rms_3d, dop_mean

    def run_spp(self):
        """station 模式：计算本文件的 SPP，并将 ENU3DRMS、DOP 均值加入 qc_value。

        ENU 参考：
          - 优先使用 --snx 提供的 SINEX 周解坐标（若可解析且能找到对应测站）；
          - 否则回退到 SPP 坐标均值 XMEAN/YMEAN/ZMEAN。
        """
        if self.nav_data is None:
            print("--> 跳过 SPP：未提供导航数据")
            return

        try:
            from spos.spp import spp as _spp, read_snx_coord
        except ImportError as e:
            print(f"--> 跳过 SPP：无法导入 spos.spp ({e})")
            return

        print("-->计算SPP")
        band_choose = self._build_spp_band_choose()
        if not band_choose:
            print("--> 跳过 SPP：无可用波段")
            return

        pos_file = self._outpath('_spp.txt')
        plot_file = self._outpath('_spp.png')

        # === 自动从文件名识别测站 & 解析 SNX ===
        station = _extract_station_from_filename(self.obs_file)
        snx_file = _resolve_snx_file(getattr(self.args, 'snx', None), obs_file = self.obs_file)

        try:
            result = _spp(
                self.obs_head, self.obs_data, self.nav_data, band_choose,
                startDatetime=self.start_dt, endDatetime=self.end_dt,
                posFile=pos_file, plot_save=plot_file,
                snx_file=snx_file, station=station,
            )
        except Exception as e:
            print(f"--> SPP 计算失败: {e}")
            return

        if result is None:
            print("--> SPP 失败：无返回")
            return
        posData, XMEAN, YMEAN, ZMEAN = result
        if posData is None or XMEAN is None:
            print("--> SPP 失败：有效历元不足")
            return

        # === 计算 ENU 3D RMS 时使用与 pos 文件相同的参考原点 ===
        refX, refY, refZ = XMEAN, YMEAN, ZMEAN
        if snx_file and station:
            try:
                sx, sy, sz = read_snx_coord(snx_file, station)
                refX, refY, refZ = sx, sy, sz
            except Exception as e:
                print(f"    警告：读取 SNX 失败（{e}），改用坐标均值作为 ENU 参考。")
                refX, refY, refZ = XMEAN, YMEAN, ZMEAN

        rms_3d, dop_mean = self._extract_spp_metrics(posData, refX, refY, refZ)

        # SPP 是整站级别的量，这里用伪 PRN 'STATION' 承载
        if rms_3d is not None:
            self.qc_value.setdefault('spp3drms', {})['STATION'] = {
                '1': {'val': rms_3d}
            }
        else:
            print("    ENU3DRMS 无法提取")

        if dop_mean is not None:
            self.qc_value.setdefault('dop', {})['STATION'] = {
                '1': {'val': dop_mean}
            }
        else:
            print("    DOP 平均值 无法提取")

    # ---------- 内部聚合：PRN 级 ----------
    def _aggregate_prn_module_val(self):
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

    # ---------- 汇总：station 级 ----------
    def get_station_qcvalue(self):
        if not self.qc_value:
            return None

        prn_module_val, modules_all = self._aggregate_prn_module_val()
        if not prn_module_val:
            return None

        station_vals = {}
        for module in modules_all:
            nums = []
            for _prn, mv in prn_module_val.items():
                v = mv.get(module, '')
                if v == '' or v is None:
                    continue
                try:
                    nums.append(float(v))
                except (ValueError, TypeError):
                    pass
            if nums:
                avg = sum(nums) / len(nums)
                if math.isinf(avg):
                    station_vals[module] = 'inf' if avg > 0 else '-inf'
                else:
                    station_vals[module] = f"{avg:.4f}"
            else:
                station_vals[module] = ''

        station = os.path.basename(self.obs_file)[:4]
        return station, station_vals, modules_all

    def write_summary(self):
        if not self.qc_value:
            return None

        prn_module_val, modules_all = self._aggregate_prn_module_val()
        if not prn_module_val:
            return None

        # ---- 1) PRN 级明细 ----
        detail_file = self._outpath('_QCValue.csv')
        prns_all = sorted(prn_module_val.keys())
        with open(detail_file, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(['prn'] + modules_all)
            for prn in prns_all:
                writer.writerow([prn] +
                                [prn_module_val[prn].get(m, '') for m in modules_all])

        # ---- 2) station 级单行 ----
        result = self.get_station_qcvalue()
        if result is None:
            return None
        station, station_vals, _ = result

        qc_value_file = self._outpath('_QCValue_station.csv')
        with open(qc_value_file, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(['station'] + modules_all)
            writer.writerow([station] +
                            [station_vals.get(m, '') for m in modules_all])

        print(f"-->PRN级统计值: {detail_file}")
        print(f"-->测站级统计值: {qc_value_file}")
        return station, station_vals, modules_all

    # ---------- 主流程 ----------
    def process(self):
        """执行完整 QC 流程，返回 (station, station_vals, modules_all) 或 None"""
        self.load_obs()
        for name in self.MODULES:
            if getattr(self.args, name, False):
                getattr(self, f'run_{name}')()
                plt.close('all')

        # ---- SPP 仅在 station/Sstation 模式 + 开关开启时计算 ----
        if (getattr(self.args, 'mode', 'prn') in ('station', 'Sstation')
                and getattr(self.args, 'spp', True)):
            self.run_spp()
            plt.close('all')

        return self.write_summary()


# ========================= 多文件对比 =========================
class MultiFileQC:
    """处理多个观测文件，汇总为一张 station 级对比表并做跨站综合评价。"""

    def __init__(self, args):
        self.args = args
        self.station_rows = []
        self.modules_all = []
        self._nav_cache = {}
        # Sstation：单站多日期文件，行标签为日期，形成时间序列
        self.is_sstation = (getattr(args, 'mode', 'prn') == 'Sstation')
        self.station_name = None

    # ---------- 判断 QC 是否全关 ----------
    def _all_modules_disabled(self):
        """所有 QC 模块关闭，且（非 station 模式 / 无导航数据 / SPP 关闭）时返回 True。"""
        for name in QCProcessor.MODULES:
            if getattr(self.args, name, False):
                return False
        # station/Sstation 模式 + 有导航数据 + SPP 开启 => SPP 会产出指标，不能走快速通道
        if getattr(self.args, 'mode', 'prn') in ('station', 'Sstation'):
            if (getattr(self.args, 'nav', None)
                    and getattr(self.args, 'spp', True)):
                return False
        return True

    # ---------- 导航数据加载 ----------
    def _load_nav_file(self, nav_file):
        if nav_file in self._nav_cache:
            return self._nav_cache[nav_file]
        try:
            data = readNav(nav_file)
        except Exception as e:
            print(f"警告：读取导航文件失败 {nav_file}: {e}")
            data = None
        self._nav_cache[nav_file] = data
        return data

    def _load_nav_for_obs(self, obs_file):
        """为给定观测文件寻找匹配的导航数据。"""
        nav_pattern = getattr(self.args, 'nav', None)
        if not nav_pattern:
            return None
        import glob as _glob
        nav_files = _glob.glob(nav_pattern)
        if not nav_files:
            print(f"警告：未找到导航文件 '{nav_pattern}'，跳过 SPP")
            return None

        # 单文件
        if len(nav_files) == 1:
            return self._load_nav_file(nav_files[0])

        # 多文件：按日期编码匹配
        obs_key = _extract_date_key(obs_file)
        if obs_key:
            for nf in nav_files:
                if _extract_date_key(nf) == obs_key:
                    return self._load_nav_file(nf)
        # 回退到第一个
        return self._load_nav_file(nav_files[0])

    # ---------- 复用已有对比表 ----------
    def _load_existing_comparison(self, csv_file):
        if not os.path.exists(csv_file):
            return False

        print(f"-->复用已有对比表: {csv_file}")
        with open(csv_file, 'r', newline='') as f:
            reader = csv.reader(f)
            try:
                header = next(reader)
            except StopIteration:
                return False
            if not header or len(header) < 2:
                return False

            # Sstation 时间序列表含 station+date 两列标签
            has_date = (len(header) >= 2
                        and header[0].strip().lower() == 'station'
                        and header[1].strip().lower() == 'date')
            modules = header[2:] if has_date else header[1:]
            rows = []
            for row in reader:
                if not row or not row[0].strip():
                    continue
                if row[0].strip().lower() in ('weight', '权重'):
                    continue
                if has_date:
                    self.station_name = row[0].strip()
                    label = row[1].strip()
                    vals = {m: v.strip() for m, v in zip(modules, row[2:])}
                else:
                    label = row[0].strip()
                    vals = {m: v.strip() for m, v in zip(modules, row[1:])}
                rows.append((label, vals))

        if not rows:
            return False

        self.modules_all = list(modules)
        self.station_rows = rows
        print(f"-->成功读取 {len(rows)} 行、{len(modules)} 个指标")
        return True

    # ---------- 结果收集 ----------
    def _collect_result(self, obs_file, station, station_vals, modules):
        """将单文件测站级结果计入对比表（Sstation 下行标签为日期）。"""
        if self.is_sstation:
            # Sstation：校验单站一致性，行标签使用日期
            if self.station_name and self.station_name != station:
                print(f"   警告：检测到不同测站 {station}（已有 {self.station_name}），"
                      "Sstation 模式建议仅处理同一测站的多日期文件")
            self.station_name = self.station_name or station
            label = _extract_date_label(obs_file) or station
        else:
            label = station
        for m in modules:
            if m not in self.modules_all:
                self.modules_all.append(m)
        self.station_rows.append((label, station_vals))

    # ---------- 断点续跑：复用单文件汇总结果 ----------
    def _reuse_file_result(self, obs_file, sub_dir):
        """若该文件的 _QCValue_station.csv 已存在，直接读取其结果并返回 True（跳过计算）。"""
        base = os.path.splitext(os.path.basename(obs_file))[0]
        summary = os.path.join(sub_dir, f'{base}_QCValue_station.csv')
        if not os.path.exists(summary):
            return False
        try:
            with open(summary, 'r', newline='') as f:
                reader = csv.reader(f)
                header = next(reader)
                row = next(reader)
        except Exception as e:
            print(f"   (已有汇总文件读取失败，将重新计算: {e})")
            return False
        if not header or not row:
            return False
        station = row[0].strip()
        modules = [m.strip() for m in header[1:]]
        station_vals = {m: v.strip() for m, v in zip(modules, row[1:])}
        self._collect_result(obs_file, station, station_vals, modules)
        print("   (检测到已有汇总结果，跳过重复计算)")
        return True

    def run(self, obs_files):
        comparison_file = os.path.join(self.args.output_dir, 'QCValue_comparison.csv')
        is_station = (getattr(self.args, 'mode', 'prn') in ('station', 'Sstation'))

        # ============ 快速通道：QC 全关 + 已有对比表 ============
        if self._all_modules_disabled():
            if self._load_existing_comparison(comparison_file):
                if getattr(self.args, 'evaluate', False):
                    self.run_evaluate()
                return
            else:
                print(f"--> QC模块关闭，但未找到有效的 {comparison_file}")
                print("   请先运行一次完整 QC，或启用至少一个 QC 模块")
                return

        # ============ 常规流程：逐文件跑 QC ============
        total = len(obs_files)
        for idx, obs_file in enumerate(obs_files, 1):
            sub_dir = os.path.join(
                self.args.output_dir,
                os.path.splitext(os.path.basename(obs_file))[0]
            )
            print(f"\n========== 处理文件 [{idx}/{total}]: {obs_file} ==========")

            # 断点续跑：已有汇总结果的文件直接复用，跳过重复计算
            if self._reuse_file_result(obs_file, sub_dir):
                continue

            # station 模式下为每个观测文件加载对应的导航数据（SPP 开启时）
            nav_data = None
            if is_station and getattr(self.args, 'spp', True):
                nav_data = self._load_nav_for_obs(obs_file)
                if nav_data is None:
                    print("   (无导航数据，SPP 将被跳过)")
            elif is_station:
                print("   (SPP 已禁用，跳过导航数据加载)")

            proc = QCProcessor(obs_file, sub_dir, self.args, nav_data=nav_data)
            try:
                result = proc.process()
            except Exception as e:
                print(f"文件 {obs_file} 处理失败，已跳过: {e}")
                continue
            if result is None:
                print(f"文件 {obs_file} 未生成有效 QCValue")
                continue
            station, station_vals, modules = result
            self._collect_result(obs_file, station, station_vals, modules)

        # Sstation：按日期排序，形成时间序列
        if self.is_sstation:
            self.station_rows.sort(key=lambda r: r[0])

        # ---- 汇总多站对比 CSV ----
        self.write_comparison()

        # ---- 跨站综合评价 ----
        if getattr(self.args, 'evaluate', False):
            self.run_evaluate()

        print(f"\n所有文件处理完成，结果保存在: {self.args.output_dir}")

    # ---------- 多站对比表 ----------
    def write_comparison(self):
        if not self.station_rows:
            print("无有效 station 数据，跳过对比表输出")
            return
        out_file = os.path.join(self.args.output_dir, 'QCValue_comparison.csv')
        ensure_dir(self.args.output_dir)
        with open(out_file, 'w', newline='') as f:
            writer = csv.writer(f)
            if self.is_sstation:
                # 单站时间序列：station + date + 指标
                writer.writerow(['station', 'date'] + self.modules_all)
                for label, vals in self.station_rows:
                    writer.writerow([self.station_name, label] +
                                    [vals.get(m, '') for m in self.modules_all])
            else:
                writer.writerow(['station'] + self.modules_all)
                for station, vals in self.station_rows:
                    writer.writerow([station] +
                                    [vals.get(m, '') for m in self.modules_all])
        if self.is_sstation:
            print(f"-->单站时间序列质量统计表已保存至: {out_file}")
        else:
            print(f"-->多站对比表已保存至: {out_file}")

    # ---------- 跨站综合评价 ----------
    def run_evaluate(self):
        if len(self.station_rows) < 2:
            unit = '日期' if self.is_sstation else 'station'
            print(f" {unit} 数量少于 2，跳过综合评价")
            return

        print("\n-->跨日期综合评价(evaluate)" if self.is_sstation
              else "\n-->跨站综合评价(evaluate)")
        weight_method = getattr(self.args, 'eval_weight', 'entropy')
        eval_method = getattr(self.args, 'eval_method', 'gra')
        rho = getattr(self.args, 'eval_rho', 0.5)
        lamb = getattr(self.args, 'eval_lambda', 0.5)
        use_pca = getattr(self.args, 'pca', False)
        pca_var = getattr(self.args, 'pca_var', 0.85)
        pca_n = getattr(self.args, 'pca_n', None)
        pca_report_file = None
        if use_pca:
            pca_report_file = os.path.join(
                self.args.output_dir,
                f'QCValue_eval_{weight_method}_{eval_method}_pca_report.csv'
            )
        stations = [s for s, _ in self.station_rows]
        modules = self.modules_all
        data = np.full((len(stations), len(modules)), np.nan, dtype=float)
        for i, (_, vals) in enumerate(self.station_rows):
            for j, m in enumerate(modules):
                v = vals.get(m, '')
                if v == '' or v is None:
                    continue
                try:
                    data[i, j] = float(v)
                except (ValueError, TypeError):
                    pass
        # ---- 有监督分支：spp3drms 作为标签，k-means + KNN ----
        if getattr(self.args, 'eval_supervised', False):
            self._run_evaluate_KKNN(stations, modules, data)
            return

        result = evaluate_matrix(
            stations, modules, data,
            weight_method=weight_method,
            eval_method=eval_method,
            rho=rho, lamb=lamb,
            use_pca=use_pca,
            pca_var_threshold=pca_var,
            pca_n_components=pca_n,
            pca_report_file = pca_report_file,
        )
        if result is None:
            return

        stations, modules, data_used, weights, scores, ranking = result
        order = np.argsort(-scores)

        # CSV 输出的列名自动反映是否使用了 PCA（PC1、PC2 ...）
        out_file = os.path.join(
            self.args.output_dir,
            f'QCValue_eval_{weight_method}_{eval_method}'
            f'{"_pca" if use_pca else ""}.csv'
        )
        with open(out_file, 'w', newline='') as f:
            writer = csv.writer(f)
            first_col = 'date' if self.is_sstation else 'station'
            writer.writerow([first_col] + list(modules) + ['Score', 'Rank'])
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
                writer.writerow([stations[i]] + raw_row +
                                [f"{scores[i]:.4f}", int(ranking[i])])
            writer.writerow([])
            writer.writerow(['Weight'] + [f"{w:.4f}" for w in weights] + ['', ''])

        eval_title = '单站时间序列评估结果' if self.is_sstation else '多站评估结果'
        row_name = 'date' if self.is_sstation else 'station'
        print(f"\n===== {eval_title} [{weight_method} + {eval_method}"
              f"{' + PCA' if use_pca else ''}] =====")
        print("指标权重:")
        for m, w in zip(modules, weights):
            print(f"    {m:<14s} {w:.4f}")
        print(f"    {row_name:<10s}{'Score':>12s}{'Rank':>8s}")
        print("    " + "-" * 30)
        for i in order:
            print(f"    {stations[i]:<10s}{scores[i]:>12.4f}{int(ranking[i]):>8d}")
        print(f"-->跨站评估结果已保存至: {out_file}")

    def _run_evaluate_KKNN(self, stations, modules, data):
        """有监督质量分级：spp3drms 作为标签，k-means + KNN。"""
        print("\n-->跨日期有监督质量分级(k-means + KNN)" if self.is_sstation
              else "\n-->跨站有监督质量分级(k-means + KNN)")

        label_module = getattr(self.args, 'eval_label_module', 'spp3drms')
        n_clusters   = getattr(self.args, 'eval_n_clusters', 10)
        knn_k        = getattr(self.args, 'eval_knn_k', None)
        use_pca      = getattr(self.args, 'pca', False)
        pca_var      = getattr(self.args, 'pca_var', 0.95)
        pca_n        = getattr(self.args, 'pca_n', None)

        pca_report_file = None
        if use_pca:
            pca_report_file = os.path.join(
                self.args.output_dir,
                'QCValue_eval_KKNN_pca_report.csv'
            )

        result = evaluate_KKNN(
            stations, list(modules), data,
            label_module=label_module,
            n_clusters=n_clusters,
            knn_k=knn_k,
            use_pca=use_pca,
            pca_var_threshold=pca_var,
            pca_n_components=pca_n,
            pca_print=True,
            pca_report_file=pca_report_file,
        )
        if result is None:
            return

        grades       = result['grades']
        cluster_ids  = result['cluster_ids']
        label_values = result['label_values']
        feature_out  = result['feature_modules']
        data_used    = result['data_used']
        knn_k_final  = result['knn_k']
        cv_acc       = result['knn_cv_acc']

        out_file = os.path.join(
            self.args.output_dir,
            f'QCValue_eval_KKNN'
            f'{"_pca" if use_pca else ""}.csv'
        )
        first_col = 'date' if self.is_sstation else 'station'
        order = np.argsort(-grades)
        with open(out_file, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow([first_col, 'grade', 'cluster', label_module]
                            + list(feature_out))
            for i in order:
                row = [stations[i], int(grades[i])]
                row.append(int(cluster_ids[i]) if cluster_ids[i] >= 0 else '')
                v = label_values[i]
                row.append('' if not np.isfinite(v) else f"{v:.4f}")
                for val in data_used[i]:
                    row.append('' if not np.isfinite(val) else f"{val:.4f}")
                writer.writerow(row)
            writer.writerow([])
            writer.writerow(['# k-means clusters', n_clusters,
                             'KNN k', knn_k_final,
                             'CV acc',
                             '' if not np.isfinite(cv_acc) else f"{cv_acc:.4f}"])

        print(f"\n===== 有监督质量分级 [k-means + KNN]"
              f"{' + PCA' if use_pca else ''} =====")
        print(f"标签列: {label_module}    聚类数 k={n_clusters}    "
              f"KNN n_neighbors={knn_k_final}")
        if np.isfinite(cv_acc):
            print(f"KNN 交叉验证平均准确率: {cv_acc * 100:.2f}%")
        print(f"    {first_col:<10s}{'Grade':>8s}{'Cluster':>10s}"
              f"{label_module:>14s}")
        print("    " + "-" * 46)
        for i in order:
            cid = cluster_ids[i] if cluster_ids[i] >= 0 else -1
            lv = label_values[i]
            lv_s = 'N/A' if not np.isfinite(lv) else f"{lv:.4f}"
            print(f"    {stations[i]:<10s}{int(grades[i]):>8d}"
                  f"{cid:>10d}{lv_s:>14s}")
        print(f"-->有监督分级结果已保存至: {out_file}")
