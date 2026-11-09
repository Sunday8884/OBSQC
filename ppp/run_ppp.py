#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
run_ppp.py — PPP 命令行入口（风格与 spp.py 一致）。

用法示例（项目根目录下）：
    python ppp/run_ppp.py --obs "./test/*.24o" --sp3 "./test/*.sp3" \
        --out_dir "./test/ppp" --config ppp/ppp.toml

参数优先级：命令行 > --config 指定的配置文件（TOML/JSON） > 默认配置（ppp/pppConfig.py）。
"""
import os
import sys
import glob
import argparse
import datetime

# 保证以 `python ppp/run_ppp.py` 方式运行时项目根目录在 sys.path 中
PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if PROJECT_ROOT not in sys.path:
    sys.path.insert(0, PROJECT_ROOT)

from com.readObs import readObsHead, readObs
from com.readSp3 import readSp3, mergeSp3
from com.readClk import readRinexClk
from ppp.pppConfig import load_config, resolve_band_choose
from ppp.ppp import ppp


# ========================= 辅助函数（与 spp.py 保持一致） =========================
def extract_station_from_filename(filename):
    """从观测文件名前 4 位提取测站名（大写），不规范返回 None。"""
    basename = os.path.basename(filename)
    if len(basename) < 4:
        return None
    station = basename[:4].upper()
    if not station.isalnum():
        return None
    return station


def resolve_snx_file(snx_pattern):
    """对 --snx 参数做通配符展开；返回第一个匹配文件或 None。"""
    if not snx_pattern:
        return None
    matches = glob.glob(snx_pattern)
    if not matches:
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
    if len(matches) > 1:
        print(f"警告：--snx 匹配到多个文件，使用第一个：{matches[0]}")
    return matches[0]


def parse_time(timestr):
    """将时间字符串（YYYY-MM-DD HH:MM:SS）解析为 datetime 对象。"""
    if not timestr or timestr.strip() == "":
        return None
    return datetime.datetime.strptime(timestr, '%Y-%m-%d %H:%M:%S')


def ensure_dir(path):
    """确保目录存在。"""
    if not os.path.exists(path):
        os.makedirs(path)


# ========================= 单文件处理 =========================
def process_one_file(obs_file, sp3Data, clkData, cfg, base_output_dir):
    """处理单个观测文件：执行 PPP，写出 *_ppp.txt，可选绘图。"""
    file_basename = os.path.splitext(os.path.basename(obs_file))[0]
    out_dir = os.path.join(base_output_dir, file_basename)
    ensure_dir(out_dir)

    start_dt = parse_time(cfg.get('start'))
    end_dt = parse_time(cfg.get('end'))

    band_choose = resolve_band_choose(cfg)
    if not band_choose:
        print("错误：未为任何系统指定波段，请检查 config 的 systems/bands。")
        return

    station = extract_station_from_filename(obs_file)
    if not station:
        # 文件名不规范时回退到 RINEX 头 MARKER NAME（海潮 BLQ 匹配依赖测站名）
        obsHead0 = readObsHead(obs_file, needSatList=False)
        station = obsHead0.get('MARKER NAME') or None
        if station:
            station = station.upper()
    print(f"{'='*60}")
    print(f"[PPP] 正在处理观测文件：{obs_file}")

    obsHead = readObsHead(obs_file, needSatList=True)
    obsData = readObs(obs_file, obsHead=obsHead)
    print(f"[PPP] 历元数量：{len(obsData)}")

    posFile = os.path.join(out_dir, f"{file_basename}_ppp.txt")
    plotSave = os.path.join(out_dir, f"{file_basename}_ppp.png") if cfg.get('plot', True) else None

    posData, refX, refY, refZ = ppp(
        obsHead, obsData, sp3Data, clkData, band_choose,
        start_dt, end_dt,
        config=cfg,
        posFile=posFile,
        plot_show=False,
        plot_save=plotSave,
        snx_file=cfg.get('snx'),
        station=station,
    )

    if posData is not None:
        print(f"[PPP] 完成，结果保存至 {posFile}")
    else:
        print("[PPP] 失败：无有效解算历元。")


# ========================= 主流程 =========================
def run_ppp(args):
    # 1. 配置：默认 -> JSON -> 命令行
    cfg = load_config(args.config)
    for key in ('start', 'end', 'snx', 'mode'):
        val = getattr(args, key, None)
        if val:
            cfg[key] = val
    # 精密产品（命令行覆盖 JSON）
    for key, attr in (('atx_file', 'atx'), ('blq_file', 'blq'),
                      ('dcb_file', 'dcb'), ('erp_file', 'erp')):
        val = getattr(args, attr, None)
        if val:
            cfg[key] = val
    if args.systems:
        cfg['systems'] = args.systems
    if args.bands:
        cfg['bands'] = args.bands
    if args.grad:
        cfg['est_trop_grad'] = True
    if args.elev is not None:
        cfg['elev_cutoff_deg'] = args.elev
    if args.no_sat_var:
        cfg['use_sat_var'] = False
    if args.no_jump:
        cfg['phase_code_jump'] = False
    if args.no_plot:
        cfg['plot'] = False
    if args.ex:
        cfg['exclude_sats'] = [s.upper() for s in args.ex]

    cfg['snx'] = resolve_snx_file(cfg.get('snx'))
    if cfg.get('snx'):
        print(f"使用SINEX文件：{cfg['snx']}")

    # 2. 观测文件
    obs_files = glob.glob(args.obs)
    if not obs_files:
        print(f"错误：未找到匹配 '{args.obs}' 的观测文件")
        return
    print(f"找到 {len(obs_files)} 个观测文件")

    # 3. 精密星历（多个文件自动合并，按时段覆盖）
    sp3_files = glob.glob(args.sp3)
    if not sp3_files:
        print(f"错误：未找到匹配 '{args.sp3}' 的精密星历文件（SP3）")
        return
    sp3List = []
    for f in sp3_files:
        print(f"读取精密星历：{f}")
        sd = readSp3(f)
        if sd['epochs']:
            sp3List.append(sd)
        else:
            print(f"警告：SP3 文件无有效数据，跳过 {f}")
    if not sp3List:
        print("错误：没有可用的精密星历数据。")
        return
    sp3Data = sp3List[0] if len(sp3List) == 1 else mergeSp3(sp3List)
    print(f"精密星历：{len(sp3Data['epochs'])} 个历元，{len(sp3Data['sat'])} 颗卫星，"
          f"时段 {sp3Data['epochs'][0]} ~ {sp3Data['epochs'][-1]}")

    # 4. 精密钟差（可选）
    clkData = None
    if args.clk:
        clk_files = glob.glob(args.clk)
        if not clk_files:
            print(f"警告：未找到匹配 '{args.clk}' 的钟差文件，改用 SP3 内嵌钟差")
        else:
            clkList = []
            for f in clk_files:
                print(f"读取精密钟差：{f}")
                cd = readRinexClk(f)
                if cd['epochs']:
                    clkList.append(cd)
            if clkList:
                clkData = {'epochs': [], 'sat': {}}
                seen = set()
                for cd in clkList:
                    for ep in cd['epochs']:
                        if ep not in seen:
                            seen.add(ep)
                            clkData['epochs'].append(ep)
                    for prn, epDict in cd['sat'].items():
                        clkData['sat'].setdefault(prn, {}).update(epDict)
                clkData['epochs'].sort()
                print(f"精密钟差：{len(clkData['epochs'])} 个历元，{len(clkData['sat'])} 颗卫星")

    # 5. 逐文件处理
    base_output_dir = args.out_dir
    ensure_dir(base_output_dir)
    for obs_file in obs_files:
        process_one_file(obs_file, sp3Data, clkData, cfg, base_output_dir)

    print(f"所有PPP处理完成，结果保存在 {base_output_dir}")


# ========================= 命令行解析器 =========================
def create_parser():
    parser = argparse.ArgumentParser(
        description="精密单点定位 PPP（SP3 精密星历 + 双频无电离层组合 + Kalman 滤波，支持通配符）")
    parser.add_argument("--obs", required=True,
                        help="输入 RINEX 观测文件（支持通配符，例如 './test/*.24o'）")
    parser.add_argument("--sp3", required=True,
                        help="精密星历文件（SP3，支持通配符，例如 './test/*.SP3'）")
    parser.add_argument("--clk", default=None,
                        help="精密钟差文件（RINEX Clock，可选；缺省用 SP3 内嵌钟差）")
    parser.add_argument("--out_dir", default="./ppp_out",
                        help="输出根目录（默认：./ppp_out）。每个观测文件一个子文件夹，自动生成 *_ppp.txt")
    parser.add_argument("--config", default=os.path.join('ppp', 'ppp.toml'),
                        help="PPP 配置文件（TOML/JSON，默认：ppp/ppp.toml）")

    # 时间过滤
    parser.add_argument("--start", help="起始时间（YYYY-MM-DD HH:MM:SS）")
    parser.add_argument("--end", help="结束时间（YYYY-MM-DD HH:MM:SS）")

    # 卫星/波段（覆盖配置文件）
    parser.add_argument("--systems", nargs="+", help="GNSS 系统，例如 G C E")
    parser.add_argument("--bands", help="波段选择，例如 'G:C1C,C2W;C:C2I,C6I'")

    # 精密产品（覆盖配置文件）
    parser.add_argument("--atx", default=None,
                        help="ANTEX 相位中心改正文件（igs20.atx 等）")
    parser.add_argument("--blq", default=None,
                        help="BLQ 海潮载荷文件（FES2004 等）")
    parser.add_argument("--dcb", default=None,
                        help="CODE DCB 文件（P1-C1 码偏差）")
    parser.add_argument("--erp", default=None,
                        help="IGS ERP 地球定向参数文件（UT1-UTC 等）")

    # 参考坐标 / 模式
    parser.add_argument("--snx", default=None, help="SINEX 文件，提供 ENU 参考周解坐标")
    parser.add_argument("--mode", choices=['static', 'kinematic'], help="定位模式（默认取配置文件）")
    # 模型开关（覆盖配置文件）
    parser.add_argument("--grad", action="store_true",
                        help="启用对流层梯度估计（MacMillan 模型，RTKLIB TROPOPT_ESTG 同源）")
    parser.add_argument("--elev", type=float, default=None,
                        help="高度角截止角（度，覆盖配置 elev_cutoff_deg）")
    parser.add_argument("--no_sat_var", action="store_true",
                        help="关闭 SP3 卫星轨道/钟差方差加权（var_rs）")
    parser.add_argument("--no_jump", action="store_true",
                        help="关闭相位-码跳变改正（毫秒级钟跳检测）")
    parser.add_argument("--no_plot", action="store_true", help="不绘制结果图")
    parser.add_argument("--ex", nargs="+", default=None,
                        help="排除指定卫星，例如 --ex C01 C02 C03 C04 C05（BDS GEO）")
    return parser


if __name__ == "__main__":
    parser = create_parser()
    args = parser.parse_args()
    run_ppp(args)
