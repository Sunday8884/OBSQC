#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import re
import glob
import argparse
import datetime
import numpy as np
from com.readObs import readObsHead, readObs
from com.readNav import readNav
from spos.spp import spp


# ========================= 辅助函数 =========================
def extract_station_from_filename(filename):
    """
    从观测文件名前 4 位提取测站名，并转为大写。
    例如 'wuh20050.24o' -> 'WUH2'。
    若文件名不规范则返回 None。
    """
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
    if len(matches) > 1:
        print(f"警告：--snx 匹配到多个文件，使用第一个：{matches[0]}")
    return matches[0]

def parse_time(timestr):
    """将时间字符串（YYYY-MM-DD HH:MM:SS）解析为datetime对象。"""
    if not timestr or timestr.strip() == "":
        return None
    return datetime.datetime.strptime(timestr, '%Y-%m-%d %H:%M:%S')

def parse_band_string(band_str, sys_list=None):
    """
    解析波段字符串，例如 "G:C1C,C2W;C:L2I,L6I"，生成适用于 spp() 的 bandChoose 字典。
    """
    if not band_str:
        return {}
    band_dict = {}
    for part in band_str.split(';'):
        if ':' not in part:
            continue
        sys, bands = part.split(':', 1)
        if sys_list and sys not in sys_list:
            continue
        band_dict[sys] = {band.strip(): {} for band in bands.split(',') if band.strip()}
    return band_dict

def ensure_dir(path):
    """确保目录存在，若不存在则创建。"""
    if not os.path.exists(path):
        os.makedirs(path)

def get_bandchoose(systems, bands):
    """结合 systems 和 bands 参数生成最终的 bandChoose 字典。"""
    default_bands = {
        'G': ['C1C', 'C2W'],
        'C': ['C2I', 'C6I'],
        'E': ['C1X', 'C5X'],
    }
    if bands:
        band_dict = parse_band_string(bands, sys_list=systems)
    else:
        band_dict = {}
        for sys in (systems or default_bands.keys()):
            if sys in default_bands:
                band_dict[sys] = {code: {} for code in default_bands[sys]}
    return band_dict

def extract_date_key(filename):
    """从文件名提取 '年份_日期编码' 组合键，用于匹配观测与导航文件。"""
    basename = os.path.basename(filename)
    match = re.search(r'(\d{4})\.(\d{2})[opn]', basename)
    if match:
        date_code = match.group(1)
        year_code = match.group(2)
        return f"{year_code}_{date_code}"
    else:
        return None

# ========================= 单文件处理 =========================
def process_one_file(obs_file, navData, args, base_output_dir):
    """处理单个观测文件：执行SPP，自动写出 *_spp.txt（含date/E/N/U），并可选绘图。"""
    file_basename = os.path.splitext(os.path.basename(obs_file))[0]
    out_dir = os.path.join(base_output_dir, file_basename)
    ensure_dir(out_dir)

    # 时间过滤
    start_dt = parse_time(args.start) if args.start else None
    end_dt = parse_time(args.end) if args.end else None

    # 构建 bandChoose
    band_choose = get_bandchoose(args.systems, args.bands)
    if not band_choose:
        print(f"错误：未为任何系统指定波段。请使用 --bands 或确保存在默认波段。")
        return

    # 从观测文件名自动提取测站名
    station = extract_station_from_filename(obs_file)
    if station is None:
        print(f"警告：无法从文件名 {obs_file} 提取测站名，将使用坐标均值作为 ENU 参考。")
    else:
        pass
    print(f"{'='*60}")
    print(f"正在处理观测文件：{obs_file}")

    # 读取观测数据
    obsHead = readObsHead(obs_file, needSatList=True)
    obsData = readObs(obs_file, obsHead=obsHead)
    print(f"历元数量：{len(obsData)}")

    # SPP 处理
    posFile = os.path.join(out_dir, f"{file_basename}_spp.txt")
    plot_save = os.path.join(out_dir, f"{file_basename}_spp.png")
    posData, Xmean, Ymean, Zmean = spp(
        obsHead, obsData, navData, band_choose,
        startDatetime=start_dt, endDatetime=end_dt,
        posFile=posFile,
        plot_save=plot_save,
        snx_file=args.snx,
        station=station,
    )

    if posData is not None:
        print(f"SPP完成，结果保存至 {posFile}")
    else:
        print("SPP 失败：有效历元不足。")

# ========================= 主函数 =========================
def run_spp(args):
    # 解析 SNX（通配符 / 大小写不敏感）
    args.snx = resolve_snx_file(args.snx)
    if args.snx:
        print(f"使用SINEX文件：{args.snx}")

    # 查找观测文件
    obs_files = glob.glob(args.obs)
    if not obs_files:
        print(f"错误：未找到匹配 '{args.obs}' 的观测文件")
        return
    print(f"找到 {len(obs_files)} 个观测文件：")
    # for f in obs_files:
    #     print(f"  - {f}")

    # 查找导航文件
    nav_files = glob.glob(args.nav)
    if not nav_files:
        print(f"错误：未找到匹配 '{args.nav}' 的导航文件")
        return
    print(f"找到 {len(nav_files)} 个导航文件：")
    # for f in nav_files:
    #     print(f"  - {f}")

    base_output_dir = args.out_dir
    ensure_dir(base_output_dir)

    # 只有一个导航文件 -> 所有观测文件共用该导航数据
    if len(nav_files) == 1:
        nav_file = nav_files[0]
        navData = readNav(nav_file)
        if not navData:
            print("错误：读取导航文件失败。")
            return
        for obs_file in obs_files:
            process_one_file(obs_file, navData, args, base_output_dir)
    else:
        # 多个导航文件 -> 根据日期编码配对
        nav_dict = {}
        for nav_file in nav_files:
            key = extract_date_key(nav_file)
            if key is None:
                print(f"警告：无法从导航文件名提取日期编码，跳过 {nav_file}")
                continue
            print(f"读取导航文件 {nav_file}，日期编码 {key}")
            navData = readNav(nav_file)
            if navData:
                nav_dict[key] = navData
            else:
                print(f"警告：读取导航文件失败 {nav_file}")

        if not nav_dict:
            print("错误：没有可用的导航数据。")
            return

        for obs_file in obs_files:
            obs_key = extract_date_key(obs_file)
            if obs_key is None:
                print(f"警告：无法从观测文件名提取日期编码，跳过 {obs_file}")
                continue
            if obs_key in nav_dict:
                print(f"匹配成功：观测文件 {obs_file} (编码 {obs_key}) -> 导航文件编码 {obs_key}")
                process_one_file(obs_file, nav_dict[obs_key], args, base_output_dir)
            else:
                print(f"警告：未找到与观测文件 {obs_file} (编码 {obs_key}) 匹配的导航文件，跳过")

    print(f"所有SPP处理完成，结果保存在 {base_output_dir}")


# ========================= 命令行解析器 =========================
def create_parser():
    parser = argparse.ArgumentParser(description="使用广播星历的标准单点定位（支持观测和导航文件通配符，自动配对）")
    parser.add_argument("--obs", required=True, help="输入 RINEX 观测文件（支持通配符，例如 './test/*.24o'）")
    parser.add_argument("--nav", required=True, help="广播导航文件（支持通配符，例如 './test/*.24p'）")
    parser.add_argument("--out_dir", default="./spp", help="输出根目录（默认：./spp）。每个文件会创建单独的子文件夹，并自动生成 *_spp.txt")

    # 时间过滤
    parser.add_argument("--start", help="起始时间（YYYY-MM-DD HH:MM:SS）")
    parser.add_argument("--end", help="结束时间（YYYY-MM-DD HH:MM:SS）")

    # 卫星/波段
    parser.add_argument("--systems", nargs="+", help="使用的 GNSS 系统，例如 G C E")
    parser.add_argument("--bands", help="波段选择，例如 'G:C1C,C2W;C:C2I,C6I'。若不提供，则使用默认波段。")

    # 参考坐标（SINEX 周解）
    parser.add_argument("--snx", default=None, help="SINEX 文件路径（.snx），用于提供 ENU 参考周解坐标")

    return parser


if __name__ == "__main__":
    parser = create_parser()
    args = parser.parse_args()
    run_spp(args)