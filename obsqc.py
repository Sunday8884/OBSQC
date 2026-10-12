#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import glob
import argparse

# ========================= 主流程 =========================
def obs_qc(args):
    obs_files = glob.glob(args.input_obs)
    if not obs_files:
        print(f"错误：未找到匹配的文件: {args.input_obs}")
        return
    mode = getattr(args, 'mode', 'prn')
    print(f"处理模式: {mode}")

    if mode in ('station', 'Sstation'):
        # ===== Station / Sstation 模式（均走多文件引擎） =====
        from qc.qc_engine_exp import MultiFileQC
        MultiFileQC(args).run(obs_files)
    else:
        # ============ PRN 模式 ============
        from qc.qc_engine import QCProcessor
        total = len(obs_files)
        for idx, obs_file in enumerate(obs_files, 1):
            sub_dir = os.path.join(
                args.output_dir,
                os.path.splitext(os.path.basename(obs_file))[0]
            )
            print(f"\n========== 处理文件 [{idx}/{total}]: {obs_file} ==========")
            # 断点续跑：已有汇总结果的文件直接跳过
            qc_csv = os.path.join(sub_dir,os.path.splitext(os.path.basename(obs_file))[0] + '_QCValue.csv')
            if os.path.exists(qc_csv):
                print("   (检测到已有汇总结果，跳过重复计算)")
                continue
            try:
                QCProcessor(obs_file, sub_dir, args).process()
            except Exception as e:
                print(f"文件 {obs_file} 处理失败，已跳过: {e}")
                continue
        print(f"所有文件处理完成，结果保存在: {args.output_dir}")

# ========================= 命令行解析 =========================
def _add_toggle(parser, name, help_text):
    parser.add_argument(f"--{name}", action="store_true", default=True, help=f"启用{help_text}")
    parser.add_argument(f"--no-{name}", dest=name, action="store_false", help=f"禁用{help_text}")

def create_parser():
    parser = argparse.ArgumentParser(description="GNSS 观测数据质量检查（QC + SPP，不含综合评价）")
    parser.add_argument("--mode", default="prn", choices=["prn", "station", "Sstation"], help="处理模式：prn=逐文件卫星评估；station=跨文件测站评估；Sstation=单站多日期文件，生成时间序列质量统计")
    parser.add_argument("--input_obs", required=True, help="输入 RINEX 观测文件路径（支持通配符）")
    parser.add_argument("--nav", help="广播星历导航文件（支持通配符），仅 station 模式下用于计算 SPP")
    parser.add_argument("--snx", default=None, help="SINEX 周解文件（支持通配符），仅 station 模式有效，提供 ENU 参考")
    parser.add_argument("--output_dir", default="./qc", help="输出根目录（默认: ./qc）")
    parser.add_argument("--start", help="起始时间（YYYY-MM-DD HH:MM:SS）")
    parser.add_argument("--end", help="结束时间（YYYY-MM-DD HH:MM:SS）")
    parser.add_argument("--systems", nargs="+", help="要处理的 GNSS 系统，例如 G C E")
    parser.add_argument("--satellites", nargs="+", help="卫星 PRN 列表，例如 G01 G02 C06")
    parser.add_argument("--bands", help="波段选择，例如 'G:1C,2W;C:2I,6I'")

    # 模块开关
    _add_toggle(parser, "satnum",     "卫星数量统计")
    _add_toggle(parser, "freq",       "频段统计")
    _add_toggle(parser, "cnr",        "载噪比分析")
    _add_toggle(parser, "slip",       "周跳检测")
    _add_toggle(parser, "phasenoise", "相位噪声计算")
    _add_toggle(parser, "multipath",  "多路径计算")
    _add_toggle(parser, "cmc",        "CMC 计算")
    _add_toggle(parser, "iod",        "IOD 计算")
    _add_toggle(parser, "highOrderDiff",  "高阶差分噪声计算")
    _add_toggle(parser, "spp", "SPP计算,仅 station 模式生效")

    return parser

if __name__ == "__main__":

    parser = create_parser()
    args = parser.parse_args()
    obs_qc(args)
