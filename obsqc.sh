#!/bin/bash
# ======================================================
#   GNSS QC 脚本
#   cd /f/OBSQC && bash obsqc.sh
#
#   说明:
#     - MODE=prn     : 逐文件 QC，不做 SPP
#     - MODE=station : 跨文件 QC，并对每个文件额外计算 SPP
#     - MODE=Sstation : 跨文件 QC，单站多日期时间序列
#     - MODE=time    : 跨文件 QC，逐历元时间序列 QCValue.csv，统计SPP选用的频段
#                     （--snx 提供 ENU 参考原点，与 station 模式一致）
#   综合评价项目输入本脚本产出的 QCValue_comparison.csv 运行
# ======================================================

#INPUT_OBS="F:/GNSS_data/O/2024/brux0050*.24o"       # 具体文件或通配符 "./test/*.24o" 或 "./test/wuh20050.24o"
#NAV_FILE="F:/GNSS_data/N/2024/*.24p"                # 导航文件（仅 station 模式用于 SPP；留空则跳过 SPP）
#SNX_FILE="F:/GNSS_data/SNX/*.SNX"
#OUTPUT_DIR="F:/GNSS_data/O/2024/brux/qc"            # 输出目录

INPUT_OBS="./test/data/brux0050.24o"   # 具体文件或通配符 "./test/*.24o" 或 "./test/wuh20050.24o"
NAV_FILE="./test/data/*.24p"           # 导航文件（仅 station 模式用于 SPP；留空则跳过 SPP）
SNX_FILE="./test/data/*.SNX"
OUTPUT_DIR="./test/qc"            # 输出目录

START_TIME=""                     # 例如："2023-01-01 00:00:00"
END_TIME=""                       # 例如："2023-01-01 23:59:59"
SYSTEMS="C"                       # 空格分隔："G C E S I J W"
SATELLITES=""                     # 空格分隔："G01 C06"
BANDS=""                          # qc频段选择 "G:L1C;C:L2I,L6I",为空则遍历obs头所有可用频段
CUTOFF=0                          # 完整率模块的理论可见截止高度角（度），需 NAV_FILE

# qc模块开关（1 = 启用，0 = 禁用）
ENABLE_SATNUM=1
ENABLE_OBSFREQ=1
ENABLE_COMPLETENESS=1
ENABLE_CNR=1
ENABLE_CYCLESLIP=1
ENABLE_PHASENOISE=1
ENABLE_MULTIPATH=1
ENABLE_CMC=1
ENABLE_IOD=1
ENABLE_HIGHORDERDIFF=1

ENABLE_SPP=1

MODE="time"                    # prn 逐文件卫星评估 | station 跨文件测站评估 | Sstation 单站多日期时间序列 | time 逐历元时间序列

# time 模式参数
TIME_ORIGIN=""                     # 时间原点（YYYY-MM-DD HH:MM:SS），留空 = 由数据中最早的实际历元自动推导
TIME_INTERVAL=""                  # 采样间隔（秒），留空 = 由实际时间戳中位差分推导
SAT_LIST="C25"                        # 选星参数（仅 time 模式），如 "C05 C06"，生成单星汇总表

# 构建命令行
CMD="python obsqc.py --input_obs \"$INPUT_OBS\" --output_dir \"$OUTPUT_DIR\" --mode $MODE"

[ -n "$START_TIME" ] && CMD="$CMD --start \"$START_TIME\""
[ -n "$END_TIME" ] && CMD="$CMD --end \"$END_TIME\""
[ -n "$SYSTEMS" ] && CMD="$CMD --systems $SYSTEMS"
[ -n "$SATELLITES" ] && CMD="$CMD --satellites $SATELLITES"
[ -n "$BANDS" ] && CMD="$CMD --bands \"$BANDS\""
[ -n "$CUTOFF" ] && CMD="$CMD --cutoff $CUTOFF"

# ---- time 模式专属参数 ----
if [ "$MODE" = "time" ]; then
    [ -n "$TIME_ORIGIN" ] && CMD="$CMD --time_origin \"$TIME_ORIGIN\""
    [ -n "$TIME_INTERVAL" ] && CMD="$CMD --time_interval $TIME_INTERVAL"
    [ -n "$SAT_LIST" ] && CMD="$CMD --sat $SAT_LIST"
fi

# ---- SPP 仅在 station/Sstation/time 模式下注入 --nav ----
if { [ "$MODE" = "station" ] || [ "$MODE" = "Sstation" ] || [ "$MODE" = "time" ]; } && [ -n "$NAV_FILE" ]; then
    CMD="$CMD --nav \"$NAV_FILE\""
fi
# ---- prn 模式下，完整率模块的"理论可见"口径也需要星历 ----
if [ "$MODE" = "prn" ] && [ "$ENABLE_COMPLETENESS" -eq 1 ] && [ -n "$NAV_FILE" ]; then
    CMD="$CMD --nav \"$NAV_FILE\""
fi
{ [ "$MODE" = "station" ] || [ "$MODE" = "Sstation" ] || [ "$MODE" = "time" ]; } && [ -n "$SNX_FILE" ] && CMD="$CMD --snx \"$SNX_FILE\""

# 模块标志
[ "$ENABLE_SATNUM" -eq 1 ] && CMD="$CMD --satnum" || CMD="$CMD --no-satnum"
[ "$ENABLE_OBSFREQ" -eq 1 ] && CMD="$CMD --freq" || CMD="$CMD --no-freq"
[ "$ENABLE_COMPLETENESS" -eq 1 ] && CMD="$CMD --completeness" || CMD="$CMD --no-completeness"
[ "$ENABLE_CNR" -eq 1 ] && CMD="$CMD --cnr" || CMD="$CMD --no-cnr"
[ "$ENABLE_CYCLESLIP" -eq 1 ] && CMD="$CMD --slip" || CMD="$CMD --no-slip"
[ "$ENABLE_PHASENOISE" -eq 1 ] && CMD="$CMD --phasenoise" || CMD="$CMD --no-phasenoise"
[ "$ENABLE_MULTIPATH" -eq 1 ] && CMD="$CMD --multipath" || CMD="$CMD --no-multipath"
[ "$ENABLE_CMC" -eq 1 ] && CMD="$CMD --cmc" || CMD="$CMD --no-cmc"
[ "$ENABLE_IOD" -eq 1 ] && CMD="$CMD --iod" || CMD="$CMD --no-iod"
[ "$ENABLE_HIGHORDERDIFF" -eq 1 ] && CMD="$CMD --highOrderDiff" || CMD="$CMD --no-highOrderDiff"
[ "$ENABLE_SPP" -eq 1 ] && CMD="$CMD --spp" || CMD="$CMD --no-spp"

echo "正在执行: $CMD"
eval $CMD
