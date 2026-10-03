#!/bin/bash
# ======================================================
#   cd /f/OBSQC    bash spp.sh
# ======================================================

OBS_FILE="./test/data/*.24o"   # obs文件,支持通配符,"./test/*.24o"
NAV_FILE="./test/data/*.24p"   # nav文件
SNX_FILE="./test/data/*.SNX"   # snx文件
OUTPUT_DIR="./test/qc"         # 输出根目录
START_TIME=""                  # 时间过滤,留空表示处理整个文件,例"2023-01-01 23:59:59"
END_TIME=""
SYSTEMS="C G"                # 卫星系统,空格分隔,"G C E"
BANDS=""                       # 波段选择,"G:C1C,C2W;C:C2I,C6I"
                               # 留空:默认波段,"GPS: C1C,C2W； BDS: C2I,C6I； GLO: C1X,C5X"

# 构建命令行
CMD="python spp_main.py --obs \"$OBS_FILE\" --nav \"$NAV_FILE\" --out_dir \"$OUTPUT_DIR\""

[ -n "$START_TIME" ] && CMD="$CMD --start \"$START_TIME\""
[ -n "$END_TIME" ]   && CMD="$CMD --end \"$END_TIME\""
[ -n "$SYSTEMS" ]    && CMD="$CMD --systems $SYSTEMS"
[ -n "$BANDS" ]      && CMD="$CMD --bands \"$BANDS\""
[ -n "$SNX_FILE" ]   && CMD="$CMD --snx \"$SNX_FILE\""

echo "正在执行: $CMD"
eval $CMD