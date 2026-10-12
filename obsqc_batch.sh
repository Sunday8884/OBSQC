#!/bin/bash
# ======================================================
#   GNSS QC 批量处理脚本
#
#   两层循环：
#     外层：逐个处理 INPUT_OBS_LIST 中的观测文件组
#     内层：对每组分别跑 SYSTEMS_LIST 中的每个系统
#   输出目录：<该组 OUTPUT_DIR>/qc/<系统>   例如 .../brux/qc/G

#   cd /f/OBSQC && bash obsqc_batch.sh
# ======================================================

# ---- 外层循环：观测文件组与各自输出根目录 ----
INPUT_OBS_LIST=(
    #"F:/GNSS_data/O/2024/brux/*.24o"
    #"F:/GNSS_data/O/2024/wuh2/*.24o"
#    "F:/GNSS_data/O/2024/ptbb/*.24o"
#    "F:/GNSS_data/O/2024/jfng/*.24o"
#    "F:/GNSS_data/O/2024/urum/*.24o"
#    "F:/GNSS_data/O/2024/abpo/*.24o"
    "F:/GNSS_data/O/2024/ffmj/*.24o"
    "F:/GNSS_data/O/2024/coco/*.24o"
)
OUTPUT_DIR_LIST=(
    #"F:/GNSS_data/O/2024/brux"
    #"F:/GNSS_data/O/2024/wuh2"
#    "F:/GNSS_data/O/2024/ptbb"
#    "F:/GNSS_data/O/2024/jfng"
#    "F:/GNSS_data/O/2024/urum"
#    "F:/GNSS_data/O/2024/abpo"
    "F:/GNSS_data/O/2024/ffmj"
    "F:/GNSS_data/O/2024/coco"
)

# ---- 内层循环：每个文件组要分别处理的系统 ----
SYSTEMS_LIST=(C G E)

# ---- 共用参数 ----
NAV_FILE="F:/GNSS_data/N/2024/*.24p"     # 导航文件（仅 station/Sstation 模式用于 SPP；留空则跳过 SPP）
SNX_FILE="F:/GNSS_data/SNX/*.SNX"
START_TIME=""                            # 例如："2023-01-01 00:00:00"
END_TIME=""                              # 例如："2023-01-01 23:59:59"
SATELLITES=""                            # 空格分隔："G01 G02 C06"
BANDS=""                                 # 波段选择 "G:L1C,L2W;C:L2I,L6I"

MODE="Sstation"                          # prn 逐文件卫星评估 | station 跨文件测站评估 | Sstation 单站多日期时间序列

# qc模块开关（1 = 启用，0 = 禁用）
ENABLE_SATNUM=1
ENABLE_OBSFREQ=1
ENABLE_CNR=1
ENABLE_CYCLESLIP=1
ENABLE_PHASENOISE=1
ENABLE_MULTIPATH=1
ENABLE_CMC=1
ENABLE_IOD=1
ENABLE_HIGHORDERDIFF=1

ENABLE_SPP=1

# ======================================================
#   以下为循环主体，一般无需修改
# ======================================================
for i in "${!INPUT_OBS_LIST[@]}"; do
    INPUT_OBS="${INPUT_OBS_LIST[$i]}"
    BASE_OUT="${OUTPUT_DIR_LIST[$i]}"

    for SYS in "${SYSTEMS_LIST[@]}"; do
        OUTPUT_DIR="$BASE_OUT/qc/$SYS"

        echo ""
        echo "############################################################"
        echo "# 输入: $INPUT_OBS"
        echo "# 系统: $SYS    输出: $OUTPUT_DIR"
        echo "############################################################"

        # 构建命令行
        CMD="python obsqc.py --input_obs \"$INPUT_OBS\" --output_dir \"$OUTPUT_DIR\" --mode $MODE"

        [ -n "$START_TIME" ] && CMD="$CMD --start \"$START_TIME\""
        [ -n "$END_TIME" ] && CMD="$CMD --end \"$END_TIME\""
        CMD="$CMD --systems $SYS"
        [ -n "$SATELLITES" ] && CMD="$CMD --satellites $SATELLITES"
        [ -n "$BANDS" ] && CMD="$CMD --bands \"$BANDS\""

        # ---- SPP 仅在 station/Sstation 模式下注入 --nav / --snx ----
        if { [ "$MODE" = "station" ] || [ "$MODE" = "Sstation" ]; } && [ -n "$NAV_FILE" ]; then
            CMD="$CMD --nav \"$NAV_FILE\""
        fi
        { [ "$MODE" = "station" ] || [ "$MODE" = "Sstation" ]; } && [ -n "$SNX_FILE" ] && CMD="$CMD --snx \"$SNX_FILE\""

        # 模块标志
        [ "$ENABLE_SATNUM" -eq 1 ] && CMD="$CMD --satnum" || CMD="$CMD --no-satnum"
        [ "$ENABLE_OBSFREQ" -eq 1 ] && CMD="$CMD --freq" || CMD="$CMD --no-freq"
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
    done
done

echo ""
echo "全部批量任务完成"
