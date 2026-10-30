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
#    "F:/GNSS_data/O/2024/ffmj/*.24o"
#    "F:/GNSS_data/O/2024/coco/*.24o"
    "./test/data/brux0050.24o"
)
OUTPUT_DIR_LIST=(
    #"F:/GNSS_data/O/2024/brux"
    #"F:/GNSS_data/O/2024/wuh2"
#    "F:/GNSS_data/O/2024/ptbb"
#    "F:/GNSS_data/O/2024/jfng"
#    "F:/GNSS_data/O/2024/urum"
#    "F:/GNSS_data/O/2024/abpo"
#    "F:/GNSS_data/O/2024/ffmj"
#    "F:/GNSS_data/O/2024/coco"
     "./test"
)

# ---- 内层循环：每个文件组要分别处理的系统 ----
SYSTEMS_LIST=(C G E)

# ---- 共用参数 ----
#NAV_FILE="F:/GNSS_data/N/2024/*.24p"     # 导航文件（station/Sstation/time 用于 SPP；prn 模式供完整率算理论可见口径；留空则跳过）
#SNX_FILE="F:/GNSS_data/SNX/*.SNX"       # SINEX 周解（station/Sstation/time 模式提供 ENU 参考原点）
NAV_FILE="./test/data/*.24p"           # 导航文件（仅 station 模式用于 SPP；留空则跳过 SPP）
SNX_FILE="./test/data/*.SNX"
START_TIME=""                            # 例如："2023-01-01 00:00:00"
END_TIME=""                              # 例如："2023-01-01 23:59:59"
SATELLITES=""                            # 空格分隔："G01 G02 C06"
BANDS=""                                 # 波段选择 "G:L1C,L2W;C:L2I,L6I"
CUTOFF=0                                # 完整率模块的理论可见截止高度角（度），需 NAV_FILE

MODE="time"                               # prn 逐文件卫星评估 | station 跨文件测站评估 | Sstation 单站多日期时间序列 | time 逐历元时间序列

# time 模式参数
TIME_ORIGIN=""                            # 时间原点，留空 = 由数据中最早的实际历元自动推导
TIME_INTERVAL=""                          # 采样间隔（秒），留空 = 由实际时间戳中位差分推导
SAT_LIST="C30 G05 E12"                               # 选星参数（仅 time 模式），如 "C05 C06"，生成单星汇总表

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
        [ -n "$CUTOFF" ] && CMD="$CMD --cutoff $CUTOFF"

        # ---- time 模式专属参数 ----
        if [ "$MODE" = "time" ]; then
            [ -n "$TIME_ORIGIN" ] && CMD="$CMD --time_origin \"$TIME_ORIGIN\""
            [ -n "$TIME_INTERVAL" ] && CMD="$CMD --time_interval $TIME_INTERVAL"
            # SAT_LIST 可能跨系统（如 "C05 G01"），按当前系统过滤，
            # 避免为不属于本系统的卫星反复打印「无数据」提示
            if [ -n "$SAT_LIST" ]; then
                SYS_SATS=""
                for S in $SAT_LIST; do
                    case "$S" in
                        "$SYS"*) SYS_SATS="$SYS_SATS $S" ;;
                    esac
                done
                [ -n "$SYS_SATS" ] && CMD="$CMD --sat$SYS_SATS"
            fi
        fi

        # ---- SPP 仅在 station/Sstation/time 模式下注入 --nav / --snx ----
        if { [ "$MODE" = "station" ] || [ "$MODE" = "Sstation" ] || [ "$MODE" = "time" ]; } && [ -n "$NAV_FILE" ]; then
            CMD="$CMD --nav \"$NAV_FILE\""
        fi
        { [ "$MODE" = "station" ] || [ "$MODE" = "Sstation" ] || [ "$MODE" = "time" ]; } && [ -n "$SNX_FILE" ] && CMD="$CMD --snx \"$SNX_FILE\""
        # ---- prn 模式下，完整率模块的"理论可见"口径也需要星历 ----
        if [ "$MODE" = "prn" ] && [ "$ENABLE_COMPLETENESS" -eq 1 ] && [ -n "$NAV_FILE" ]; then
            CMD="$CMD --nav \"$NAV_FILE\""
        fi

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
    done
done

echo ""
echo "全部批量任务完成"
