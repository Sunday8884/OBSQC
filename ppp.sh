#!/bin/bash
# ======================================================
#   PPP 精密单点定位封装脚本，留空用配置文件./ppp/ppp.toml
#   cd /f/OBSQC    bash ppp.sh
# ======================================================

OBS_FILE="./ppptest/data/wuh22670.22o"       # obs文件,支持通配符
SP3_FILE="./ppptest/data/WUM0MGXFIN_20222670000_01D_15M_ORB.SP3" # SP3精密星历,支持通配符(多个文件自动合并)
CLK_FILE="./ppptest/data/WUM0MGXFIN_20222670000_01D_30S_CLK.CLK" # 精密钟差文件(RINEX Clock),留空则用SP3内嵌钟差
SNX_FILE="./ppptest/data/*.snx"       # snx文件,提供ENU参考周解坐标
ATX_FILE="./ppptest/data/igs20.atx"   # ANTEX天线改正文件
BLQ_FILE="./ppptest/data/FES2004_DXS.BLQ"                         # 海潮载荷文件
DCB_FILE="./ppptest/data/CAS0MGXRAP_20222670000_01D_01D_DCB.BSX"  # 码偏差文件
ERP_FILE="./ppptest/data/igs22287.erp"                            # 地球定向参数
OUTPUT_DIR="./ppptest/ppp"            # 输出根目录
CONFIG_FILE="./ppp/ppp.toml"          # PPP配置文件
START_TIME=""                         # 时间过滤,留空表示处理整个文件,例"2024-01-05 00:00:00"
END_TIME=""
SYSTEMS="G"                           # 卫星系统,空格分隔
BANDS=""                              # 双频波段选择,"G:C1C,C2W;C:C2I,C6I;E:C1X,C5X"
EX_SATS="C01 C02 C03 C04 C05 C06 C07 C08 C09 C10 C11 C12 C13 C14" # 排除的卫星,空格分隔,BDS GEO"C01 C02 C03 C04 C05"; BDS-2"C01 C02 C03 C04 C05 C06 C07 C08 C09 C10 C11 C12 C13 C14"
MODE="static"                         # static(静态) | kinematic(动态)
NO_PLOT=""                            # 置"1"则不绘图

# 构建命令行
CMD="python ppp/run_ppp.py --obs \"$OBS_FILE\" --sp3 \"$SP3_FILE\" --out_dir \"$OUTPUT_DIR\" --config \"$CONFIG_FILE\""

[ -n "$CLK_FILE" ]   && CMD="$CMD --clk \"$CLK_FILE\""
[ -n "$SNX_FILE" ]   && CMD="$CMD --snx \"$SNX_FILE\""
[ -n "$ATX_FILE" ]   && CMD="$CMD --atx \"$ATX_FILE\""
[ -n "$BLQ_FILE" ]   && CMD="$CMD --blq \"$BLQ_FILE\""
[ -n "$DCB_FILE" ]   && CMD="$CMD --dcb \"$DCB_FILE\""
[ -n "$ERP_FILE" ]   && CMD="$CMD --erp \"$ERP_FILE\""
[ -n "$START_TIME" ] && CMD="$CMD --start \"$START_TIME\""
[ -n "$END_TIME" ]   && CMD="$CMD --end \"$END_TIME\""
[ -n "$SYSTEMS" ]    && CMD="$CMD --systems $SYSTEMS"
[ -n "$BANDS" ]      && CMD="$CMD --bands \"$BANDS\""
[ -n "$EX_SATS" ]    && CMD="$CMD --ex $EX_SATS"
[ -n "$MODE" ]       && CMD="$CMD --mode $MODE"
[ "$NO_PLOT" = "1" ] && CMD="$CMD --no_plot"

echo "正在执行: $CMD"
eval $CMD
