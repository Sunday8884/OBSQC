#!/bin/bash
# ======================================================
#   GNSS QC 脚本
#   python obsqc.py --input_obs "./test/*.24o" --output_dir ./test --systems C
#   cd /f/OBSQC    bash obsqc.sh
#
#   说明:
#     - MODE=prn     : 逐文件 QC，不做 SPP
#     - MODE=station : 跨文件 QC，并对每个文件额外计算 SPP
#     - MODE=Sstation : 跨文件 QC，单站多日期时间序列
#     - ENABLE_SUPERVISED=1：k-means + KNN 有监督质量分级，spp3drms 作为标签
# ======================================================

INPUT_OBS="F:/GNSS_data/O/2024/brux/*.24o"          # 具体文件或通配符 "./test/*.24o" 或 "./test/wuh20050.24o"
NAV_FILE="F:/GNSS_data/N/2024/*.24p"           # 导航文件（仅 station 模式用于 SPP；留空则跳过 SPP）
SNX_FILE="F:/GNSS_data/SNX/*.SNX"
OUTPUT_DIR="F:/GNSS_data/O/2024/brux/qc"            # 输出目录（默认：./qc_results）
START_TIME=""                     # 例如："2023-01-01 00:00:00"
END_TIME=""                       # 例如："2023-01-01 23:59:59"
SYSTEMS="C"                     # 空格分隔："G C E"
SATELLITES=""                     # 空格分隔："G01 G02 C06"
BANDS=""                          # 波段选择 "G:L1C,L2W;C:L2I,L6I"

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

# 评估开关
ENABLE_EVALUATE=1                  # 1 = 启用评估，0 = 禁用
MODE="Sstation"                    # prn 逐文件卫星评估 | station 跨文件测站评估 | Sstation 单站多日期时间序列

# PCA 降维
ENABLE_PCA=1                      # 1 = 评估前做 PCA 降维
PCA_VAR=0.95                      # 累计贡献率阈值，ENABLE_PCA=1 且 PCA_N 为空时生效
PCA_N=""                          # 指定主成分个数，覆盖 PCA_VAR

EVAL_WEIGHT="ahp"                 # entropy | ahp | bwm | equal（有监督模式下不生效）
EVAL_METHOD="income"              # gra | topsis | waspas | income（有监督模式下不生效）
EVAL_RHO=0.5
EVAL_LAMBDA=0.5
INCOME_LEVELS=3                   # 每个指标的离散等级数
INCOME_K=3                        # kNN 邻居数
INCOME_BINNING="equal_width"      # 离散化方式: equal_width | equal_freq

# 有监督质量分级（k-means + KNN，spp3drms 作为标签）
ENABLE_SUPERVISED=1               # 1 = 启用有监督分级
EVAL_LABEL_MODULE="spp3drms"      # 监督标签列名
EVAL_N_CLUSTERS=5                 # k-means 聚类数
EVAL_KNN_K=""                     # KNN 邻居数，留空 = 交叉验证自动选取



# 构建命令行
CMD="python obsqc.py --input_obs \"$INPUT_OBS\" --output_dir \"$OUTPUT_DIR\" --mode $MODE"

[ -n "$START_TIME" ] && CMD="$CMD --start \"$START_TIME\""
[ -n "$END_TIME" ] && CMD="$CMD --end \"$END_TIME\""
[ -n "$SYSTEMS" ] && CMD="$CMD --systems $SYSTEMS"
[ -n "$SATELLITES" ] && CMD="$CMD --satellites $SATELLITES"
[ -n "$BANDS" ] && CMD="$CMD --bands \"$BANDS\""

# ---- SPP 仅在 station/Sstation 模式下注入 --nav ----
if { [ "$MODE" = "station" ] || [ "$MODE" = "Sstation" ]; } && [ -n "$NAV_FILE" ]; then
    CMD="$CMD --nav \"$NAV_FILE\""
fi
{ [ "$MODE" = "station" ] || [ "$MODE" = "Sstation" ]; } && [ -n "$SNX_FILE" ] && CMD="$CMD --snx \"$SNX_FILE\""

# PCA 降维
if [ "$ENABLE_PCA" -eq 1 ]; then
    CMD="$CMD --pca --pca_var $PCA_VAR"
    [ -n "$PCA_N" ] && CMD="$CMD --pca_n $PCA_N"
fi

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

# 评估参数
[ "$ENABLE_EVALUATE" -eq 1 ] && CMD="$CMD --evaluate" || CMD="$CMD --no-evaluate"
CMD="$CMD --eval_weight $EVAL_WEIGHT --eval_method $EVAL_METHOD"
CMD="$CMD --eval_rho $EVAL_RHO --eval_lambda $EVAL_LAMBDA"

# ---- INCOME 超参----
if [ "$EVAL_METHOD" = "income" ]; then
    CMD="$CMD --income_levels $INCOME_LEVELS"
    CMD="$CMD --income_k $INCOME_K"
    CMD="$CMD --income_binning $INCOME_BINNING"
fi

# ---- 有监督分级参数（k-means + KNN，spp3drms 作为标签）----
if [ "$ENABLE_SUPERVISED" -eq 1 ]; then
    CMD="$CMD --eval_supervised"
    CMD="$CMD --eval_label_module $EVAL_LABEL_MODULE"
    CMD="$CMD --eval_n_clusters $EVAL_N_CLUSTERS"
    [ -n "$EVAL_KNN_K" ] && CMD="$CMD --eval_knn_k $EVAL_KNN_K"
fi

echo "正在执行: $CMD"
eval $CMD
