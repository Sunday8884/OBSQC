#!/bin/bash
# ======================================================
#   带观测质量门控的 SPP（sppgate.py）参数配置
#   cd /f/OBSQC    bash sppgate.sh
#
#   门控层次（可叠加）:
#     [1] 阈值极端剔除  MP>MP_THRESH 或 CNR<CNR_THRESH 剔星,
#         剔后保底 MIN_KEEP 星（默认开启; 实验证实有效 +3.0%~+5.4%）
#     [2] 卫星级 S 逐历元剔除  GATE_SAT_S=worst1|fixed
#     [3] 历元级 S 筛选  EPOCH_S_PCT, 解后统计筛选, 解不变
#         S 评分器/赋权: S_SCORE_METHOD / S_WEIGHT_MODE
#   开启 [2]或[3] 时 QC 指标自动升级为全模块（9类）, 耗时增加
# ======================================================

# ---------- 输入/输出 ----------
OBS_FILE="./test/data/wuh20050.24o"   # obs文件, 支持通配符 "./test/data/*.24o"
NAV_FILE="./test/data/brdc0050.24p"   # nav文件, 支持通配符（多文件按日期配对）
SNX_FILE="./test/data/*.SNX"          # snx文件（提供 ENU 参考; 留空则无 RMS 摘要）
OUTPUT_DIR="./test/sppgate"           # 输出根目录

START_TIME=""                  # "2024-01-05 00:00:00", 解算时间窗口，留空处理整个文件
END_TIME=""                    # "2024-01-05 12:00:00"
SYSTEMS="C"                    # 卫星系统, 空格分隔 "G C E"
BANDS=""                       # 波段选择 "G:C1C,C2W;C:C2I,C6I;E:C1C,C5Q"  留空: 默认波段

# ---------- [1] 阈值极端剔除 ----------
NO_GATE=0                      # 1 = 关闭阈值门控，正常SPP解算
MP_THRESH=3.0                  # MP 极端阈值 (m)
CNR_THRESH=25.0                # CNR 极端下限 (dB-Hz)
MIN_KEEP=5                     # 每历元最少保留卫星数

# ---------- [2] 卫星级 S 逐历元剔除（剔星） ----------
GATE_SAT_S="off"               # off | worst1 | fixed
                               # ⚠ 全天实验全负(critic T=0.75 -0.20% /
                               #   T=0.70 -3.51% / T=0.60 -16.96%,
                               #   worst1 -1.70%), 仅作诊断对照, 默认 off
SAT_S_THRESH=0.75              # fixed 规则 S_sat 剔除阈值（0~1, 剔 S_sat>=阈值的星;
                               #  sat_score 双重归一化 bug 修复后 S_sat∈[0,1],
                               #  >=0.75 仅占 1.8%: 阈值需按此重新标定）
SAT_S_MIN_SAT=7                # worst1 规则要求的可见星数下限

# ---------- [3] 历元级 S 筛选（整历元剔除） ----------
EPOCH_S_PCT="20"               # 剔除比例(%), 空格分隔, 例"10 20"; 留空关闭
S_SCORE_METHOD="gra"           # S 评分器: wavg/gra/topsis/waspas/ecdf
S_WEIGHT_MODE="critic"         # 赋权模式: fixed/critic/entropy/equal/ahp/game
S_WEIGHT_CSV=""                # 权重表, 留空用默认 sppgate/weights_QC9.csv
# ---------- 运行模式 ----------
COMPARE=1                      # 1 = 基线与门控双跑并输出 RMS 对比（需 SNX_FILE）

# ==================== 构建命令行 ====================
CMD="python sppgate.py --obs \"$OBS_FILE\" --nav \"$NAV_FILE\" --out_dir \"$OUTPUT_DIR\""

[ -n "$SNX_FILE" ]   && CMD="$CMD --snx \"$SNX_FILE\""
[ -n "$START_TIME" ] && CMD="$CMD --start \"$START_TIME\""
[ -n "$END_TIME" ]   && CMD="$CMD --end \"$END_TIME\""
[ -n "$SYSTEMS" ]    && CMD="$CMD --systems $SYSTEMS"
[ -n "$BANDS" ]      && CMD="$CMD --bands \"$BANDS\""

[ "$NO_GATE" = "1" ] && CMD="$CMD --no-gate"
CMD="$CMD --mp_thresh $MP_THRESH --cnr_thresh $CNR_THRESH --min_keep $MIN_KEEP"

[ "$GATE_SAT_S" != "off" ] && CMD="$CMD --gate_sat_s $GATE_SAT_S"
[ "$GATE_SAT_S" = "fixed" ] && CMD="$CMD --sat_s_thresh $SAT_S_THRESH"
[ "$GATE_SAT_S" = "worst1" ] && CMD="$CMD --sat_s_min_sat $SAT_S_MIN_SAT"

[ -n "$EPOCH_S_PCT" ] && CMD="$CMD --epoch_s_pct $EPOCH_S_PCT"
[ -n "$S_WEIGHT_CSV" ] && CMD="$CMD --s_weight_csv \"$S_WEIGHT_CSV\""
[ "$S_WEIGHT_MODE" != "fixed" ] && CMD="$CMD --s_weight_mode $S_WEIGHT_MODE"
[ "$S_SCORE_METHOD" != "wavg" ] && CMD="$CMD --s_score_method $S_SCORE_METHOD"

[ "$COMPARE" = "1" ] && CMD="$CMD --compare"

echo "正在执行: $CMD"
eval $CMD
