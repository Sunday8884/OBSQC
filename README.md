# OBSQC — GNSS 观测数据质量检查与综合评价工具

基于 RINEX 观测文件的质量检查（QC）与多指标综合评价（MADM）工具集。
支持**逐文件卫星级评估（prn 模式）**与**跨文件测站级评估（station 模式）**两种工作方式，
输出每个指标的图片、统计文本与 CSV 明细，并可对 QC 指标矩阵做综合评价与排名。

---

## 1. 目录结构与模块职责

```
OBSQC/
├── obsqc.py                # 主入口：命令行解析 + 模式路由
├── obsqc.sh                # Bash 封装脚本：集中配置参数，拼装命令行调用 obsqc.py
├── qc/
│   ├── qc_engine.py        # PRN 模式引擎：单文件 QC → PRN 级汇总 → 综合评价
│   ├── qc_engine_exp.py    # Station 模式引擎：多文件 QC → 测站级汇总 → SPP → 跨站综合评价
│   ├── satNum.py / freq.py / CNR.py / slip.py / phasenoise.py /
│   │   multipath.py / CMC.py / IOD.py / highOrderDiff.py   # 各 QC 指标计算与写文件
├── plot/                   # 各指标的绘图模块（plotSatNum / plotCnr / plotCycleSlip ...）
├── com/                    # 基础库：readObs / readNav / xyz2neu 等
├── spos/spp.py             # SPP 单点定位（station 模式专用）
└── madm/
    ├── EQCValue.py         # 综合评价调度层（正向化→预处理→PCA→权重→评价）
    ├── GRA.py              # 权重（entropy/ahp/bwm）+ GRA/TOPSIS/WASPAS 算法
    ├── INCOME.py           # INCOME（离散化 + kNN）评价方法
    └── PCA.py              # PCA 降维
```

---

## 2. 整体流程

```
obsqc.py（解析参数）
   │
   ├── mode = prn  ──► qc_engine.QCProcessor（逐文件循环）
   │                     读 RINEX → 逐 QC 模块 → _QCValue.csv（PRN 级）
   │                     → evaluate_qcvalue() → _eval_<权重>_<方法>.csv
   │
   └── mode = station / Sstation ──► qc_engine_exp.MultiFileQC
                         逐文件 QCProcessor（含 SPP）→ 每文件聚合成一行
                         → QCValue_comparison.csv
                           · station ：多站对比表（行=测站）
                           · Sstation：单站时间序列表（行=日期，station+date 两列）
                         → evaluate_matrix() → QCValue_eval_<权重>_<方法>.csv
```

### 2.1 三种模式对比

| 项目 | prn 模式（默认） | station 模式 | Sstation 模式 |
|---|---|---|---|
| 引擎 | `qc/qc_engine.py` | `qc/qc_engine_exp.py` | `qc/qc_engine_exp.py` |
| 评估对象 | 单文件内各卫星 PRN | 多个测站 | 同一测站的多个日期 |
| SPP 定位 | 不支持 | 支持（需 `--nav`，可选 `--snx`） | 支持（同 station） |
| 汇总输出 | `<name>_QCValue.csv`（PRN×指标） | `QCValue_comparison.csv`（行=测站） | `QCValue_comparison.csv`（行=日期，`station,date` 两列，按日期排序） |
| 综合评价 | 卫星间排名 | 测站间排名 | 日期间排名（时间序列质量评估） |

### 2.2 单文件 QC 流程（QCProcessor.process）

1. **读取观测文件**：`readObsHead` + `readObs`，记录历元数；
2. **按开关顺序执行 9 个 QC 模块**（每模块结束 `plt.close('all')` 释放绘图内存）；
3. station 模式且开启 SPP：执行 SPP 解算；
4. **汇总**：把各模块返回的 `{模块: {PRN: {波段: {val}}}}` 结构
   - 同 PRN 跨波段取均值 → PRN 级 `_QCValue.csv`；
   - 同模块跨 PRN 取均值 → 测站级 `_QCValue_station.csv`（仅 station 模式）；
5. **综合评价**（可选）。

### 2.3 QC 指标模块

| 模块开关 | 含义 | 输出文件后缀 |
|---|---|---|
| `--satnum` | 卫星数量统计 | `_sat.png/.txt` |
| `--freq` | 频段利用率统计 | `_freq.png/.txt` |
| `--cnr` | 载噪比分析 | `_cnr.png/.txt/.csv` |
| `--slip` | 周跳检测（TurboEdit） | `_slip.png/.txt/.csv` |
| `--phasenoise` | 相位噪声 | `_phasenoise.png/.txt/.csv` |
| `--multipath` | 多路径效应 | `_multipath.png/.txt/.csv` |
| `--cmc` | 历元间差分码减相位（CMC） | `_cmc.png/.txt/.csv` |
| `--iod` | 电离层延迟变化率 | `_iod.png/.txt/.csv` |
| `--highOrderDiff` | 载波相位高阶差分（4 阶） | `_highOrderDiff.png/.txt/.csv` |
| `--spp` | SPP 定位（仅 station 模式） | `_spp.txt/.png`，产出 `spp3drms`、`dop` 指标 |

所有开关默认为**启用**，可用 `--no-<name>` 单独关闭。

### 2.4 SPP 指标提取（station 模式）

- 波段：优先 `--bands`，否则使用默认值（G: C1C/C2W，C: C2I/C6I，E: C1X/C5X）；
- 导航文件匹配：多个 nav 文件时按文件名中的"年份_年积日"编码与 obs 文件配对；
- ENU 参考原点：优先 `--snx` SINEX 周解坐标（测站名取自 obs 文件名前 4 位），失败则回退到 SPP 坐标均值；
- 提取指标：**ENU 3D RMS**（相对参考原点的三维均方根误差）与 **DOP 平均值**（优先 PDOP，回退 GDOP，再回退 HDOP+VDOP），以伪 PRN `STATION` 写入 `qc_value`。

### 2.5 综合评价流水线（madm/EQCValue.py）

```
原始矩阵 (labels × modules)
  1. 正向化：按 DEFAULT_IMPACTS 把"越小越好"指标做 max-x 变换
  2. 预处理：NaN→列最小值、±inf→列极值、负值列整体平移
  3. PCA 降维（可选）：标准化 + 符号对齐，按累计贡献率阈值或指定主成分数
     输出 pca_report（贡献率 + 载荷矩阵 CSV），列名变为 PC1..PCk
  4. 权重计算：entropy / ahp / bwm / equal
  5. 评价打分：GRA（ρ）/ TOPSIS / WASPAS（λ）/ INCOME（离散化+kNN）
  6. 排名：按得分降序，写 <标签, 指标, Score, Rank> CSV + Weight 行
```

指标方向（`DEFAULT_IMPACTS`）：CNR、freq 为**越大越好**；CSR、multipath、phaseNoise、CMC、IOD、highOrder、spp3drms、dop 为**越小越好**。未匹配的列默认按"越大越好"。

---

## 3. 使用方法

### 3.1 命令行

```bash
# prn 模式：对单文件做卫星级 QC 与评估
python obsqc.py --input_obs "./test/data/*.24o" --output_dir ./test/qc --systems C

# station 模式：多站对比 + SPP + 跨站综合评价
python obsqc.py --mode station \
    --input_obs "./test/data/*.24o" \
    --nav "./test/data/*.24p" \
    --snx "./test/data/*.SNX" \
    --output_dir ./test/qc \
    --eval_weight ahp --eval_method income --pca --pca_var 0.95
```

### 3.2 主要参数

| 参数 | 说明 |
|---|---|
| `--mode` | `prn`（默认）/ `station` / `Sstation`（单站多日期时间序列） |
| `--input_obs` | RINEX 观测文件路径，**必填**，支持通配符 |
| `--nav` | 广播星历文件（通配符），仅 station 模式用于 SPP |
| `--snx` | SINEX 周解文件，提供 ENU 参考坐标（仅 station 模式） |
| `--output_dir` | 输出根目录，默认 `./qc` |
| `--start / --end` | 时间窗过滤（`YYYY-MM-DD HH:MM:SS`） |
| `--systems` | GNSS 系统筛选，如 `G C E` |
| `--satellites` | 卫星筛选，如 `G01 G02 C06` |
| `--bands` | 波段筛选，格式 `"G:1C,2W;C:2I,6I"` |
| `--evaluate / --no-evaluate` | 综合评价开关（默认开） |
| `--pca / --pca_var / --pca_n` | 评估前 PCA 降维（默认开，累计贡献率阈值 0.95） |
| `--eval_weight` | 权重方法：`entropy`（默认）/ `ahp` / `bwm` / `equal` |
| `--eval_method` | 评价方法：`gra`（默认）/ `topsis` / `waspas` / `income` |
| `--eval_rho` | GRA 分辨系数（默认 0.5） |
| `--eval_lambda` | WASPAS 系数（默认 0.5） |
| `--income_levels / --income_k / --income_binning` | INCOME 离散等级数 / kNN 邻居数 / 分箱方式 |

### 3.3 一键脚本

`obsqc.sh` 把所有参数集中在脚本头部变量中（输入路径、模块开关 0/1、评估参数），
修改后直接运行：

```bash
cd /f/OBSQC && bash obsqc.sh
```

脚本会自动拼装命令行：`--nav`/`--snx` 仅在 station 模式注入；INCOME 超参仅在 `--eval_method income` 时追加。

---

## 4. 输出说明

输出根目录下每个观测文件一个子目录（以文件名命名）：

```
<output_dir>/
├── <obs_name>/
│   ├── <obs_name>_sat.png / _freq.png / _cnr.png ...     # 各模块图
│   ├── <obs_name>_cnr.txt / _cnr.csv ...                 # 统计文本与明细
│   ├── <obs_name>_QCValue.csv                            # PRN×指标汇总
│   ├── <obs_name>_QCValue_station.csv                    # 测站级单行（station 模式）
│   ├── <obs_name>_spp.txt / _spp.png                     # SPP 结果（station 模式）
│   └── <obs_name>_eval_<w>_<m>[_pca].csv                 # 评估得分与排名（prn 模式）
├── QCValue_comparison.csv                                # 多站对比表（station 模式）
├── QCValue_eval_<w>_<m>[_pca].csv                        # 跨站评估结果（station 模式）
└── QCValue_eval_<w>_<m>_pca_report.csv                   # PCA 报告（开启 PCA 时）
```

评估 CSV 末尾附带 `Weight` 行记录各指标（或主成分）权重。

---

## 5. 依赖

- Python 3 + `numpy`、`matplotlib`
- 内部包：`com`（RINEX 读取、坐标转换）、`qc`（指标计算）、`plot`（绘图）、`spos`（SPP）、`madm`（综合评价）
- 绘图使用 `Agg` 后端，无需图形界面，可在服务器/批处理环境运行

## 6. 备注

- 断点续跑：station/Sstation 模式下，若文件子目录中已有 `<name>_QCValue_station.csv`，重跑时直接读取该结果计入对比表并跳过计算（prn 模式检查 `<name>_QCValue.csv`）；坏文件只警告跳过，不中断整批。如需强制重算，删除对应子目录或汇总文件即可。
- station 模式的快速通道：当全部 QC 模块关闭且 `QCValue_comparison.csv` 已存在时，直接读表做综合评价，便于只调整评估参数反复实验。
- 综合评价需要至少 2 个评估对象（测站/PRN），否则自动跳过。
- 文件匹配导航星历依赖文件名中的年积日编码（如 `wuh20050.24o` ↔ `wuh20050.24p`）。
