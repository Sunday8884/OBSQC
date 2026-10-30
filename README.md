# OBSQC — GNSS 观测数据质量检查工具

基于 RINEX 观测文件的质量检查（QC）工具,输出每个指标的图片、统计文本与 CSV 明细。

---

## 1. 目录结构与模块职责

```
OBSQC/
├── obsqc.py                # 主入口：命令行解析 
├── obsqc.sh                # Bash 封装脚本：集中配置参数，拼装命令行调用 obsqc.py
├── obsqc_batch.sh          # Bash 封装脚本: 批处理脚本
├── qc/ 
│   ├── qc_common.py        # 公共基础设施：QCProcessorBase 基类 + 波段归一化
│   ├── qc_engine.py        # PRN 模式引擎：单文件 QC → PRN 级汇总
│   ├── qc_engine_exp.py    # Station 模式引擎：多文件 QC → 测站级汇总 → SPP → 多站对比表
│   ├── qc_time.py          # Time 模式引擎：多文件 QC → 逐历元时间序列 QCValue.csv
│   ├── satNum.py / freq.py / completeness.py / CNR.py / slip.py / phasenoise.py /
│   │   multipath.py / CMC.py / IOD.py / highOrderDiff.py   # 各 QC 指标计算与写文件
├── plot/                   # 各指标的绘图模块（plotSatNum / plotCnr / plotCycleSlip ...）
├── com/                    # 基础库：readObs / readNav / xyz2neu 等
├── spos/spp.py             # SPP 单点定位
├── spp_main.py             # spp计算入口：命令行解析 
├── spp.sh                  # spp封装脚本：集中配置参数，拼装命令行调用 spp_main.py 

```

### 1.1 引擎分层

`qc_engine.py` 与 `qc_engine_exp.py` 中的 `QCProcessor` 共享同一套执行骨架，
该骨架抽取在 `qc/qc_common.py` 的 `QCProcessorBase` 中：

```
qc_common.QCProcessorBase          ← 共享骨架
├── load_obs()                     读 RINEX 头 + 数据
├── run_modules()                  按开关顺序执行 10 个 QC 模块
├── run_satnum / run_freq / run_completeness / run_cnr / run_slip / run_phasenoise
├── run_multipath / run_cmc / run_iod / run_highOrderDiff
├── _aggregate_prn_module_val()    PRN 跨波段取均值
├── _write_prn_summary()           写 _QCValue.csv
└── write_summary() / process()    ← 子类覆写

├── qc_engine.QCProcessor          仅覆写输出：PRN 级 _QCValue.csv
└── qc_engine_exp.QCProcessor      覆写输出：PRN 级 + 测站级 _QCValue_station.csv
    └── MultiFileQC                多文件编排 + QCValue_comparison.csv

qc_time.TimeSeriesQC               复用 qc_engine_exp.QCProcessor 跑完单文件 QC，
                                   再读各模块逐历元明细 CSV → 按 SPP 频段过滤
                                   → 逐历元汇总 → QCValue.csv + QCValue_sppband.txt
```

两个引擎的差异仅在输出组织；`qc_engine_exp.py` 额外承担 SPP（ENU 3D RMS、DOP 均值）
与测站级聚合、多站对比表。

---

## 2. 整体流程

```
obsqc.py（解析参数）
   │
   ├── mode = prn  ──► qc_engine.QCProcessor（逐文件循环）
   │                     读 RINEX → 逐 QC 模块 → _QCValue.csv（PRN 级）
   │
   ├── mode = station / Sstation ──► qc_engine_exp.MultiFileQC
   │                     逐文件 QCProcessor（含 SPP）→ 每文件聚合成一行
   │                     → QCValue_comparison.csv
   │                       · station ：多站对比表（行=测站）
   │                       · Sstation：单站时间序列表（行=日期，station+date 两列）
   │
   └── mode = time ──► qc.qc_time.TimeSeriesQC
                         逐文件 QCProcessor（含 SPP）→ 读取各模块逐历元明细 CSV
                         → 按「SPP 实际选用频段」过滤 → 逐历元汇总
                         → QCValue_<标签>.csv（行=历元）+ QCValue_sppband.txt
                         → --sat 时另出 QCValue_<sat>_<标签>.csv（按频段分列）
```

### 2.1 四种模式对比

| 项目 | prn 模式（默认） | station 模式 | Sstation 模式 | time 模式 |
|---|---|---|---|---|
| 引擎 | `qc/qc_engine.py` | `qc/qc_engine_exp.py` | `qc/qc_engine_exp.py` | `qc/qc_time.py` |
| 评估对象 | 单文件内各卫星 PRN | 多个测站 | 同一测站的多个日期 | 单站多日期的**逐历元**时序 |
| SPP 定位 | 不支持 | 支持（需 `--nav`，可选 `--snx`） | 支持（同 station） | 支持（同 station） |
| 汇总输出 | `<name>_QCValue.csv`（PRN×指标） | `QCValue_comparison.csv`（行=测站） | `QCValue_comparison.csv`（行=日期，`station,date` 两列，按日期排序） | `QCValue_<标签>.csv`（行=历元）+ `QCValue_sppband.txt` |
| 频段范围 | 全部 | 全部 | 全部 | **仅 SPP 实际选用的频段** |
| 选星输出 | 无 | 无 | 无 | `--sat` → `QCValue_<sat>_<标签>.csv`（按频段分列） |

> **关于综合评价**：MADM 综合评价（正向化 → 预处理 → PCA → 权重 → 评价）已从本仓库拆出，
> 由独立项目读取本工具产出的 `_QCValue.csv` / `QCValue_comparison.csv` 完成。
> 本仓库只负责 QC + SPP。

### 2.2 单文件 QC 流程（QCProcessor.process）

1. **读取观测文件**：`readObsHead` + `readObs`，记录历元数；
2. **按开关顺序执行 9 个 QC 模块**（每模块结束 `plt.close('all')` 释放绘图内存）；
3. station 模式且开启 SPP：执行 SPP 解算；
4. **汇总**：把各模块返回的 `{模块: {PRN: {波段: {val}}}}` 结构
   - 同 PRN 跨波段取均值 → PRN 级 `_QCValue.csv`；
   - 同模块跨 PRN 取均值 → 测站级 `_QCValue_station.csv`（仅 station 模式）；
5. **综合评价**（可选）。

### 2.2.1 time 模式：逐历元时间序列

`time` 模式在上述 QC 流程之后追加一步「明细归并」：

```
各模块逐历元明细 CSV                SPP 选段（_sppband.csv）
_cnr.csv  _slip.csv  _multipath.csv      C: C2I+C6I
_cmc.csv  _iod.csv  _phasenoise.csv      E: C1C+C5Q
_highOrderDiff.csv  _completeness.csv    G: C1C+C2W
        │                                      │
        └──────── 频段过滤（只保留 SPP 频段）───┘
                          │
                          ▼
              逐历元汇总为一行（均值/RMS/计数）
                          │
                          ▼
        QCValue.csv（行=历元） + QCValue_sppband.txt
```

要点：

- **频段过滤**：各模块明细 CSV 的列名形如 `C19_C2I` / `C19_L2I` / `C05_S2I`，
  去掉观测类型前缀后得到频段键（`2I`），只保留 SPP 选段涉及的频段键；
  非 SPP 频点的列整列丢弃。
- **时间轴**：`Epoch` 取各模块 CSV 中的实际时间戳；`Idx` 为相对时间原点的
  采样序号，原点默认取**数据中最早的历元**（`--time_origin` 可指定），
  采样间隔由实际时间戳的中位差分推导（`--time_interval` 可覆盖）。
- **统计量**：CNR 取均值/最小值；多路径/CMC/IOD/相位噪声/高阶差分取 RMS/峰值；
  周跳与完整率取非零计数/有效样本数。每列都带 `_n` 有效样本数。
- **SPP 频段统计汇总**（`QCValue_sppband.txt`）：各系统 SPP 选段的观测码对、
  中心频率、Δf、λW、MW 放大倍数，以及逐 PRN 的该频段观测次数与有效率。

#### 2.2.2 time 模式的文件命名与选星参数

**标签（`<标签>`）**：由各系统的「系统字母 + SPP 选段各码的频段键」组成，
系统内**双频用 `_` 连接**，各系统之间也用 `_` 连接，系统按字母序排列。
系统字母只与第一个频段键直接相连（保持 `C2I` 这样的完整形式）：

| 系统 | SPP 选段 | 标签 |
|---|---|---|
| BDS | `C2I` + `C6I` | `C2I_6I` |
| Galileo | `C1C` + `C5Q` | `E1C_5Q` |
| GPS | `C1C` + `C2W` | `G1C_2W` |

三种系统共存时，全局标签为 `C2I_6I_E1C_5Q_G1C_2W`：

| 文件 | 内容 |
|---|---|
| `QCValue_C2I_6I_E1C_5Q_G1C_2W.csv` | 全卫星逐历元时间序列（各频段合并统计） |
| `QCValue_C19_C2I_6I.csv` | 单星（C19，BDS 系统）逐历元时间序列 |
| `QCValue_G01_G1C_2W.csv` | 单星（G01，GPS 系统）逐历元时间序列 |
| `QCValue_sppband.txt` | SPP 频段统计汇总（不带标签，全局唯一） |

单星文件的标签只取**该星所属系统**的标识，因此不同系统的同名卫星不会冲突。

**选星参数 `--sat`**（仅 time 模式）：为指定的每颗卫星单独输出一份时间序列表。

```bash
python obsqc.py --mode time \
    --input_obs "./test/data/*.24o" \
    --nav "./test/data/*.24p" \
    --output_dir ./test/qc_time \
    --sat C19 G01
```

单星表的统计口径与全星表**不同**：单星表**按频段分列**，双频各出一组统计量，
不做跨频段平均；与频段无关的指标（周跳 CSR 等）只出一列，不带频段后缀：

```
CNR_2I      CNR_6I                   ← B1 / B3 频段（B2I / B3I）
CSR                                ← 周跳按卫星整体判定，与频段无关
MP_2I        MP_6I              ← 多路径 RMS
CMC_2I       CMC_6I             ← CMC RMS
COMP                             ← 弧段完整率(0~1)，不按频段分列，各行取值相同
SAT                             ← 统计范围内各系统卫星数之和
SPP_E   SPP_N   SPP_U   SPP_3D   ← 整站量，置于最右
```

每个模块只出**一列主统计量**，列名不带统计量后缀：

| 模块 | 统计量 | 列名示例 |
|---|---|---|
| CNR | 均值 | `CNR` / `CNR_2I` |
| MP / CMC / IOD / PN / HOD | RMS | `MP` / `MP_2I` |
| CSR | 计数 | `CSR` |
| COMP | 弧段完整率 (0~1) | `COMP`（不按频段分列） |
| SAT | 卫星数之和 | `SAT`（不按频段分列） |

时间列为 `Epoch`（实际时间戳）与 `Idx`（采样序号），不含相对秒数。

`--sat` 与 `--satellites` 的区别：`--satellites` 会**过滤各 QC 模块本身的计算
范围**（影响所有输出文件），`--sat` 只影响**单星汇总表**的生成，不改变全星表
与各模块明细。若某颗星在 SPP 频段内无数据，则跳过其单星表并给出提示。

### 2.3 QC 指标模块

| 模块开关 | 含义 | 输出文件后缀 |
|---|---|---|
| `--satnum` | 卫星数量统计 | `_sat.png/.txt/.csv` |
| `--freq` | 频段利用率统计 | `_freq.png/.txt` |
| `--completeness` | 卫星数据完整率统计 | `_completeness.png/.txt/.csv`、`_gaps.png` |
| `--cnr` | 载噪比分析 | `_cnr.png/.txt/.csv` |
| `--slip` | 周跳检测（TurboEdit） | `_slip.png/.txt/.csv` |
| `--phasenoise` | 相位噪声 | `_phasenoise.png/.txt/.csv` |
| `--multipath` | 多路径效应 | `_multipath.png/.txt/.csv` |
| `--cmc` | 历元间差分码减相位（CMC） | `_cmc.png/.txt/.csv` |
| `--iod` | 电离层延迟变化率 | `_iod.png/.txt/.csv` |
| `--highOrderDiff` | 载波相位高阶差分（4 阶） | `_highOrderDiff.png/.txt/.csv` |
| `--spp` | SPP 定位（仅 station 模式） | `_spp.txt/.png`，产出 `spp3drms`、`dop` 指标 |

所有开关默认为**启用**，可用 `--no-<name>` 单独关闭。

### 2.4 各模块的频段（波段）是如何确定的

各模块对"用哪个频段"的处理方式并不相同，归纳如下。

#### 2.4.1 三类频段匹配方式

| 匹配方式 | 依据 | 使用该方式的模块 |
|---|---|---|
| **观测头全量** | 遍历观测头 `OBS TYPES`，不使用任何筛选 | `satnum` |
| **载波观测码** | 按 `L` + 频段键匹配（如 `L2I`），需该系统存在对应载波观测码 | `freq` |
| **载噪比观测码** | 按 `S` + 频段键匹配（如 `S2I`） | `cnr` |
| **频段键** | 按 `band[1:]` 匹配（如伪距码 `C2I` → 键 `2I`） | `multipath`、`cmc`、`iod`、`phasenoise`、`highOrderDiff` |
| **频段键 + 载波码命名** | 按频段键分组，但 Band 列以载波观测码全称输出（如 `L2I`） | `completeness` |
| **固定优先级组合表** | 见下方 `bandComb` | `slip` |

#### 2.4.2 `--bands` 参数的语义

`--bands` 使用**频段键**格式（不含观测类型前缀）：

```bash
--bands "G:1C,2W;C:2I,6I"
```

其中 `1C` 表示 L1 频点 + C 跟踪模式，与 RINEX 观测码 `L1C`（载波）、`C1C`（伪距）、`S1C`（载噪比）一一对应。
为兼容习惯，也允许直接写完整观测码 `--bands "G:L1C,L2W"`，程序会自动去掉前缀归一化为频段键。

归一化在 `qc/qc_common.py` 的 `parse_band_string()` 中完成，随后由 `band_choose_variants()`
分发出三种形式发给不同模块，保证 `--bands` 对所有模块一致生效：

| 形式 | 示例 | 用途 |
|---|---|---|
| 频段键 | `{'C': ['2I','6I']}` | multipath / cmc / iod / phasenoise / highOrderDiff / completeness |
| 载波观测码 | `{'C': ['L2I','L6I']}` | freq |
| 载噪比观测码 | `{'C': ['S2I','S6I']}` | cnr |

> 早期版本中 `--bands` 直接透传给各模块，导致 `freq` 匹配不到任何 `L` 波段、
> `cnr` 匹配不到任何 `S` 观测码，二者实际失效；现已统一归一化。

#### 2.4.3 双频模块（multipath / cmc / iod / phasenoise / highOrderDiff）的配对规则

这 5 个模块都基于**双频电离层无关（GF）组合**，流程一致：

1. **候选伪距码**：扫描观测头，取所有以 `C` 或 `P` 开头的观测码，按**频点号**（`band[1]`，如 `1`/`2`/`5`/`6`）分组；
2. **双频要求**：该系统至少要有 2 个不同频点，否则该系统整体跳过；
3. **配对**：以观测头中**最先出现的两个频点**为一组——
   - 若当前伪距码属于第 1 个频点，则搭配第 2 个频点的第一个伪距码；
   - 否则搭配第 1 个频点的第一个伪距码；
4. **载波配对**：载波码由伪距码直接推导（`C2I` → `L2I`），`L` 码或第 2 频点码不存在时跳过该组合；
5. **频率与波长**：由 `com/gnssParameter.py` 的 `getBandFreq()` 按系统查表得到（如 `C` 系统的 `2I` → B2a = 1207.14 MHz）；
6. **组合系数**：`COEFL1 = -(f1²+f2²)/(f1²-f2²)`，`COEFL2 = 2f2²/(f1²-f2²)`，另有 `lamdaW = c/(f1-f2)`。

**注意**：`--bands` 只筛选第 1 频点（`band`），**不限制配对的第 2 频点**——第 2 频点始终取自观测头顺序。
因此 `--bands "C:2I,6I"` 表示"以 2I/6I 作为主频点"，其搭配频点仍是观测头中最先出现的那个。

#### 2.4.4 slip（TurboEdit）的频段选择

`slip` 不使用 `--bands`，而是采用 `com/gnssParameter.py` 中的固定优先级组合表 `bandComb`：

```python
bandComb = {'G': [[1,2],[1,5]], 'C': [[2,6],[1,5]], 'E': [[1,5]],
            'S': [[1,5]], 'I': [[5,9]], 'J': [[1,2],[1,5]], 'W': [[1,2],[1,5]]}
```

选取逻辑：

1. 仅取载波观测码 `L*`，且要求该频点存在**同频点的伪距码**（`C` 或 `P`）才纳入候选；
2. 按 `bandComb` 的优先级逐个尝试，**第一个两频点都存在且频率可查到的组合即被采用**（`break` 后不再回退）；
3. 例如 GPS 优先取 L1/L2，不可用时退到 L1/L5；BDS 优先取 B2/B3（即 `2I`/`6I`），其次 B1/B2a。

#### 2.4.5 SPP 的频段选择

SPP 需要**完整伪距观测码**（如 `C1C`/`C2W`），与 QC 模块的频段键不同：

- **默认**：按 `DEFAULT_SPP_BANDS` 的候选配对优先级检测，当前配对在观测头中不可用时自动换用下一优先级
  （G: C1C/C2W → C1C/C2L → C1W/C2W；C: C2I/C6I → C2I/C7I → C1X/C5X；E: C1X/C5X → C1C/C5Q → C1X/C7X）；
- **指定 `--bands`**：按频段键在观测头中还原为实际存在的伪距观测码（`C` 优先于 `P`，
  同频点号可放宽匹配如 `2I`→`2X`），每系统取前两个可用码组成双频；不足两个则跳过该系统。

#### 2.4.6 `BANDS` 为空（未指定 `--bands`）时的自动选段逻辑

`--bands` 留空时**不是"不做处理"**，而是各模块退回各自的自动选段策略。
下表以实测文件 `wuh20050.24o`（RINEX3 MIXED，含 G/C/E/I/J/R/S 七系统）为例。

| 模块 | 空 BANDS 时的行为 | 实测结果（wuh20050.24o） |
|---|---|---|
| `satnum` | 统计全部系统全部卫星，不涉及频段 | 与频段无关 |
| `freq` | **自动使用观测头中所有 `L` 波段** | G: L1C/L1W/L1X/L2W/L2X/L5X；C: 8 个；E: 5 个… |
| `completeness` | **按载波观测码全称逐个统计**（同频点不同跟踪模式不合并） | G: L1C/L1W/L1X/L2W/L2X/L5X；C: L1X/L2I/L5X/L6I/L7D/L7I/L7Z/L8X；E: 5 个… |
| `cnr` | **不加任何过滤**，所有 `S` 观测码全部统计 | G: S1C/S1W/S1X/S2W/S2X/S5X；C: 8 个… |
| `multipath` `cmc` `iod` `phasenoise` `highOrderDiff` | 见下方"双频自动配对" | G: C1C+C2W；C: C1X+C2I；E: C1X+C5X |
| `slip` | 走 `bandComb` 固定优先级表 | G: L1C+L2W；C: **L2I+L6I**；E: L1X+L5X |
| `spp` | 走 `DEFAULT_SPP_BANDS` 优先级表 | G: C1C/C2W；C: C2I/C6I；E: C1X/C5X |

**双频模块的自动配对（BANDS 为空）**：

1. 扫描该系统观测头中所有 `C`/`P` 开头的伪距码，按频点号 `band[1]` 分组，**保持观测头中的出现顺序**；
2. 取**前两个频点组**作为配对基准（第 1 频点的首个码、第 2 频点的首个码）；
3. 遍历该系统每个伪距码：若它属于第 1 频点组，则配第 2 频点的首个码；否则一律配第 1 频点的首个码；
4. 载波码由伪距码直接推导（`C2I` → `L2I`），缺任一 `L` 码则跳过该组合；
5. 该系统伪距频点 < 2 个，或系统不在 `obsFreq` 频率表中，则该系统整体跳过。

实测（G 系统，伪距分组 `1:[C1C,C1W,C1X] 2:[C2W,C2X] 5:[C5X]`）：

| 主频点 | 配对频点 | 载波组合 | 频率 (MHz) |
|---|---|---|---|
| C1C | C2W | L1C / L2W | 1575.42 / 1227.60 |
| C1W | C2W | L1W / L2W | 1575.42 / 1227.60 |
| C1X | C2W | L1X / L2W | 1575.42 / 1227.60 |
| C2W | C1C | L2W / L1C | 1227.60 / 1575.42 |
| C2X | C1C | L2X / L1C | 1227.60 / 1575.42 |
| C5X | C1C | L5X / L1C | 1176.45 / 1575.42 |

可见第 1 频点组内的 3 个码（C1C/C1W/C1X）都各自作为"主频点"参与计算，
即**同一频点的不同跟踪模式会被分别统计**。

**空 BANDS 时的三个要点**：

1. **配对基准依赖观测头的码序**，不同文件可能不同。实测 `brux0050.24o` 的 GPS 为
   `C1C + C2L`（而非 `wuh20050.24o` 的 `C1C + C2W`），
   因为其观测头中 L2 频点的首个码是 `L2L`。
2. **`slip` 与双频模块的选段结果不同**：BDS 在 `slip` 中用 `L2I+L6I`（B2+B3），
   而在双频模块中用 `C1X+C2I`（B1C+B1），因为二者的候选码来源与优先级表不同。
3. **R（GLONASS）系统被跳过**：`obsFreq` 频率表中没有 `R`（GLONASS 为 FDMA 单频，
   无双频组合），故双频模块不处理；但 `freq`/`cnr`/`satnum` 仍会统计 R 的
   `L1C/L2C/L3X`、`S1C/S2C/S3X` 等。`I`（BeiDou-3 惯性解）/ `S`（SBAS）因只有
   单频，同样只在 `freq`/`cnr`/`satnum` 中出现。

### 2.5 卫星数据完整率（completeness 模块）

#### 2.5.1 为什么不能用 `freq` 代替

`freq` 模块给出的是「某频段观测次数 / 文件总历元数」。这个值会被
**卫星本身不可见的时间严重稀释**，无法反映真实的数据缺口。

实测 `wuh20050.24o` 的 G01：

| 指标 | 数值 | 含义 |
|---|---|---|
| `freq` 的出现频率 | 670 / 2880 = **23.26%** | 含大量"卫星本就在地平线以下"的时间 |
| 完整率（弧段口径） | 670 / 698 = **95.99%** | 卫星可见期内，数据实际完好的比例 |

前者低不代表数据有问题，后者才是数据质量。

#### 2.5.2 核心概念：可见弧段（presence arc）

```
presence arc = 该 PRN 首个「有值历元」~ 末个「有值历元」之间的所有历元
```

- **有值历元**：该历元中该 PRN 至少有一个观测码的数值非 `None`；
- 弧段**之外**的历元视为"卫星不可见"，不计入缺测；
- 弧段以 PRN 为单位（而非频段），因为一颗卫星是否可见与频段无关。

这样分母才真正对应"应该有数据的那段时间"。

#### 2.5.3 输出的指标

| 指标 | 定义 | 说明 |
|---|---|---|
| `R_arc` | `N_band / N_arc` | **主指标**，弧段内该频段有值历元数 / 弧段历元数 |
| `R_span` | `N_band / N_span` | 该频段自身首末有值历元之间的完好率（与 `R_arc` 互补） |
| `R_present` | `N_band / N_present` | 频段层面缺测率 |
| `R_total` | `N_band / total_epochs` | 相对整个观测时段的数据覆盖率（≈ `freq` 的口径） |
| `R_visible` | `N_band / N_visible` | **理论可见口径**，需 `--nav`，见 2.5.5 |
| `gap_count` | 连续缺测段数 | 连续 ≥2 历元无值记为一段缺口 |
| `max_gap` / `max_gap_sec` | 最长缺测历元数 / 时长 | 定位最严重的数据中断 |

#### 2.5.4 频段粒度与 Band 列的含义

输出的 `Band` 列是**载波观测码全称**（如 `L2I`、`L1C`、`L7D`），与 RINEX
观测码书写一致，可直接与观测头 `OBS TYPES` 对照。

- **分组粒度 = 频点号 + 跟踪模式**（band key，如 `2I`），**不跨跟踪模式合并**。
  例如 GPS 的 `L1C` / `L1W` / `L1X` 分列为三行，各自独立统计完整率。
- **同组内的观测码**：同一 band key 下的全部观测码归为一组
  （伪距 `C` / 伪距 `P` / 载波 `L` / 载噪比 `S` / 多普勒 `D`），
  判定"有值"时按 `C` > `P` > `L` > `S` > `D` 优先级取第一个非 `None` 值。
  例如 `L2I` 组含 `C2I` / `L2I` / `S2I` / `D2I`。
- **指定 `--bands`**：严格按用户给定的频段键分组，Band 名取 `L` + 频段键
  （`--bands "C:2I,6I"` → `L2I` / `L6I`）。
- **自动剔除空频段**：整个时段无任何数据的频段不参与统计。
  RINEX3 MIXED 的观测头会声明全部可能类型，但接收机实际只提供一部分
  （如 BDS 各星的 `L1X`/`C5X` 全空），若保留会产生大量 0% 的假告警。
- **`Freq(MHz)` 列**：由 `obsCodes`/`obsFreq` 静态表反查该观测码的中心频率，
  便于识别同名不同频的情况——BDS 尤其需要留意：
  `L1X` → B1C 1575.42、`L2I` → **B1** 1561.098、`L5X` → B2a 1176.45、
  `L6I` → B3 1268.52、`L7I` → B2 1207.14、`L7D`/`L7Z` → B2b 1207.14、
  `L8X` → B2a+B2b 1191.795。也就是说 BDS 的 `L1X` 与 `L2I` 属同一代
  B1 频段族，频差仅 14.3 MHz。

实测输出示例：

```
   PRN     Band  Freq(MHz)   ObsNum   ArcNum  SpanNum   Arc(%)  Span(%)   GapCnt   MaxGap   MaxGapSec
  G01      L1C    1575.420      670      698      670   95.99   100.00        1       28       840.0
  G01      L1W    1575.420      670      698      670   95.99   100.00        1       28       840.0
  G01      L2X    1227.600      690      698      690   98.85   100.00        1        8       240.0
  C01      L2I    1561.098     2880     2880     2880  100.00   100.00        0        0         0.0
```

#### 2.5.5 理论可见口径（需 `--nav`）

提供广播星历时，额外按卫星高度角（`--cutoff`，默认 10°）计算理论可见历元：

- 星历不可用（无对应星历 / `SVhealth` 异常）的历元**排除在分母之外**，
  否则会把"星历问题"误判成"数据缺失"；
- `R_visible` 允许略微超过 100%——广播星历间隔大（常为 2 小时）且部分卫星
  `SVhealth` 异常，插值精度有限，星历口径偏保守；
- 实测 `wuh20050.24o` 有 165 项 `R_visible > 100%`，最大偏差约 3%。

#### 2.5.6 输出文件

| 文件 | 内容 |
|---|---|
| `_completeness.txt` | 系统级 + 卫星级汇总表（`ObsNum`/`ArcNum`/`SpanNum`/`Arc(%)`/`Span(%)`/`GapCnt`/`MaxGap`）。`+SYS` 段的 `BandNum` 是该系统下的 **(PRN, Band) 组合条目数**，不是频段数 |
| `_completeness.csv` | 逐历元的 0/1 标记宽表，用于二次分析与绘图 |
| `_completeness.png` | 完整率柱状图（按系统分面，低于 90% 标红，虚线为 90% 参考线，不标注柱顶数值） |
| `_gaps.png` | 缺口时间线图（横轴时间，红线为缺口持续时段） |

`completeness` 指标以 `R_arc` 参与 `_QCValue.csv` 的 PRN 级汇总与
`_QCValue_station.csv` 的测站级汇总（跨 PRN 取均值）。

### 2.6 SPP 指标提取（station 模式）

- 波段：优先 `--bands`，否则使用默认值（G: C1C/C2W，C: C2I/C6I，E: C1X/C5X）；
- 导航文件匹配：多个 nav 文件时按文件名中的"年份_年积日"编码与 obs 文件配对；
- ENU 参考原点：优先 `--snx` SINEX 周解坐标（测站名取自 obs 文件名前 4 位），失败则回退到 SPP 坐标均值；
- 提取指标：**ENU 3D RMS**（相对参考原点的三维均方根误差）与 **DOP 平均值**（优先 PDOP，回退 GDOP，再回退 HDOP+VDOP），以伪 PRN `STATION` 写入 `qc_value`。


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
    --output_dir ./test/qc 

# 仅统计卫星数据完整率（配合 --nav 可得"理论可见"口径）
python obsqc.py --mode prn \
    --input_obs "./test/data/wuh20050.24o" \
    --nav "./test/data/brdc0050.24p" \
    --cutoff 10 \
    --output_dir ./test/qc \
    --no-satnum --no-freq --no-cnr --no-slip --no-phasenoise \
    --no-multipath --no-cmc --no-iod --no-highOrderDiff

# time 模式：逐历元时间序列（只统计 SPP 实际选用频段的相关数据）
python obsqc.py --mode time \
    --input_obs "./test/data/*.24o" \
    --nav "./test/data/*.24p" \
    --output_dir ./test/qc_time

# time 模式 + 选星：额外输出单星汇总表 _QCValue_C19_C2I.csv
python obsqc.py --mode time \
    --input_obs "./test/data/*.24o" \
    --nav "./test/data/*.24p" \
    --output_dir ./test/qc_time \
    --sat C19 G01
```

### 3.2 主要参数

| 参数 | 说明 |
|---|---|
| `--mode` | `prn`（默认）/ `station` / `Sstation`（单站多日期时间序列）/ `time`（逐历元时间序列） |
| `--input_obs` | RINEX 观测文件路径，**必填**，支持通配符 |
| `--nav` | 广播星历文件（通配符）。station/time 模式用于 SPP；三种模式均用于完整率的"理论可见"口径 |
| `--snx` | SINEX 周解文件，提供 ENU 参考坐标（仅 station 模式） |
| `--cutoff` | 理论可见性判定的截止高度角（度），默认 `10`；仅在提供 `--nav` 时生效 |
| `--output_dir` | 输出根目录，默认 `./qc` |
| `--start / --end` | 时间窗过滤（`YYYY-MM-DD HH:MM:SS`） |
| `--systems` | GNSS 系统筛选，如 `G C E` |
| `--satellites` | 卫星筛选，如 `G01 G02 C06` |
| `--bands` | 波段筛选，格式 `"G:1C,2W;C:2I,6I"` |
| `--sat` | 选星参数（**仅 time 模式**），如 `C19 G01`。为指定的每颗卫星单独输出 `_QCValue_<sat>_<标签>.csv`，不改变全星表与各模块明细 |
| `--time_origin` | time 模式时间原点；留空则由数据中最早的**实际历元**自动推导 |
| `--time_interval` | time 模式采样间隔（秒）；留空则由实际时间戳的中位差分推导 |

### 3.3 一键脚本

`obsqc.sh` 把所有参数集中在脚本头部变量中（输入路径、模块开关 0/1、评估参数），
修改后直接运行：

```bash
cd /f/OBSQC && bash obsqc.sh
```

脚本会自动拼装命令行：`--nav`/`--snx` 仅在 station/Sstation/time 模式注入。

---

## 4. 输出说明

输出根目录下每个观测文件一个子目录（以文件名命名）：

```
<output_dir>/
├── <obs_name>/
│   ├── <obs_name>_sat.png / _freq.png / _cnr.png ...     # 各模块图
│   ├── <obs_name>_cnr.txt / _cnr.csv ...                 # 统计文本与明细
│   ├── <obs_name>_QCValue.csv                            # PRN×指标汇总
│   ├── <obs_name>_QCValue_station.csv                    # 测站级单行（station/Sstation/time 模式）
│   ├── <obs_name>_spp.txt / _spp.png                     # SPP 结果（station/Sstation/time 模式）
│   ├── <obs_name>_sppband.csv                            # SPP 选段记录（time 模式）
│   ├── <obs_name>_sppband_freq.csv                       # SPP 频段逐 PRN 观测计数（time 模式）
├── QCValue_comparison.csv                                # 多站对比表（station/Sstation 模式）
├── QCValue_<标签>.csv                                    # 逐历元时间序列（time 模式，全卫星合并统计）
├── QCValue_<sat>_<标签>.csv                              # 单星逐历元时间序列（time 模式，--sat，按频段分列）
├── QCValue_sppband.txt                                   # SPP 频段统计汇总（time 模式）
```

time 模式的 `QCValue_<标签>.csv`（全卫星）表头：

| 列 | 含义 |
|---|---|
| `STATION` | 测站名，取观测文件名的**前 4 位**（大写），逐历元标记数据来源 |
| `Epoch` | 实际时间戳（取自各模块明细 CSV） |
| `Idx` | 采样序号 = 相对时间原点秒数 / 采样间隔 |
| `CNR` | 载噪比均值 (dB-Hz) |
| `CSR` | 周跳次数 |
| `MP` | 多路径 RMS (m) |
| `CMC` | 历元间差分码减相位 RMS (m) |
| `IOD` | 电离层延迟变化率 RMS |
| `PN` | 相位噪声 RMS |
| `HOD` | 载波相位高阶差分 RMS |
| `COMP` | 弧段完整率 (0~1)，整段弧段统计量，各行取值相同；单星表为该星各频段均值 |
| `SAT` | 逐历元统计范围内各系统的可见卫星数之和（`--systems` 限定单系统时即该系统数量） |
| `SPP_E` / `SPP_N` / `SPP_U` | 逐历元单点定位的东/北/天残差 (m)，相对参考原点 |
| `SPP_3D` | 逐历元三维误差 sqrt(E²+N²+U²) (m) |

SPP 与卫星数量都是**整站量**，与卫星无关，因此全星表与单星表的这些列取值相同。
SPP 结果取自各文件的 `_spp.txt`，卫星数量取自 `_sat.csv`（未提供 `--nav`
或对应模块未启用时无相应列）。

> **完整率的口径**：分母为该星**实际可见**的历元数，且按**连续可见段**切分——
> 卫星中途长时间消失（换星座、遮挡、信号中断）后又重新出现时，中间不可见期
> 不计入分母。否则该段会被误判为缺测，把完整率严重拉低（实测 C19 实际可见
> 1165 历元，若跨越 878 历元的中断只算出 57%，实际应为 99.7%）。

> **卫星数量的系统口径**：`_sat.csv` 的各系统列统计范围**与 `--systems` 一致**——
> 只统计被选中的系统；而其 `ALL` 列**统计所有 PRN**（不过滤系统），
> 与 `_sat.txt` / `_sat.png` 中 `plotSatNum` 的 ALL 口径一致。
> 默认系统列表为 `G C R E L S J W`（不含 BeiDou-3 惯性解 `I`），
> 故默认下各系统列之和小于 `ALL`。
>
> 汇总表 `QCValue_*.csv` 只出**一列 `SAT`**，取值为统计范围内各系统之和
> （不含 `ALL`）；`--systems` 限定单个系统时，即为该系统的卫星数。

单星表 `QCValue_<sat>_<标签>.csv` 的列名在此基础上**按频段分列**（后缀为频段键），
与频段无关的指标不带后缀，详见 2.2.2 节。

> **例外：`COMP` 与 `SAT` 不按频段分列**
> - `COMP`（完整率）是**整段弧段**的统计量，双频取值几乎相同
>   （例：C19 的 L2I 56.79% / L6I 56.94%），拆列无意义，故只出一列，
>   值为该星各频段完整率的均值
> - `SAT`（卫星数）同理，只出一列，取统计范围内各系统之和

## 5. 依赖

- Python 3 + `numpy`、`matplotlib`
- 内部包：`com`（RINEX 读取、坐标转换）、`qc`（指标计算）、`plot`（绘图）、`spos`（SPP）
- 绘图使用 `Agg` 后端，无需图形界面，可在服务器/批处理环境运行

## 6. 备注

- 断点续跑：station/Sstation 模式下，若文件子目录中已有 `<name>_QCValue_station.csv`，重跑时直接读取该结果计入对比表并跳过计算（prn 模式检查 `<name>_QCValue.csv`）；time 模式检查 `<name>_QCValue_station.csv` + `<name>_sppband.csv`，两者都在才跳过；只有前者时仅重算 SPP 选段（只读观测头，不重跑 QC 模块）。坏文件只警告跳过，不中断整批。如需强制重算，删除对应子目录或汇总文件即可。
- 综合评价需要至少 2 个评估对象（测站/PRN），否则自动跳过。
- 文件匹配导航星历依赖文件名中的年积日编码（如 `wuh20050.24o` ↔ `wuh20050.24p`）。
