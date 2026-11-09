#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
pppConfig.py — PPP 配置管理：默认参数、配置文件加载（TOML/JSON）、命令行参数合并。

优先级：命令行参数 > TOML 配置文件 > DEFAULT_CONFIG。
"""
import json
import os

# 与 ppp/ppp.toml 保持一致的默认参数
DEFAULT_CONFIG = {
    # ---- 系统/波段 ----
    'systems': ['G', 'C', 'E'],
    'bands': 'G:C1C,C2W;C:C2I,C6I;E:C1X,C5X',
    # ---- 定位模式与高度角 ----
    'elev_cutoff_deg': 7.0,
    'mode': 'static',            # static | kinematic
    # ---- 观测噪声（米）----
    'sigma_p': 3.0,              # 无电离层组合伪距
    'sigma_l': 0.03,             # 无电离层组合载波
    'sigma_elev_a': 0.3,         # 高度角定权 σ = σ0 * (a + b/sin(el))
    'sigma_elev_b': 0.7,
    # ---- 天顶湿延迟（随机游走）----
    'est_zwd': True,
    'sigma_zwd_init': 0.1,       # 初值方差 (m)^2
    'q_zwd': 1e-8,               # 过程噪声方差率 (m^2/s)
    # ---- 对流层梯度（MacMillan 模型，需 est_zwd）----
    'est_trop_grad': False,      # 估计北/东梯度 Gn/Ge（RTKLIB TROPOPT_ESTG 同源）
    'sigma_grad_init': 0.01,     # 梯度初值方差 (状态量)^2（RTKLIB VAR_GRA）
    'q_grad': 1e-8,              # 梯度过程噪声方差率（RTKLIB SQR(prn_trop*0.1)）
    # ---- 滤波初值与过程噪声 ----
    'sigma_pos_init': 100.0,     # 坐标初值方差 (m)^2
    'sigma_clk': 100.0,          # 接收机钟差：白噪声标准差 (m)，方差=平方（RTKLIB VAR_CLK=1e4 同源）
    'sigma_amb': 100.0,          # 模糊度初值方差（距离域，(m)^2）
    'q_pos_kinematic': 100.0,    # 动态模式位置随机游走方差率 (m^2/s)
    # ---- 误差模型 ----
    'trop_prior': 'saas_nmf',    # saas_nmf: Saastamoinen+NMF(RTKLIB 同源) | collins: 旧经验模型
    'windup': True,
    'solid_tide': True,
    'sat_pco': True,             # 卫星天线 PCO（ANTEX，含在改正 corr 中）
    'rec_pco': True,             # 接收机天线 PCO/PCV（ANTEX）
    'rec_ant_delta': False,      # RINEX 头 DELTA H/E/N：仅解算基准为标石坐标时开启
                                 # （IGS SINX 参考坐标指 ARP，默认关闭避免 U 偏高 antdel_U）
    'pole_tide': True,           # 极潮改正（IERS 2010 eq.7.26，RTKLIB tide_pole 同源）
    'windup_flip': False,        # 消融开关：缠绕改正取反（诊断符号用）
    # 载波码-相一致性平移阈值（米）：连续弧段载波新息中位数超过该值时，
    # 整体平移模糊度以恢复码-相相干（接收机钟 Stepping 码跟随、相位连续所致）
    'phase_jump_thresh_m': 10.0,
    # 载波拒绝门限下限（米）：模糊度方差收敛后 3σ 门限收缩到厘米级，
    # 而接收机钟漂（码漂移数百米/历元、载波不含）使载波新息达数百米，
    # 新息由白噪声钟差吸收（RTKLIB 载波无拒绝门同源行为）；门限下限
    # 防止载波被永久拒绝导致模糊度方差不收缩的死锁
    'carrier_reject_thresh_m': 1000.0,
    # ---- 精密产品（可选）----
    'atx_file': None,            # ANTEX 相位中心改正文件（igs20.atx 等）
    'blq_file': None,            # BLQ 海潮载荷文件（FES2004 等）
    'dcb_file': None,            # CODE P1-C1 DCB 文件（GPS C/A 码基准归算）
    'erp_file': None,            # IGS ERP（UT1-UTC 修正 GMST 等）
    'rec_antenna': None,         # 接收机天线类型（缺省取 RINEX 头 ANT # / TYPE）
    # ---- 精密星历内插 ----
    'sp3_pos_order': 10,
    'sp3_clk_interp': 'linear',  # linear | lagrange
    'sp3_apply_rel': True,       # 补周期相对论钟差（IGS 精密钟差不含该项；若钟差产品已含则关闭）
    'use_sat_var': True,         # 量测 R 叠加 SP3 卫星轨道/钟差方差（RTKLIB var_rs 同源）
    # ---- 相位-码跳变改正 ----
    'phase_code_jump': True,     # 毫秒级钟跳检测：模糊度整体平移保持相位-码相干（RTKLIB udbias 同源）
    # ---- 周跳探测 ----
    'slip_mw_sigma': 4.0,
    'slip_gf_thresh': 0.10,      # (m)
    'gap_reset_sec': 120.0,
    # ---- 质量控制 ----
    'outlier_nsigma': 3.0,
    'min_sat': 4,
    'exclude_sats': [],          # 排除的卫星（如 BDS GEO: ['C01','C02','C03','C04','C05']）
    'ref_mean_skip_min': 0,      # ENU 参考均值跳过前 N 分钟收敛段
    'ppp_conv_thresh': 0.3,      # 绘图收敛判据：3D 误差阈值（m）
    'ppp_conv_hold': 5,          # 收敛判据需持续保持的时长（分钟）
    # ---- 输出 ----
    'plot': True,
}


def _stripJsonComments(text):
    """
    去除 JSON 文本中的 // 行注释（字符串内部的 // 不受影响）。
    使 ppp.json 支持注释（JSONC 风格）。
    """
    out = []
    inStr = False
    i, n = 0, len(text)
    while i < n:
        c = text[i]
        if inStr:
            out.append(c)
            if c == '\\':          # 转义字符，跳过下一个
                i += 2
                continue
            if c == '"':
                inStr = False
            i += 1
        else:
            if c == '"':
                inStr = True
                out.append(c)
                i += 1
            elif c == '/' and i + 1 < n and text[i + 1] == '/':
                while i < n and text[i] != '\n':   # 跳到行尾
                    i += 1
            else:
                out.append(c)
                i += 1
    return ''.join(out)


def _parseSimpleToml(text):
    """
    简易 TOML 解析（回退用，覆盖本配置文件的子集）：
        - # 注释（整行或行尾，字符串内 # 不受影响）
        - 扁平 key = value（无 [table]）
        - 值类型：字符串（"..." 或 '...'）、整数、浮点、true/false、字符串数组 ["A", "B"]
    仅在 tomllib/tomli 均不可用时使用。
    """
    def stripComment(s):
        out, inStr, quote = [], False, None
        for c in s:
            if inStr:
                out.append(c)
                if c == quote:
                    inStr = False
            elif c in ('"', "'"):
                inStr, quote = True, c
                out.append(c)
            elif c == '#':
                break
            else:
                out.append(c)
        return ''.join(out)

    def parseVal(s):
        s = s.strip()
        if not s:
            return None
        if s.startswith('['):
            body = s[1:s.rindex(']')] if ']' in s else s[1:]
            items = []
            cur, inStr, quote = [], False, None
            for c in body:
                if inStr:
                    if c == quote:
                        inStr = False      # 闭合引号不并入内容
                    else:
                        cur.append(c)
                elif c in ('"', "'"):
                    inStr, quote = True, c
                elif c == ',':
                    items.append(''.join(cur).strip())
                    cur = []
                else:
                    cur.append(c)
            tail = ''.join(cur).strip()
            if tail:
                items.append(tail)
            return [parseVal(it) for it in items if it != '']
        if (s.startswith('"') and s.endswith('"')) or \
           (s.startswith("'") and s.endswith("'")):
            return s[1:-1]
        if s == 'true':
            return True
        if s == 'false':
            return False
        try:
            return int(s)
        except ValueError:
            pass
        try:
            return float(s)
        except ValueError:
            pass
        return s

    out = {}
    for raw in text.splitlines():
        line = stripComment(raw).strip()
        if not line or line.startswith('['):
            continue
        eq = line.find('=')
        if eq <= 0:
            continue
        key = line[:eq].strip().strip('"\'')
        if key:
            out[key] = parseVal(line[eq + 1:])
    return out


def _loadToml(configPath):
    """加载 TOML：优先 tomllib（Py>=3.11），其次 tomli，最后内置简易解析器。"""
    try:
        import tomllib
    except ImportError:
        try:
            import tomli as tomllib
        except ImportError:
            with open(configPath, 'r', encoding='utf-8') as f:
                return _parseSimpleToml(f.read())
    with open(configPath, 'rb') as f:
        return tomllib.load(f)


def load_config(configPath=None):
    """
    加载配置文件并与默认配置合并（文件值覆盖默认值）。
    按扩展名选择解析器：
        .toml — TOML 格式（原生支持 # 注释，推荐，见 ppp/ppp.toml；
                tomllib(Py>=3.11) / tomli / 内置简易解析器三级回退）
        .json — JSON 格式（兼容 // 注释）
    configPath 为 None 时返回默认配置副本。
    """
    cfg = dict(DEFAULT_CONFIG)
    if configPath is None:
        return cfg
    if not os.path.isfile(configPath):
        raise FileNotFoundError('PPP 配置文件不存在: %s' % configPath)
    if configPath.lower().endswith('.toml'):
        userCfg = _loadToml(configPath)
    else:
        with open(configPath, 'r', encoding='utf-8') as f:
            userCfg = json.loads(_stripJsonComments(f.read()))
    for k, v in userCfg.items():
        if k in cfg:
            cfg[k] = v
        else:
            print('警告：配置文件中出现未知键 "%s"，已忽略' % k)
    # 空字符串归一化为 None（TOML 无 null，如 rec_antenna）
    if cfg.get('rec_antenna') == '':
        cfg['rec_antenna'] = None
    return cfg


def resolve_band_choose(cfg):
    """
    根据配置中的 systems/bands 生成 bandChoose 字典：
        {'G': {'C1C': {}, 'C2W': {}}, ...}
    bands 为空字符串时使用各系统默认波段。
    """
    from spp import parse_band_string   # 复用 spp.py 的解析逻辑

    systems = cfg.get('systems') or ['G', 'C', 'E']
    bands = cfg.get('bands') or ''
    if bands:
        bandChoose = parse_band_string(bands, sys_list=set(systems))
        if bandChoose:
            return bandChoose
    # 默认波段
    defaultBands = {
        'G': ['C1C', 'C2W'],
        'C': ['C2I', 'C6I'],
        'E': ['C1X', 'C5X'],
    }
    bandChoose = {}
    for s in systems:
        if s in defaultBands:
            bandChoose[s] = {code: {} for code in defaultBands[s]}
    return bandChoose
