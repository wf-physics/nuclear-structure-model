# Bond-Network Nuclear Structure Model (v5.0 frozen)

四键饱和键网络核结构模型（v5.0 定档冻结版）的完整代码、数据与结果存档。

**理论核心**：核子间唯一相互作用为线性键阱 `V(r) = −5·(2.2−r) fm⁻¹·MeV`（窗口 [0.8, 2.2] fm）+ 泡利硬核；每核子最多 4 键（k_max=4，互相同意制）。该规则在无限核物质极限下收敛到**理想四面体四配位网络的简并族**（立方金刚石为代表构型，非唯一基态；族内能量差 ≤0.03 MeV/核子）——其四面体配位正是 α 簇团结构，与独立几何验证（电四极矩/α 簇团/电荷半径）同源自洽。

## 主要结果（全部见 `results/`，逐核素清单随附）

| 模块 | 结果 | 基线对照 | 存档 |
|---|---|---|---|
| 衰变判定层（A≤40 全量） | **99.0%（311/314）** | 纯 Q 符号基线 92.6%；零结构管线 97.5% | `results/chain_results.json` |
| 衰变判定层盲测（30 seeds, P1 协议） | mean 97.4±1.1% / 剔补丁核素 99.2% | 置换对照跌至 30% 以下 | `results/blind_results.json` |
| 核物质饱和点 | ρ0=0.160 fm⁻³ ✓, E/A=−15.8 MeV ✓ | 经验 0.16 / −15.8 | `nm_formula_v4.json` |
| K（不可压缩） | 547 MeV（不达标，不预言 GMR——limitation） | 经验 240±30 | 同上 |
| 金刚石团簇求和 B(A,Z) | 2811 核素 MAE 12.6 MeV（零逐核拟合） | 拟合液滴 RMSE 2–4 | `q_formula_v5_meta.json` |
| F1 S_n 闭式（A≥40, 85 核扩样） | MAE 1.58，零"不过" | — | `f1_expanded_85.json` |
| F3 Q_α（A≥180） | MAE 1.32 MeV | — | `q_formula_v4.json` |
| F4 Q_fs 裂变释放能 | MAE 3.83 MeV，5/5 过（判定带 20%·Q≈35 MeV） | — | `q_formula_v5.json` |
| F2 Q_β⁻ | MAE 3.80（Lane 基座），−3~−9 系统负偏（壳效应，limitation） | — | `q_formula_v6_f2_patch.json` |
| 电四极矩/α 簇团/电荷半径（几何大验证 H→Ca） | 见汇总报告 | — | `big_validation_table.csv` |
| 阀门验证（S_n 换闭式接衰变层） | **33.9%/66.7% < 90% → 不宣称"独立预言"** | 判据预注册 | `valve_test.json` |

## 立场声明（防审稿雷区，主动披露）

1. **能量层系数"是长出来的，不是拟合的"**：液滴系数 a_v=16.06 / a_s=19.01 / a_s/a_v=1.18（经验 1.14）来自冻结线性井的金刚石团簇求和 × 井深重标，零逐核拟合；v4 闭式四系数全部由模型网格能量最小二乘导出，零实验借用。唯一显式借用：F3 的 B_α 锚点 28.296 MeV（⁴He 实验结合能，文中标注）。
2. **99% 衰变判定中模型的净贡献 = 1.5 个百分点**（97.5% 零结构管线 → 99.0%），能量 Q 值来自 AME2020 质量表——不宣称能量独立（阀门验证实测否决，见 `valve_test.json`）。
3. **自由参数清单**（预注册口径，逐层声明）：

| 层 | 参数 | 数值/来源 |
|---|---|---|
| v5.0 冻结法则 | 井深 5、窗口 [0.8,2.2]、k_max=4、泡利硬核 | 理论定档，不可调 |
| 几何引擎 | C_REP=12, R_EQ=1.608 fm | 键标定脚本（`bond_calibration.py`） |
| 核物质标定 | k=56.97, C=699.1, δ=0.2（+alpha=3, p=6） | 280 组解析解 → 120 有效解标定（`nm_formula_v4.py`） |
| v4 闭式 | a_v=13.851, a_s=1.592, a_c=0.7265, a_sym=19.50 | 模型网格 99 点最小二乘，零实验借用 |
| 重标闭式 | a_v=16.060, a_s=19.010, a_c=0.756, a_sym=19.50 | 团簇求和×井深重标（2.61），零逐核拟合 |
| 衰变判定层 | sig_k, NU_TRY, S_ALPHA, EC_F, FT_N, T_STABLE（6 个） | 盲测协议：仅训练集拟合，冻结后预测盲测集 |
| F3 锚点 | B_α=28.296 MeV | 实验值借用，唯一例外，文中标注 |

4. **Limitations（八条，全部实测落定，与论文 Paper I §11 一致）**：
   ① K=547 不预言 GMR 能量（窗口冻结解析墙）；② A≤40 能量侧闭式失效需显式团簇（阀门验证实测）；③ F2 有 −3~−9 MeV 系统负偏（残差被壳修正差完全吸收，闭式无壳结构）；④ F3 Q_α 仅 A≥180 区适用；⑤ 滴线区（需 0.1 MeV）精度不足；⑥ 几何形状观测量对退火种子敏感——能量侧散布 mean 3.4% 稳健，但 A≥20 的 Q₀/α 计数随种子散布（std 最高 1154 / 计数 1–5），A≤14 完全稳定；⑦ 12 个滴线补丁核素的判定正确性不来自模型（名单随文公开）；⑧ 电荷半径中位比值 0.794（系统性偏小 ~21%，缺零点振动修正，结构性缺陷）。

## 目录结构

```
repo/
├── src/                  # 全部脚本 (25 个, 按复现顺序见 run_all.sh)
├── data/                 # 官方实验数据 (未改动)
│   ├── mass_1.mas20      # AME2020 原子质量表
│   ├── nubase_4.mas20    # NUBASE2020
│   ├── charge_radii_iaea.csv   # IAEA 电荷半径 (Angeli/Marinova 2013)
│   └── elec_mom_iaea.csv       # IAEA 电四极矩推荐值
├── results/              # 关键结果存档 (18 个)
├── theory/               # v5.0 法则文档 + 能量公式体系审计链 + 验证汇总
├── run_all.sh            # 一键复现 (分两层: 核心分钟级 / 全量小时级)
└── REPRODUCE.md          # 逐脚本复现说明
```

## 快速复现

```bash
pip install -r requirements.txt
bash run_all.sh core      # 核心结果 (~10 分钟)
bash run_all.sh full      # 全量几何/盲测 (~数小时, 结果已存档于 results/)
```

## 数据来源

- AME2020 / NUBASE2020: IAEA (https://www-nds.iaea.org/amdc/)
- 电荷半径: Angeli & Marinova 2013, IAEA 镜像
- 电四极矩: IAEA 电矩推荐值 (Stone 2016)

## License

MIT（见 LICENSE）

## 论文数字 → 脚本 → 输出文件映射表（可复现性命门）

| 论文数字 | 生成脚本 | 输出存档 |
|---|---|---|
| 核物质 ρ0=0.160 / E/A=−15.8 / K=547（§3） | `src/nm_formula_v4.py` | `results/nm_formula_v4.json` |
| 金刚石四键网络 120 组有效解（§3 表） | `src/nm_formula_v4.py` | 同上 |
| 2811 核素 B(A,Z) MAE 12.64 / 分区 6.48/9.30/14.33（§4） | DeepSeek 重标闭式对照（`src/q_formula_v5.py` 内嵌系数） | `results/q_formula_v5_meta.json` |
| 路线二四系数 a_v=13.851/a_s=1.592/a_c=0.7265/a_sym=19.50（§4） | `src/q_formula_v4.py` | `results/model_drop_coef.json` |
| F1 S_n 85 核 MAE 1.58（59 过/25 半/1 不过）（Paper B §3） | `src/f1_expanded.py` | `results/f1_expanded_42.json` |
| F2 Q_β MAE 3.80（Lane 采纳）（Paper B §4） | `src/q_formula_v6_f2patch.py` | `results/q_formula_v6_f2_patch.json` |
| F3 Q_α MAE 1.32 / F4 Q_fs 5/5 MAE 3.83（Paper B §5–6） | `src/q_formula_v4.py` / `src/q_formula_v5.py` | `results/q_formula_v4.json` / `results/q_formula_v5.json` |
| T4 四口径矩阵 33.93/57.19/66.74/66.74（Paper B §8） | `src/t4_err_sensitivity.py` + `src/t4_cal03_gap.py` | `results/t4_err_sensitivity.json` |
| T5 壳修正回归 R²=0.978 / 残差 0.41 MeV（Paper B §4） | `src/t5_shell_correlation.py` + `src/t5_expanded.py` | `results/t5_shell_correlation_v3.json` |
| 衰变判定 99.0%（311/314）与盲测 P1 97.4（Paper C） | `src/blind_test.py` | `results/blind_results.json` |
| 反应链 99.0%（H→Ca 389 核素） | `src/chain_predict.py` | `results/chain_results.json` |
| 几何四项（B 69.2% / Q_s 符号 9/13 / r_ch 66% / α 8/10）（§5） | `src/big_validation.py` + `src/geometry_validation.py` | `results/big_validation_table.csv` / `results/geometry_validation.json` |
| 种子散布（E 3.4% / Q0 std 1154 / α 计数 1–5）（§2.4、§7） | `src/seed_dispersion.py` | `results/seed_dispersion.json` |
| 阀门验证 33.9/66.7（A≤40 能量侧负结果）（§7 / Paper B §8） | `src/valve_test.py` | `results/valve_test.json` |
| 配对两路证伪（噪声/信号 ~50×、引擎奇偶缺额 −3.9）（Paper B §4） | `src/q_formula_v6_pairing_v2.py` / `src/q_engine_pairing.py` | `results/pairing_probe_v2.json` / `results/engine_pairing_probe_v2.json` |
| 附录 A 汇总表（315 核素 E_model/E_bond/k 分布） | `src/appendixB_export.py` | `results/appendixA_summary.csv` |
| 附录 B 逐核子表（⁴He/¹⁶O/⁴⁸Ca 正文；全量 203,662 键行补充材料） | `src/appendixB_export.py` | `results/appendixB_{He-4,O-16,Ca-48}.json` / `results/appendixB_all_nucleons.csv` |
| 全部 9 张图（A-图1…B-图4） | `src/make_figures.py` | `figures/*.png` |

**口径说明（两引擎差异）**：判定层/结构生成引擎（`structures.json`，315 核素）与几何验证引擎（`geometry_validation.py`，C=12, r_eq=1.608）为两个独立标定的实现；附录 B 正文三表统一用几何引擎现场弛豫（seed=13）保证表内自洽，全量长表沿用结构存档坐标。表注已在论文中声明。

## 数据再分发说明

- `mass_1.mas20` / `nubase_4.mas20`：IAEA/AMDC 供科研自由使用，转载保留原始表头与出处声明
- IAEA 电荷半径/电四极矩 CSV：来自 IAEA 结构图数据服务，引用 Angeli & Marinova 2013 与 Stone 2016
