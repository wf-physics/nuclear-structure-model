# 逐脚本复现说明

所有脚本在 `src/` 下运行（`run_all.sh` 已处理 cwd 与 results 软链）。数据依赖 `data/` 四个官方文件，无需下载。

## core 层（分钟级）

| 顺序 | 脚本 | 作用 | 产出 | 预期 |
|---|---|---|---|---|
| 1 | `nm_formula_v4.py` | 金刚石网络核物质标定（解析直解 280 组→120 有效） | `results/nm_formula_v4.json` | ρ0=0.160, E/A=−15.8, K=547 |
| 2 | `q_formula_v4.py` | v4 闭式（模型网格 99 点最小二乘）+ F1–F4 首测 | `q_formula_v4.json`, `model_drop_coef.*` | F1 MAE 0.83 (10/10) |
| 3 | `q_formula_v5.py` | 重标闭式（团簇求和×2.61）重测 F1–F4 | `q_formula_v5.json`, `q_formula_v5_meta.json` | **F4 MAE 3.83 (5/5)**, F3 退化(保留 v4) |
| 4 | `q_formula_v6_f2patch.py` | F2 收口：Lane/Lane+配对矩阵（δ 独立拟合） | `q_formula_v6_f2_patch.json` | Lane 采纳 3.80；δ<0 否决 |
| 5 | `q_formula_v6_pairing_v2.py` | 配对项零借用提取探测（198 点扩充网格） | `pairing_probe_v2.json` | 未检测到（噪声/信号≈50×） |
| 6 | `q_engine_pairing.py` | 引擎级同 A 差分预研 | `engine_pairing_probe_v2.json` | 不可提取，缺额 ~4 MeV 定量证据 |
| 7 | `f1_expanded.py` | F1 扩样 85 核 | `f1_expanded_85.json` | A≥40 零不过, MAE 1.58 |
| 8 | `valve_test.py` | 阀门验证（S_n 换闭式接衰变层, 30 seeds） | `valve_test.json` | **33.9%/66.7% <90% → 阀门关闭** |
| 9 | `finite_nuclei_v1.py` | 有限核 23 核液滴近似 | `finite_nuclei_v1.json` | 中核 A40–65 最佳工作区 |

## full 层（小时级）

| 顺序 | 脚本 | 作用 | 产出 | 预期 |
|---|---|---|---|---|
| 10 | `bond_calibration.py` | 几何引擎键标定 | C_REP=12, R_EQ=1.608 | — |
| 11 | `relax_all.py` | A≤40 全量结构弛豫（最慢，约数小时） | `results/structures.json`（已随附存档） | — |
| 12 | `chain_predict.py` | 反应链推演 + 判定准确率（依赖 11 的存档） | `chain_results.json` | 模型 92.6%→见存档 |
| 13 | `big_validation.py` | 几何大验证 H→Ca 全量（依赖 11） | `big_validation_table.csv` | B/r_ch/Q_s/α 四项 |
| 14 | `geometry_validation.py` | Q₀ 解析验证 + α 探测器 + L-BFGS 引擎 | `geometry_validation.json` | Q₀ 解析解精确 |
| 15 | `blind_test.py` | 随机盲测 P1–P4 协议（30 seeds） | `blind_results.json` | P1 mean 97.44% / free 99.2% |
| 16 | `criterion_test.py` | 判据单元测试 | — | — |

## 依赖关系

- `big_validation.py` / `blind_test.py` 内部 `exec` 前缀自 `relax_v2_test.py`（势与力函数）；
- `chain_predict.py` 需 `results/structures.json`（脚本 11 产物，**已随附存档**——不重跑 11 也能复现 12–15）；
- `blind_test.py` 需 `results/baseline_judge.json`（已随附）。

## 环境说明

- Python ≥3.10, numpy, scipy；无 GPU 依赖；
- 交叉平台注意：Metropolis/L-BFGS 有固定随机种子（seed=13/42/43/44 等），结果可精确复现；
- AME2020/NUBASE2020/IAEA 数据为官方原文件，未做任何修改。

## 三层分工框架预印本（英文版 v1.1）新增映射

| 顺序 | 脚本 | 作用 | 产出 | 预期 |
|---|---|---|---|---|
| 17 | `holdout_refit.py` | §3.3 六参数全局拟合复现 + §3.4 留出验证（规则冻结：Z≥88 全剔） | `results/holdout_Zge88_result.json` | 样本内 MAE 2.296（论文 2.297）；留出后 U-238: B +5.44 / S_n +0.47 / Q_α +0.95 / Q_fs(Pd+Pd) −4.57 MeV |
| 18 | `radius_scaling.py` | §2.4 电荷半径标度因子诊断 | `results/radius_scaling_result.json` | 中位比 0.815（存档 95 核；论文口径 0.794/97 核）；标度后 ±15% 带内 30.5%→52.6% |
| 23 | `frdm_same_set.py` | FRDM2012 同集对比（Paper I v2.2 §4.4 脚注；亦用于预印本 §3.3 表） | `results/frdm_same_set.json` | A=40–240∩FRDM n=2910 MAE 0.449/RMS 0.631；A≥16 全集 MAE 0.603/RMS 0.988 |

## 预言卡（2026-11-24 预注册）

| 顺序 | 脚本 | 作用 | 产出 | 预期 |
|---|---|---|---|---|
| 19 | `prediction_ca34.py` | ³⁴Ca 预言卡一条龙（仅实验拟合→局部误差棒→预言数字） | `results/prediction_ca34_numbers.json` | ⁹C+²⁵Si→³⁴Ca Q=+18.03±1.55；S_2p 衰变能 +2.33±1.55 MeV |
| 20 | `results/fusion_channels_bound.json` | 束缚核过滤后的全部 A+B→C 聚变道枚举（52 个 C 核素） | 同左 | 头牌选择依据 |
| 21 | `results/ca34_structure_seed13.json` | ³⁴Ca 几何引擎现场弛豫汇总（seed=13） | 同左 | α 簇 3/3（理想上限），Q0=+170.33 fm² |
| 22 | `results/nubase_bound.json` | NUBASE2020 基态衰变模式分类（束缚/粒子发射/未知） | 同左 | 束缚核过滤器 |
| 24 | `results/alpha_cluster_spearman.json` | §5.3 α 簇团识别比例对 A 的 Spearman 秩相关（由 `results/big_validation_rows.json` 实算） | 同左 | r=−0.603、p=0.065（n=10） |
| 25 | `results/radius_scaling_result.json` | §5.2 电荷半径标度因子诊断（由 `results/big_validation_rows.json` 实算） | 同左 | ±15% 带内 30.5%→52.6% |

## 九张图重制修订（2026-11-24，配合排版审图）

`src/make_figures.py` 重写后全部 9 张图重出，图值与存档/正文逐位对应：
图1 修复面板索引错位（ρ0/E/A 各归其位、FCC 解析直解画"无有效解"注记）；
图2 图例区间统一 100–240、数据源移入 `data/mass_dev_2811_recalibrated.csv`（消除外部硬编码路径）；
图3 数据源由 `f1_expanded_42.json` 改为正式 85 核档 `f1_expanded_85.json`（过档全部落回 ±2 带内）并标注 ¹⁸F；
图4 P4 柱明确为盲测侧 86.6%、30.1% 训练侧画为置换对照线（消除柱/注口径混用）；
图5 (a) 标题更正为引擎层探针（原误标 closed-form），RMS 121.7/δ_model=+1.46 与存档一致，闭式层 δ=−3.03/RMS 14.35 只入正文（`q_formula_v6_f2patch.py` 复跑验证）；
图6 五道裂变道标签改为"母核→碎片道"逐一区分；
图7 解析 v3 嵌套存档补上失效轮，F1 曲线标注 n=10 选型集并加 85 核正式成绩 1.58 星标；
图8 分区标签统一 Z 升序（Z<=8 / Z 9–14 / Z 15–20），柱内标注各观测量分母 n。
