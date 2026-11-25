# -*- coding: utf-8 -*-
"""Paper I 全部 9 张图。数据源全部来自 repo/results 与 repo/data 存档，
图值与论文图注/正文表格逐位对应（修订版 v2：修复 A_fig1 面板索引错位、
A_fig3 换 85 核存档、P4 双口径、B_fig1 补 v3、B_fig3 标题改用本实验存档数、
B_fig2 裂变道标签逐一区分、A_fig5 分区标签统一 Z 升序并标 n）。"""
import os, json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams["font.family"] = ["DejaVu Sans"]
plt.rcParams["axes.unicode_minus"] = False

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, "..", "results")
DAT = os.path.join(HERE, "..", "data")
FIG = os.path.join(HERE, "..", "figures")
os.makedirs(FIG, exist_ok=True)

def save(fig, name):
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, name), dpi=200)
    plt.close(fig)
    print("  OK", name)

# ============ 图1 (A_fig1): 参数化×晶格扫描病态对照 ============
# 数据源 = 论文 §3.2 表（金刚石行存档于 nm_formula_v4.json, FCC 行为 §3.2 网格扫描记录）
dia = json.load(open(os.path.join(RES, "nm_formula_v4.json")))
assert abs(dia["rho0"] - 0.16) < 1e-9 and abs(dia["eA"] + 15.8) < 1e-9
data = {  # (rho0 fm^-3, E/A MeV, K MeV)
    "FCC + F(rho)\n(7p DE)": (0.30, -15.0, 846),
    "FCC + Pauli core\n(5p DE)": (0.30, -11.4, 1446),
    "FCC analytic\n(280-grid)": (None, None, 6.7e4),
    "Diamond 4-bond\nnetwork": (dia["rho0"], dia["eA"], dia["K"]),
}
fig, ax = plt.subplots(1, 2, figsize=(11, 4.4))
xs = np.arange(len(data))
cols = ["#c44", "#c66", "#999", "#267"]
names = list(data.keys())
# (a) 饱和密度
vals = [data[k][0] for k in names]
ok = [i for i, v in enumerate(vals) if v is not None]
ax[0].bar([xs[i] for i in ok], [vals[i] for i in ok], color=[cols[i] for i in ok])
ax[0].axhline(0.16, color="k", ls="--", lw=1, label="empirical 0.16")
ax[0].set_ylim(0, 0.36)
ax[0].annotate("+88%", (0, 0.305), ha="center", fontsize=9)
ax[0].annotate("+88%", (1, 0.305), ha="center", fontsize=9)
ax[0].text(2, 0.05, "no valid solution\n(all 280 grid:\nK > 350)", ha="center", va="center",
           fontsize=8, color="#c44")
ax[0].set_xticks(xs); ax[0].set_xticklabels(names, fontsize=8)
ax[0].set_ylabel(r"$\rho_0$ (fm$^{-3}$)")
ax[0].legend(fontsize=8); ax[0].set_title("(a) saturation density")
# (b) 每核子结合能 (MeV)
vals = [data[k][1] for k in names]
ok = [i for i, v in enumerate(vals) if v is not None]
ax[1].bar([xs[i] for i in ok], [vals[i] for i in ok], color=[cols[i] for i in ok])
ax[1].axhline(-15.8, color="k", ls="--", lw=1, label="empirical -15.8")
ax[1].set_ylim(-17.5, 0.5)
ax[1].annotate("-5%", (0, -14.3), ha="center", fontsize=9)
ax[1].annotate("-28%", (1, -10.7), ha="center", fontsize=9)
ax[1].text(2, -8.5, "no valid solution", ha="center", fontsize=8, color="#c44")
ax[1].set_xticks(xs); ax[1].set_xticklabels(names, fontsize=8)
ax[1].set_ylabel("E/A (MeV)")
ax[1].legend(fontsize=8); ax[1].set_title("(b) energy per nucleon")
save(fig, "A_fig1_lattice_scan.png")

# ============ 图2 (A_fig2): 2811 核素 B 残差散点 ============
df = pd.read_csv(os.path.join(DAT, "mass_dev_2811_recalibrated.csv"))
fig, ax = plt.subplots(figsize=(7.5, 5))
for lo, hi, c in [(40, 65, "#26a"), (65, 100, "#2a7"), (100, 240, "#a62")]:
    m = (df["A"] >= lo) & (df["A"] < hi)
    ax.scatter(df["A"][m], df["偏差(MeV)"][m], s=3, alpha=0.5, c=c,
               label=f"A {lo}-{hi}: MAE {abs(df['偏差(MeV)'][m]).mean():.1f}")
ax.axhline(0, color="k", lw=0.8)
ax.set_xlabel("A"); ax.set_ylabel("B_model - B_exp (MeV)")
ax.legend(fontsize=9)
ax.set_title("2811 nuclei, zero per-nucleus fitting (MAE 12.6 MeV)")
save(fig, "A_fig2_mass_table.png")

# ============ 图3 (A_fig3): 85 核 F1 偏差分布 ============
rows = json.load(open(os.path.join(RES, "f1_expanded_85.json")))
assert len(rows) == 85
fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
VN = {"过": "pass", "半": "half", "不过": "fail"}
vc = {"过": "#267", "半": "#fa0", "不过": "#c44"}
for v in ["过", "半", "不过"]:
    m = [r for r in rows if r["v"] == v]
    ax[0].scatter([r["Z"] + r["N"] for r in m], [r["dev"] for r in m], c=vc[v], s=25,
                  label=f"{VN[v]} (n={len(m)})", alpha=0.8)
f18 = [r for r in rows if r["v"] == "不过"][0]
ax[0].annotate(f"{f18['nucl']} ({f18['dev']:+.1f})",
               (f18["Z"] + f18["N"], f18["dev"]), fontsize=8,
               xytext=(8, -2), textcoords="offset points", color="#c44")
ax[0].axhspan(-2, 2, color="#267", alpha=0.08)
ax[0].axhspan(-5, 5, color="#fa0", alpha=0.06)
ax[0].axhline(0, color="k", lw=0.8)
ax[0].set_xlabel("A"); ax[0].set_ylabel("S_n dev (MeV)")
ax[0].legend(fontsize=9); ax[0].set_title("(a) F1 residuals vs A (n=85)")
mae = float(np.mean([abs(r["dev"]) for r in rows]))
ax[1].hist([abs(r["dev"]) for r in rows], bins=20, color="#26a", alpha=0.8)
ax[1].axvline(mae, color="k", ls="--", label=f"MAE {mae:.2f}")
ax[1].set_xlabel("|S_n dev| (MeV)"); ax[1].legend(fontsize=9)
ax[1].set_title("(b) |dev| distribution")
save(fig, "A_fig3_f1_expanded.png")

# ============ 图4 (A_fig4): 盲测 P1-P4 + T4 四口径矩阵 ============
bl = json.load(open(os.path.join(RES, "blind_results.json")))
t4 = json.load(open(os.path.join(RES, "t4_err_sensitivity.json")))
fig, ax = plt.subplots(figsize=(10, 4.5))
labels, vals, errs, cols4 = [], [], [], []
for p in ["P1", "P2", "P3", "P4"]:
    r = bl[p]
    mv = r.get("mean", r.get("mean_mid", r.get("mean_blind")))
    sd = r.get("std", None)
    labels.append(p); vals.append(mv * 100)
    errs.append(sd * 100 if sd else 0)
    cols4.append("#1f77b4")
p4_train = bl["P4"]["mean_train_shuffled"] * 100   # 30.1% 训练侧(标签置换)
p4_blind = bl["P4"]["mean_blind"] * 100            # 86.6% 盲测侧(柱值)
mc = t4.get("matrix_complete", {})
quad = [("raw @gate=AMEerr", mc.get("raw_ame")), ("raw @gate=0.3", mc.get("raw_03")),
        ("cal @gate=AMEerr", mc.get("cal_ame")), ("cal @gate=0.3", mc.get("cal_03"))]
for lab, mv in quad:
    if mv is None:
        continue
    labels.append(lab); vals.append(mv * 100); errs.append(0)
    cols4.append({"raw @gate=AMEerr": "#1f77b4", "raw @gate=0.3": "#7fb3d8",
                  "cal @gate=AMEerr": "#e1812c", "cal @gate=0.3": "#f0b27a"}[lab])
xs = np.arange(len(labels))
ax.bar(xs, vals, yerr=errs, color=cols4, capsize=3)
ax.axhline(p4_train, color="#888", ls=":", lw=1.2)
ax.text(6.5, p4_train - 5, f"P4 train-side after\nlabel permutation {p4_train:.1f}%",
        fontsize=8, color="#555", ha="center")
ax.annotate("blind side\n86.6%", (3, 78), ha="center", fontsize=8,
            color="white", fontweight="bold")
ax.axhline(90, color="k", ls="--", lw=0.8)
ax.text(0.1, 91, "preregistered 90% line", fontsize=8)
ax.set_xticks(xs)
ax.set_xticklabels(labels, fontsize=7, rotation=20, ha="right")
ax.set_ylabel("accuracy (%)\nP error bar = 1 std over seeds")
ax.set_title("Blind protocols P1-P4 (bar = blind side)  +  T4 four-quadrant gate matrix (all < 90%)")
save(fig, "A_fig4_blind_and_T4.png")

# ============ 图5 (A_fig5): 几何大验证分区 (Z 分区, 升序, 标 n) ============
t = pd.read_csv(os.path.join(RES, "big_validation_table.csv"))
zmax = {"轻(Z≤8)": 8, "中(9-14)": 14, "重(15-20)": 20}
ren = {"轻(Z≤8)": "Z<=8", "中(9-14)": "Z 9-14", "重(15-20)": "Z 15-20"}
regions = sorted(t["region"].unique(), key=lambda r: zmax.get(r, 99))
def rlabel(r):
    n = (t["region"] == r).sum()
    return f"{ren.get(r, r)}\n(n={n})"
fig, axes = plt.subplots(1, 4, figsize=(14, 3.8))
for axi, (col, name) in zip(axes, [("B_v", "Binding energy"), ("Q_v", "Quadrupole sign"),
                                   ("A_v", "alpha cluster"), ("R_v", "charge radius")]):
    s = t[col].dropna()
    ok = (s == "过").mean() * 100
    sub = t.dropna(subset=[col])
    by_reg = sub.groupby("region")[col].apply(lambda x: (x == "过").mean() * 100)
    n_reg = sub.groupby("region")[col].size()
    order = [r for r in regions if r in by_reg.index]
    axi.bar(range(len(order)), [by_reg[r] for r in order], color="#267")
    for i, r in enumerate(order):
        nv = int(n_reg[r]); yv = by_reg[r]
        if yv >= 12:
            axi.annotate(f"n={nv}", (i, yv - 4), ha="center", fontsize=8,
                         color="white", fontweight="bold")
        else:
            axi.annotate(f"n={nv}", (i, yv + 3), ha="center", fontsize=8, color="#333")
    axi.axhline(ok, color="k", ls="--", lw=0.8)
    axi.set_xticks(range(len(order)))
    axi.set_xticklabels([ren.get(r, r) for r in order], fontsize=8)
    axi.set_ylim(0, 100); axi.set_title(f"{name}\n(overall {ok:.0f}%)", fontsize=9)
    axi.set_ylabel("pass rate (%)")
save(fig, "A_fig5_geometry_regions.png")

# ============ 图6 (B_fig1): 差分公式族五轮迭代 ============
vers = ["v1", "v2", "v3", "v4", "v5"]
keys = ["F1_Sn", "F2_Qbm", "F3_Qa", "F4_Qfs"]
def get_mae(f, key):
    d = json.load(open(os.path.join(RES, f)))
    rows = d.get(key)
    if rows is None:  # v3 为嵌套结构 {lit: {...}, net: {...}}, 取 lit 路(延续到 v4/v5 的路线)
        rows = d.get("lit", {}).get(key)
    if isinstance(rows, list) and rows and "dev" in rows[0]:
        return float(np.mean([abs(r["dev"]) for r in rows]))
    return None
fig, ax = plt.subplots(figsize=(8, 4.5))
key_off = {"F1_Sn": (0, 5), "F2_Qbm": (0, -14), "F3_Qa": (0, 5), "F4_Qfs": (0, 5)}
for key, c in zip(keys, ["#26a", "#fa0", "#267", "#c44"]):
    pts = [(v, get_mae(f"q_formula_{v}.json", key)) for v in vers]
    pts = [(v, y) for v, y in pts if y is not None]
    lab = key + (" (n=10 selection set)" if key == "F1_Sn" else "")
    ax.plot([p[0] for p in pts], [p[1] for p in pts], "o-", color=c, label=lab)
    for x, y in pts:
        ax.annotate(f"{y:.2g}", (x, y), fontsize=7, xytext=key_off[key],
                    textcoords="offset points")
# F1 正式成绩 (85 核, f1_expanded_85.json) 星标于 v5 侧
f1_85 = float(np.mean([abs(r["dev"]) for r in json.load(open(os.path.join(RES, "f1_expanded_85.json")))]))
ax.plot([4.35], [f1_85], "*", color="#26a", ms=14)
ax.annotate(f"F1 formal score {f1_85:.2f} (n=85)", (4.35, f1_85), fontsize=8,
            color="#26a", ha="center", va="top", xytext=(0, -14),
            textcoords="offset points")
ax.set_yscale("log")
ax.set_ylabel("MAE (MeV, log)")
ax.legend(fontsize=8)
ax.set_title("Differential formula family: 5-round iteration (v3 = failed round)")
save(fig, "B_fig1_evolution.png")

# ============ 图7 (B_fig2): F4 五裂变道 v4 vs v5 ============
f4o = json.load(open(os.path.join(RES, "q_formula_v4.json")))["F4_Qfs"]
f4n = json.load(open(os.path.join(RES, "q_formula_v5.json")))["F4_Qfs"]
assert len(f4o) == len(f4n) == 5
fig, ax = plt.subplots(figsize=(8.6, 4.8))
xs = np.arange(len(f4o))
ax.bar(xs - 0.2, [abs(r["dev"]) for r in f4o], 0.38, color="#c44",
       label=f"v4 (no surface) MAE {np.mean([abs(r['dev']) for r in f4o]):.1f}")
ax.bar(xs + 0.2, [abs(r["dev"]) for r in f4n], 0.38, color="#267",
       label=f"v5 (rescaled) MAE {np.mean([abs(r['dev']) for r in f4n]):.2f}")
ax.set_yscale("log"); ax.set_xticks(xs)
ax.set_xticklabels([f"{r['parent']}\n{r['frags']}" for r in f4o], fontsize=7)
ax.axhline(35, color="k", ls="--", lw=0.8, label="pass band 20%*Q ~ 35 MeV")
ax.set_ylabel("|Q_fs dev| (MeV, log)"); ax.legend(fontsize=8)
ax.set_title("Fission Q: surface term revives F4 (5/5 pass)")
save(fig, "B_fig2_f4_revival.png")

# ============ 图8 (B_fig3): F2 配对两路证伪 ============
pp = json.load(open(os.path.join(RES, "pairing_probe_v2.json")))
ep = json.load(open(os.path.join(RES, "engine_pairing_probe_v2.json")))
fig, ax = plt.subplots(1, 2, figsize=(11, 4.4))
pm = pp["parity_means"]
VN2 = {"偶偶": "even-even", "奇A": "odd-A", "奇奇": "odd-odd"}
xs_ = list(pm.keys())
ys_ = [pm[k][0] for k in xs_]
es_ = [pm[k][1] / np.sqrt(pm[k][2]) for k in xs_]
ax[0].bar([VN2[k] for k in xs_], ys_, yerr=es_, color="#26a", alpha=0.8, capsize=4)
ax[0].set_ylabel("mean residual (MeV)")
ax[0].axhline(0, color="k", lw=0.8)
ax[0].set_title(f"(a) engine route (model-own energy, n={pp['n']}):\n"
                f"parity classes indistinguishable within RMS {pp['rms5']:.1f} MeV\n"
                f"(delta_model = {pp['delta_model']:+.2f} MeV, negligible)")
axd = ep["rows"]
devs_e = [r["dev"] for r in axd]
ax[1].bar(np.arange(len(axd)), devs_e, color="#a62")
ax[1].axhline(0, color="k", lw=0.8)
ax[1].axhline(float(np.mean(devs_e)), color="#c44", ls="--", lw=1,
              label=f"mean {np.mean(devs_e):+.1f} MeV ({sum(1 for d in devs_e if d < 0)}/{len(devs_e)} negative)")
ax[1].set_xticks(np.arange(len(axd)))
ax[1].set_xticklabels([f"A={r['A']}" for r in axd], fontsize=8)
ax[1].legend(fontsize=8); ax[1].set_ylabel("B_ee - B_oo dev (MeV)")
ax[1].set_title("(b) engine route: same-A odd-even deficit")
save(fig, "B_fig3_pairing_falsify.png")

# ============ 图9 (B_fig4): 阀门验证 + T5 回归前残差 ============
valve = json.load(open(os.path.join(RES, "valve_test.json")))
t5 = json.load(open(os.path.join(RES, "t5_shell_correlation.json")))
fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
labels = ["A_raw\n(zero calib.)", "B_calibrated\n(train-region)", "AME ref (P1)"]
vals = [valve["A_raw(原始闭式, 零校准)"]["mean"] * 100,
        valve["B_calibrated(训练区线性校准)"]["mean"] * 100, 97.44]
ax[0].bar(np.arange(len(labels)), vals, color=["#c44", "#c66", "#26a"])
ax[0].axhline(90, color="k", ls="--", lw=0.8)
ax[0].text(0.1, 91, "preregistered 90% line", fontsize=8)
ax[0].set_xticks(np.arange(len(labels)))
ax[0].set_xticklabels(labels, fontsize=8)
ax[0].set_ylabel("blind accuracy (%)")
ax[0].set_title("(a) valve test: energy-side independence fails")
rows5 = t5["rows"]
dmin = [r["dmin_NZ"] for r in rows5]; dev = [r["dev"] for r in rows5]
r_pear = t5["pearson_dmin"]
ax[1].scatter(dmin, dev, c="#26a", s=40)
for rr in rows5:
    ax[1].annotate(rr["nucl"], (rr["dmin_NZ"], rr["dev"]), fontsize=7,
                   xytext=(3, 3), textcoords="offset points")
z = np.polyfit(dmin, dev, 1)
xr = np.linspace(0, max(dmin), 10)
ax[1].plot(xr, np.polyval(z, xr), "k--", lw=1,
           label=f"Pearson r={r_pear:.2f}, Spearman={t5['spearman_dmin']:.2f}")
ax[1].axhline(0, color="k", lw=0.8)
ax[1].legend(fontsize=8); ax[1].set_xlabel("distance to nearest magic number")
ax[1].set_ylabel("Q_beta dev (MeV)")
ax[1].set_title("(b) T5 pre-regression probe: shell attribution (n=10)")
save(fig, "B_fig4_valve_and_T5.png")

print("\n全部图片输出:", FIG)
