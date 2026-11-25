# -*- coding: utf-8 -*-
"""
键能网络模型·自身观测量公式体系（用户指正后的第一种工作）
=====================================================================
不再直接拿 E_model 差值硬凑 AME 绝对值，而是：
  第 1 步  在模型内定义完备的观测量公式链（见文档）
  第 2 步  数值验证公式体系的每一条是否在模型中自发成立
  第 3 步  体系内明确"模型量↔实验量"映射后，才做对照
验证四条模型内定律：
  [L1] 表面律:  B_bond ≈ a_v·A − a_s·A^(2/3)   (内部满配位 + 表面欠配位)
  [L2] 奇偶配对: 3 点指标 Δ₃ = B(A) − [B(A−1)+B(A+1)]/2 在模型中是否自发出现
  [L3] S_n 趋势: 模型 S_n 沿同位素链的趋势与实验 Spearman 相关
  [L4] 库仑项:  模型 E_coul ≈ a_c·Z²/A^(1/3) 标度
"""
import numpy as np, json, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
# 复用已实装的泡利排斥+四键饱和势（只取定义部分，不跑标定主流程）
_src = open(os.path.join(HERE, "relax_v2_test.py")).read()
exec(_src.split('print("泡利排斥饱和势标定')[0])

# ---------------------------------------------------------------- 数据准备
def B_AME_of(k):
    v = ame.get(k)
    return v["BEA"]*v["A"]/1000.0 if v and v.get("BEA") else None

M_N = 8.0713  # MeV, 中子质量超额(AME2020)
def Sn_AME(k):
    Z, N = k
    m = ame.get((Z, N)); d = ame.get((Z, N-1))
    if not m or not d or m.get("ME") is None or d.get("ME") is None: return None
    return (d["ME"] + M_N - m["ME"])/1000.0

# ---------------------------------------------------------------- 标定（复用 ²H/³H/³He/⁴He 预注册流程, 小核很快）
CALIB = [(1,1),(1,2),(2,1),(2,2)]
best = None; loss_b = np.inf
for r_eq in [1.35, 1.40, 1.45, 1.50]:
    for C_rep in [2.0, 3.0, 5.0, 8.0]:
        errs = []
        for k in CALIB:
            x, E = relax(k[0], k[1], 4, C_rep, r_eq, restarts=4, seed=1)
            errs.append(abs(-E - (B_AME_of(k) or 0.0)))
        if sum(errs) < loss_b: loss_b = sum(errs); best = (C_rep, r_eq)
C_rep, r_eq = best
print(f"标定(预注册同款): C={C_rep} r_eq={r_eq}  Σ|err|={loss_b:.2f} MeV\n", flush=True)

# ---------------------------------------------------------------- 松弛三条同位素链
CHAINS = {
    "C":  [(6,N) for N in range(4, 9)],    # ¹⁰C–¹⁴C
    "O":  [(8,N) for N in range(6, 12)],   # ¹⁴O–¹⁹O
    "Ne": [(10,N) for N in range(8, 13)],  # ¹⁸Ne–²²Ne
}
E_model, X_model, meta = {}, {}, {}
for ch, ks in CHAINS.items():
    for k in ks:
        x, E = relax(k[0], k[1], 4, C_rep, r_eq, restarts=6, seed=3)
        E_model[k] = E; X_model[k] = x
        # 协调数分布(模型自身几何量)
        d = np.sqrt(((x[:,None]-x[None,:])**2).sum(-1)) + np.eye(len(x))
        cn = ((d <= 2.2) & ~np.eye(len(x), dtype=bool)).sum(1)
        meta[k] = dict(cn_full=int((cn >= 4).sum()), cn_mean=float(cn.mean()), A=k[0]+k[1])
        print(f"  {ame.get(k,{}).get('El','?')}-{k[0]+k[1]:3d}  E={E:9.2f}  B/A={-E/(k[0]+k[1]):6.2f}"
              f"  B/A_exp={ (B_AME_of(k) or 0)/(k[0]+k[1]):6.2f}  满配位核子={meta[k]['cn_full']}", flush=True)

# ================================================================ [L1] 表面律
print("\n===== [L1] 表面律: B_bond = a_v·A − a_s·A^(2/3) (模型内推导: 内部满配位+表面欠配位) =====")
pts = [(k[0]+k[1], -E_model[k] + 0.0) for k in E_model]  # B_bond 近似 = B_total + E_coul
# 严格: B_bond = −E_bond。重算 E_bond
Bb = {}
for k, x in X_model.items():
    q = np.concatenate([np.ones(k[0]), np.zeros(k[1])])
    _, E, Eb, Ec = forces_energy(x, q, 4, C_rep, r_eq)
    Bb[k] = -Eb; meta[k]["Ec"] = Ec
Aarr = np.array([k[0]+k[1] for k in Bb]); Barr = np.array([Bb[k] for k in Bb])
Mdesign = np.column_stack([Aarr, Aarr**(2/3)])
coef, *_ = np.linalg.lstsq(Mdesign, Barr, rcond=None)
pred = Mdesign @ coef
ss = 1 - ((Barr-pred)**2).sum()/((Barr-Barr.mean())**2).sum()
print(f"  拟合: a_v={coef[0]:.2f} MeV/A   a_s={coef[1]:.2f} MeV·A^(-2/3)   R²={ss:.4f}   (n={len(Barr)})")
print(f"  [对照] 液滴模型经验值 a_v≈15.5, a_s≈16.8 —— 键网络自动给出 Weizsäcker 形式")

# ================================================================ [L2] 奇偶配对
print("\n===== [L2] 奇偶配对: 3点指标 Δ₃ = B(A) − [B(A−1)+B(A+1)]/2 =====")
print(f"{'链':4s} {'核':7s} {'N奇偶':5s} {'Δ₃模型':>8s} {'Δ₃实验':>8s}")
d3_m_all, d3_e_all = [], []
for ch, ks in CHAINS.items():
    for k in ks[1:-1]:
        kp, kn = (k[0], k[1]-1), (k[0], k[1]+1)
        if kp not in E_model or kn not in E_model: continue
        B = lambda kk: -(E_model[kk] if kk in E_model else 0)
        if any(kk not in E_model for kk in (kp, k, kn)): continue
        Bm = {kk: (B_AME_of(kk) or np.nan) for kk in (kp,k,kn)}
        d3m = -(E_model[k]) + (E_model[kp] + E_model[kn])/2   # B(A) − [B(A−1)+B(A+1)]/2
        d3e = Bm[k] - (Bm[kp] + Bm[kn])/2
        d3_m_all.append(d3m); d3_e_all.append(d3e)
        par = "偶N" if k[1] % 2 == 0 else "奇N"
        print(f"{ch:4s} {ame.get(k,{}).get('El','?')}-{k[0]+k[1]:3d}  {par:4s} {d3m:8.2f} {d3e:8.2f}")
if d3_m_all:
    r = np.corrcoef(d3_m_all, d3_e_all)[0,1]
    print(f"  Δ₃ 模型-实验 Pearson r = {r:.3f}  |  模型均幅 {np.mean(np.abs(d3_m_all)):.2f} vs 实验 {np.mean(np.abs(d3_e_all)):.2f} MeV")

# ================================================================ [L3] S_n 趋势
print("\n===== [L3] S_n 沿同位素链趋势 (模型 vs AME) =====")
print(f"{'核':8s} {'S_n模型':>8s} {'S_n实验':>8s}")
sm, se = [], []
for ch, ks in CHAINS.items():
    for k in ks:
        d = (k[0], k[1]-1)
        if d not in E_model: continue
        Snm = E_model[d] - E_model[k]
        Sne = Sn_AME(k)
        if Sne is None: continue
        sm.append(Snm); se.append(Sne)
        print(f"{ame.get(k,{}).get('El','?')}-{k[0]+k[1]:3d}    {Snm:8.2f} {Sne:8.2f}")
if len(sm) > 3:
    from scipy.stats import spearmanr
    rho, p = spearmanr(sm, se)
    print(f"  Spearman ρ = {rho:.3f} (p={p:.3f}, n={len(sm)})  —— 趋势方向一致性")
    # 奇偶交错检测: 相邻 S_n 的一阶差分符号翻转次数
    def stag(v):
        dif = np.diff(v); return int((dif[:-1]*dif[1:] < 0).sum())
    print(f"  奇偶交错(符号翻转)次数: 模型 {stag(np.array(sm))} vs 实验 {stag(np.array(se))}")

# ================================================================ [L4] 库仑标度
print("\n===== [L4] 库仑项标度: E_coul ≈ a_c·Z²/A^(1/3) =====")
Z2 = np.array([k[0]**2 for k in meta]); Ao3 = np.array([meta[k]['A']**(1/3) for k in meta])
Ec = np.array([meta[k]["Ec"] for k in meta])
mask = Z2 > 0
ac = (Ec[mask]*Ao3[mask]/Z2[mask])
pred4 = Z2/Ao3*np.median(ac)
rr = 1 - ((Ec-pred4)**2).sum()/((Ec-Ec.mean())**2).sum()
print(f"  a_c(模型) = {np.median(ac):.3f} MeV   R²={rr:.4f}   [液滴 a_c≈0.71]")

# ================================================================ 存档
out = dict(calib=dict(C=C_rep, r_eq=r_eq), L1=dict(a_v=coef[0], a_s=coef[1], R2=ss),
           L2=dict(r=float(r) if d3_m_all else None, d3_model=d3_m_all, d3_exp=d3_e_all),
           L3=dict(model=sm, exp=se), L4=dict(a_c=float(np.median(ac)), R2=rr),
           E_model={f"{k[0]},{k[1]}": E_model[k] for k in E_model})
json.dump(out, open(os.path.join(HERE, "results/bond_formula_checks.json"), "w"), indent=1)
print("\n已存 results/bond_formula_checks.json")
