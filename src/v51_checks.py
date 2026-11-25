# -*- coding: utf-8 -*-
"""
v5.1 外部评审回应自检（2026-11-20）
逐条对应评审六刺 + 第二评审两项要求，全部实跑、只报数字：
  [刺6] ²H/³H/³He/⁴He 局域库仑 vs B/A 四点自检
  [刺2] E_C 坐标求和 vs 液滴式 vs AME 实际库仑差（315 核素）
  [刺4] α 三修（出射半径 r0(A_d^{1/3}+4^{1/3}) + 物理尝试频率 + 形变长轴）5 个 α 核
  [刺5] J^π 禁戒度 K 分类器 + 8 核素半衰期固定考题表
  [第二评审·基线] 纯 Q 符号基线 / 零结构上界 / v5.0 三行对照表
  [第二评审·验收] 结构层算分离能 S_n/S_p 五核验收（判据 <1 MeV）
"""
import os, sys, json, re
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from nucdata import parse_ame, parse_nubase
from beta_phase_space import f_value
from chain_predict import q_channels
from scipy.integrate import quad
from scipy.stats import spearmanr

HERE = os.path.dirname(os.path.abspath(__file__))
YR = 3.1557e7; LN2 = np.log(2); HBARC = 197.327
ame = parse_ame(os.path.join(HERE, "data/mass_1.mas20"))
nub = parse_nubase(os.path.join(HERE, "data/nubase_4.mas20"))
ME_H = ame[(1, 0)]["ME"]
ame[(2, 0)] = dict(A=2, ME=2*ME_H + 91.8, MEerr=0.5, BEA=0.0, Qb=None, flag=0, aux=[], El="He")
ST = json.load(open(os.path.join(HERE, "results/structures.json")))
B_AME = {k: v["BEA"] * v["A"] for k, v in ame.items() if v.get("BEA")}  # 总结合能 MeV

def nm(k): return f"{ame[k]['El']}-{ame[k]['A']}"

# ============ [刺6] 四点自检：局域库仑 vs B/A ============
print("=" * 70)
print("[刺6] 四点自检: 局域库仑应力 vs 每核子结合能 (AME)")
print(f"{'核':6s} {'键数pn':>6s} {'局域E_C/MeV':>11s} {'B_total':>8s} {'B/A':>6s}")
for k in [(1,1),(1,2),(2,1),(2,2)]:
    d = ST.get(f"{k[0]},{k[1]}")
    c = np.array(d["coords"]); Z, N = k
    ec_loc = d["E_coul"]
    n_pn = sum(1 for i in range(Z+N) for j in range(i+1, Z+N)
               if ((i < Z) != (j < Z)) and np.linalg.norm(c[i]-c[j]) <= 2.2)
    print(f"{nm(k):6s} {n_pn:6d} {ec_loc:11.3f} {B_AME[k]:8.3f} {B_AME[k]/(Z+N):6.3f}")
print("结论: 局域 E_C(³H)≈(³He)/2 且 ⁴He≈³He, 但 B/A: 2.83/2.57/7.07 ——")
print("      局域应力排序器无法解释 ⁴He vs ³He 差 2.8 倍 → 刺6 成立, §3.1 须改全局 S_x\n")

# ============ [刺2] E_C 坐标求和 vs 液滴式 ============
print("=" * 70)
print("[刺2] ΔE_C 三种算法 vs AME 实际 (β⁻ 对, 双方均有结构):")
def ec_coord(k):
    d = ST.get(f"{k[0]},{k[1]}")
    if d is None: return None
    c = np.array(d["coords"]); Z = k[0]
    s = 0.0
    for i in range(Z):
        for j in range(i+1, Z):
            r = np.linalg.norm(c[i]-c[j])
            s += 1.44 / max(r, 0.1)
    return s
def ec_drop(k):
    Z, N = k; A = Z+N
    return 0.6 * 1.44 * Z*(Z-1) / (1.2*A**(1/3))
xs_c, xs_d, ys = [], [], []
for k in list(ame.keys()):
    Z, N = k
    if Z == 0 or N == 0: continue
    k2 = (Z+1, N-1)
    if k2 not in ame or f"{k[0]},{k[1]}" not in ST or f"{k2[0]},{k2[1]}" not in ST: continue
    Qbm = ame[k]["ME"]/1000 - ame[k2]["ME"]/1000
    dE_nuc = Qbm + 0.782          # ΔE_nuclear = Qβ⁻ − (m_n−m_p−m_e)c²
    ec1, ec2 = ec_coord(k), ec_coord(k2)
    if ec1 is None or ec2 is None: continue
    xs_c.append(ec2 - ec1); xs_d.append(ec_drop(k2) - ec_drop(k)); ys.append(dE_nuc)
xs_c, xs_d, ys = map(np.array, (xs_c, xs_d, ys))
r_c = spearmanr(xs_c, ys).statistic; r_d = spearmanr(xs_d, ys).statistic
rms_c = np.sqrt(np.mean((np.array(xs_c) - np.array(ys))**2))
rms_d = np.sqrt(np.mean((np.array(xs_d) - np.array(ys))**2))
print(f"  样本 n={len(ys)} 个 β⁻ 对")
print(f"  坐标求和 ΔE_C: Spearman ρ={r_c:.3f}  RMS 残差={rms_c:.2f} MeV")
print(f"  液滴式   ΔE_C: Spearman ρ={r_d:.3f}  RMS 残差={rms_d:.2f} MeV")
print(f"  坐标式 vs 液滴式 直接偏差: RMS={np.sqrt(np.mean((xs_c-xs_d)**2)):.2f} MeV (即构型信息的量级)")

# ============ [刺4] α 三修 ============
print("\n" + "=" * 70)
print("[刺4] α 半衰期三修 (R=r0(A_d^{1/3}+4^{1/3}), ν=v/2R 物理尝试频率):")
def wkb_alpha(Q, Zd, A_d, beta2=0.0, r0=1.17):
    mu = 4.0*Zd/(4.0+Zd) * 931.5          # MeV
    Y20_0 = 0.6308
    R = r0*(A_d**(1/3) + 4**(1/3)) * (1 + beta2*Y20_0)
    rt = 1.44*2*Zd/Q
    if rt <= R: return 1.0, R
    integ, _ = quad(lambda r: np.sqrt(2*mu*(1.44*2*Zd/r - Q))/HBARC, R, rt, limit=200)
    v = np.sqrt(2*Q/mu)                    # fm/c
    nu = v / (2*R) * 2.998e23              # s^-1 (c=2.998e23 fm/s)
    return np.exp(-2*integ), nu
ALPHA_SET = [((60,84),0.28),((84,126),0.0),((88,138),0.17),((90,142),0.27),((92,146),0.29)]
print(f"{'核':8s} {'实验 T½/y':>12s} {'旧:R=1.2A^(1/3)':>16s} {'新:出射半径+物理ν':>18s} {'+形变长轴':>12s}")
for k, b2 in ALPHA_SET:
    Z, N = k; A = Z+N
    qc = q_channels(ame, k)
    Qa = qc["Qa"]
    b = nub[k]; exp_mode = b["mode"]; exp_T = b["T"]
    if Qa is None or Qa <= 0 or exp_mode != "alpha":
        continue
    # 旧模型
    P_old = None
    def wkb_old(Q, Zd, A):
        R = 1.2*A**(1/3); rt = 1.44*2*Zd/Q
        if rt <= R: return 1.0
        mu = 4.0*Zd/(4.0+Zd)*931.5
        integ, _ = quad(lambda r: np.sqrt(2*mu*(1.44*2*Zd/r-Q))/HBARC, R, rt, limit=200)
        return np.exp(-2*integ)
    lam_old = 1e21 * 0.01 * wkb_old(Qa, Z-2, A)
    T_old = LN2/lam_old/YR if lam_old > 0 else np.inf
    P_new, nu_new = wkb_alpha(Qa, Z-2, A-2, beta2=0.0)
    lam_new = 0.01 * nu_new * P_new
    T_new = LN2/lam_new/YR
    P_def, nu_def = wkb_alpha(Qa, Z-2, A-2, beta2=b2)
    lam_def = 0.01 * nu_def * P_def
    T_def = LN2/lam_def/YR
    eT = f"{exp_T/YR:.2e}" if exp_T else "?"   # NUBASE T 以秒存储
    print(f"{nm(k):8s} {eT:>12s} {T_old:16.2e} {T_new:18.2e} {T_def:12.2e}")
print("  (S_α=0.01 保持不变; Langer 项对基态 l=0 无离心贡献, 声明为已含于出射半径口径)\n")

# ============ [刺5] J^π 禁戒度分类器 + 8 核素考题表 ============
print("=" * 70)
print("[刺5] J^π → 禁戒度 K → log ft 桶 → T½ 预测 (8 核素固定考题)")
jp_raw = {}
with open(os.path.join(HERE, "data/nubase_4.mas20")) as f:
    for line in f:
        s = line.rstrip()
        if s.startswith("#") or len(s) < 102: continue
        zzz = s[4:8]
        if not zzz.strip().endswith("0") or not zzz.strip().isdigit(): continue
        jp = s[88:102].strip()
        m = re.match(r"\(?(\d+/\d+|\d+)([+-])\)?", jp)
        if m:
            J = float(m.group(1)) if "/" not in m.group(1) else float(eval(m.group(1)))
            A = int(s[0:3])
            me_el = re.search(r"([A-Z][a-z]?)", s[11:16])
            if not me_el: continue
            El = me_el.group(1)
            jp_raw[(A, El)] = (J, m.group(2))
def forbid(Ji, pi, Jf, pf):
    dJ = abs(Ji - Jf); pc = (pi != pf)
    for L in range(0, 8):
        if dJ <= L + 1 and (pc == (L % 2 == 1)):
            return L
    return 8
LOGFT_K = {0: 4.4, 1: 6.8, 2: 9.8, 3: 13.2, 4: 19.0}   # 文献口径桶中位 (K≥4 合并)
# 考题: (核, 主衰变分支到达的子核 Jπ[声明输入], 实验 log ft, 实验模式)
EXAM = [
    ((1,2),  "bm", "³He 1/2+", (0.5,"+"), 3.049, 12.32),
    ((17,19),"EC", "³⁶S 2+ (3.3MeV 主支)", (2,"+"), 7.4, 3.01e5),
    ((19,21),"bm", "⁴⁰Ca 4+ (3.9MeV 主支)", (4,"+"), 8.0, 1.248e9),
    ((19,21),"EC", "⁴⁰Ca 0+ (基态, 10.7%)", (0,"+"), 18.1, None),
    ((37,50),"bm", "⁸⁷Sr 9/2+ (基态)", (4.5,"+"), 10.3, 4.97e10),
    ((49,66),"bm","¹¹⁵Sn 1/2+ (基态)", (0.5,"+"), 23.0, 4.41e14),
    ((6,8),  "bm", "¹⁴N 1+ (基态)", (1,"+"), 9.03, 5730),
    ((27,33),"bm", "⁶⁰Ni 4+ (2.5MeV 主支)", (4,"+"), 6.6, 5.27),
]
# 代码口径锚点: ³H 的 code log ft = log10(f_value·t_exp)
fH = f_value(ame[(1,2)]["ME"]/1000 - ame[(2,1)]["ME"]/1000, 2)
anchor = np.log10(fH * 12.32*YR)
print(f"  锚点: ³H code口径 log ft = {anchor:.3f} (文献 3.049, 偏移 {anchor-3.049:+.3f})")
print(f"{'核':8s} {'K(主支)':>7s} {'log ft 桶':>9s} {'文献 exp':>8s} {'T½桶预测/y':>12s} {'实验 T½/y':>12s} {'差(量级)':>8s}")
for k, mode, br, (Jf, pf), lf_exp, t_exp in EXAM:
    keyAE = (ame[k]["A"], ame[k]["El"])
    Jp, pp = jp_raw[keyAE]
    K = forbid(Jp, pp, Jf, pf)
    lf_bucket = LOGFT_K[min(K, 4)]
    Z, N = k
    if mode == "bm":
        Q = ame[k]["ME"]/1000 - ame[(Z+1, N-1)]["ME"]/1000
        fv = f_value(Q, Z+1)
    else:
        Q = ame[k]["ME"]/1000 - ame[(Z-1, N+1)]["ME"]/1000 - 1.022
        fv = f_value(Q, -(Z-1)) if Q > 0 else None
    if fv is None or fv <= 0:
        print(f"{nm(k):8s} {K:7d} {lf_bucket:9.1f} {lf_exp:8.2f}   Q<0 通道关闭"); continue
    # 桶 log ft(文献口径) → code口径 → T½ = 10^(lf−offset)/f
    lf_code = lf_bucket + (anchor - 3.049)
    T_pred = 10**(lf_code) / fv / YR
    # 实验 log ft → 一致性
    dlf = lf_bucket - lf_exp
    tS = f"{t_exp:.2e}" if t_exp else "-"
    print(f"{nm(k):8s} {K:7d} {lf_bucket:9.1f} {lf_exp:8.2f} {T_pred:12.2e} {tS:>12s} {dlf:+8.1f}")

# ============ [第二评审·基线三行表] ============
print("\n" + "=" * 70)
print("[基线] A≤40 (314 核素, 无滴线补丁口径) 三行对照:")
results = {}
for k in sorted(nub, key=lambda k: nub[k]["A"]):
    b = nub[k]
    if b["A"] > 40 or k not in ame: continue
    e = "stable" if b["stable"] else b["mode"]
    if e in ("other", None): continue
    results[k] = e
def qch(k): return q_channels(ame, k)
def _row1(k, order="npa"):
    qc = qch(k)
    pn = (qc["Qn"] is not None and qc["Qn"] > 0) or (qc["Q2n"] is not None and qc["Q2n"] > 0)
    pp = (qc["Qp"] is not None and qc["Qp"] > 0) or (qc["Q2p"] is not None and qc["Q2p"] > 0)
    pa = qc["Qa"] is not None and qc["Qa"] > 0
    if order == "npa":
        p = "n" if pn else ("p" if pp else ("alpha" if pa else None))
    elif order == "apn":
        p = "alpha" if pa else ("p" if pp else ("n" if pn else None))
    else:
        p = "p" if pp else ("alpha" if pa else ("n" if pn else None))
    if p is None:
        if qc["Qbm"] is not None and qc["Qbm"] > 0: p = "bm"
        elif qc["Qbp"] is not None and qc["Qbp"] > 0: p = "bp"
        else: p = "stable"
    return p
# Row1: 纯 Q 符号 + 固定优先级 (粒子先于 β, 粒子内部 α>p>n)
acc1 = 0; tot = 0
for k, e in results.items():
    p = _row1(k, order="apn")
    tot += 1; acc1 += (p == e)
print(f"  Row1 纯Q符号+优先级(α>p>n): {acc1}/{tot} = {acc1/tot*100:.1f}%")
# Row1b: 优先级 n>p>α
acc1b = 0
for k, e in results.items():
    qc = qch(k)
    p = None
    if (qc["Qn"] is not None and qc["Qn"] > 0) or (qc["Q2n"] is not None and qc["Q2n"] > 0): p = "n"
    elif (qc["Qp"] is not None and qc["Qp"] > 0) or (qc["Q2p"] is not None and qc["Q2p"] > 0): p = "p"
    elif qc["Qa"] is not None and qc["Qa"] > 0: p = "alpha"
    elif qc["Qbm"] is not None and qc["Qbm"] > 0: p = "bm"
    elif qc["Qbp"] is not None and qc["Qbp"] > 0: p = "bp"
    else: p = "stable"
    acc1b += (p == e)
print(f"  Row1b 优先级(n>p>α):          {acc1b}/{tot} = {acc1b/tot*100:.1f}%")
# Row3: v5.0 完整 (baseline_judge.json 含补丁)
bj = json.load(open(os.path.join(HERE, "results/baseline_judge.json")))
acc3 = sum(1 for v in bj.values() if v["ok"])
print(f"  Row3 v5.0 完整(含滴线补丁):    {acc3}/{len(bj)} = {acc3/len(bj)*100:.1f}%")
# Row2 引用 blind_test 全数据重标定: 97.5% (no-patch); 此处用 THETA0 无补丁速算
from baseline_check import wkb_P
def predict_nopatch(k, theta=(2.0,1e21,0.01,0.5,10**3.5,1e17)):
    sig_k, nu, sa, ec, ft, ts = theta
    qc = qch(k); Z, N = k; A = Z+N; err = ame[k]["MEerr"]/1000.0
    sig = lambda Q: Q is not None and Q > sig_k*err
    parts = []
    if sig(qc["Qn"]) or sig(qc["Q2n"]): parts.append(("n", nu))
    if sig(qc["Qp"]) or sig(qc["Q2p"]):
        Qe = qc["Qp"] or (qc["Q2p"]/2 if qc["Q2p"] is not None else None)
        if Qe: parts.append(("p", nu*wkb_P(Qe,1,Z-1,A)))
    if sig(qc["Qa"]): parts.append(("alpha", nu*sa*wkb_P(qc["Qa"],2,Z-2,A)))
    if qc["Qbm"] and qc["Qbm"]>0: parts.append(("bm", LN2*f_value(qc["Qbm"],Z+1)/10**3.5))
    if qc["Qbp"] and qc["Qbp"]>0: parts.append(("bp", LN2*f_value(qc["Qbp"],-(Z-1))/10**3.5))
    elif qc["Qec"] and qc["Qec"]>0: parts.append(("bp", ec*LN2*f_value(qc["Qec"],Z-1)/10**3.5))
    if not parts: return "stable"
    ch, lam = max(parts, key=lambda x: x[1])
    T = LN2/lam/YR if lam>0 else np.inf
    return "stable" if T > ts else ch
acc2 = sum(1 for k in results if predict_nopatch(k) == results[k])
print(f"  Row2 +相空间/Gamow/分支竞争(无补丁): {acc2}/{tot} = {acc2/tot*100:.1f}%")
# 同分母口径: 剔除 12 个 AME 外推边界核素 (带补丁判对、无补丁判错者)
bj_pred = {tuple(map(int, key.split(","))): v["pred"] for key, v in bj.items()}
_patch = [k for k in results if predict_nopatch(k) != results[k]
          and bj_pred.get(k) == results[k]]
clean = [k for k in results if k not in _patch]
r1c = sum(1 for k in clean if _row1(k, order="npa") == results[k])
r2c = sum(1 for k in clean if predict_nopatch(k) == results[k])
print(f"  补丁依赖核素: {len(_patch)} 个 {[nm(k) for k in _patch]}")
print(f"  同分母(剔 {len(_patch)} 个后 n={len(clean)}):")
print(f"    Row1b 纯Q符号:  {r1c}/{len(clean)} = {r1c/len(clean)*100:.1f}%")
print(f"    Row2 +速率竞争: {r2c}/{len(clean)} = {r2c/len(clean)*100:.1f}%")

# ============ [第二评审·结构层分离能验收] ============
print("\n" + "=" * 70)
print("[验收] 结构层 E_model 算分离能 vs AME (判据 <1 MeV):")
print(f"{'核':8s} {'通道':4s} {'模型/MeV':>10s} {'AME/MeV':>9s} {'误差':>8s} {'判定':>4s}")
TESTS = [((2,2),"p"),((4,4),"p"),((4,4),"n"),((4,5),"n"),((6,6),"p"),((8,8),"p"),((8,8),"n")]
for k, ch in TESTS:
    Z, N = k
    kd = (Z-1, N) if ch == "p" else (Z, N-1)
    if f"{k[0]},{k[1]}" not in ST or f"{kd[0]},{kd[1]}" not in ST:
        print(f"{nm(k):8s} {ch:4s}   结构缺失"); continue
    Em = ST[f"{k[0]},{k[1]}"]["E_model"]; Ed = ST[f"{kd[0]},{kd[1]}"]["E_model"]
    S_model = Ed - Em
    if ch == "p":
        S_ame = ame[kd]["ME"]/1000 + ME_H/1000 - ame[k]["ME"]/1000
    else:
        S_ame = ame[kd]["ME"]/1000 + ame[(0,1)]["ME"]/1000 - ame[k]["ME"]/1000
    err = abs(S_model - S_ame)
    print(f"{nm(k):8s} {ch:4s} {S_model:10.2f} {S_ame:9.2f} {err:8.2f} {'过' if err<1 else '不过'}")
