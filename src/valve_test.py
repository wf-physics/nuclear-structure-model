# -*- coding: utf-8 -*-
"""
阀门验证 (DeepSeek 判据写死): 把衰变判定层的能量侧 Q 值从 AME 换成模型闭式
判据: ≥95% 论文价值翻倍 / 90–95% 可发表带注 / <90% 放弃"独立预言"表述
基准: AME 版 P1 盲测 mean 97.44% (blind_results.json, mean_free 99.2%)

口径设计 (留出集纪律):
  A_raw        : v4 闭式原样, 零校准 (最严格口径)
  B_calibrated : B_cal = B_model − λ·A − c, (λ,c) 仅用 Z≤10 元素 (A≤40, 束缚核) 拟合,
                 测试区 Z 11–20 不参与任何拟合
不换的部分: 判定管线 (blind_test.py 逐字)、结构闸门 Wb_min/EC_max (本就是模型结构输出)、
           NUBASE 实验标签 (判定层的"考卷", 当然来自实验)、MEerr (仅作闸门尺度参数, 声明)
"""
import os, sys, json
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from nucdata import parse_ame, parse_nubase
from beta_phase_space import f_value
from chain_predict import q_channels

HERE = os.path.dirname(os.path.abspath(__file__))
YR = 3.1557e7
LN2 = np.log(2)
ame = parse_ame(os.path.join(HERE, "data/mass_1.mas20"))
nub = parse_nubase(os.path.join(HERE, "data/nubase_4.mas20"))
ME_H, ME_N = ame[(1, 0)]["ME"]/1000, ame[(0, 1)]["ME"]/1000
ame[(2, 0)] = dict(A=2, ME=2*ME_H + 91.8, MEerr=0.5, BEA=0.0, Qb=None, flag=0, aux=[], El="He")

AV, AS, AC, ASYM = 13.850813715082808, 1.5924382285401872, 0.7264962440025764, 19.50

def B_model(Z, N):
    A = Z + N
    return AV*A - AS*A**(2/3) - AC*Z**2/A**(1/3) - ASYM*(N-Z)**2/A

# ---- 口径 B: 训练区 (Z<=10) 拟合线性校准 ----
tr = [(Z, N) for (Z, N), v in ame.items() if Z <= 10 and Z >= 1 and N >= 1
      and v.get("BEA") is not None and v["BEA"]*(Z+N)/1000. > 0.2]
dB = np.array([B_model(Z, N) - ame[(Z, N)]["BEA"]*(Z+N)/1000. for Z, N in tr])
AA = np.array([Z+N for Z, N in tr], float)
lam, c0 = np.polyfit(AA, dB, 1)
print(f"口径 B 校准 (训练区 Z<=10, n={len(tr)}): B_cal = B_model − ({lam:.4f}·A + {c0:.2f})")
print(f"  训练区扣校准后 RMS = {np.sqrt(((dB - np.polyval([lam, c0], AA))**2).mean()):.2f} MeV")
te = [(Z, N) for (Z, N), v in ame.items() if 10 < Z <= 20 and N >= 1
      and v.get("BEA") is not None and v["BEA"]*(Z+N)/1000. > 0.2]
dte = np.array([B_model(Z, N) - ame[(Z, N)]["BEA"]*(Z+N)/1000. for Z, N in te])
Ate = np.array([Z+N for Z, N in te], float)
print(f"  测试区 (11<=Z<=20, n={len(te)}) 扣校准后 RMS = {np.sqrt(((dte - np.polyval([lam, c0], Ate))**2).mean()):.2f} MeV (泛化检验)")

def make_ame_model(calibrated):
    m = {}
    for (Z, N), v in ame.items():
        if v.get("ME") is None: continue
        B = B_model(Z, N) - (lam*(Z+N) + c0 if calibrated else 0.0)
        me_model = (Z*ME_H + N*ME_N - B)*1000.0     # keV
        d = dict(v); d["ME"] = me_model; d["MEerr"] = v.get("MEerr", 10.0)
        m[(Z, N)] = d
    m[(2, 0)] = dict(ame.get((2, 0)) or dict(A=2))
    m[(2, 0)] = dict(A=2, ME=2*ME_H + 91.8, MEerr=0.5, BEA=0.0, Qb=None, flag=0, aux=[], El="He")
    return m

# ---------------- 判定管线 (复刻 blind_test.py, 逐字协议) ----------------
def build(keys, ame_src):
    out = []
    for k in keys:
        b = nub[k]
        e = "stable" if b["stable"] else b["mode"]
        if e in ("other", None):
            continue
        qc = q_channels(ame_src, k)
        Z, N = k; A = Z + N
        err = ame_src[k]["MEerr"] / 1000.0
        Qe1 = qc["Qp"]
        Qe2 = qc["Q2p"] / 2 if qc["Q2p"] is not None else None
        w_p1 = wkb_P(Qe1, 1, Z-1, A) if Qe1 is not None else np.nan
        w_p2 = wkb_P(Qe2, 1, Z-1, A) if Qe2 is not None else np.nan
        w_a = wkb_P(qc["Qa"], 2, Z-2, A) if qc["Qa"] is not None else np.nan
        f_bm = f_bp = np.nan
        if qc["Qbm"] is not None and qc["Qbm"] > 0:
            f_bm = float(f_value(qc["Qbm"], Z+1))
        if qc["Qbp"] is not None and qc["Qbp"] > 0:
            f_bp = float(f_value(qc["Qbp"], -(Z-1)))
        elif qc["Qec"] is not None and qc["Qec"] > 0:
            f_bp = float(f_value(qc["Qec"], Z-1))
        out.append(dict(key=k, A=A, Z=Z, exp=e, err=err,
                        Qn=qc["Qn"], Q2n=qc["Q2n"], Qp=qc["Qp"], Q2p=qc["Q2p"], Qa=qc["Qa"],
                        w_p1=w_p1, w_p2=w_p2, w_a=w_a, f_bm=f_bm, f_bp=f_bp))
    return out

def wkb_P(Q, l, Zp, A):
    sys.path.insert(0, HERE)
    from baseline_check import wkb_P as _w
    return _w(Q, l, Zp, A)

POOL_KEYS = sorted([k for k in nub if k in ame and nub[k]["A"] <= 40], key=lambda k: nub[k]["A"])

def run_protocol(ame_src, label, n_seeds=30, n_try=2500):
    POOL = build(POOL_KEYS, ame_src)
    nN = len(POOL)
    ERR  = np.array([d["err"] for d in POOL])
    EXP  = np.array([d["exp"] for d in POOL])
    Qn   = np.array([np.nan if d["Qn"] is None else d["Qn"] for d in POOL])
    Q2n  = np.array([np.nan if d["Q2n"] is None else d["Q2n"] for d in POOL])
    Qp   = np.array([np.nan if d["Qp"] is None else d["Qp"] for d in POOL])
    Q2p  = np.array([np.nan if d["Q2p"] is None else d["Q2p"] for d in POOL])
    Qa   = np.array([np.nan if d["Qa"] is None else d["Qa"] for d in POOL])
    Wp1  = np.array([d["w_p1"] for d in POOL]); Wp2 = np.array([d["w_p2"] for d in POOL])
    Wa   = np.array([d["w_a"] for d in POOL])
    Fbm  = np.array([d["f_bm"] for d in POOL]); Fbp = np.array([d["f_bp"] for d in POOL])
    ISB_BM, ISB_BP = ~np.isnan(Fbm), ~np.isnan(Fbp)

    def predict(theta, idx=None):
        sig_k, nu, sa, ec, ft, ts = theta
        ii = np.arange(nN) if idx is None else np.asarray(idx)
        thr = sig_k * ERR[ii]
        open_n = np.nanmax(np.c_[Qn[ii], Q2n[ii]], axis=1) > thr
        open_a = Qa[ii] > thr
        cond_p1 = Qp[ii] > thr
        cond_p2 = Q2p[ii] > thr
        open_p = cond_p1 | cond_p2
        w_p = np.where(cond_p1, Wp1[ii], Wp2[ii])
        r_n  = np.where(open_n, nu, -np.inf)
        r_p  = np.where(open_p, nu * w_p, -np.inf)
        r_a  = np.where(open_a, nu * sa * Wa[ii], -np.inf)
        bsc  = LN2 / ft
        r_bm = np.where(ISB_BM[ii], bsc * Fbm[ii], -np.inf)
        r_bp = np.where(ISB_BP[ii], ec * bsc * Fbp[ii], -np.inf)
        R = np.c_[r_n, r_p, r_a, r_bm, r_bp]
        modes = np.array(["n", "p", "alpha", "bm", "bp"])
        j = np.argmax(R, axis=1)
        lam_ = R[np.arange(len(ii)), j]
        T = np.where(lam_ > 0, LN2 / lam_ / YR, np.inf)
        pred = np.where(T > ts, "stable", modes[j])
        pred[np.isnan(R).all(axis=1) | (lam_ == -np.inf)] = "stable"
        return pred

    def acc(pred, idx):
        return float(np.mean(pred == EXP[idx]))

    RANGES = [(0.0, 1.3), (18, 24), (-4, 0), (-1, 0), (2.5, 4.5), (8, 25)]
    def fit(idx, rng):
        best = (-1, None)
        ii = np.asarray(idx)
        for _ in range(n_try):
            th = tuple(10**rng.uniform(*r) for r in RANGES)
            a = acc(predict(th, ii), ii)
            if a > best[0]:
                best = (a, th)
        return best

    scores, scores_free = [], []
    patch_keys = set()
    # 补丁依赖核素沿用原存档标记 (诚实条款)
    try:
        br = json.load(open(os.path.join(HERE, "results/blind_results.json")))
        patch_keys = {tuple(int(x) for x in k.split(",")) for k in br.get("patch_dependent", [])}
    except Exception:
        pass
    for seed in range(n_seeds):
        rng = np.random.default_rng(seed)
        idx = rng.permutation(nN)
        ntr = int(0.7*nN)
        tr_i, te_i = idx[:ntr], idx[ntr:]
        _, th = fit(tr_i, rng)
        pred = predict(th, te_i)
        scores.append(acc(pred, te_i))
        free = [j for j in te_i if POOL[j]["key"] not in patch_keys]
        scores_free.append(acc(pred, np.array(free)) if free else acc(pred, te_i))
        if (seed+1) % 10 == 0:
            print(f"  [{label}] seed {seed+1}/{n_seeds}  累计 mean={np.mean(scores):.4f}", flush=True)
    return dict(n=n_seeds, mean=float(np.mean(scores)), std=float(np.std(scores)),
                min=float(np.min(scores)), max=float(np.max(scores)),
                mean_free=float(np.mean(scores_free)))

print()
out = {}
for cal, label in [(False, "A_raw(原始闭式, 零校准)"), (True, "B_calibrated(训练区线性校准)")]:
    print("="*84)
    print(f"阀门验证 {label} (P1 协议: 30 seeds × 70/30, 训练集拟合参数)")
    print("="*84)
    ame_m = make_ame_model(cal)
    out[label] = run_protocol(ame_m, label)
    r = out[label]
    print(f"  → 盲测 mean={r['mean']*100:.2f}%  std={r['std']*100:.2f}  min={r['min']*100:.2f}%  "
          f"mean_free(剔补丁核素)={r['mean_free']*100:.2f}%")
    v = ("≥95% → 论文价值翻倍" if r['mean'] >= 0.95 else
         "90–95% → 可发表, 注明闸门需 AME 兜底" if r['mean'] >= 0.90 else
         "<90% → 放弃'独立预言'表述")
    print(f"  → DeepSeek 判据裁决: {v}")

print()
print("="*84)
print("阀门验证汇总 (基准: AME 版 P1 mean=97.44%, mean_free=99.2%)")
print("="*84)
for k, v in out.items():
    print(f"  {k:34s} mean={v['mean']*100:6.2f}%  min={v['min']*100:6.2f}%  free={v['mean_free']*100:6.2f}%")
json.dump(out, open(os.path.join(HERE, "results/valve_test.json"), "w"), indent=1, ensure_ascii=False)
print("存档: results/valve_test.json")
