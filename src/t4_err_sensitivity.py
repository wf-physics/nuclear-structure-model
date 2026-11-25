# -*- coding: utf-8 -*-
"""
T4 敏感性实验 (DeepSeek 建议): 判定阈值 thr = sig_k·ERR 中的 ERR 从 AME MEerr 换成固定 0.3 MeV,
检验衰变判定结论是否对闸门尺度参数敏感。
对照设计:
  ① AME 能量侧: ERR=AME (原版基准 97.44%) vs ERR=0.3 (敏感性版)
  ② 闭式 A_raw: ERR=AME (原版 33.93%) vs ERR=0.3
参数搜索范围 RANGES 冻结不变 (协议一致性)。
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

def make_ame_model_raw():
    m = {}
    for (Z, N), v in ame.items():
        if v.get("ME") is None: continue
        B = B_model(Z, N)
        me_model = (Z*ME_H + N*ME_N - B)*1000.0
        d = dict(v); d["ME"] = me_model; d["MEerr"] = v.get("MEerr", 10.0)
        m[(Z, N)] = d
    m[(2, 0)] = dict(A=2, ME=2*ME_H + 91.8, MEerr=0.5, BEA=0.0, Qb=None, flag=0, aux=[], El="He")
    return m

def wkb_P(Q, l, Zp, A):
    from baseline_check import wkb_P as _w
    return _w(Q, l, Zp, A)

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
        Qe2 = qc["Q2p"] / 2 if qc["Q2p"] is not None else None
        out.append(dict(key=k, A=A, exp=e, err=err,
                        Qn=qc["Qn"], Q2n=qc["Q2n"], Qp=qc["Qp"], Q2p=qc["Q2p"], Qa=qc["Qa"],
                        w_p1=wkb_P(qc["Qp"], 1, Z-1, A) if qc["Qp"] is not None else np.nan,
                        w_p2=wkb_P(Qe2, 1, Z-1, A) if Qe2 is not None else np.nan,
                        w_a=wkb_P(qc["Qa"], 2, Z-2, A) if qc["Qa"] is not None else np.nan,
                        f_bm=float(f_value(qc["Qbm"], Z+1)) if (qc["Qbm"] is not None and qc["Qbm"] > 0) else np.nan,
                        f_bp=(float(f_value(qc["Qbp"], -(Z-1))) if (qc["Qbp"] is not None and qc["Qbp"] > 0)
                              else (float(f_value(qc["Qec"], Z-1)) if (qc["Qec"] is not None and qc["Qec"] > 0) else np.nan))))
    return out

POOL_KEYS = sorted([k for k in nub if k in ame and nub[k]["A"] <= 40], key=lambda k: nub[k]["A"])

def run_protocol(POOL, err_mode, label, n_seeds=30, n_try=2500):
    nN = len(POOL)
    ERR  = np.array([d["err"] for d in POOL])
    if err_mode == "const":
        ERR = np.full(nN, 0.3)
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

print("T4 敏感性: 闸门尺度 ERR=AME(MEerr) vs ERR=0.3 MeV 固定")
print("参数搜索范围冻结 (RANGES 同原协议), 30 seeds × 70/30")
print()
out = {}
jobs = [("AME能量侧", ame, "ame_err", "AME+ERR=AME"), ("AME能量侧", ame, "const", "AME+ERR=0.3"),
        ("闭式A_raw", make_ame_model_raw(), "ame_err", "闭式+ERR=AME"), ("闭式A_raw", make_ame_model_raw(), "const", "闭式+ERR=0.3")]
for tag, src, mode, label in jobs:
    print("="*84); print(f"{label}  ({tag}, err_mode={mode})"); print("="*84)
    POOL = build(POOL_KEYS, src)
    out[label] = run_protocol(POOL, mode, label)
    r = out[label]
    print(f"  → mean={r['mean']*100:.2f}%  std={r['std']*100:.2f}  min={r['min']*100:.2f}%  free={r['mean_free']*100:.2f}%")

print()
print("="*84)
print("T4 敏感性汇总")
print("="*84)
for k, v in out.items():
    print(f"  {k:20s} mean={v['mean']*100:6.2f}%  free={v['mean_free']*100:6.2f}%")
ame_d = abs(out["AME+ERR=AME"]["mean"] - out["AME+ERR=0.3"]["mean"])
raw_d = abs(out["闭式+ERR=AME"]["mean"] - out["闭式+ERR=0.3"]["mean"])
print(f"\n敏感性差: AME 侧 Δ={ame_d*100:.2f} pp   闭式侧 Δ={raw_d*100:.2f} pp")
verdict = ("结论对闸门尺度不敏感 (Δ<1 pp) → 加分项成立"
           if max(ame_d, raw_d) < 0.01 else
           "闸门尺度影响有限但不为零 → 需在正文声明两种口径" if max(ame_d, raw_d) < 0.03 else
           "结论对闸门尺度敏感 → 判定层表述必须依赖 ERR 来源声明")
print(f"裁决 → {verdict}")
json.dump({k: dict(v, verdict=verdict) for k, v in out.items()},
          open(os.path.join(HERE, "results/t4_err_sensitivity.json"), "w"), indent=1, ensure_ascii=False)
print("存档: results/t4_err_sensitivity.json")
