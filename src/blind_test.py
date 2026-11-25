# -*- coding: utf-8 -*-
"""
随机盲点测试 —— 检验"法则一·方案D"是学到了规律还是背住了答案

作弊嫌疑堵漏:
  A) 六个可调参数 (sig_k, NU_TRY, S_ALPHA, EC_F, FT_N, T_STABLE_Y) 原是在全数据集上
     调定的 → 盲测中每次随机分割后仅在训练集上重新拟合, 冻结后预测盲测集
  B) 滴线区补丁读取 NUBASE 实测 T<1μs (属于"看答案") → 盲测全程禁用,
     受补丁影响的核素单独统计 (patch-dependent)

协议:
  P1 个体盲测: 30 seeds × 随机 70/30 分割, 训练集拟合 → 盲测集预测
  P2 整元素封锁: 30 seeds × 随机封锁 30% 元素族 (同位素链整族隐藏, 防近邻泄题)
  P3 外推盲测: 仅用 A≤20 轻核拟合 → 预测 A 21–40 全部 + A>40 重核 22 个
  P4 置换对照: 打乱训练集答案后拟合 → 盲测集真实标签上评分 (应跌至随机水平)
"""
import os, sys, json
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from nucdata import parse_ame, parse_nubase
from beta_phase_space import f_value
from chain_predict import q_channels
from baseline_check import wkb_P

HERE = os.path.dirname(os.path.abspath(__file__))
YR = 3.1557e7
LN2 = np.log(2)

ame = parse_ame(os.path.join(HERE, "data/mass_1.mas20"))
nub = parse_nubase(os.path.join(HERE, "data/nubase_4.mas20"))
ME_H = ame[(1, 0)]["ME"]
ame[(2, 0)] = dict(A=2, ME=2*ME_H + 91.8, MEerr=0.5, BEA=0.0, Qb=None, flag=0, aux=[], El="He")

HEAVY = sorted(set([(26,30),(38,48),(47,60),(79,118),(26,33),(38,52),(55,82),(82,126),
                    (84,126),(88,138),(90,142),(90,144),(92,143),(92,146),(20,21),
                    (27,33),(28,31),(43,55),(53,76),(62,82),(44,60),(62,90)]))

# ---------------- 数据集构建（通道基速率与参数无关，只算一次） ----------------
def build(keys):
    """每核素: 通道 [n, p, alpha, bm, bp]; 返回 dict 数组"""
    out = []
    for k in keys:
        b = nub[k]
        e = "stable" if b["stable"] else b["mode"]
        if e in ("other", None):
            continue
        qc = q_channels(ame, k)
        Z, N = k; A = Z + N
        err = ame[k]["MEerr"] / 1000.0  # MeV
        Qn, Q2n = qc["Qn"], qc["Q2n"]
        Qp, Q2p = qc["Qp"], qc["Q2p"]
        Qa = qc["Qa"]
        # p 通道两候选: 经 1p (Qe=Qp) 或经 2p (Qe=Q2p/2)
        Qe1 = Qp
        Qe2 = Q2p / 2 if Q2p is not None else None
        w_p1 = wkb_P(Qe1, 1, Z-1, A) if Qe1 is not None else np.nan
        w_p2 = wkb_P(Qe2, 1, Z-1, A) if Qe2 is not None else np.nan
        w_a = wkb_P(Qa, 2, Z-2, A) if Qa is not None else np.nan
        f_bm = f_bp = np.nan
        if qc["Qbm"] is not None and qc["Qbm"] > 0:
            f_bm = float(f_value(qc["Qbm"], Z+1))
        if qc["Qbp"] is not None and qc["Qbp"] > 0:
            f_bp = float(f_value(qc["Qbp"], -(Z-1)))
        elif qc["Qec"] is not None and qc["Qec"] > 0:
            f_bp = float(f_value(qc["Qec"], Z-1))  # EC: 与 bp 同槽, 乘 EC_F
        out.append(dict(key=k, A=A, Z=Z, exp=e, err=err,
                        Qn=Qn, Q2n=Q2n, Qp=Qp, Q2p=Q2p, Qa=Qa,
                        w_p1=w_p1, w_p2=w_p2, w_a=w_a, f_bm=f_bm, f_bp=f_bp))
    return out

POOL = build(sorted([k for k in nub if k in ame and nub[k]["A"] <= 40],
                    key=lambda k: nub[k]["A"]))
HEAVY_D = build([k for k in HEAVY if k in ame and k in nub])
print(f"数据集: A≤40 共 {len(POOL)} 核素; 重核抽检 {len(HEAVY_D)} 个")

# 向量化数组
nN = len(POOL)
ERR  = np.array([d["err"] for d in POOL])
EXP  = np.array([d["exp"] for d in POOL])
Qn   = np.array([np.nan if d["Qn"] is None else d["Qn"] for d in POOL])
Q2n  = np.array([np.nan if d["Q2n"] is None else d["Q2n"] for d in POOL])
Qp   = np.array([np.nan if d["Qp"] is None else d["Qp"] for d in POOL])
Q2p  = np.array([np.nan if d["Q2p"] is None else d["Q2p"] for d in POOL])
Qa   = np.array([np.nan if d["Qa"] is None else d["Qa"] for d in POOL])
Wp1  = np.array([d["w_p1"] for d in POOL])
Wp2  = np.array([d["w_p2"] for d in POOL])
Wa   = np.array([d["w_a"] for d in POOL])
Fbm  = np.array([d["f_bm"] for d in POOL])
Fbp  = np.array([d["f_bp"] for d in POOL])
ISBETA_BM = ~np.isnan(Fbm)
ISBETA_BP = ~np.isnan(Fbp)

def predict(theta, idx=None):
    """theta = (sig_k, NU_TRY, S_ALPHA, EC_F, FT_N, T_STABLE); 返回预测模式数组"""
    sig_k, nu, sa, ec, ft, ts = theta
    ii = np.arange(nN) if idx is None else np.asarray(idx)
    thr = sig_k * ERR[ii]
    open_n = np.nanmax(np.c_[Qn[ii], Q2n[ii]], axis=1) > thr        # n: Qn 或 Q2n 显著
    open_a = Qa[ii] > thr
    cond_p1 = Qp[ii] > thr
    cond_p2 = Q2p[ii] > thr
    open_p = cond_p1 | cond_p2
    w_p = np.where(cond_p1, Wp1[ii], Wp2[ii])                       # Qp 优先, 否则 Q2p/2
    # 速率 (s^-1)
    r_n  = np.where(open_n, nu, -np.inf)
    r_p  = np.where(open_p, nu * w_p, -np.inf)
    r_a  = np.where(open_a, nu * sa * Wa[ii], -np.inf)
    bsc  = LN2 / ft
    r_bm = np.where(ISBETA_BM[ii], bsc * Fbm[ii], -np.inf)
    r_bp = np.where(ISBETA_BP[ii], ec * bsc * Fbp[ii], -np.inf)
    R = np.c_[r_n, r_p, r_a, r_bm, r_bp]
    modes = np.array(["n", "p", "alpha", "bm", "bp"])
    j = np.argmax(R, axis=1)
    lam = R[np.arange(len(ii)), j]
    T = np.where(lam > 0, LN2 / lam / YR, np.inf)
    pred = np.where(T > ts, "stable", modes[j])
    pred[np.isnan(R).all(axis=1) | (lam == -np.inf)] = "stable"
    return pred

def acc(pred, idx, exp_labels=None):
    e = EXP[idx] if exp_labels is None else exp_labels
    return float(np.mean(pred == e))

# ---------------- 参数拟合（仅训练集, 随机搜索） ----------------
RANGES = [(0.0, 1.3), (18, 24), (-4, 0), (-1, 0), (2.5, 4.5), (8, 25)]  # log10

def fit(idx, n_try=2500, rng=None, exp_labels=None):
    rng = rng or np.random.default_rng()
    best = (-1, None)
    ii = np.asarray(idx)
    e = EXP[ii] if exp_labels is None else np.asarray(exp_labels)
    for _ in range(n_try):
        th = tuple(10**rng.uniform(*r) for r in RANGES)
        a = acc(predict(th, ii), ii, EXP[ii] if exp_labels is None else e)
        if a > best[0]:
            best = (a, th)
    return best  # (train_acc, theta)

# ---------------- 滴线补丁依赖核素（诚实条款统计） ----------------
THETA0 = (2.0, 1e21, 0.01, 0.5, 10**3.5, 1e17)
pred_nopatch0 = predict(THETA0)
# 补丁版预测（读 NUBASE T —— 仅用于标记依赖核素, 不参与盲测）
import re
lt_T = {}
with open(os.path.join(HERE, "data/nubase_4.mas20")) as f:
    for line in f:
        s = line.rstrip()
        m = re.search(r"<\s*([0-9.]+)\s*(ys|zs|ps|ns|us|ms|s)\b", s[:100])
        if m:
            m2 = re.search(r"([1-9]\d{0,2})([A-Z][a-z]?)\b", s[6:26])
            if m2:
                lt_T[(int(m2.group(1)), m2.group(2))] = float(m.group(1)) * \
                    {"ys": 1e-24, "zs": 1e-21, "ps": 1e-12, "ns": 1e-9,
                     "us": 1e-6, "ms": 1e-3, "s": 1}[m.group(2)]
key_of = {(v["A"], v["El"]): k for k, v in ame.items()}
# 正确口径: 带补丁(原模型)判对 且 无补丁判错 → 补丁依赖核素
_base = json.load(open(os.path.join(HERE, "results/baseline_judge.json")))
PATCH_KEYS = []
for i, d in enumerate(POOL):
    k = d["key"]
    rec = _base.get(f"{k[0]},{k[1]}")
    if rec and rec["pred"] == rec["exp"] and pred_nopatch0[i] != d["exp"]:
        PATCH_KEYS.append(k)
PATCH_IDX = np.array([i for i, d in enumerate(POOL) if d["key"] in set(PATCH_KEYS)])
FREE_IDX = np.setdiff1d(np.arange(nN), PATCH_IDX)
patch_names = ["%s-%d" % (ame[k]["El"], ame[k]["A"]) for k in PATCH_KEYS]
print(f"滴线补丁依赖核素: {len(PATCH_KEYS)} 个 {patch_names}")

RESULTS = {"patch_dependent": [f"{ame[k]['El']}-{ame[k]['A']}" for k in PATCH_KEYS]}

# ---------------- P1: 个体盲测 ----------------
print("\n=== P1 个体盲测 (30 seeds, 70/30 随机) ===")
p1 = []
for s in range(30):
    rng = np.random.default_rng(1000 + s)
    perm = rng.permutation(nN)
    ntr = int(0.7 * nN)
    tr, bl = perm[:ntr], perm[ntr:]
    a_tr, th = fit(tr, rng=rng)
    a_bl = acc(predict(th, bl), bl)
    a_bl_free = acc(predict(th, np.intersect1d(bl, FREE_IDX)), np.intersect1d(bl, FREE_IDX))
    p1.append(dict(seed=s, train=a_tr, blind=a_bl, blind_free=a_bl_free, theta=th))
accs = np.array([x["blind"] for x in p1])
accf = np.array([x["blind_free"] for x in p1])
print(f"盲测准确率: mean={accs.mean()*100:.1f}%  std={accs.std()*100:.1f}%  "
      f"min={accs.min()*100:.1f}%  max={accs.max()*100:.1f}%")
print(f"剔除补丁依赖核素后: mean={accf.mean()*100:.1f}%  min={accf.min()*100:.1f}%")
th_med = [float(np.median([x["theta"][j] for x in p1])) for j in range(6)]
print("参数中位数: sig_k=%.2f NU_TRY=%.2e S_ALPHA=%.2e EC_F=%.2f FT_N=%.2f T_S=%.1e" %
      (th_med[0], th_med[1], th_med[2], th_med[3], th_med[4], th_med[5]))
RESULTS["P1"] = dict(n=30, mean=float(accs.mean()), std=float(accs.std()),
                     min=float(accs.min()), max=float(accs.max()),
                     mean_free=float(accf.mean()), min_free=float(accf.min()),
                     theta_median=th_med, runs=p1)

# ---------------- P2: 整元素封锁 ----------------
print("\n=== P2 整元素封锁盲测 (30 seeds, 封锁 30% 元素族) ===")
Zs = np.array([d["Z"] for d in POOL])
uniqZ = sorted(set(Zs.tolist()))
p2 = []
for s in range(30):
    rng = np.random.default_rng(2000 + s)
    hid = set(rng.choice(uniqZ, size=max(1, int(0.3 * len(uniqZ))), replace=False).tolist())
    bl = np.where(np.isin(Zs, list(hid)))[0]
    tr = np.where(~np.isin(Zs, list(hid)))[0]
    a_tr, th = fit(tr, rng=rng)
    a_bl = acc(predict(th, bl), bl)
    bl_f = np.intersect1d(bl, FREE_IDX)
    a_bl_free = acc(predict(th, bl_f), bl_f) if len(bl_f) else float("nan")
    p2.append(dict(seed=s, hidden=sorted(hid), train=a_tr, blind=a_bl, blind_free=a_bl_free, theta=th))
accs2 = np.array([x["blind"] for x in p2])
accf2 = np.array([x["blind_free"] for x in p2])
print(f"盲测准确率: mean={accs2.mean()*100:.1f}%  std={accs2.std()*100:.1f}%  "
      f"min={accs2.min()*100:.1f}%  max={accs2.max()*100:.1f}%")
print(f"剔除补丁依赖核素后: mean={np.nanmean(accf2)*100:.1f}%  min={np.nanmin(accf2)*100:.1f}%")
RESULTS["P2"] = dict(n=30, mean=float(accs2.mean()), std=float(accs2.std()),
                     min=float(accs2.min()), max=float(accs2.max()),
                     mean_free=float(np.nanmean(accf2)), min_free=float(np.nanmin(accf2)),
                     runs=p2)

# ---------------- P3: 外推盲测 ----------------
print("\n=== P3 外推盲测: 仅 A≤20 标定 → 预测 A 21–40 + 重核 ===")
light = np.array([i for i, d in enumerate(POOL) if d["A"] <= 20])
heavy_in_pool = np.array([i for i, d in enumerate(POOL) if d["A"] > 20])
p3 = []
for s in range(10):
    rng = np.random.default_rng(3000 + s)
    a_tr, th = fit(light, rng=rng)
    a_mid = acc(predict(th, heavy_in_pool), heavy_in_pool)
    p3.append(dict(seed=s, train_light=a_tr, blind_mid=a_mid, theta=th))
accs3 = np.array([x["blind_mid"] for x in p3])
print(f"A 21–40 盲测准确率: mean={accs3.mean()*100:.1f}%  min={accs3.min()*100:.1f}%")
th_med3 = [float(np.median([x["theta"][j] for x in p3])) for j in range(6)]
print("轻核拟合参数中位数:", [f"{v:.3g}" for v in th_med3])
RESULTS["P3"] = dict(n=10, mean_mid=float(accs3.mean()), min_mid=float(accs3.min()),
                     theta_median=th_med3, runs=p3)

# 重核逐个预测（P3 最优 seed 的参数中位数）
thH = tuple(th_med3)
print(f"\n重核抽检 (参数={['%.3g'%v for v in thH]}):")
hok = 0
ok = 0
hdet = []
for d in HEAVY_D:
    k = d["key"]; Z, N = k; A = Z + N
    err = d["err"]; thr = thH[0] * err
    parts = []
    if (not np.isnan(d["Qn"]) and d["Qn"] > thr) or (not np.isnan(d["Q2n"]) and d["Q2n"] > thr):
        parts.append(("n", thH[1]))
    if (not np.isnan(d["Qp"]) and d["Qp"] > thr) or (not np.isnan(d["Q2p"]) and d["Q2p"] > thr):
        Qe = d["Qp"] if not np.isnan(d["w_p1"]) else d["Q2p"]/2
        parts.append(("p", thH[1] * (d["w_p1"] if not np.isnan(d["w_p1"]) else d["w_p2"])))
    if not np.isnan(d["Qa"]) and d["Qa"] > thr:
        parts.append(("alpha", thH[1] * thH[2] * d["w_a"]))
    if not np.isnan(d["f_bm"]):
        parts.append(("bm", LN2/thH[4] * d["f_bm"]))
    if not np.isnan(d["f_bp"]):
        parts.append(("bp", thH[3]*LN2/thH[4] * d["f_bp"]))
    if not parts:
        p = "stable"
    else:
        ch, lam = max(parts, key=lambda x: x[1])
        T = LN2/lam/YR if lam > 0 else np.inf
        p = "stable" if T > thH[5] else ch
    ok += p == d["exp"]
    hdet.append((f"{ame[k]['El']}-{A}", d["exp"], p))
    print(f"   {'✓' if p == d['exp'] else '✗'} {ame[k]['El']}-{A}: 实验={d['exp']}  盲测模型={p}")
print(f"重核外推盲测: {ok}/{len(HEAVY_D)}")
RESULTS["P3_heavy"] = dict(ok=ok, n=len(HEAVY_D), detail=hdet)

# ---------------- P4: 置换对照 ----------------
print("\n=== P4 置换对照 (10 seeds, 打乱训练答案) ===")
p4 = []
for s in range(10):
    rng = np.random.default_rng(4000 + s)
    perm = rng.permutation(nN)
    ntr = int(0.7 * nN)
    tr, bl = perm[:ntr], perm[ntr:]
    shuffled = rng.permutation(EXP[tr])
    a_tr_sh, th = fit(tr, rng=rng, exp_labels=shuffled)
    a_bl = acc(predict(th, bl), bl)  # 真实标签
    p4.append(dict(seed=s, train_shuffled=a_tr_sh, blind_true=a_bl))
accs4 = np.array([x["blind_true"] for x in p4])
trs4 = np.array([x["train_shuffled"] for x in p4])
print(f"打乱答案后的训练拟合精度: mean={trs4.mean()*100:.1f}% (最高类占比 {max(np.unique(EXP, return_counts=True)[1])/nN*100:.1f}%)")
print(f"打乱答案后模型在盲测集(真实标签): mean={accs4.mean()*100:.1f}%  max={accs4.max()*100:.1f}%")
RESULTS["P4"] = dict(n=10, mean_train_shuffled=float(trs4.mean()),
                     mean_blind=float(accs4.mean()), max_blind=float(accs4.max()), runs=p4)

# ---------------- 全数据重标定对照（说明参数敏感度） ----------------
a_all, th_all = fit(np.arange(nN), n_try=4000, rng=np.random.default_rng(7))
print(f"\n参考: 全数据重标定(不含补丁) 准确率 = {a_all*100:.1f}%  参数={['%.3g'%v for v in th_all]}")
RESULTS["full_refit"] = dict(acc=float(a_all), theta=list(th_all))

with open(os.path.join(HERE, "results/blind_results.json"), "w") as f:
    json.dump(RESULTS, f, indent=1, ensure_ascii=False, default=str)
print("\n已存 results/blind_results.json")
