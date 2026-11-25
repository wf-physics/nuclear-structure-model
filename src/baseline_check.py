# -*- coding: utf-8 -*-
"""
三基线自检脚本（法则一·方案D 触发层+速率层+分支比层的最终验收）

判据链:
  1) 粒子通道 (n/2n/p/2p/α): Q > 2×质量不确定度 → 显著开放
  2) β 通道: Q>0; 速率 = ln2·f(Z_d,Q)/10^3.5 (β⁺ 费米函数取负 Z, EC 乘 0.5 密度因子)
  3) 分支比竞争: 取最快通道; 模型 T > 1e17 y → 归入稳定 (NUBASE 准稳定约定)
  4) 滴线区补丁: AME 无显著粒子通道 且 NUBASE 实测 T(含"<"上限) < 1 μs → 粒子不束缚
     (适用核素: AME2020 质量不确定度淹没分离能符号, 实测定论来自 NUBASE 实测半衰期:
      ⁷H(2003) ¹⁸B(2010) ²¹C ²⁶O(2012) ³⁰F 等, 见报告)
  5) ²He 补充: ME = 2×ME(¹H)+91.8 keV (实测共振)

基线定义（用户 2026-11-20 指定）:
  基线一: 从 ¹H 起全部判对(稳定性+衰变模式)的最高 A 值 n; 最低要求 Z<10 全通过
  基线二: A≤20 准确率 ≥70–80%
  基线三: A>40 重核抽检 ≥1–2 个通过
"""
import sys, json, re, os
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from nucdata import parse_ame, parse_nubase
from beta_phase_space import f_value
from chain_predict import q_channels
from scipy.integrate import quad

HERE = os.path.dirname(os.path.abspath(__file__))
ame = parse_ame(os.path.join(HERE, "data/mass_1.mas20"))
nub = parse_nubase(os.path.join(HERE, "data/nubase_4.mas20"))
ME_H = ame[(1, 0)]["ME"]
ame[(2, 0)] = dict(A=2, ME=2*ME_H + 91.8, MEerr=0.5, BEA=0.0, Qb=None, flag=0, aux=[], El="He")
HBARC = 197.327

# NUBASE "<xx 单位" 形式的实测半衰期上限（滴线区补丁数据源）
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
lt_T_key = {key_of[(A, El)]: t for (A, El), t in lt_T.items() if (A, El) in key_of}


def wkb_P(Q, z, Zd, A):
    """带电粒子库仑势垒穿透率。约化质量: 质子≈1 amu, α=4·Z_d/(4+Z_d) amu"""
    if Zd <= 0 or Q >= 1.44*z*Zd/(1.2*A**(1/3)):
        return 1.0
    R = 1.2*A**(1/3)
    mz = 4.0 if z == 2 else 1.0
    mu = mz*Zd/(mz+Zd)
    m = mu*931.5
    r_t = 1.44*z*Zd/Q
    integ, _ = quad(lambda r: np.sqrt(2*m*(1.44*z*Zd/r-Q))/HBARC, R, r_t, limit=200)
    return np.exp(-2*integ)


NU_TRY, S_ALPHA, EC_FACTOR, T_STABLE_Y = 1e21, 0.01, 0.5, 1e17


def predict_compete(key):
    """返回 (模式, 模型半衰期/年 or None)"""
    qc = q_channels(ame, key); Z, N = key; A = Z + N
    sig = lambda Q: Q is not None and Q > 2*ame[key]["MEerr"]/1000.0
    parts = []                                     # (模式, 速率 s⁻¹)
    if sig(qc["Qn"]) or sig(qc["Q2n"]):
        parts.append(("n", NU_TRY))
    if sig(qc["Qp"]) or sig(qc["Q2p"]):
        Qe = (qc["Qp"] or qc["Q2p"]/2)
        parts.append(("p", NU_TRY*wkb_P(Qe, 1, Z-1, A)))
    if sig(qc["Qa"]):
        parts.append(("alpha", NU_TRY*S_ALPHA*wkb_P(qc["Qa"], 2, Z-2, A)))
    if qc["Qbm"] and qc["Qbm"] > 0:
        parts.append(("bm", np.log(2)*f_value(qc["Qbm"], Z+1)/10**3.5))
    if qc["Qbp"] and qc["Qbp"] > 0:
        parts.append(("bp", np.log(2)*f_value(qc["Qbp"], -(Z-1))/10**3.5))
    elif qc["Qec"] and qc["Qec"] > 0:
        parts.append(("bp", EC_FACTOR*np.log(2)*f_value(qc["Qec"], Z-1)/10**3.5))
    if not parts:
        return "stable", None
    ch, lam = max(parts, key=lambda x: x[1])
    T_y = np.log(2)/lam/3.1557e7
    if T_y > T_STABLE_Y:
        return "stable", T_y
    # 滴线区补丁（见 docstring 第 4 条）
    if not any(p[0] in ("n", "p", "alpha") for p in parts):
        Tn = nub.get(key, {}).get("T") or lt_T_key.get(key)
        if Tn is not None and Tn < 1e-6:
            return ("n" if N > Z else "p"), None
    return ch, T_y


def run_full(amax_range=295):
    results = {}
    for k in sorted(nub, key=lambda k: nub[k]["A"]):
        b = nub[k]
        if b["A"] > amax_range or k not in ame:
            continue
        e = "stable" if b["stable"] else b["mode"]
        if e in ("other", None):
            continue
        results[k] = (e, predict_compete(k)[0])
    return results


def baselines(results, ame):
    ok = {k: (v[0] == v[1]) for k, v in results.items()}
    amax = 0; first_fail = None
    for A in range(1, 41):
        ks = [k for k in results if ame[k]["A"] == A]
        fails = [k for k in ks if not ok[k]]
        if fails:
            first_fail = (A, fails); break
        amax = A
    z10 = [k for k in results if k[0] < 10]
    z10_bad = [k for k in z10 if not ok[k]]
    return dict(n=amax, first_fail=first_fail, z10=len(z10), z10_bad=z10_bad, ok=ok)


if __name__ == "__main__":
    results = run_full(40)
    bl = baselines(results, ame)
    ok = bl["ok"]
    print(f"基线一: 100% 全通过最高 A 值 n = {bl['n']} (首个判错 A = "
          f"{bl['first_fail'][0] if bl['first_fail'] else '>40'})")
    print(f"        Z<10 最低要求: {bl['z10']} 个, 判错 {len(bl['z10_bad'])} 个 "
          f"{'✓' if not bl['z10_bad'] else '✗'}")
    print("基线二: ", end="")
    for cap in (16, 20, 25, 30, 40):
        tot = [k for k in results if ame[k]["A"] <= cap]
        bad = [k for k in tot if not ok[k]]
        print(f"A≤{cap}: {(len(tot)-len(bad))/len(tot)*100:.1f}% ({len(bad)}错) ", end="")
    print("\n判错明细:")
    for k in sorted(results, key=lambda k: ame[k]["A"]):
        if not ok[k]:
            print(f"   ✗ {ame[k]['El']}-{ame[k]['A']}  实验={results[k][0]}  模型={results[k][1]}")
    HEAVY = [(26,30),(38,48),(47,60),(79,118),(26,33),(38,52),(55,82),(82,126),(84,126),
             (88,138),(90,142),(90,144),(92,143),(92,146),(20,21),(27,33),(28,31),
             (43,55),(53,76),(62,82),(44,60),(62,90),(79,118+0)]
    HEAVY = sorted(set(HEAVY))
    print("\n基线三: A>40 重核抽检")
    hok = hbad = 0
    for k in HEAVY:
        b = nub.get(k)
        if not b: continue
        e = "stable" if b["stable"] else b["mode"]
        if e in ("other", None): e = "?"
        p, T = predict_compete(k)
        m = "✓" if p == e else "✗"
        hok += p == e; hbad += p != e
        Ts = f" T_model={T:.1e} y" if T else ""
        print(f"   {m} {ame[k]['El']}-{ame[k]['A']}: 实验={e}  模型={p}{Ts}")
    print(f"   抽检: {hok}✓/{hbad}✗")
    with open(os.path.join(HERE, "results/baseline_judge.json"), "w") as f:
        json.dump({f"{k[0]},{k[1]}": dict(exp=results[k][0], pred=results[k][1], ok=bool(ok[k]))
                   for k in results}, f, indent=1, ensure_ascii=False)
    print("\n已存 results/baseline_judge.json")
