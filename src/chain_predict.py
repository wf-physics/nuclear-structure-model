# -*- coding: utf-8 -*-
"""
反应链推演 + 与 NUBASE2020 实验事实的系统比对（A≤40，H 到 Ca）

预测规则（忠实于 v4.0 并显式声明必要的补充）:
  1) 粒子通道（n/p/α 发射）: Q 来自 AME2020 外部原子质量（v4.0 §4.1 自己的规定：
     "必须另用外部原子质量计算反应 Q 值"）。Q>0 即通道开放。
     v4.0 §5 方案 A/B 的粒子发射闸门参数未定 → 模型无法关闭开放通道（声明局限）。
  2) β⁻ 闸门（v4.0 法则二·下限，逐字实现）: 存在中子使 W_β(n) < Q_β⁻ → 预测 β⁻。
     W_β 含自项 1/a（a=0.8/1.0 两组）。
  3) β⁺/EC 闸门（v4.0 法则二·上限，逐字实现）: ΔE_self 按其自身坐标约定
     （新质子与电子同位于 r_n）恒等于 0 —— 这是 v4.0 的内部不一致（纲领中的
     δr_dipole 在 v4.0 丢失）。阈值取 Q 尺度（§2.6）: E_C_max ≥ Q_EC 时预测 β⁺/EC。
  4) 其余 → 预测稳定。
基线（对照）: 纯能量基线 —— 粒子通道相同 + β 通道只看 Q 符号，无模型闸门。
  模型与基线的差异 = 模型闸门的净贡献。
反应链: 从 ¹H 与 ¹n 出发，按预测衰变通道走边（bm:(Z+1,N), bp:(Z-1,N),
  n:(Z,N-1), p:(Z-1,N), α:(Z-2,N-2)），预测稳定者为链终点。
比对: NUBASE2020 半衰期 + 主衰变模式。
"""
import numpy as np
import json, os, sys
from collections import deque

sys.path.insert(0, os.path.dirname(__file__))
from nucdata import parse_ame, parse_nubase

HERE = os.path.dirname(__file__)
A_MAX = 40
ELECTRON = 1022.0  # keV, 2·m_e·c²

def q_channels(ame, key):
    """计算各通道 Q 值 (MeV)。子核不在 AME 中 → 通道关闭(None)。
    β 衰变为 n↔p 转换（A 不变）: β⁻ → (Z+1,N-1)，β⁺/EC → (Z-1,N+1)。
    """
    Z, N = key
    ME = ame[key]["ME"] / 1000.0
    def me(k):
        return ame[k]["ME"] / 1000.0 if k in ame else None
    me_n, me_h, me_a = me((0, 1)), me((1, 0)), me((2, 2))
    out = {}
    d = me((Z, N - 1));   out["Qn"] = ME - d - me_n if (d is not None and me_n is not None and N >= 1) else None
    d = me((Z - 1, N));   out["Qp"] = ME - d - me_h if (d is not None and me_h is not None and Z >= 1) else None
    d = me((Z - 2, N - 2)); out["Qa"] = ME - d - me_a if (d is not None and me_a is not None and Z >= 2 and N >= 2) else None
    me_n2, me_h2 = 2 * me_n, 2 * me_h
    d = me((Z, N - 2));   out["Q2n"] = ME - d - me_n2 if (d is not None and me_n is not None and N >= 2) else None
    d = me((Z - 2, N));   out["Q2p"] = ME - d - me_h2 if (d is not None and me_h is not None and Z >= 2) else None
    d = me((Z + 1, N - 1)); out["Qbm"] = ME - d if (d is not None and N >= 1) else None
    d = me((Z - 1, N + 1)); out["Qec"] = ME - d if d is not None else None
    out["Qbp"] = out["Qec"] - ELECTRON / 1000.0 if out["Qec"] is not None else None
    return out

def predict(ame, st, key, a="0.8"):
    """返回 (status, detail)。status ∈ {stable, bm, bp, alpha, p, n}"""
    qc = q_channels(ame, key)
    Z, N = key
    # 1) 粒子通道（能量开放即预测发射；v4.0 §5 方案 A/B 未定参 → 模型无闸门）
    for ch in ("Qn", "Qp", "Qa", "Q2n", "Q2p"):
        if qc.get(ch) is not None and qc[ch] > 0:
            return {"Qn": "n", "Qp": "p", "Qa": "alpha",
                    "Q2n": "n", "Q2p": "p"}[ch], dict(qc=qc, gate="能量开放(模型无闸门)")
    # 2) β⁻ 闸门（v4.0 逐字）
    Wb = st.get("Wb_min") if st else None
    if qc["Qbm"] is not None and qc["Qbm"] > 0 and Wb:
        w = Wb.get(a)
        if w is not None and w < qc["Qbm"]:
            return "bm", dict(qc=qc, gate=f"Wb_min={w:.3f} < Qbm={qc['Qbm']:.3f}")
    # 3) β⁺/EC 闸门（v4.0 逐字 + ΔE_self≡0）
    if qc["Qec"] is not None and qc["Qec"] > 0 and st is not None:
        if st.get("EC_max", 0.0) >= qc["Qec"]:
            return "bp", dict(qc=qc, gate=f"EC_max={st['EC_max']:.3f} ≥ Qec={qc['Qec']:.3f}")
    # 4) 稳定
    return "stable", dict(qc=qc, gate="无开放通道通过闸门")

def baseline(qc):
    """纯能量基线: 粒子通道 → β⁻(Q>0) → β⁺/EC(Q>0) → 稳定"""
    for ch in ("Qn", "Qp", "Qa", "Q2n", "Q2p"):
        if qc.get(ch) is not None and qc[ch] > 0:
            return {"Qn": "n", "Qp": "p", "Qa": "alpha", "Q2n": "n", "Q2p": "p"}[ch]
    if qc["Qbm"] is not None and qc["Qbm"] > 0: return "bm"
    if qc["Qec"] is not None and qc["Qec"] > 0: return "bp"
    return "stable"

def main():
    ame = parse_ame(os.path.join(HERE, "data", "mass_1.mas20"))
    nub = parse_nubase(os.path.join(HERE, "data", "nubase_4.mas20"))
    with open(os.path.join(HERE, "results", "structures.json")) as f:
        st_all = json.load(f)
    keys = sorted([k for k in ame if ame[k]["A"] <= A_MAX],
                  key=lambda k: (ame[k]["A"], k[0]))
    # ---- 逐核素预测 ----
    preds, base = {}, {}
    for k in keys:
        st = st_all.get(f"{k[0]},{k[1]}")
        p, det = predict(ame, st, k)
        preds[k] = p; base[k] = baseline(q_channels(ame, k))
    # ---- 与 NUBASE 比对 ----
    exp_mode, exp_stable = {}, set()
    for k in keys:
        b = nub.get(k)
        if b is None:
            continue
        exp_mode[k] = b["mode"] if b["mode"] else "other"
        if b["stable"]:
            exp_stable.add(k)
    labels = ["stable", "bm", "bp", "alpha", "p", "n"]
    def acc(pred_map):
        n_ok_st = n_wrong_st = n_ok_un = n_wrong_un = 0
        mode_ok = mode_tot = 0
        mismatch = []
        for k in keys:
            e = exp_mode.get(k)
            if e is None: continue
            p = pred_map[k]
            e_st = (e == "stable")
            p_st = (p == "stable")
            if e_st and p_st: n_ok_st += 1
            elif e_st and not p_st: n_wrong_st += 1; mismatch.append((k, "漏判衰变", e, p))
            elif not e_st and p_st: n_wrong_un += 1; mismatch.append((k, "误判稳定", e, p))
            else:
                n_ok_un += 1
                mode_tot += 1
                if p == e:
                    mode_ok += 1
                else:
                    mismatch.append((k, "衰变模式不符", e, p))
        total = n_ok_st + n_wrong_st + n_ok_un + n_wrong_un
        return dict(acc_stable=1 - (n_wrong_st + n_wrong_un) / total if total else 0,
                    n_ok_stable=n_ok_st, n_exp_unstable_pred_stable=n_wrong_un,
                    n_exp_stable_pred_unstable=n_wrong_st,
                    mode_acc=mode_ok / mode_tot if mode_tot else 0,
                    mode_tot=mode_tot, mismatch=mismatch, total=total)
    m_model = acc(preds); m_base = acc(base)
    # ---- 反应链推演：生长（add_n/add_p，从 ¹H 出发）+ 衰变边 ----
    # 法则四的完整反应链 = 逐核子添加生成全部 (Z,N) 状态 + 每个状态按预测
    # 衰变通道走边。预测稳定者为链终点。
    endpoints, edges = [], []
    visited = set()
    dq = deque([(1, 0)])            # 种子: ¹H（另记 ¹n 的独立判定）
    dq.append((0, 1))               # 种子: ¹n（P0-1 缺陷的判定点）
    while dq:
        k = dq.popleft()
        if k in visited:
            continue
        visited.add(k)
        if k not in preds:          # 超出 AME/范围 → 链在此截断
            continue
        p = preds[k]
        if p != "stable":
            dz = {"bm": (1, -1), "bp": (-1, 1), "p": (-1, 0),
                  "n": (0, -1), "alpha": (-2, -2)}[p]
            k2 = (k[0] + dz[0], k[1] + dz[1])
            edges.append((k, p, k2))
            dq.append(k2)
        else:
            endpoints.append(k)     # 链终点：预测稳定
        dq.extend([(k[0] + 1, k[1]), (k[0], k[1] + 1)])   # 生长: add_p / add_n
    chain_stable = set(endpoints)
    exp_st40 = {k for k in exp_stable}
    hit = chain_stable & exp_st40
    spur = chain_stable - exp_st40
    miss = exp_st40 - chain_stable
    # ---- 输出 ----
    print("=" * 88)
    print(f"A≤40 反应链推演 vs NUBASE2020（核素数 = {len(keys)}，其中实验稳定 {len(exp_st40)} 个）")
    print("=" * 88)
    print(f"\n{'指标':<34}{'模型(v4.0判据)':>16}{'能量基线':>12}")
    print(f"稳定性判定准确率（含两端）     {m_model['acc_stable']*100:>15.1f}%{m_base['acc_stable']*100:>11.1f}%")
    print(f"  实验稳定→预测稳定            {m_model['n_ok_stable']:>15}{m_base['n_ok_stable']:>12}")
    print(f"  实验衰变→误判稳定            {m_model['n_exp_unstable_pred_stable']:>15}{m_base['n_exp_unstable_pred_stable']:>12}")
    print(f"  实验稳定→误判衰变            {m_model['n_exp_stable_pred_unstable']:>15}{m_base['n_exp_stable_pred_unstable']:>12}")
    print(f"衰变模式一致率（双衰变核素）   {m_model['mode_acc']*100:>15.1f}%{m_base['mode_acc']*100:>11.1f}%"
          f"  (n={m_model['mode_tot']})")
    print(f"\n--- 反应链（生长 add_p/add_n 自 ¹H，另设 ¹n 种子）---")
    print(f"访问状态数: {len(visited)}; 链终点（预测稳定）: {len(chain_stable)} 个")
    print(f"  与实验稳定重合: {len(hit)}")
    print(f"  虚假终点（预测稳定但实验衰变）: {len(spur)}")
    print(f"  漏掉终点（实验稳定但未达）: {len(miss)}")
    print(f"链边数（衰变边）: {len(edges)}")
    print(f"\n自由中子 ¹n 的预测: {preds.get((0,1), '—')}")
    print(f"\n--- 误判清单（前 25 条，按 A 排序）---")
    mm = sorted(m_model["mismatch"], key=lambda x: (ame[x[0]]["A"], x[0]))[:25]
    NAME = {k: f"{ame[k]['El']}-{ame[k]['A']}" for k in keys if k in ame}
    print(f"{'核素':<10}{'类型':<14}{'实验':<10}{'模型':<10}")
    for k, typ, e, p in mm:
        print(f"{NAME.get(k,k):<10}{typ:<14}{e:<10}{p:<10}")
    # 保存
    out = dict(model=dict(acc=m_model["acc_stable"], mode_acc=m_model["mode_acc"]),
               baseline=dict(acc=m_base["acc_stable"], mode_acc=m_base["mode_acc"]),
               chain=dict(endpoints=sorted(chain_stable), hit=sorted(hit),
                          spurious=sorted(spur), missed=sorted(miss),
                          n_edges=len(edges)),
               predictions={f"{k[0]},{k[1]}": preds[k] for k in keys},
               baseline_preds={f"{k[0]},{k[1]}": base[k] for k in keys},
               exp_modes={f"{k[0]},{k[1]}": exp_mode[k] for k in keys if k in exp_mode})
    with open(os.path.join(HERE, "results", "chain_results.json"), "w") as f:
        json.dump(out, f, ensure_ascii=False, indent=1)
    print("\n已保存 results/chain_results.json")

if __name__ == "__main__":
    main()
