# -*- coding: utf-8 -*-
"""
F2 收口 (方案2): 实验标定配对项 + Lane 项, 重测 Q_β⁻
用户批准边界: 公式层可自由修改, 但 (1) 逻辑/数学/物理必须严谨 (2) 不得越改越差
            (3) 不动 v5.0 法则与已成功的 F1/F3/F4 结果 (旧 JSON 一律不覆盖)

严谨性设计 (审稿人视角):
  R1. 配对形式 = 标准液滴配对项 Δ_pair = δ·A^(-α)·p,  p=+1(偶偶)/0(奇A)/-1(奇奇)
      α ∈ {1/2, 3/4} 两种经典形式都测, 由独立拟合集选定 (不偷看测试点)
  R2. δ 从 AME 独立样本拟合: 拟合集 = 稳定区 Z∈[20,82], A∈[40,240] 的同 A 偶偶/奇奇核,
      并显式排除全部 10 个 F2 测试核及其子核
  R3. δ 提取法: 同 A 的 B_ee − B_oo 实验差, 扣除模型液滴差分(含 Lane)后, 剩余即配对间隙
      δ_A = [ΔB_exp − ΔB_model] · A^α / 2,  对拟合集最小二乘
  R4. 结果矩阵: v4 基座 × {a_sym 19.5, 23} × {无配对, +配对}, 判定标准不变
      → 若任何新版本 MAE/判定不如旧版, 该版本标记"不采纳"
  R5. 身份标注: 半经验补丁 (δ 来自实验), 适用边界单列
"""
import numpy as np, json, sys
sys.path.insert(0, '.')
from nucdata import parse_ame
AME = parse_ame('data/mass_1.mas20')

AV, AS, AC = 13.850813715082808, 1.5924382285401872, 0.7264962440025764  # v4 模型自有闭式 (不动)

def Bf(Z, N, a_sym=19.50):
    A = Z + N
    return AV*A - AS*A**(2/3) - AC*Z**2/A**(1/3) - a_sym*(N-Z)**2/A

def Bexp(Z, N):
    v = AME.get((Z, N))
    return v['BEA']*v['A']/1000.0 if v and v.get('BEA') else None

def parity(Z, N):
    A = Z + N
    if A % 2 == 1:
        return 0.0
    return 1.0 if Z % 2 == 0 else -1.0

def judge(dev, Qe):
    return ('过' if abs(dev) <= max(0.4, 0.2*abs(Qe)) else
            '半' if abs(dev) <= max(1.0, 0.5*abs(Qe)) else '不过')

TESTS = [(6,8,'C-14'),(11,13,'Na-24'),(15,17,'P-32'),(16,19,'S-35'),
         (19,23,'K-42'),(25,31,'Mn-56'),(26,33,'Fe-59'),
         (38,52,'Sr-90'),(55,82,'Cs-137'),(61,86,'Pm-147')]
TEST_KEYS = set()
for Z, N, _ in TESTS:
    TEST_KEYS |= {(Z, N), (Z+1, N-1)}

# ============ R2/R3: δ 独立拟合 ============
print("="*84)
print("第一步: 配对系数 δ 独立拟合 (AME 稳定区, 排除全部 F2 测试核)")
print("="*84)
samples = {0.5: [], 0.75: []}
for A in range(40, 241, 2):
    for NZ in range(0, 15, 2):          # 偶 A 偶 NZ → ee(Z 偶) 与 oo(Z 奇) 对
        Zee = (A-NZ)//2
        if Zee < 20 or Zee > 80 or A-Zee > 160:
            continue
        Zoo = Zee + 1                   # 奇奇: Z+1, N-1, 同 A 同 NZ-2? 
        # ee: (Zee, A-Zee); oo: (Zee+1, A-Zee-1)
        k_ee, k_oo = (Zee, A-Zee), (Zee+1, A-Zee-1)
        if k_oo in TEST_KEYS or k_ee in TEST_KEYS:
            continue
        Bee, Boo = Bexp(*k_ee), Bexp(*k_oo)
        if not (Bee and Boo):
            continue
        # 模型液滴差分 (同 A, NZ 相差 2): 含 Lane 的模型液滴对 ee−oo 的能量差
        a_sym = 23.0
        NZee, NZoo = (A-2*Zee), (A-2*Zoo)
        dmodel = (Bf(*k_ee, a_sym) - Bf(*k_oo, a_sym))
        dexp = Bee - Boo
        # δ_A = (dexp − dmodel_非配对部分)·A^α/2 ; 模型差分里无配对项, 全额归配对+残差
        for al in (0.5, 0.75):
            samples[al].append((dexp - dmodel)*A**al/2.0)

delta_fit = {}
for al, s in samples.items():
    s = np.array(s)
    delta_fit[al] = float(np.median(s))
    print(f"  α={al}: n={len(s)}  δ 中位={delta_fit[al]:.3f} MeV  IQR=[{np.percentile(s,25):.2f},{np.percentile(s,75):.2f}]")

# 用拟合集自身 RMSE 选 α (不偷看测试点)
for al in (0.5, 0.75):
    s = np.array(samples[al])
    # 一致性检验: 用 δ 拟合值重构每个样本的配对间隙, 残差 RMS
    resid = s - delta_fit[al]
    print(f"  α={al}: 拟合集重构 RMS = {np.sqrt((resid**2).mean()):.2f} MeV")
ALPHA = min((0.5, 0.75), key=lambda al: np.sqrt(((np.array(samples[al])-delta_fit[al])**2).mean()))
DELTA = delta_fit[ALPHA]
print(f"  → 选定 α={ALPHA}, δ={DELTA:.3f} MeV (对照经验 ~12: {'' if abs(DELTA-12)<4 else '偏离经验值, 需在报告中说明'})")

def pair_term(Z, N):
    return DELTA*(Z+N)**(-ALPHA)*parity(Z, N)

# ============ R4: F2 重测矩阵 ============
print()
print("="*84)
print("第二步: F2 Q_β⁻ 重测矩阵 (判定标准不变: 过≤max(0.4,20%Q) 半≤max(1.0,50%Q))")
print("="*84)
CONFIGS = [
    ('v4基座 a_sym19.5 (原版, 基准)', 19.50, False),
    ('v4基座 a_sym23 (Lane)',        23.00, False),
    ('v4基座 a_sym23+配对 (半经验)',  23.00, True),
    ('v4基座 a_sym19.5+配对 (半经验)', 19.50, True),
]
results = {}
for name, asym, use_pair in CONFIGS:
    rows = []
    for Z, N, nm in TESTS:
        e = [AME.get(kk) for kk in [(Z, N), (Z+1, N-1)]]
        if any(v is None or v.get('ME') is None for v in e):
            continue
        Qm = Bf(Z+1, N-1, asym) - Bf(Z, N, asym) - 0.782
        if use_pair:
            Qm += pair_term(Z+1, N-1) - pair_term(Z, N)
        Qe = (e[0]['ME']-e[1]['ME'])/1000.
        dev = Qm - Qe
        rows.append(dict(nucl=nm, Qm=float(Qm), Qe=float(Qe), dev=float(dev), v=judge(dev, Qe)))
    ok = sum(1 for r in rows if r['v'] == '过')
    hf = sum(1 for r in rows if r['v'] == '半')
    mae = float(np.mean([abs(r['dev']) for r in rows]))
    results[name] = dict(rows=rows, ok=ok, half=hf, fail=len(rows)-ok-hf, mae=mae)

base_mae = results[CONFIGS[0][0]]['mae']
print(f"\n  {'配置':38s} {'过':>3s} {'半':>3s} {'不过':>4s} {'MAE':>8s}  vs基准")
for name, _, _ in CONFIGS:
    r = results[name]
    verdict = '基准' if '原版' in name else ('✅ 改善' if r['mae'] < base_mae else '❌ 变差, 不采纳')
    print(f"  {name:38s} {r['ok']:3d} {r['half']:3d} {r['fail']:4d} {r['mae']:8.2f}  {verdict}")

print()
best = max(CONFIGS, key=lambda c: -results[c[0]]['mae'] if results[c[0]]['mae'] < base_mae else 0)
bname = best[0]
print(f"逐核素明细 (最优可行配置: {bname}):")
for r in results[bname]['rows']:
    print(f"  {r['nucl']:>8s} 实验={r['Qe']:8.3f} 模型={r['Qm']:8.3f} 偏差={r['dev']:+8.3f} {r['v']}")

json.dump(dict(delta=DELTA, alpha=ALPHA, delta_fit_all={str(k): v for k, v in delta_fit.items()},
               n_fitset={str(k): len(v) for k, v in samples.items()},
               configs={k: {kk: vv for kk, vv in v.items() if kk != 'rows'} for k, v in results.items()},
               rows={k: v['rows'] for k, v in results.items()},
               identity='半经验补丁: δ 来自 AME 独立拟合集(排除全部测试核), 非零借用'),
          open('results/q_formula_v6_f2_patch.json', 'w'), indent=1, ensure_ascii=False)
print("\n存档: results/q_formula_v6_f2_patch.json (未覆盖任何旧结果)")
