# -*- coding: utf-8 -*-
"""
差分公式族 v5 — 采用 DeepSeek 重标井深 B 闭式 (F4 表面项复活后重测)
路线: DeepSeek 用冻结线性井做金刚石团簇求和 + 井深重标(×2.61), 给出
     B_DS(A,Z) = 16.060·A − 19.010·A^(2/3) − 0.7560·Z²/A^(1/3) − 19.50·(N−Z)²/A
     零逐核拟合, 其对照表 (A=40..240, 2811 核素) MAE=12.64 已验收.
本脚本: 同一组测试点/判定标准 (与 q_formula_v4.py 完全一致), 只换 B 函数 → 重跑 F1-F4.
     附 Lane 变体: a_sym 19.50→23.0, 只测 F2 (对称能缺额归因检验).
     A<40 与 A>240 为闭式外推区, 结果单独标注.
"""
import numpy as np, json, sys
sys.path.insert(0, '.')
from nucdata import parse_ame
AME = parse_ame('data/mass_1.mas20')

# DeepSeek 重标闭式 (由其对照表反推, RMS=0.000, 精确公式而非拟合)
AV, AS, AC, ASYM = 16.060, 19.010, 0.7560, 19.50

def Bf(Z, N, a_sym=ASYM):
    A = Z + N
    return AV*A - AS*A**(2/3) - AC*Z**2/A**(1/3) - a_sym*(N-Z)**2/A

def Bexp(Z, N):
    v = AME.get((Z, N))
    return v['BEA']*v['A']/1000.0 if v and v.get('BEA') else None

def judge(dev, Qe):
    return ('过' if abs(dev) <= max(0.4, 0.2*abs(Qe)) else
            '半' if abs(dev) <= max(1.0, 0.5*abs(Qe)) else '不过')

rows = {'F1_Sn': [], 'F2_Qbm': [], 'F3_Qa': [], 'F4_Qfs': [], 'F2_Qbm_Lane': []}

print("="*84)
print("v5: DeepSeek 重标井深闭式 B 重跑 F1-F4 (同 v4 测试点/判定标准)")
print(f"B = {AV}·A − {AS}·A^(2/3) − {AC}·Z²/A^(1/3) − {ASYM}·(N−Z)²/A")
print("="*84)

print("--- F1  S_n = B(A,Z) − B(A−1,Z) ---")
for Z, N, nm in [(20,20,'Ca-40'),(20,28,'Ca-48'),(26,30,'Fe-56'),(28,30,'Ni-58'),
                 (30,34,'Zn-64'),(36,44,'Kr-80'),(38,50,'Sr-88'),(40,50,'Zr-90'),
                 (50,62,'Sn-112'),(82,126,'Pb-208')]:
    v1, v2 = AME.get((Z,N)), AME.get((Z,N-1))
    if not (v1 and v2 and v1.get('BEA') and v2.get('BEA')):
        continue
    Sn_m = Bf(Z, N) - Bf(Z, N-1)
    Sn_e = v1['BEA']*v1['A']/1000. - v2['BEA']*v2['A']/1000.
    dev = Sn_m - Sn_e
    vd = '过' if abs(dev) <= 2.0 else ('半' if abs(dev) <= 5.0 else '不过')
    print(f"  {nm:>8s} 模型={Sn_m:8.2f} 实验={Sn_e:8.2f} 偏差={dev:+8.2f} {vd}")
    rows['F1_Sn'].append(dict(nucl=nm, Sn_m=Sn_m, Sn_e=Sn_e, dev=dev, v=vd))

print("--- F2  Q_β⁻ = B(A,Z+1) − B(A,Z) − 0.782   [双列: 原版 a_sym=19.50 | Lane 变体 23.0] ---")
for Z, N, nm in [(6,8,'C-14'),(11,13,'Na-24'),(15,17,'P-32'),(16,19,'S-35'),
                 (19,23,'K-42'),(25,31,'Mn-56'),(26,33,'Fe-59'),
                 (38,52,'Sr-90'),(55,82,'Cs-137'),(61,86,'Pm-147')]:
    e = [AME.get(kk) for kk in [(Z,N),(Z+1,N-1)]]
    if any(v is None or v.get('ME') is None for v in e):
        continue
    Qm  = Bf(Z+1, N-1) - Bf(Z, N) - 0.782
    QmL = Bf(Z+1, N-1, 23.0) - Bf(Z, N, 23.0) - 0.782
    Qe = (e[0]['ME']-e[1]['ME'])/1000.
    dev, devL = Qm-Qe, QmL-Qe
    vd, vdL = judge(dev, Qe), judge(devL, Qe)
    print(f"  {nm:>8s} 实验={Qe:8.3f}  原19.5: {Qm:8.3f} ({dev:+7.3f}) {vd}   Lane23: {QmL:8.3f} ({devL:+7.3f}) {vdL}")
    rows['F2_Qbm'].append(dict(nucl=nm, Qm=Qm, Qe=Qe, dev=dev, v=vd))
    rows['F2_Qbm_Lane'].append(dict(nucl=nm, Qm=QmL, Qe=Qe, dev=devL, v=vdL))

print("--- F3  Q_α = B(A−4,Z−2) + B_α(锚28.296) − B(A,Z) ---")
for Z, N, nm in [(62,85,'Sm-147'),(64,88,'Gd-152'),(74,106,'W-180'),(84,126,'Po-210'),
                 (86,136,'Rn-222'),(88,138,'Ra-226'),(90,142,'Th-232'),
                 (92,143,'U-235'),(92,146,'U-238'),(94,145,'Pu-239'),
                 (94,146,'Pu-240'),(95,146,'Am-241'),(96,148,'Cm-244'),(98,154,'Cf-252')]:
    v1, v2 = AME.get((Z,N)), AME.get((Z-2,N-2))
    if not (v1 and v2 and v1.get('ME') is not None and v2.get('ME') is not None):
        continue
    Qm = Bf(Z-2, N-2) + 28.2957 - Bf(Z, N)
    Qe = (v1['ME']-v2['ME']-2424.916)/1000.
    dev = Qm - Qe
    vd = judge(dev, Qe)
    print(f"  {nm:>8s} 模型={Qm:8.2f} 实验={Qe:8.2f} 偏差={dev:+8.2f} {vd}")
    rows['F3_Qa'].append(dict(nucl=nm, Qm=Qm, Qe=Qe, dev=dev, v=vd))

print("--- F4  Q_fs = B(f1)+B(f2) − B(母核) [重测: 表面项已复活] ---")
FISSIONS = [
    (92,146,'U-238', [(54,86,'Xe-140'),(38,57,'Sr-95')], 3),
    (92,146,'U-238', [(46,73,'Pd-119'),(46,73,'Pd-119')], 0),
    (92,143,'U-235', [(56,85,'Ba-141'),(36,53,'Kr-89')], 3),
    (94,146,'Pu-240', [(57,88,'La-145'),(37,55,'Rb-92')], 3),
    (98,154,'Cf-252', [(56,86,'Ba-142'),(42,65,'Mo-107')], 3),
]
for Z, N, nm, frags, xn in FISSIONS:
    vp = AME.get((Z, N))
    if not (vp and vp.get('BEA')):
        continue
    Qm = -Bf(Z, N) + sum(Bf(fz, fn) for fz, fn, _ in frags)
    Qe = -vp['BEA']*vp['A']/1000.
    ok = True
    for fz, fn, _ in frags:
        vf = AME.get((fz, fn))
        if vf is None or vf.get('BEA') is None:
            ok = False; break
        Qe += vf['BEA']*vf['A']/1000.
    if not ok:
        print(f"  {nm:>8s} 碎片实验质量缺"); continue
    dev = Qm - Qe
    vd = judge(dev, Qe)
    frag_s = '+'.join(f[2] for f in frags) + (f'+{xn}n' if xn else '')
    print(f"  {nm:>8s} {frag_s:>18s} 模型={Qm:8.1f} 实验={Qe:8.1f} 偏差={dev:+8.1f} {vd}")
    rows['F4_Qfs'].append(dict(parent=nm, frags=frag_s, Qm=Qm, Qe=Qe, dev=dev, v=vd))

print()
print("="*84)
print("汇总: v5 (重标闭式) vs v4 (模型自有闭式)")
print("="*84)
v4 = json.load(open('results/q_formula_v4.json'))
print(f"  {'公式':10s} {'n':>3s} {'过':>3s} {'半':>3s} {'不过':>4s} {'MAE_v5':>8s} {'MAE_v4':>8s}")
for fid in ['F1_Sn', 'F2_Qbm', 'F2_Qbm_Lane', 'F3_Qa', 'F4_Qfs']:
    t = rows[fid]
    if not t:
        continue
    ok = sum(1 for r in t if r['v'] == '过')
    hf = sum(1 for r in t if r['v'] == '半')
    mae = np.mean([abs(r['dev']) for r in t])
    f4key = 'F2_Qbm' if fid == 'F2_Qbm_Lane' else fid
    old = v4.get(f4key, [])
    mae4 = np.mean([abs(r['dev']) for r in old]) if old else float('nan')
    tag = ' (Lane a_sym→23)' if fid == 'F2_Qbm_Lane' else ''
    print(f"  {fid:10s}{tag:18s} {len(t):3d} {ok:3d} {hf:3d} {len(t)-ok-hf:4d} {mae:8.2f} {mae4:8.2f}")

json.dump(rows, open('results/q_formula_v5.json', 'w'), indent=1, ensure_ascii=False)
json.dump(dict(source='DeepSeek重标闭式(团簇求和×2.61井深重标, 零逐核拟合)',
               a_v=AV, a_s=AS, a_c=AC, a_sym=ASYM,
               accept={'MAE_all':12.64, 'MAE_40_65':6.48, 'MAE_65_100':9.30, 'MAE_100_240':14.33},
               extrapolate='A<40 与 A>240 为闭式外推区'),
          open('results/q_formula_v5_meta.json', 'w'), indent=1, ensure_ascii=False)
print("存档: results/q_formula_v5.json  q_formula_v5_meta.json")
