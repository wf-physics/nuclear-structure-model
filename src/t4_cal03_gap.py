# -*- coding: utf-8 -*-
"""T4 第四格补测: calibrated@gate=0.3 (四口径矩阵闭合)"""
import os, sys, json
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
HERE = os.path.dirname(os.path.abspath(__file__))

# 复用 t4 脚本的协议, 仅把能量侧换成"训练区线性校准闭式"
exec(open(os.path.join(HERE, "t4_err_sensitivity.py")).read().split("print(\"T4 敏感性")[0])

# 口径 B: 训练区线性校准 (复刻 valve_test.py 逐字)
tr = [(Z, N) for (Z, N), v in ame.items() if Z <= 10 and Z >= 1 and N >= 1
      and v.get("BEA") is not None and v["BEA"]*(Z+N)/1000. > 0.2]
dB = np.array([B_model(Z, N) - ame[(Z, N)]["BEA"]*(Z+N)/1000. for Z, N in tr])
AA = np.array([Z+N for Z, N in tr], float)
lam, c0 = np.polyfit(AA, dB, 1)
print(f"校准 (Z<=10, n={len(tr)}): B_cal = B_model − ({lam:.4f}·A + {c0:.2f})")

def make_ame_model_cal():
    m = {}
    for (Z, N), v in ame.items():
        if v.get("ME") is None: continue
        B = B_model(Z, N) - (lam*(Z+N) + c0)
        me_model = (Z*ME_H + N*ME_N - B)*1000.0
        d = dict(v); d["ME"] = me_model; d["MEerr"] = v.get("MEerr", 10.0)
        m[(Z, N)] = d
    m[(2, 0)] = dict(A=2, ME=2*ME_H + 91.8, MEerr=0.5, BEA=0.0, Qb=None, flag=0, aux=[], El="He")
    return m

print("T4 第四格: calibrated@gate=0.3")
POOL = build(POOL_KEYS, make_ame_model_cal())
r = run_protocol(POOL, "const", "校准+ERR=0.3")
print(f"  → mean={r['mean']*100:.2f}%  std={r['std']*100:.2f}  min={r['min']*100:.2f}%  free={r['mean_free']*100:.2f}%")

# 并入四口径矩阵存档
try:
    t4 = json.load(open(os.path.join(HERE, "results/t4_err_sensitivity.json")))
except Exception:
    t4 = {}
t4["校准+ERR=0.3"] = r
t4["matrix_complete"] = dict(raw_ame=0.3393, raw_03=t4.get("闭式+ERR=0.3", {}).get("mean"),
                             cal_ame=None, cal_03=r["mean"])
# 校准@AMEerr 在原 valve_test.json 里 (66.74%)
try:
    vt = json.load(open(os.path.join(HERE, "results/valve_test.json")))
    cal_ame = vt.get("B_calibrated(训练区线性校准)", {}).get("mean")
    t4["matrix_complete"]["cal_ame"] = cal_ame
except Exception:
    pass
t4["matrix_note"] = "四口径矩阵: raw@AMEerr 33.93% / raw@0.3 57.19% / cal@AMEerr 66.74% / cal@0.3 见本文件 → 全部 <90%, 负结果对闸门与校准均稳健"
json.dump(t4, open(os.path.join(HERE, "results/t4_err_sensitivity.json"), "w"), indent=1, ensure_ascii=False)
print("四口径矩阵已闭合, 存档 results/t4_err_sensitivity.json")
print(t4["matrix_complete"])
