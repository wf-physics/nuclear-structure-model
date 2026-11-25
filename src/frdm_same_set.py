"""FRDM2012 vs 本文闭式的同批核素对比（Paper I §4.4 口径脚注）。
本文闭式评估集：AME2020 A=40-240 共 2811 核（路线一，MAE 12.6）。
对比集：该 2811 核与 FRDM2012 表的交集；误差定义 = B_model - B_exp（MeV）。
另附 A>=16 全集口径（3487 核）供参考。"""
import csv, json, os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))

ame = {}
with open(os.path.join(HERE, "..", "data", "mass_1.mas20"), encoding="utf-8", errors="replace") as f:
    lines = f.readlines()[15:]
for ln in lines:
    if len(ln) < 67:
        continue
    try:
        n = int(ln[4:9]); z = int(ln[9:14]); a = int(ln[14:19])
    except ValueError:
        continue
    if a == 0 or z == 0:
        continue
    try:
        bea = float(ln[54:67].strip().replace("#", ""))
    except ValueError:
        continue
    ame[(z, n)] = bea * a / 1000.0

frdm = {}
with open(os.path.join(HERE, "..", "data", "frdm2012_masses.csv")) as fh:
    for row in csv.DictReader(fh):
        try:
            z = int(row["Z"]); n = int(row["N"]); b = float(row["E_bind"])
        except (ValueError, KeyError, TypeError):
            continue
        frdm[(z, n)] = b

out = {"frdm_nuclides": len(frdm)}
for label, cond in [("paper_set_A40_240", lambda z, n: 40 <= z + n <= 240),
                    ("full_set_Age16", lambda z, n: z + n >= 16)]:
    keys = [k for k in ame if k in frdm and cond(*k)]
    dev = np.array([frdm[k] - ame[k] for k in keys])
    out[label] = {"n_common": len(keys),
                  "MAE": round(float(np.mean(np.abs(dev))), 4),
                  "RMS": round(float(np.sqrt(np.mean(dev ** 2))), 4)}
    print(f"{label}: n={len(keys)} MAE={out[label]['MAE']} RMS={out[label]['RMS']}")
json.dump(out, open(os.path.join(HERE, "..", "results", "frdm_same_set.json"), "w"), indent=1)
print("已存 results/frdm_same_set.json")
