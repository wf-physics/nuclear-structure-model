# -*- coding: utf-8 -*-
"""
AME2020 / NUBASE2020 官方 ASCII 表解析模块
来源: IAEA AMDC (https://www-nds.iaea.org/amdc/)
  - mass_1.mas20  : 原子质量 [质量过剩, BE/A, Q(beta-), 原子质量]
  - nubase_4.mas20: 半衰期, 自旋宇称, 衰变模式与强度
引用: Huang et al., Chin. Phys. C 45, 030002 (2021);
      Kondev et al., Chin. Phys. C 45, 030001 (2021)
"""
import re

def _is_int(s):
    try:
        int(s); return True
    except ValueError:
        return False

def _is_num(s):
    try:
        float(s.rstrip("#")); return True
    except ValueError:
        return False

def _is_aux(s):
    """O 辅助列 / 粒子不稳定标志: 'x','-','--','+3n','p4n','2p-n','-pp','IT'(同核异能素) 等（排除数字）"""
    if re.fullmatch(r"-?\d+(\.\d+)?[eE]?-?\d*", s):
        return False
    return re.fullmatch(r"[a-z0-9+\-]{1,5}|[A-Z]{2}", s) is not None

_EL_RE = re.compile(r"([A-Z][a-z]?|n|p)$")

def parse_ame(path):
    """解析 mass_1.mas20 → dict[(Z,N)] = {A, ME(keV), MEerr, BEA, Qb(keV), flag, aux, El}
    行结构: [flag?] NZ N Z A El (aux...) ME MEerr BEA BEAerr B- (Qb err | *) A mass err
    注意: 首列 flag 与 N-Z 字段可能合并（flag 为空格时只有 4 个前导整数）。
    """
    nuclei = {}
    with open(path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            tok = line.split()
            if "B-" not in tok:
                continue
            i = tok.index("B-")
            if i < 5 or i + 1 >= len(tok):
                continue
            try:
                ME = float(tok[i-4].rstrip("#")); MEerr = float(tok[i-3].rstrip("#"))
                BEA = float(tok[i-2].rstrip("#")); BEAerr = float(tok[i-1].rstrip("#"))
            except ValueError:
                continue
            # 元素符号: 从 i-5 向左跳过 aux/O 标志（'n'/'p' 是元素符号，优先判定）
            k = i - 5
            while k >= 0 and not _EL_RE.fullmatch(tok[k]) and _is_aux(tok[k]):
                k -= 1
            if k < 0 or not _EL_RE.fullmatch(tok[k]):
                continue
            El = tok[k]
            aux = [tok[m] for m in range(k+1, i-4)]
            # 前导整数: 5 个 = flag,NZ,N,Z,A；4 个 = NZ,N,Z,A（flag 空缺）
            lead = []
            m = k - 1
            while m >= 0 and _is_int(tok[m]):
                lead.append(int(tok[m])); m -= 1
            if len(lead) == 5:
                lead = lead[::-1]     # 从右往左收集，需反转
                flag, NZ, N, Z, A = lead[0], lead[1], lead[2], lead[3], lead[4]
            elif len(lead) == 4:
                lead = lead[::-1]
                flag, NZ, N, Z, A = 0, lead[0], lead[1], lead[2], lead[3]
            else:
                continue
            Qb = None if tok[i+1] == "*" else float(tok[i+1].rstrip("#"))
            nuclei[(Z, N)] = dict(A=A, ME=ME, MEerr=MEerr, BEA=BEA, Qb=Qb,
                                  flag=flag, aux=aux, El=El)
    return nuclei

# 半衰期单位 → 秒
_T_UNITS = {"ys": 1e-24, "zs": 1e-21, "as": 1e-18, "fs": 1e-15, "ps": 1e-12,
            "ns": 1e-9, "us": 1e-6, "ms": 1e-3, "s": 1.0, "m": 60.0,
            "h": 3600.0, "d": 86400.0, "y": 3.1557e7, "ky": 3.1557e10,
            "my": 3.1557e13, "gy": 3.1557e16, "ty": 3.1557e19}

def _parse_decay_mode(br):
    """解析 BR 字段的主衰变模式 → 'bm'/'bp'/'alpha'/'p'/'n'/'other'/None
    按衰变分支强度取主模式；IS=（同位素丰度）不算衰变，仅当无任何衰变分支时返回 'stable'。"""
    if not br.strip():
        return None
    has_is = False
    best = (None, -1.0)          # (mode, intensity)
    for part in br.split(";"):
        part = part.strip()
        m = re.match(r"^\s*(B-|B\+|EC|IS|2A|2P|2N|2p|2n|[Aptn])", part)
        if not m:
            continue
        code = m.group(1)
        rest = part[m.end():].strip()
        inten = None
        rm = re.match(r"=([\d.]+)", rest)
        if rm:
            inten = float(rm.group(1))
        if code == "IS":
            has_is = True
            continue
        mode = {"B-": "bm", "B+": "bp", "EC": "bp", "A": "alpha", "2A": "alpha",
                "P": "p", "2P": "p", "p": "p", "2p": "p",
                "N": "n", "2N": "n", "n": "n", "2n": "n"}.get(code, "other")
        iv = inten if inten is not None else 0.0
        if iv > best[1]:
            best = (mode, iv)
    if best[0] is not None:
        return best[0]
    if has_is:
        return "stable"
    return None

def parse_nubase(path):
    """解析 nubase_4.mas20（仅基态 i=0）→ dict[(Z,N)] = {A, T(s|None), stable, mode}"""
    nuclei = {}
    with open(path) as f:
        for line in f:
            if line.startswith("#") or len(line.rstrip()) < 100:
                continue
            try:
                A = int(line[0:3])
                zzi = int(line[4:8])
            except ValueError:
                continue
            if zzi % 10 != 0:          # 跳过同核异能态/共振/能级行
                continue
            Z = zzi // 10
            N = A - Z
            T_str = line[69:78].strip()
            unit = line[78:80].strip()
            br = line[119:209] if len(line) > 119 else ""
            rec = dict(A=A, T=None, stable=False, mode=None, raw_T=T_str)
            if T_str == "stbl":
                rec["stable"] = True
                rec["T"] = float("inf")
                rec["mode"] = "stable"
            elif T_str == "p-unst":
                rec["T"] = 0.0
                rec["mode"] = _parse_decay_mode(br) or "n"
            else:
                try:
                    val = float(T_str.split("#")[0])
                    rec["T"] = val * _T_UNITS.get(unit.lower(), None)
                    rec["stable"] = rec["T"] > 1e12 * 3.1557e7   # > 1e12 y 视作实际稳定
                    rec["mode"] = _parse_decay_mode(br)
                except (ValueError, TypeError):
                    rec["T"] = None
                    rec["mode"] = _parse_decay_mode(br)
            if rec["stable"] and rec["mode"] is None:
                rec["mode"] = "stable"
            nuclei[(Z, N)] = rec
    return nuclei

if __name__ == "__main__":
    import sys, os
    d = os.path.join(os.path.dirname(__file__), "data")
    ame = parse_ame(os.path.join(d, "mass_1.mas20"))
    nub = parse_nubase(os.path.join(d, "nubase_4.mas20"))
    print(f"AME2020 基态核素: {len(ame)};  NUBASE2020 基态核素: {len(nub)}")
    for key in [(0,1), (1,0), (1,1), (1,2), (2,2), (6,8), (20,20)]:
        a = ame.get(key); b = nub.get(key)
        print(key, "AME:", {k: a[k] for k in ("A","ME","BEA","Qb","aux") if a and k in a},
              "| NUBASE:", b)
    n_le40 = [k for k in ame if ame[k]["A"] <= 40]
    print(f"A<=40 的 AME 核素数: {len(n_le40)}")
