#!/bin/bash
# 一键复现: 分两层
#   core : 分钟级, 覆盖能量公式族全部成绩 + 阀门验证 + F1 扩样 + 核物质标定复核
#   full : 小时级, 几何大验证 + 盲测协议 + 反应链 (结果已存档, 重跑用于验证存档)
set -e
cd "$(dirname "$0")/src"

CORE="
nm_formula_v4.py
q_formula_v4.py
q_formula_v5.py
q_formula_v6_f2patch.py
q_formula_v6_pairing_v2.py
q_engine_pairing.py
f1_expanded.py
valve_test.py
finite_nuclei_v1.py
"

FULL="
bond_calibration.py
relax_all.py
chain_predict.py
big_validation.py
geometry_validation.py
blind_test.py
criterion_test.py
"

mkdir -p ../results
if [ "$1" = "core" ]; then
  echo "=== 核心复现 (分钟级) ==="
  for s in $CORE; do echo "--- $s"; python3 "$s"; done
elif [ "$1" = "full" ]; then
  echo "=== 全量复现 (小时级) ==="
  for s in $FULL; do echo "--- $s"; python3 "$s"; done
else
  echo "用法: bash run_all.sh [core|full]"
fi
