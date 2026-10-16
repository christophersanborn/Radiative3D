#!/bin/sh
set -eu

ROOT=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
cd "$ROOT"
make >/dev/null

./main --num-phonons 8 --threads 2 --seed 9 --timetolive 0.001 \
  >tests/.accumulation.log 2>&1
grep "__SIMULATION_COMPLETE__" tests/.accumulation.log >/dev/null
grep "100% of 8" tests/.accumulation.log >/dev/null

rm -f tests/.accumulation.log
echo "multi-worker accumulation checks passed"
