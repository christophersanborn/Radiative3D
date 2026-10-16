#!/bin/sh
set -eu

ROOT=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
cd "$ROOT"
make >/dev/null

./main --num-phonons 8 --threads 2 --seed 77 --timetolive 0.001 \
  --reports GEN,SCT --report-file tests/.reports.dat >tests/.output.log 2>&1
grep "Workers: 2" tests/.output.log >/dev/null
grep "Random Seed: 77" tests/.output.log >/dev/null
grep "__SIMULATION_COMPLETE__" tests/.output.log >/dev/null
test -s tests/.reports.dat

rm -f tests/.output.log tests/.reports.dat
echo "multi-worker output checks passed"
