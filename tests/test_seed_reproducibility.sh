#!/bin/sh
set -eu

ROOT=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
cd "$ROOT"
make >/dev/null

run_once() {
  output=$1
  ./main --num-phonons 2 --threads 1 --seed 123 --timetolive 0.001 \
    --reports GEN --report-file "$output.report" >"$output.log" 2>&1
  grep "__SIMULATION_COMPLETE__" "$output.log" >/dev/null
}

run_once tests/.seed-a
run_once tests/.seed-b
test -s tests/.seed-a.report
sed -E 's/cell: 0x[0-9a-f]+/cell: 0xADDR/g' tests/.seed-a.report > tests/.seed-a.normalized
sed -E 's/cell: 0x[0-9a-f]+/cell: 0xADDR/g' tests/.seed-b.report > tests/.seed-b.normalized
cmp tests/.seed-a.normalized tests/.seed-b.normalized

rm -f tests/.seed-a.* tests/.seed-b.*
echo "seed reproducibility checks passed"
