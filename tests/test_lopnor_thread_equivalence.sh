#!/bin/sh
set -eu

ROOT=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
cd "$ROOT"

make >/dev/null

work=$(mktemp -d "${TMPDIR:-/tmp}/radiative3d-lopnor.XXXXXX")
trap 'rm -rf "$work"' EXIT HUP INT TERM

run_case() {
  "$ROOT/tests/run_lopnor_thread_case.sh" "$1" "$work/$2"
}

run_case 1 one
run_case 12 twelve
run_case auto automatic

for case in one twelve automatic; do
  grep '__SIMULATION_COMPLETE__' "$work/$case/stdout.txt" >/dev/null
  grep '^Workers:' "$work/$case/stdout.txt" >/dev/null
  test -s "$work/$case/seis_traces_asc.dat"
  grep -E '^[[:space:]]*[0-9.eE+-]+[[:space:]]+' "$work/$case/seis_traces_asc.dat" >/dev/null
done

normalize_trace() {
  sed -E 's/cell: 0x[0-9a-f]+/cell: 0xADDR/g' "$1"
}

normalize_trace "$work/one/seis_traces_asc.dat" > "$work/one.trace"
normalize_trace "$work/twelve/seis_traces_asc.dat" > "$work/twelve.trace"
normalize_trace "$work/automatic/seis_traces_asc.dat" > "$work/automatic.trace"
cmp "$work/one.trace" "$work/twelve.trace"
cmp "$work/one.trace" "$work/automatic.trace"

echo "Lop Nor thread-equivalence checks passed without Octave"
