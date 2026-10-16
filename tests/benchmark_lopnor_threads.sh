#!/bin/sh
set -eu

ROOT=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
cd "$ROOT"

phonons=${LOPNOR_BENCH_PHONONS:-1000000}
target_serial_seconds=${LOPNOR_BENCH_TARGET_SECONDS:-600}
output_root=${LOPNOR_BENCH_OUTPUT_DIR:-}

make >/dev/null

if [ -z "$output_root" ]; then
  output_root=$(mktemp -d "${TMPDIR:-/tmp}/radiative3d-lopnor-bench.XXXXXX")
  cleanup=1
else
  mkdir -p "$output_root"
  cleanup=0
fi
trap '[ "$cleanup" -eq 1 ] && rm -rf "$output_root"' EXIT HUP INT TERM

run_benchmark_case() {
  mode=$1
  name=$2
  start=$(python3 -c 'import time; print(time.time_ns())')
  "$ROOT/tests/run_lopnor_thread_case.sh" "$mode" "$output_root/$name" \
    "$phonons" "$output_root/$name" >/dev/null
  end=$(python3 -c 'import time; print(time.time_ns())')
  elapsed_ns=$((end - start))
  [ "$elapsed_ns" -gt 0 ] || elapsed_ns=1
  elapsed=$(awk -v ns="$elapsed_ns" 'BEGIN {printf "%.3f", ns/1000000000}')
  workers=$(sed -n 's/^Workers: //p' "$output_root/$name/stdout.txt" | tail -1)
  completed=$(sed -n 's/^100% of \([0-9][0-9]*\) have been cast.$/\1/p' \
    "$output_root/$name/stderr.txt" | tail -1)
  rate=$(awk -v n="$completed" -v t="$elapsed" 'BEGIN {printf "%.2f", n/t}')
  printf '%s\t%s\t%s\t%s\t%s\n' "$name" "$workers" "$completed" "$elapsed" "$rate" \
    >> "$output_root/benchmark.tsv"
}

: > "$output_root/benchmark.tsv"
run_benchmark_case 1 serial
run_benchmark_case 12 threads12
run_benchmark_case auto automatic

serial_seconds=$(awk -F '\t' '$1 == "serial" {print $4}' "$output_root/benchmark.tsv")
if awk -v actual="$serial_seconds" -v target="$target_serial_seconds" 'BEGIN {exit !(actual < target)}'; then
  echo "warning: serial run took ${serial_seconds}s; target was ${target_serial_seconds}s" >&2
fi

awk -F '\t' -v serial="$serial_seconds" '
  BEGIN {print "Mode\tWorkers\tPhonons\tSeconds\tPhonons/s\tSpeedup"}
  {printf "%s\t%s\t%s\t%s\t%s\t%.2fx\n", $1,$2,$3,$4,$5,serial/$4}
' "$output_root/benchmark.tsv" | tee "$output_root/benchmark-table.tsv"

echo "Benchmark output: $output_root"
