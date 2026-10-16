#!/bin/sh
set -eu

ROOT=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
cd "$ROOT"
make >/dev/null

run_and_check() {
  expected=$1
  shift
  ./main --num-phonons 1 --threads "$@" --timetolive 0.001 2>&1 | grep "Workers: $expected" >/dev/null
}

run_and_check 1 1
run_and_check 7 7

auto_workers=$(./main --num-phonons 32 --timetolive 0.001 2>&1 |
  sed -n 's/^Workers: //p' | tail -1)
performance_workers=$(sysctl -n hw.perflevel0.logicalcpu 2>/dev/null || true)
if [ -n "$performance_workers" ] && [ "$performance_workers" -gt 0 ]; then
  expected=$performance_workers
  if [ "$expected" -gt 2 ]; then expected=$((expected - 1)); fi
  test "$auto_workers" -eq "$expected"
fi

echo "worker policy checks passed"
