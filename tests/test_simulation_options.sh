#!/bin/sh
set -eu

ROOT=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
cd "$ROOT"

make >/dev/null

expect_failure() {
  if ./main "$@" >"$ROOT/tests/.options.out" 2>&1; then
    echo "expected failure for: $*" >&2
    exit 1
  fi
  grep -E "Error processing command-line option|positive|non-negative|integer" "$ROOT/tests/.options.out" >/dev/null
}

./main --threads 1 --seed 123 --help >/dev/null
expect_failure --threads 0 --help
expect_failure --threads -1 --help
expect_failure --threads nope --help
expect_failure --threads 999999999999999999999 --help
expect_failure --seed -1 --help
expect_failure --seed nope --help
expect_failure --seed 999999999999999999999999999999 --help

rm -f "$ROOT/tests/.options.out"
echo "simulation option checks passed"
