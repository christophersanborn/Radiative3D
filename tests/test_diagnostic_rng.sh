#!/bin/sh
set -eu

ROOT=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
cd "$ROOT"
make >/dev/null
./main --rtcoef-test >tests/.rtcoef.log 2>&1
test -s tests/.rtcoef.log
rm -f tests/.rtcoef.log
echo "diagnostic RNG checks passed"
