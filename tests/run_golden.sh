#!/usr/bin/env bash
# Golden tests: diff mytools against real bedtools on the fixtures in data/.
# bedtools is the oracle -- if we differ on valid input, we are wrong, not it.
#
# Usage: ./tests/run_golden.sh
#        MYTOOLS=/path/to/mytools ./tests/run_golden.sh    (test another build)
#
# Cases live in tests/golden/<subcommand>.sh, one file per subcommand, sourced
# below. Add a subcommand by adding a file -- never by editing this one. That is
# what keeps parallel branches from colliding here.
set -uo pipefail

HERE=$(cd "$(dirname "$0")" && pwd)

# Default to THIS checkout's mytools, not whatever `mytools` resolves to on PATH:
# the symlink in ~/.local/bin points at one particular clone, so `${MYTOOLS:-mytools}`
# would silently test a different working tree and report someone else's bugs.
MYTOOLS=${MYTOOLS:-$HERE/../mytools}
DATA=$HERE/../data

if ! command -v bedtools >/dev/null; then
  echo "bedtools is not installed -- the golden tests need the oracle" >&2
  exit 2
fi
if [[ ! -x $MYTOOLS ]]; then
  echo "not executable: $MYTOOLS" >&2
  exit 2
fi

tmp=$(mktemp -d); trap 'rm -rf "$tmp"' EXIT
pass=0; fail=0; skip=0

# check <name> -- <args...>
#   Runs "$MYTOOLS <args>" and "bedtools <args>" and compares exit code first,
#   then stdout byte for byte.
#
#   Set STDIN=<file> for one call to feed both on stdin.
#
#   Two things deliberately not done, both from tests/README.md:
#     - stderr is not compared. bedtools' wording is its own.
#     - output is never sorted before diffing. Row order is part of the answer.
check() {
  local name=$1; shift; shift        # drop the literal --
  local stdin=${STDIN:-/dev/null}

  "$MYTOOLS" "$@" > "$tmp/got"  2>"$tmp/got.err"  < "$stdin"
  local got_rc=$?
  bedtools   "$@" > "$tmp/want" 2>"$tmp/want.err" < "$stdin"
  local want_rc=$?

  if [[ $got_rc -ne $want_rc ]]; then
    echo "FAIL $name (exit $got_rc, bedtools gave $want_rc)"
    sed 's/^/      /' "$tmp/got.err" | head -3
    (( fail++ )); return
  fi
  if diff -q "$tmp/want" "$tmp/got" >/dev/null; then
    echo "ok   $name"; (( pass++ ))
  else
    echo "FAIL $name"
    diff -u "$tmp/want" "$tmp/got" | sed 's/^/      /' | head -20
    (( fail++ ))
  fi
}

# skip <name> <why> -- for a case that cannot run here. A scale test that skips
# with a reason is a test; one that silently passes is a lie (tests/README.md).
skip() { echo "skip $1 ($2)"; (( skip++ )); }

# ---------------------------------------------------------------------------
# Shared fixtures. a.bed and b.bed are deliberately unsorted, so merge and
# closest need sorted copies -- prepared once here rather than in each case file.
# ---------------------------------------------------------------------------
bedtools sort -i "$DATA/a.bed"     > "$tmp/sorted_a.bed"
bedtools sort -i "$DATA/b.bed"     > "$tmp/sorted_b.bed"
bedtools sort -i "$DATA/genes.bed" > "$tmp/sorted_genes.bed"
export SORTED_A=$tmp/sorted_a.bed
export SORTED_B=$tmp/sorted_b.bed
export SORTED_GENES=$tmp/sorted_genes.bed

gzip -c "$DATA/a.bed"      > "$tmp/a.bed.gz"
gzip -c "$tmp/sorted_a.bed" > "$tmp/sorted_a.bed.gz"
export GZ_A=$tmp/a.bed.gz
export GZ_SORTED_A=$tmp/sorted_a.bed.gz

: > "$tmp/empty.bed"
export EMPTY=$tmp/empty.bed

# Headers plus sorted data, so both `sort -header` and `merge -header` can use it.
{ printf '# a comment\ntrack name=workshop\nbrowser position chr1\n'
  cat "$tmp/sorted_a.bed"; } > "$tmp/headed.bed"
export HEADED=$tmp/headed.bed

# ---------------------------------------------------------------------------
# Cases, one file per subcommand.
# ---------------------------------------------------------------------------
shopt -s nullglob
cases=("$HERE"/golden/*.sh)
if [[ ${#cases[@]} -eq 0 ]]; then
  echo "no case files in $HERE/golden/ -- nothing to check yet"
fi
for case_file in "${cases[@]}"; do
  # shellcheck source=/dev/null
  source "$case_file"
done

echo "---"
echo "$pass passed, $fail failed, $skip skipped"
[[ $fail -eq 0 ]]
