# intersect: the off-by-one subcommand. Normal intervals overlap on strict <, so
# bookended features miss; zero-length features use an inclusive predicate
# (SPEC.md §3).
check "intersect default" -- intersect -a "$DATA/a.bed" -b "$DATA/b.bed"
check "intersect -u"      -- intersect -u  -a "$DATA/a.bed" -b "$DATA/b.bed"
check "intersect -v"      -- intersect -v  -a "$DATA/a.bed" -b "$DATA/b.bed"
check "intersect -wa"     -- intersect -wa -a "$DATA/a.bed" -b "$DATA/b.bed"
check "intersect -wb"     -- intersect -wb -a "$DATA/a.bed" -b "$DATA/b.bed"
check "intersect -c"      -- intersect -c  -a "$DATA/a.bed" -b "$DATA/b.bed"
check "intersect -s"      -- intersect -s  -a "$DATA/a.bed" -b "$DATA/b.bed"
check "intersect -S"      -- intersect -S  -a "$DATA/a.bed" -b "$DATA/b.bed"
check "intersect -f 0.5"  -- intersect -f 0.5 -a "$DATA/a.bed" -b "$DATA/b.bed"

# Combinations and the reversed argument order, exercising the same paths from a
# different angle: -wa with -wb, and a fraction landing exactly on a boundary.
check "intersect -wa -wb" -- intersect -wa -wb -a "$DATA/a.bed" -b "$DATA/b.bed"
check "intersect -c -s"   -- intersect -c -s  -a "$DATA/a.bed" -b "$DATA/b.bed"
check "intersect -v -s"   -- intersect -v -s  -a "$DATA/a.bed" -b "$DATA/b.bed"
check "intersect -f 1.0"  -- intersect -f 1.0 -a "$DATA/a.bed" -b "$DATA/b.bed"

# Swapping the files makes b.bed's zero-length features the -a side and a.bed's the
# -b side. a12 has to come out first: a zero-length feature at position 0 in -b
# widens to [-1, 1) and bedtools' binning index rejects the negative coordinate,
# aborting with "ERROR: Received illegal bin number -1 from getBin call." and exit 1
# before printing anything. Measured: it is specific to the -b side, since the same
# record as -a exits 0. A bedtools limitation, not an answer to diff against; we
# handle the record fine. See SPEC.md §3.
grep -v $'^chr2\t0\t0\t' "$DATA/a.bed" > "$tmp/a_no_zero_at_origin.bed"
check "intersect b vs a"  -- intersect -a "$DATA/b.bed" -b "$tmp/a_no_zero_at_origin.bed"

check "intersect genes"   -- intersect -wb -a "$DATA/genes.bed" -b "$DATA/hg002.highconf.bed"
