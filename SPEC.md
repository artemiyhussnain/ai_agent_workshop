# SPEC.md — mytools

A small reimplementation of a subset of bedtools. Real `bedtools` (v2.31.1, installed on
this VM) is the oracle: where this document and bedtools disagree on valid input,
bedtools is right and this document is a bug.

Every behaviour below marked **(measured)** was established by running bedtools v2.31.1
on the fixtures in `data/`, not by reading its documentation. Where bedtools' behaviour
is deliberately *not* copied, the section says so and says why.

---

## 1. Scope

v1 ships five subcommands: `sort`, `merge`, `intersect`, `subtract`, `closest`.

Explicitly **not** in v1 — say no to all of it:

- `-split` (BED12 block-aware mode). Deferred; see §4.
- `-c` / `-o` column aggregation on `merge`.
- `-wo`, `-loj`, and the other `intersect` output variants not listed in §4.
- `genomecov`, `flank`, `slop`, `getfasta`, or any subcommand not in the five above.
- GFF3, VCF, and BAM input. BED only (§2).
- `-g` genome files and custom chromosome ordering.

## 2. Input formats

- **BED only.** BED3 through BED12, tab-separated. No GFF3, VCF or BAM.
- Column count may vary between lines. Missing trailing columns are absent, not empty.
- **Both plain and gzipped input**, detected by content (gzip magic bytes `1f 8b`), not
  by filename. bedtools accepts `.gz` transparently — exit 0 on a gzipped `data/a.bed`
  **(measured)**.
- Read from a file argument or stdin. **`-` means stdin** **(measured:
  `cat data/a.bed | bedtools sort -i -` works)**. At most one input may be `-`.
- `#`, `track` and `browser` lines are **skipped silently** by default, and **echoed
  ahead of the output when `-header` is given** **(measured: `bedtools sort -i` drops
  all three; `bedtools sort -header -i` emits all three, in input order, before the
  data)**.
- Blank lines are skipped.

## 3. Interval semantics

BED is **0-based, half-open**. `chr1 100 200` covers bases 100..199. This is not a
decision; BED says so.

**Normal intervals overlap iff `a.start < b.end AND b.start < a.end`.** Strict `<` on
both sides.

- **Bookended intervals do not overlap** **(measured: `chr1 100 200` vs `chr1 200 300`
  → `intersect -u` gives 0 hits)**.
- **Bookended intervals do merge at `-d 0`** **(measured: the same pair merges to
  `chr1 100 300`)**. Overlap and mergeability are different questions: `merge` joins
  features whose gap is `<= d`, and the gap here is 0.

**Zero-length intervals (`start == end`) are legal, and they do not use the predicate
above.** A zero-length feature at point *p* overlaps any feature where
`b.start <= p <= b.end` — inclusive on both ends. **(measured)**, for a zero-length `-a`
at `chr1 500 500`:

| `-b` feature | overlaps? | note |
|---|---|---|
| `400..600` (contains *p*) | yes | |
| `500..500` (identical) | yes | |
| `500..700` (starts at *p*) | **yes** | a normal interval would miss this |
| `300..500` (ends at *p*) | **yes** | a normal interval would miss this |
| `499..500` | yes | |
| `498..499` | no | |
| `501..502` | no | |

**Zero-length intervals also distort output coordinates, not just hit/miss.**
**(measured on the fixtures)**:

    bedtools merge -i <sorted a.bed>        ->  chr1 499 600
    bedtools subtract -a a.bed -b b.bed     ->  chr1 50 99   (a01 was 0..100)
                                            ->  chr1 101 180 (a02 was 100..200)

`merge` expands the zero-length `a07` at 500 **leftwards** to 499. `subtract` against the
zero-length `b02` at 100 removes **both** bases adjacent to the point: `a01` loses base
99, and `a02` loses base 100. A zero-length feature therefore behaves as `[p-1, p+1)`
when it is subtracted. Nobody predicts this correctly; encode what bedtools prints, add a
comment, move on.

**A zero-length feature at position 0 makes bedtools abort.** **(measured)** — as the
`-b` side it exits 1 having printed nothing:

    $ bedtools intersect -a data/b.bed -b data/a.bed
    ERROR: Received illegal bin number -1 from getBin call.
    Maximum values is: 2396745

`data/a.bed`'s `a12` is `chr2 0 0`, which widens to `[-1, 1)` under the rule above, and
bedtools' binning index rejects the negative coordinate. It is the same `[p-1, p+1)`
widening that `subtract` shows, surfacing as a crash instead of a coordinate.

It is **specific to the `-b` side**: the identical record as `-a` exits 0. So there is no
correct output to diff against, and a golden case that puts `a.bed` on `-b` has to exclude
`a12` — `tests/golden/intersect.sh` does exactly that, with the reason in a comment. This
is a bedtools limitation, not a behaviour to reproduce: `mytools` handles the record.

So there are **two overlap predicates in the tool**: strict for normal intervals,
inclusive for zero-length ones. Do not unify them and do not reason about zero-length
intervals from first principles — run bedtools and match it. `data/a.bed` contains three
zero-length features (`a07`, `a12`, `a16`) and `data/b.bed` two (`b02`, `b07`)
specifically so this is discovered early.

**Minimum overlap:** 1 bp by default. `intersect -f <frac>` requires that fraction of
the `-a` feature to be covered; see §4.

## 4. Flags per subcommand

**This table is the one place in this spec that is a choice rather than a measurement.**
"Everything bedtools has" is roughly 30 flags on `intersect` alone and is not finishable
in the time available, and scope creep is the main failure mode here. Flag *names and
meanings* match bedtools exactly; the *subset* is ours.

| Subcommand  | Flags in v1 | Notes |
|-------------|-------------|-------|
| `sort`      | `-i`, `-header` | No `-faidx`, no `-g`. |
| `merge`     | `-i`, `-d N`, `-s`, `-header` | Requires sorted input (§7). No `-c`/`-o`. |
| `intersect` | `-a`, `-b`, `-u`, `-v`, `-wa`, `-wb`, `-c`, `-s`, `-S`, `-f`, `-header` | `-u`, `-v`, `-c` and `-wb` are mutually exclusive. |
| `subtract`  | `-a`, `-b`, `-s`, `-S`, `-A`, `-header` | `-A` removes the whole `-a` feature on any overlap. |
| `closest`   | `-a`, `-b`, `-d`, `-t <first\|last\|all>`, `-s`, `-S`, `-header` | Requires both inputs sorted (§7). |

- **Strand flags**: `-s` = same strand only, `-S` = opposite strand only. Available
  because input is BED6+ (§2); on a record with fewer than 6 columns there is no strand,
  and such a record never matches `-s` or `-S`.
- **`-split` is deferred to v1.1.** BED12 blocks are accepted and preserved as opaque
  columns in v1; every subcommand treats a 12-column record as its single
  `start..end` span. This is bedtools' default behaviour — `-split` is what changes it —
  so v1 is still oracle-correct, just narrower.
- **`merge` does not sort for you**, matching bedtools (§7).
- **`closest` does not sort for you**, matching bedtools (§7).

## 5. Output

- Tab-separated, LF line endings, trailing newline on the final line.
- **Empty result: print nothing, exit 0** **(measured: `intersect` with no overlaps
  gives 0 bytes on stdout and exit 0)**.
- `sort`: all input columns preserved. Order is **chrom lexicographically, then start,
  then end** — so `chr1`, `chr17`, `chr2`, `chr7`, `chrX` **(measured; note `chr17`
  sorts before `chr2`, which is not natural order and is a common wrong guess)**.
- `merge`: **BED3 only** (`chrom start end`); input columns are dropped **(measured on
  sorted BED6 input)**. With `-s`, strand is retained as a fourth column.
- `intersect`: default prints **the intersected region carrying `-a`'s trailing
  columns**, once per overlapping `-b` feature **(measured: `a05` at `320..350` against
  two overlapping `-b` features yields `chr1 320 350 a05 40 +` and
  `chr1 340 360 a05 40 +`)**. `-wa` prints `-a`'s original interval instead; `-wb`
  appends the `-b` feature; `-u` prints each `-a` feature at most once; `-v` prints `-a`
  features with no overlap; `-c` appends an overlap count.
- `subtract`: `-a` features with `-b` regions removed. One feature may become two, or
  vanish entirely.
- `closest`: the `-a` feature, the nearest `-b` feature, and with `-d` a final distance
  column (0 when they overlap).
- **Input order is preserved for `intersect`, `subtract` and `closest`.** bedtools does
  not sort `-a` for you and neither do we.
- `-header` echoes the input's header lines before the data (§2).

## 6. Memory model

- `sort` holds the whole input in memory. Acceptable at our sizes.
- `merge` **streams**, assuming sorted input.
- `intersect`, `subtract` and `closest` load `-b` into memory, grouped by chrom and
  sorted by start, then **stream `-a`** and binary-search `-b`.
- **Target: inputs up to ~10^6 intervals.** This is sized for the scale stretch goal,
  where `bedtools bamtobed` over the HG002 neighbourhoods produces ~500,000 intervals.
- No mmap, no index files, no threads, no temporary files.
- **`intersect` must not be quadratic.** Being slower than bedtools is expected — it is
  C and we are not. Being O(n·m) is a bug.

## 7. Errors and exit codes

Errors go to **stderr**. stdout carries data only, because stdout gets piped.

**This section deliberately does not match bedtools.** What bedtools actually does
**(measured)**:

| Situation | bedtools | note |
|---|---|---|
| non-integer coordinate | **exit 134, SIGABRT, core dumped** | uncaught `std::invalid_argument` |
| `start > end` | exit 1 | clean message, names the line |
| `start > end` via `intersect` | exit 1 | but the message is `unable to determine types for file` — file-type sniffing fails first, so the error misreports the cause |
| missing input file | exit 1 | |
| unknown flag | exit 1 | `*****ERROR: Unrecognized parameter: --bogus *****` plus usage, all on stderr |
| no arguments at all | **exit 0, no output** | `-i` defaults to stdin, so this is an empty input rather than a usage error |
| unsorted input to `merge` / `closest` | exit 1 | `Error: Sorted input specified, but the file ... ` |

Copying that faithfully would mean shipping a crash on malformed input, exiting 0 on a
bare `mytools`, and having no distinction between a bad file and a bad flag. It also
contradicts `CLAUDE.md` and the already-closed issue #1, which requires `mytools` with no
arguments to print usage and exit 2. So:

| Situation | stderr message | exit |
|---|---|---|
| Success, including empty output | — | `0` |
| `--version`, `--help` | (to **stdout**) | `0` |
| Malformed line / non-integer coordinate | `a.bed:14: malformed coordinate 'notanumber'` | `1` |
| `start > end` | `a.bed:14: start > end (500 > 400)` | `1` |
| Unsorted input to `merge` or `closest` | `a.bed:9: input is not sorted (chr1:100 after chr1:900)` | `1` |
| Unknown flag / missing required argument | usage, to stderr | `2` |
| No arguments at all | usage, to stderr | `2` |
| Input file does not exist | `mytools: no such file: nope.bed` | `2` |
| Mutually exclusive flags (e.g. `-u` with `-v`) | names both flags | `2` |

**Data problems are 1. Caller problems are 2.** Every data message names the file and
the line number. Never abort, never core-dump: a malformed coordinate is a diagnosable
error, not a crash.

This costs nothing against the oracle, because the golden tests compare **stdout on
valid input** — where we match bedtools byte for byte — and never compare exit codes on
invalid input.

## 8. Correctness

Oracle: real `bedtools` on the files in `data/`. Non-negotiable.

Each of these must produce **byte-identical stdout** to its bedtools equivalent:

    mytools sort -i data/a.bed
    mytools sort -i data/b.bed
    mytools merge -i <sorted a.bed>
    mytools merge -d 10 -i <sorted a.bed>
    mytools merge -s -i <sorted a.bed>
    mytools intersect -a data/a.bed -b data/b.bed
    mytools intersect -u  -a data/a.bed -b data/b.bed
    mytools intersect -v  -a data/a.bed -b data/b.bed
    mytools intersect -wa -a data/a.bed -b data/b.bed
    mytools intersect -wb -a data/a.bed -b data/b.bed
    mytools intersect -c  -a data/a.bed -b data/b.bed
    mytools intersect -s  -a data/a.bed -b data/b.bed
    mytools intersect -S  -a data/a.bed -b data/b.bed
    mytools subtract -a data/a.bed -b data/b.bed
    mytools subtract -A -a data/a.bed -b data/b.bed
    mytools closest -a <sorted a.bed> -b <sorted b.bed>
    mytools closest -d -a <sorted a.bed> -b <sorted b.bed>
    mytools sort -i <gzipped a.bed>
    cat data/a.bed | mytools sort -i -

Also required: `mytools --version` prints a version and exits 0 (issue #1, done).

**Unit tests must run without bedtools**, one per edge case that required thought:
bookended, zero-length (all seven rows of the table in §3), nested, identical,
position 0, unsorted input, both strands, and the overlap predicate itself.

**Accepted deviations from bedtools:** error handling and exit codes only (§7). If you
find any deviation in stdout on valid input, it is a bug in `mytools`, not here.

## 9. Language and layout

- **Implementation language: Python** — see `CLAUDE.md`. Every subcommand and every
  test. Standard library only, no third-party runtime dependencies.
- **Entry point:** `mytools`, an executable file at the repo root with a
  `#!/usr/bin/env python3` shebang, symlinked onto `PATH`. Invoked as
  `mytools <subcommand> [options]`.
- **Tests live in `tests/`**, driven by `tests/run_golden.sh` for the golden diffs and a
  stdlib `unittest` suite for the unit tests. Both run on every push; see issue #3.
