"""`mytools intersect` -- report features shared between -a and -b.

    mytools intersect -a <file|-> -b <file> [-u|-v|-c|-wb] [-wa] [-s|-S] [-f FRAC] [-header]

`-b` is loaded into memory, grouped by chrom and indexed by bin; `-a` is streamed
and the index is binary-searched (SPEC.md §6). Nothing here is quadratic.

Three pieces of bedtools behaviour are measured, not reasoned about. All three were
established by running bedtools v2.31.1 against `data/`:

1. Normal features overlap on strict `<` (SPEC.md §3), so bookended features miss.
2. A zero-length feature at point *p* is widened to `[p-1, p+1)` before any overlap
   question is asked -- see `_span` below.
3. When one `-a` feature hits several `-b` features, the *order* of the output rows
   is the order bedtools' bin index happens to visit them in -- see `BinIndex`.
"""

from __future__ import annotations

import struct
import sys
from bisect import bisect_left

from . import cli
from .bed import UsageError, read_bed

FLAGS = ("-u", "-v", "-wa", "-wb", "-c", "-s", "-S", "-header")
OPTIONS = ("-a", "-b", "-f")

# bedtools' own default: "1E-9 (i.e., 1bp)". Keeping the literal default rather than
# a "no filter" short circuit means one code path, and it is what bedtools compares.
DEFAULT_FRACTION = 1e-9


# ---------------------------------------------------------------------------
# Interval semantics
# ---------------------------------------------------------------------------


def _span(record):
    """The half-open span bedtools actually compares, widening zero-length features.

    A feature with `start == end` is a point *p*. bedtools widens it to `[p-1, p+1)`
    as it reads the file and every later decision -- overlap, overlap length, which
    bin it lands in -- uses the widened span. That single rule reproduces all seven
    rows of the measured table in SPEC.md §3, including the two that a strict
    half-open predicate would miss: a point at 500 hits both `500..700` and
    `300..500`.

    Normal features are returned unchanged, so they keep the strict `<` predicate
    and bookended features still do not overlap.
    """
    if record.start == record.end:
        return record.start - 1, record.end + 1
    return record.start, record.end


def overlaps(a, b):
    """True if BED records `a` and `b` overlap, by bedtools' rules.

    Strict `<` on both sides (SPEC.md §3): `chr1 100 200` and `chr1 200 300` are
    bookended and do **not** overlap. Zero-length features are widened first by
    `_span`, which is what makes their predicate effectively inclusive.
    """
    if a.chrom != b.chrom:
        return False
    a_start, a_end = _span(a)
    b_start, b_end = _span(b)
    return a_start < b_end and b_start < a_end


def _f32(value):
    """Round a Python float to C `float` precision.

    bedtools holds `-f` and the computed overlap fraction in 32-bit floats, and the
    comparison is close enough to the boundary for that to show. Measured:
    `-f 0.70000001` still accepts an overlap of 7 bases out of 10, because both
    sides collapse to the same `float`. Comparing in double precision rejects it and
    we would differ from the oracle.
    """
    return struct.unpack("f", struct.pack("f", value))[0]


def overlap_fraction(a, b):
    """Fraction of `-a` feature `a` covered by `b`, as bedtools computes it.

    Both spans are widened first, so a zero-length `-a` has a denominator of 2 and a
    one-base overlap reads as 0.5 -- measured: `a12` at `chr2 0 0` against
    `chr2 0 10` passes `-f 0.5` and fails `-f 1.0`.
    """
    a_start, a_end = _span(a)
    b_start, b_end = _span(b)
    covered = min(a_end, b_end) - max(a_start, b_start)
    return _f32(_f32(covered) / _f32(a_end - a_start))


def intersected_span(a, b):
    """The `(start, end)` that default output prints for a hit of `b` on `a`.

    The widened `-b` span is what gets clipped -- measured: `a01` at `chr1 0 100`
    against the zero-length `b02` at 100 prints `chr1 99 100`, one base that a
    literal reading of `b02` does not contain.

    A zero-length `-a` is the exception: it prints its own point back, never the
    widened span. Measured on every zero-length `-a` in the fixtures -- `a07`
    against `b07` prints `chr1 500 500`, not `chr1 499 501`.
    """
    if a.start == a.end:
        return a.start, a.end
    b_start, b_end = _span(b)
    return max(a.start, b_start), min(a.end, b_end)


def _strand_matches(a, b, same, opposite):
    """Apply `-s` / `-S`.

    A record narrower than BED6 has no strand, and so does one whose strand column
    is neither `+` nor `-`. Either way it matches neither `-s` nor `-S` -- measured:
    a BED3 `-b` gives no output at all under `-s`, and two `.` strands do not count
    as the same strand.
    """
    if not (same or opposite):
        return True
    a_strand, b_strand = a.strand, b.strand
    if a_strand not in ("+", "-") or b_strand not in ("+", "-"):
        return False
    return (a_strand == b_strand) if same else (a_strand != b_strand)


# ---------------------------------------------------------------------------
# The -b index
# ---------------------------------------------------------------------------

# bedtools' BinTree constants, copied verbatim: a hierarchy of bins, 16kb at the
# finest level and 8x wider at each of the six levels above it. A feature is filed
# in the smallest single bin that contains it.
#
# We copy the scheme rather than binary-searching a start-sorted list because the
# bin layout *is* the output order. When one -a feature hits several -b features,
# bedtools emits them in the order its index visits them: finest level first, bins
# in increasing order, file order within a bin. That is not the same as ordering by
# -b start -- measured: `a02` reports `b03` (at 180) before `b02` (at 100), because
# both share a bin and `b03` comes first in b.bed. Sorting hits by start there would
# be a one-line golden-test failure.
BIN_OFFSETS = (37359, 4681, 585, 73, 9, 1, 0)
BIN_FIRST_SHIFT = 14
BIN_NEXT_SHIFT = 3


def bin_of(start, end):
    """The single bin holding `[start, end)`, smallest level that contains it."""
    low = start >> BIN_FIRST_SHIFT
    high = (end - 1) >> BIN_FIRST_SHIFT
    for offset in BIN_OFFSETS:
        if low == high:
            return offset + low
        low >>= BIN_NEXT_SHIFT
        high >>= BIN_NEXT_SHIFT
    return 0


class BinIndex:
    """`-b` in memory: chrom -> bin -> the features filed in that bin, in file order.

    Lookup binary-searches the sorted bin numbers of the chromosome, so empty
    regions of a sparse `-b` cost nothing and a query touches only populated bins.
    """

    def __init__(self):
        self._chroms = {}
        self._sorted_bins = {}

    def add(self, record):
        start, end = _span(record)
        bins = self._chroms.setdefault(record.chrom, {})
        bins.setdefault(bin_of(start, end), []).append(record)

    def finish(self):
        """Freeze the index: one sorted bin-number list per chrom, to bisect."""
        self._sorted_bins = {
            chrom: sorted(bins) for chrom, bins in self._chroms.items()
        }

    def candidates(self, record):
        """Yield the `-b` features whose bin could contain a hit for `record`.

        Order is bedtools' index order, which the caller must not disturb.
        """
        bins = self._chroms.get(record.chrom)
        if bins is None:
            # A chrom in -a but not -b. The reverse -- chr3/b19 in the fixtures --
            # costs nothing: nothing ever asks for it.
            return
        numbers = self._sorted_bins[record.chrom]
        start, end = _span(record)
        low = start >> BIN_FIRST_SHIFT
        high = (end - 1) >> BIN_FIRST_SHIFT
        for offset in BIN_OFFSETS:
            first, last = low + offset, high + offset
            index = bisect_left(numbers, first)
            while index < len(numbers) and numbers[index] <= last:
                yield from bins[numbers[index]]
                index += 1
            low >>= BIN_NEXT_SHIFT
            high >>= BIN_NEXT_SHIFT


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------


def run(args):
    parsed, rest = cli.parse_flags(
        args, flags=FLAGS, options=OPTIONS, required=("-a", "-b")
    )
    if rest:
        raise UsageError(f"mytools intersect: unexpected argument: {rest[0]}")
    cli.exclusive(parsed, "-u", "-v", "-c", "-wb")
    cli.exclusive(parsed, "-s", "-S")
    if parsed["-a"] == "-" and parsed["-b"] == "-":
        raise UsageError("mytools intersect: only one of -a and -b may be stdin")

    fraction = cli.float_option(parsed, "-f", DEFAULT_FRACTION)
    if not 0.0 < fraction <= 1.0:
        raise UsageError(
            f"mytools intersect: -f must be in the range (0.0, 1.0], got {fraction}"
        )
    threshold = _f32(fraction)
    same, opposite = parsed["-s"], parsed["-S"]

    index = BinIndex()
    for record in read_bed(parsed["-b"]):
        index.add(record)
    index.finish()

    def hits(record):
        """The `-b` features that really hit `record`, in bedtools' output order."""
        for other in index.candidates(record):
            if (
                overlaps(record, other)
                and _strand_matches(record, other, same, opposite)
                and overlap_fraction(record, other) >= threshold
            ):
                yield other

    out = sys.stdout
    # -header echoes -a's header lines ahead of the data (SPEC.md §2). -b's are
    # dropped, which is why only this read_bed gets a callback.
    on_header = cli.header_echo(parsed["-header"], out)
    write = out.write

    count_only, any_only, invert = parsed["-c"], parsed["-u"], parsed["-v"]
    write_a, write_b = parsed["-wa"], parsed["-wb"]

    for record in read_bed(parsed["-a"], on_header):
        if count_only:
            write(f"{record.line()}\t{sum(1 for _ in hits(record))}\n")
        elif any_only or invert:
            # Both only care whether there was a hit, so stop at the first one.
            found = next(hits(record), None) is not None
            if found != invert:
                write(record.line() + "\n")
        else:
            for other in hits(record):
                if write_a:
                    left = record.line()
                else:
                    start, end = intersected_span(record, other)
                    left = record.with_span(start, end).line()
                write(f"{left}\t{other.line()}\n" if write_b else left + "\n")
    return 0
