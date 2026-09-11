"""Unit tests for `mytools intersect` (issue #7).

These need no bedtools -- that is the point of having them alongside the golden
tests. Run: python3 -m unittest discover -s tests

Where a test looks surprising it is because bedtools is surprising there. Each one
says which measurement it encodes; none of them encode what the behaviour ought to
be. SPEC.md §3 is the reference for the zero-length rows.
"""

import os
import subprocess
import sys
import unittest

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, REPO)

from mytoolslib import bed  # noqa: E402
from mytoolslib.cmd_intersect import (  # noqa: E402
    BinIndex,
    _f32,
    intersected_span,
    overlap_fraction,
    overlaps,
)

MYTOOLS = os.path.join(REPO, "mytools")
DATA = os.path.join(REPO, "data")


def rec(chrom, start, end, name="x", score="0", strand=None):
    extra = [name, score] + ([strand] if strand is not None else [])
    return bed.BedRecord(chrom, start, end, extra)


def run(args, stdin=b""):
    return subprocess.run(
        [sys.executable, MYTOOLS, *args], input=stdin, capture_output=True
    )


def rows(args):
    result = run(args)
    assert result.returncode == 0, result.stderr
    return result.stdout.decode().splitlines()


class TestOverlapPredicate(unittest.TestCase):
    """The comparison every off-by-one in this project lives in."""

    def test_partial_overlap(self):
        self.assertTrue(overlaps(rec("chr1", 100, 200), rec("chr1", 150, 250)))

    def test_bookended_does_not_overlap(self):
        # The off-by-one: strict `<`, not `<=`. a01/a02 meet at 100 in the fixtures.
        self.assertFalse(overlaps(rec("chr1", 100, 200), rec("chr1", 200, 300)))

    def test_bookended_the_other_way_round_does_not_overlap(self):
        self.assertFalse(overlaps(rec("chr1", 200, 300), rec("chr1", 100, 200)))

    def test_one_base_of_overlap_counts(self):
        self.assertTrue(overlaps(rec("chr1", 100, 200), rec("chr1", 199, 300)))

    def test_nested(self):
        # a06 inside a05; b08 inside a09.
        self.assertTrue(overlaps(rec("chr1", 300, 400), rec("chr1", 320, 350)))
        self.assertTrue(overlaps(rec("chr1", 320, 350), rec("chr1", 300, 400)))

    def test_identical(self):
        self.assertTrue(overlaps(rec("chr1", 700, 800), rec("chr1", 700, 800)))

    def test_disjoint(self):
        self.assertFalse(overlaps(rec("chr1", 100, 200), rec("chr1", 900, 1000)))

    def test_different_chroms_never_overlap(self):
        self.assertFalse(overlaps(rec("chr1", 100, 200), rec("chr2", 100, 200)))

    def test_position_zero(self):
        # a01 and b01 both start at 0; an interval there is ordinary, not a special
        # case, and 0..100 vs 0..50 overlaps like anything else.
        self.assertTrue(overlaps(rec("chr1", 0, 100), rec("chr1", 0, 50)))
        # b16 is chrX 0..5; the bookend rule still holds at the origin.
        self.assertFalse(overlaps(rec("chrX", 0, 5), rec("chrX", 5, 10)))


class TestZeroLengthTable(unittest.TestCase):
    """All seven measured rows of SPEC.md §3, for a zero-length -a at chr1 500.

    A zero-length feature does not use the strict predicate: it behaves as
    [p-1, p+1), which makes the test inclusive at both ends. Rows 3 and 4 are the
    ones a normal interval would miss.
    """

    def setUp(self):
        self.point = rec("chr1", 500, 500)

    def test_contained_in_b(self):
        self.assertTrue(overlaps(self.point, rec("chr1", 400, 600)))

    def test_identical_zero_length(self):
        self.assertTrue(overlaps(self.point, rec("chr1", 500, 500)))

    def test_b_starts_at_the_point(self):
        self.assertTrue(overlaps(self.point, rec("chr1", 500, 700)))

    def test_b_ends_at_the_point(self):
        self.assertTrue(overlaps(self.point, rec("chr1", 300, 500)))

    def test_b_ends_one_past_the_point(self):
        self.assertTrue(overlaps(self.point, rec("chr1", 499, 500)))

    def test_b_just_left_of_the_point(self):
        self.assertFalse(overlaps(self.point, rec("chr1", 498, 499)))

    def test_b_just_right_of_the_point(self):
        self.assertFalse(overlaps(self.point, rec("chr1", 501, 502)))

    def test_zero_length_b_against_a_normal_a(self):
        # Same rule seen from the other side: b02 at chr1 100 is a hit for a01
        # (0..100) and a02 (100..200), both of which end or start exactly at it.
        zero = rec("chr1", 100, 100)
        self.assertTrue(overlaps(rec("chr1", 0, 100), zero))
        self.assertTrue(overlaps(rec("chr1", 100, 200), zero))
        self.assertFalse(overlaps(rec("chr1", 0, 99), zero))
        self.assertFalse(overlaps(rec("chr1", 101, 200), zero))

    def test_zero_length_at_position_zero(self):
        # a12 is chr2 0 0. Widening runs it to -1, which is fine for us -- real
        # bedtools rejects the negative bin if such a record is in -b, so the
        # golden suite only exercises this direction.
        self.assertTrue(overlaps(rec("chr2", 0, 0), rec("chr2", 0, 10)))
        self.assertFalse(overlaps(rec("chr2", 0, 0), rec("chr2", 1, 10)))


class TestIntersectedSpan(unittest.TestCase):
    """Zero-length features distort the printed coordinates too (SPEC.md §3)."""

    def test_normal_case_is_the_plain_intersection(self):
        span = intersected_span(rec("chr1", 300, 400), rec("chr1", 320, 350))
        self.assertEqual(span, (320, 350))

    def test_zero_length_b_widens_the_printed_region(self):
        # Measured: a01 (0..100) against b02 (100 100) prints chr1 99 100.
        b02 = rec("chr1", 100, 100)
        self.assertEqual(intersected_span(rec("chr1", 0, 100), b02), (99, 100))
        # And a02 (100..200) against the same b02 prints chr1 100 101.
        self.assertEqual(intersected_span(rec("chr1", 100, 200), b02), (100, 101))

    def test_zero_length_b_in_the_middle_of_a_prints_both_sides(self):
        # Measured: 100..200 against a point at 150 prints chr1 149 151.
        a = rec("chr1", 100, 200)
        self.assertEqual(intersected_span(a, rec("chr1", 150, 150)), (149, 151))

    def test_zero_length_a_prints_its_own_point(self):
        # Measured on every zero-length -a in the fixtures: a07 against b07 prints
        # chr1 500 500, not the widened 499 501, whatever the -b feature is.
        point = rec("chr1", 500, 500)
        for other in (rec("chr1", 400, 600), rec("chr1", 500, 500),
                      rec("chr1", 300, 500)):
            self.assertEqual(intersected_span(point, other), (500, 500))


class TestFraction(unittest.TestCase):
    def test_exactly_at_the_threshold_passes(self):
        # a21 (chrX 30..40) against b17 (25..35): 5 bases of 10. -f 0.5 keeps it.
        fraction = overlap_fraction(rec("chrX", 30, 40), rec("chrX", 25, 35))
        self.assertGreaterEqual(fraction, 0.5)
        self.assertLess(fraction, 0.500001)

    def test_just_under_the_threshold_fails(self):
        # a05 (300..400) against b05 (320..350): 30 of 100.
        fraction = overlap_fraction(rec("chr1", 300, 400), rec("chr1", 320, 350))
        self.assertLess(fraction, 0.5)

    def test_full_containment_is_one(self):
        nested = overlap_fraction(rec("chr1", 320, 350), rec("chr1", 300, 400))
        self.assertEqual(nested, 1.0)

    def test_zero_length_a_has_a_denominator_of_two(self):
        # Widening applies here too, so a point's fraction is over 2 bases, not 0 --
        # measured: a12 (chr2 0 0) against b10 (0..10) passes -f 0.5, fails -f 1.0,
        # because only one of its two widened bases is covered.
        self.assertEqual(overlap_fraction(rec("chr2", 0, 0), rec("chr2", 0, 10)), 0.5)
        # Both bases covered when -b spans the point: a07 against b07, -f 1.0 keeps it.
        both = overlap_fraction(rec("chr1", 500, 500), rec("chr1", 500, 500))
        self.assertEqual(both, 1.0)

    def test_compared_in_c_float_precision(self):
        # bedtools holds -f in a 32-bit float, and 7/10 collapses onto 0.70000001
        # there. Measured: -f 0.70000001 still keeps an overlap of 7 bases in 10.
        # Both sides of the comparison have to be narrowed, which is what run()
        # does with the -f value; in plain double precision the same comparison is
        # False and we would drop a row bedtools keeps.
        seven_of_ten = overlap_fraction(rec("chr1", 0, 10), rec("chr1", 0, 7))
        self.assertGreaterEqual(seven_of_ten, _f32(0.70000001))
        self.assertLess(seven_of_ten, 0.70000001)  # the double-precision answer


class TestHitOrder(unittest.TestCase):
    """Row order is part of the answer, and it is not -b's start order."""

    def test_hits_come_back_in_bedtools_bin_order(self):
        # Measured: for a02 (chr1 100..200), bedtools prints b03 (at 180) before
        # b02 (the point at 100). Both land in the same finest-level bin, where
        # file order wins -- b03 is the first line of b.bed. Ordering hits by -b
        # start would put b02 first and fail the golden diff.
        index = BinIndex()
        b03 = rec("chr1", 180, 220, "b03")
        b02 = rec("chr1", 100, 100, "b02")
        for record in (b03, rec("chr1", 0, 50, "b01"), b02):
            index.add(record)
        index.finish()
        a02 = rec("chr1", 100, 200)
        hits = [r for r in index.candidates(a02) if overlaps(a02, r)]
        self.assertEqual([r.extra[0] for r in hits], ["b03", "b02"])

    def test_finer_bins_come_before_coarser_ones(self):
        # Measured: a 10-base -a inside three nested -b features reports the
        # smallest first, because the bin index visits the finest level first --
        # again regardless of file order or start coordinate.
        index = BinIndex()
        for record in (
            rec("chr1", 0, 100000, "big"),
            rec("chr1", 500, 600, "small"),
            rec("chr1", 0, 1000000, "huge"),
        ):
            index.add(record)
        index.finish()
        hits = list(index.candidates(rec("chr1", 550, 560)))
        self.assertEqual([r.extra[0] for r in hits], ["small", "big", "huge"])


class TestIndexIsNotAScan(unittest.TestCase):
    """SPEC.md §6: -b is indexed, not scanned. A nested scan is the bug to catch."""

    def test_lookup_examines_a_handful_of_records_not_all_of_them(self):
        index = BinIndex()
        for i in range(20000):
            index.add(rec("chr1", i * 1000, i * 1000 + 100, f"b{i}"))
        index.finish()
        query = rec("chr1", 5_000_000, 5_000_100)
        examined = list(index.candidates(query))
        # The finest bin level is 16kb wide, so it holds ~16 of these 1kb-spaced
        # features and a lookup can never see many more than that. 20000 would mean
        # every -b record was examined for every -a record, which is the quadratic
        # implementation this test exists to reject.
        self.assertLessEqual(len(examined), 40, "looks like a linear scan over -b")
        hits = [r.extra[0] for r in examined if overlaps(query, r)]
        self.assertEqual(hits, ["b5000"])

    def test_chrom_missing_from_b_costs_nothing(self):
        index = BinIndex()
        index.add(rec("chr1", 100, 200))
        index.finish()
        self.assertEqual(list(index.candidates(rec("chr9", 100, 200))), [])


class TestEndToEnd(unittest.TestCase):
    """Driven through the real fixtures, still without bedtools."""

    def setUp(self):
        self.a = os.path.join(DATA, "a.bed")
        self.b = os.path.join(DATA, "b.bed")

    def test_default_prints_intersected_region_once_per_b_feature(self):
        out = rows(["intersect", "-a", self.a, "-b", self.b])
        self.assertIn("chr1\t320\t350\ta05\t40\t+", out)
        self.assertIn("chr1\t340\t360\ta05\t40\t+", out)

    def test_input_order_of_a_is_preserved(self):
        # a.bed starts with a05, then a01 -- we do not sort, and neither does bedtools.
        out = rows(["intersect", "-wa", "-a", self.a, "-b", self.b])
        self.assertEqual(out[0].split("\t")[3], "a05")
        self.assertEqual(out[2].split("\t")[3], "a01")

    def test_u_prints_each_a_feature_at_most_once(self):
        out = rows(["intersect", "-u", "-a", self.a, "-b", self.b])
        names = [line.split("\t")[3] for line in out]
        self.assertEqual(len(names), len(set(names)))
        self.assertIn("a05", names)

    def test_v_is_the_complement_of_u(self):
        def names(flag):
            out = rows(["intersect", flag, "-a", self.a, "-b", self.b])
            return {line.split("\t")[3] for line in out}

        hit, missed = names("-u"), names("-v")
        self.assertEqual(hit & missed, set())
        self.assertEqual(len(hit | missed), 22)  # every feature in a.bed, once

    def test_c_reports_zero_rather_than_omitting_the_row(self):
        out = rows(["intersect", "-c", "-a", self.a, "-b", self.b])
        self.assertEqual(len(out), 22)
        counts = {line.split("\t")[3]: line.split("\t")[6] for line in out}
        self.assertEqual(counts["a05"], "2")
        self.assertEqual(counts["a14"], "0")

    def test_wb_appends_the_b_feature_unwidened(self):
        out = rows(["intersect", "-wb", "-a", self.a, "-b", self.b])
        # b02 is the zero-length point at 100. The printed region is widened; the
        # -b columns are not.
        self.assertIn("chr1\t99\t100\ta01\t10\t+\tchr1\t100\t100\tb02\t0\t-", out)

    def test_strand_flags_split_a03_from_a04(self):
        # a03 and a04 are identical but for strand; b03 (180..220) is '+'.
        same = rows(["intersect", "-s", "-a", self.a, "-b", self.b])
        opposite = rows(["intersect", "-S", "-a", self.a, "-b", self.b])
        self.assertIn("chr1\t180\t220\ta03\t30\t+", same)
        self.assertNotIn("chr1\t180\t220\ta04\t30\t-", same)
        self.assertIn("chr1\t180\t220\ta04\t30\t-", opposite)
        self.assertNotIn("chr1\t180\t220\ta03\t30\t+", opposite)

    def test_record_without_a_strand_matches_neither_s_nor_S(self):
        bed3 = "chr1\t150\t250\n"
        for flag in ("-s", "-S"):
            result = run(["intersect", flag, "-a", "-", "-b", self.b], stdin=bed3.encode())
            self.assertEqual(result.returncode, 0)
            self.assertEqual(result.stdout, b"", f"{flag} matched a record with no strand")

    def test_f_boundary_is_inclusive(self):
        # a21 (chrX 30..40) keeps exactly 5 of its 10 bases against b17.
        kept = rows(["intersect", "-f", "0.5", "-a", self.a, "-b", self.b])
        self.assertIn("chrX\t30\t35\ta21\t90\t-", kept)
        dropped = rows(["intersect", "-f", "0.500001", "-a", self.a, "-b", self.b])
        self.assertNotIn("chrX\t30\t35\ta21\t90\t-", dropped)

    def test_chrom_in_b_but_not_a_is_ignored(self):
        # chr3/b19 exists only in b.bed; it must not appear or disturb anything.
        out = rows(["intersect", "-wb", "-a", self.a, "-b", self.b])
        self.assertFalse([line for line in out if "b19" in line])

    def test_header_echoed_only_when_asked(self):
        text = "# comment\ntrack name=x\nchr1\t150\t250\tz\t0\t+\n".encode()
        with_header = run(["intersect", "-header", "-a", "-", "-b", self.b], stdin=text)
        self.assertTrue(with_header.stdout.decode().startswith("# comment\ntrack name=x\n"))
        without = run(["intersect", "-a", "-", "-b", self.b], stdin=text)
        self.assertNotIn(b"#", without.stdout)

    def test_stdin_matches_the_file(self):
        with open(self.a, "rb") as handle:
            piped = run(["intersect", "-a", "-", "-b", self.b], stdin=handle.read())
        from_file = rows(["intersect", "-a", self.a, "-b", self.b])
        self.assertEqual(piped.stdout.decode().splitlines(), from_file)


class TestUsageErrors(unittest.TestCase):
    """Exit 2 for caller problems (SPEC.md §7). bedtools says 1; we deviate."""

    def test_mutually_exclusive_flags_name_both(self):
        for pair in (("-u", "-v"), ("-u", "-c"), ("-v", "-wb"), ("-c", "-wb")):
            result = run(["intersect", *pair, "-a", os.path.join(DATA, "a.bed"),
                          "-b", os.path.join(DATA, "b.bed")])
            self.assertEqual(result.returncode, 2, f"{pair} was accepted")
            self.assertIn(pair[0].encode(), result.stderr)
            self.assertIn(pair[1].encode(), result.stderr)
            self.assertEqual(result.stdout, b"")

    def test_s_and_S_are_exclusive(self):
        result = run(["intersect", "-s", "-S", "-a", os.path.join(DATA, "a.bed"),
                      "-b", os.path.join(DATA, "b.bed")])
        self.assertEqual(result.returncode, 2)

    def test_missing_required_input(self):
        result = run(["intersect", "-a", os.path.join(DATA, "a.bed")])
        self.assertEqual(result.returncode, 2)
        self.assertIn(b"-b", result.stderr)

    def test_f_out_of_range(self):
        for value in ("0", "0.0", "1.5", "-0.2", "notanumber"):
            result = run(["intersect", "-f", value, "-a", os.path.join(DATA, "a.bed"),
                          "-b", os.path.join(DATA, "b.bed")])
            self.assertEqual(result.returncode, 2, f"-f {value} was accepted")

    def test_missing_file_is_exit_2(self):
        result = run(["intersect", "-a", os.path.join(DATA, "nope.bed"),
                      "-b", os.path.join(DATA, "b.bed")])
        self.assertEqual(result.returncode, 2)
        self.assertIn(b"no such file", result.stderr)

    def test_malformed_data_is_exit_1_on_stderr(self):
        result = run(["intersect", "-a", "-", "-b", os.path.join(DATA, "b.bed")],
                     stdin=b"chr1\tnotanumber\t200\n")
        self.assertEqual(result.returncode, 1)
        self.assertEqual(result.stdout, b"")
        self.assertIn(b"malformed coordinate", result.stderr)
        self.assertNotIn(b"Traceback", result.stderr)


if __name__ == "__main__":
    unittest.main()
