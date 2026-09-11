"""Unit tests for the shared BED I/O and dispatch layer (issue #4).

These need no bedtools -- that is the point of having them alongside the golden
tests. Run: python3 -m unittest discover -s tests
"""

import gzip
import io
import os
import subprocess
import sys
import tempfile
import unittest

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, REPO)

from mytoolslib import bed  # noqa: E402
from mytoolslib.bed import DataError, UsageError  # noqa: E402

MYTOOLS = os.path.join(REPO, "mytools")


def run(args, stdin=b""):
    return subprocess.run(
        [sys.executable, MYTOOLS, *args], input=stdin, capture_output=True
    )


class TestParsing(unittest.TestCase):
    def test_bed3(self):
        rec = bed.parse_line("chr1\t100\t200", "t", 1)
        self.assertEqual((rec.chrom, rec.start, rec.end), ("chr1", 100, 200))
        self.assertEqual(rec.extra, [])
        self.assertIsNone(rec.strand)

    def test_bed6_columns_preserved(self):
        rec = bed.parse_line("chr1\t100\t200\tname\t40\t+", "t", 1)
        self.assertEqual(rec.extra, ["name", "40", "+"])
        self.assertEqual(rec.strand, "+")
        self.assertEqual(rec.line(), "chr1\t100\t200\tname\t40\t+")

    def test_bed12_columns_preserved_opaquely(self):
        line = "chr1\t100\t200\tn\t0\t+\t100\t200\t0,0,0\t2\t10,10\t0,90"
        rec = bed.parse_line(line, "t", 1)
        self.assertEqual(len(rec.fields()), 12)
        self.assertEqual(rec.line(), line)

    def test_zero_length_is_legal(self):
        rec = bed.parse_line("chr1\t500\t500", "t", 1)
        self.assertEqual((rec.start, rec.end), (500, 500))

    def test_position_zero(self):
        rec = bed.parse_line("chr1\t0\t100", "t", 1)
        self.assertEqual(rec.start, 0)

    def test_malformed_coordinate_names_file_and_line(self):
        with self.assertRaises(DataError) as ctx:
            bed.parse_line("chr1\tnotanumber\t200", "a.bed", 14)
        self.assertIn("a.bed:14", str(ctx.exception))
        self.assertEqual(ctx.exception.exit_code, 1)

    def test_start_greater_than_end(self):
        with self.assertRaises(DataError) as ctx:
            bed.parse_line("chr1\t500\t400", "a.bed", 9)
        self.assertIn("start > end (500 > 400)", str(ctx.exception))

    def test_too_few_columns(self):
        with self.assertRaises(DataError):
            bed.parse_line("chr1\t100", "t", 1)


class TestHeaders(unittest.TestCase):
    def test_prefixes_recognised(self):
        for line in ("# c", "track name=x", "browser position chr1:1-2"):
            self.assertTrue(bed.is_header(line))
        self.assertFalse(bed.is_header("chr1\t1\t2"))


class TestReading(unittest.TestCase):
    def setUp(self):
        self.dir = tempfile.mkdtemp()
        self.plain = os.path.join(self.dir, "a.bed")
        with open(self.plain, "w") as handle:
            handle.write("# comment\ntrack name=x\n\nchr1\t100\t200\tkeep\n")
        self.gz = os.path.join(self.dir, "a.bed.gz")
        with gzip.open(self.gz, "wt") as handle:
            handle.write("# comment\ntrack name=x\n\nchr1\t100\t200\tkeep\n")

    def test_headers_and_blanks_skipped(self):
        recs = list(bed.read_bed(self.plain))
        self.assertEqual(len(recs), 1)
        self.assertEqual(recs[0].extra, ["keep"])

    def test_on_header_sees_them_in_order(self):
        seen = []
        list(bed.read_bed(self.plain, on_header=seen.append))
        self.assertEqual(seen, ["# comment", "track name=x"])

    def test_gzip_and_plain_agree(self):
        plain = [r.line() for r in bed.read_bed(self.plain)]
        gzipped = [r.line() for r in bed.read_bed(self.gz)]
        self.assertEqual(plain, gzipped)

    def test_missing_file_is_usage_error(self):
        with self.assertRaises(UsageError) as ctx:
            list(bed.read_bed(os.path.join(self.dir, "nope.bed")))
        self.assertEqual(ctx.exception.exit_code, 2)


class TestSortKeyAndSortedness(unittest.TestCase):
    def test_chrom_order_is_lexicographic(self):
        chroms = ["chr7", "chr17", "chr2", "chrX", "chr1"]
        recs = [bed.BedRecord(c, 1, 2) for c in chroms]
        ordered = [r.chrom for r in sorted(recs, key=bed.sort_key)]
        self.assertEqual(ordered, ["chr1", "chr17", "chr2", "chr7", "chrX"])

    def test_chr17_before_chr2(self):
        self.assertLess(
            bed.sort_key(bed.BedRecord("chr17", 1, 2)),
            bed.sort_key(bed.BedRecord("chr2", 1, 2)),
        )

    def test_ties_broken_by_end(self):
        a = bed.BedRecord("chr1", 100, 150)
        b = bed.BedRecord("chr1", 100, 200)
        self.assertLess(bed.sort_key(a), bed.sort_key(b))

    def test_require_sorted_passes_sorted_input(self):
        recs = [bed.BedRecord("chr1", 0, 10), bed.BedRecord("chr1", 10, 20)]
        self.assertEqual(len(list(bed.require_sorted(recs))), 2)

    def test_require_sorted_rejects_unsorted(self):
        recs = [
            bed.BedRecord("chr1", 900, 1000, source="a.bed", lineno=1),
            bed.BedRecord("chr1", 100, 200, source="a.bed", lineno=2),
        ]
        with self.assertRaises(DataError) as ctx:
            list(bed.require_sorted(recs))
        self.assertIn("not sorted", str(ctx.exception))
        self.assertEqual(ctx.exception.exit_code, 1)


class TestCli(unittest.TestCase):
    def test_version(self):
        result = run(["--version"])
        self.assertEqual(result.returncode, 0)
        self.assertTrue(result.stdout.decode().startswith("mytools "))

    def test_no_args_usage_on_stderr_exit_2(self):
        result = run([])
        self.assertEqual(result.returncode, 2)
        self.assertEqual(result.stdout, b"")
        self.assertIn(b"usage:", result.stderr)

    def test_help_on_stdout_exit_0(self):
        result = run(["--help"])
        self.assertEqual(result.returncode, 0)
        self.assertIn(b"usage:", result.stdout)

    def test_unknown_command_exit_2(self):
        result = run(["frobnicate"])
        self.assertEqual(result.returncode, 2)

    def test_unimplemented_subcommand_exit_2_not_silent(self):
        result = run(["sort", "-i", "-"], stdin=b"chr1\t1\t2\n")
        self.assertEqual(result.returncode, 2)
        self.assertIn(b"not implemented yet", result.stderr)
        self.assertEqual(result.stdout, b"")

    def test_no_traceback_leaks(self):
        for args in ([], ["frobnicate"], ["sort"]):
            result = run(args)
            self.assertNotIn(b"Traceback", result.stderr, f"traceback from {args}")


class TestFlagParsing(unittest.TestCase):
    def test_bool_and_value_flags(self):
        from mytoolslib import cli

        parsed, rest = cli.parse_flags(
            ["-u", "-d", "10", "file.bed"], flags=("-u",), options=("-d",)
        )
        self.assertTrue(parsed["-u"])
        self.assertEqual(parsed["-d"], "10")
        self.assertEqual(rest, ["file.bed"])

    def test_unknown_option_is_usage_error(self):
        from mytoolslib import cli

        with self.assertRaises(UsageError):
            cli.parse_flags(["--bogus"], flags=("-u",))

    def test_missing_value_is_usage_error(self):
        from mytoolslib import cli

        with self.assertRaises(UsageError):
            cli.parse_flags(["-d"], options=("-d",))

    def test_dash_is_not_an_option(self):
        from mytoolslib import cli

        _, rest = cli.parse_flags(["-"], flags=())
        self.assertEqual(rest, ["-"])

    def test_required_missing(self):
        from mytoolslib import cli

        with self.assertRaises(UsageError):
            cli.parse_flags([], options=("-i",), required=("-i",))

    def test_exclusive(self):
        from mytoolslib import cli

        with self.assertRaises(UsageError):
            cli.exclusive({"-u": True, "-v": True}, "-u", "-v")
        cli.exclusive({"-u": True, "-v": False}, "-u", "-v")

    def test_int_option_rejects_junk(self):
        from mytoolslib import cli

        self.assertEqual(cli.int_option({"-d": "10"}, "-d", 0), 10)
        self.assertEqual(cli.int_option({"-d": None}, "-d", 0), 0)
        with self.assertRaises(UsageError):
            cli.int_option({"-d": "ten"}, "-d", 0)


if __name__ == "__main__":
    unittest.main()
