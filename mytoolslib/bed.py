"""BED records and input handling.

Deliberately contains no interval logic -- no overlap, no merging, no distance.
Those live in the cmd_* modules, one per subcommand.

BED is 0-based half-open: `chr1 100 200` covers bases 100..199. See SPEC.md §3.
"""

from __future__ import annotations

import gzip
import io
import sys
from dataclasses import dataclass, field

GZIP_MAGIC = b"\x1f\x8b"

# bedtools treats lines starting with any of these as header lines. It is a prefix
# test, not a word test, so a chrom literally named "trackwise" would be misread --
# bedtools has the same flaw and we match it rather than improve on it.
HEADER_PREFIXES = ("#", "track", "browser")


class MyToolsError(Exception):
    """Base for errors we report cleanly. `exit_code` is what main() returns."""

    exit_code = 1


class DataError(MyToolsError):
    """Something wrong with the input data. Exit 1 (SPEC.md §7)."""

    exit_code = 1


class UsageError(MyToolsError):
    """Something wrong with how mytools was called. Exit 2 (SPEC.md §7)."""

    exit_code = 2


@dataclass
class BedRecord:
    """One BED feature. Columns beyond the first three are carried opaquely.

    BED12 block columns are just `extra` entries; v1 never interprets them
    (SPEC.md §4 -- `-split` is deferred).
    """

    chrom: str
    start: int
    end: int
    extra: list = field(default_factory=list)
    source: str = ""
    lineno: int = 0

    @property
    def strand(self):
        """The strand column, or None for records narrower than BED6.

        A record with no strand never matches `-s` or `-S` (SPEC.md §4).
        """
        return self.extra[2] if len(self.extra) > 2 else None

    def fields(self):
        return [self.chrom, str(self.start), str(self.end), *self.extra]

    def line(self):
        return "\t".join(self.fields())

    def with_span(self, start, end):
        """A copy at new coordinates, carrying the same trailing columns."""
        return BedRecord(self.chrom, start, end, list(self.extra), self.source, self.lineno)


def is_header(line):
    return line.startswith(HEADER_PREFIXES)


def parse_line(line, source, lineno):
    """Parse one data line into a BedRecord, or raise DataError naming file:line."""
    fields = line.split("\t")
    if len(fields) < 3:
        raise DataError(
            f"{source}:{lineno}: need at least 3 tab-separated columns, got {len(fields)}"
        )
    coords = []
    for value in (fields[1], fields[2]):
        try:
            coords.append(int(value))
        except ValueError:
            raise DataError(f"{source}:{lineno}: malformed coordinate {value!r}") from None
    start, end = coords
    if start < 0 or end < 0:
        raise DataError(f"{source}:{lineno}: negative coordinate ({start}, {end})")
    if start > end:
        raise DataError(f"{source}:{lineno}: start > end ({start} > {end})")
    return BedRecord(fields[0], start, end, fields[3:], source, lineno)


def _open_text(path):
    """Open a path or '-' as text, transparently un-gzipping.

    Compression is detected by magic bytes, not by filename (SPEC.md §2), so a
    gzipped file called `a.bed` works and a plain file called `a.bed.gz` also works.

    Returns (text_stream, underlying_binary_handle, label).
    """
    if path == "-":
        raw = sys.stdin.buffer
        label = "stdin"
    else:
        try:
            raw = open(path, "rb")
        except FileNotFoundError:
            raise UsageError(f"mytools: no such file: {path}") from None
        except IsADirectoryError:
            raise UsageError(f"mytools: not a file: {path}") from None
        except PermissionError:
            raise UsageError(f"mytools: cannot read: {path}") from None
        label = path

    if not hasattr(raw, "peek"):
        raw = io.BufferedReader(raw)
    if raw.peek(2)[:2] == GZIP_MAGIC:
        stream = io.TextIOWrapper(gzip.GzipFile(fileobj=raw), encoding="utf-8")
    else:
        stream = io.TextIOWrapper(raw, encoding="utf-8")
    # GzipFile does not own a fileobj handed to it, so closing `stream` would leave
    # `raw` open. Hand both back and let the caller close both.
    return stream, raw, label


def read_bed(path, on_header=None):
    """Yield BedRecords from `path` ('-' for stdin).

    Header lines are skipped; if `on_header` is given it is called with each one,
    which is how `-header` echoes them (SPEC.md §2).
    """
    stream, raw, label = _open_text(path)
    try:
        for lineno, raw_line in enumerate(stream, 1):
            line = raw_line.rstrip("\r\n")
            if not line.strip():
                continue
            if is_header(line):
                if on_header is not None:
                    on_header(line)
                continue
            yield parse_line(line, label, lineno)
    finally:
        if path != "-":
            stream.close()
            raw.close()
        else:
            # Don't close stdin; just drop the text wrapper's hold on it.
            stream.detach()


def sort_key(record):
    """Chrom lexicographically, then start, then end (SPEC.md §5).

    Lexicographic means `chr17` sorts before `chr2`. That is bedtools' behaviour,
    measured, and natural-order sorting is the common wrong answer.
    """
    return (record.chrom, record.start, record.end)


def require_sorted(records):
    """Pass records through, raising DataError on the first out-of-order one.

    `merge` and `closest` need this. Neither sorts for the user -- bedtools doesn't
    (SPEC.md §7).
    """
    previous = None
    for record in records:
        if previous is not None and sort_key(record) < sort_key(previous):
            raise DataError(
                f"{record.source}:{record.lineno}: input is not sorted "
                f"({record.chrom}:{record.start} after {previous.chrom}:{previous.start})"
            )
        previous = record
        yield record
