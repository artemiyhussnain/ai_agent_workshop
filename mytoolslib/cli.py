"""Subcommand dispatch, flag parsing and exit codes.

Exit codes, per SPEC.md §7: 0 success (including empty output), 1 data problems,
2 caller problems. Errors go to stderr; stdout carries data only, because stdout
gets piped.

This deliberately differs from bedtools, which exits 1 for everything, exits 0 on
no arguments, and aborts with SIGABRT (exit 134) on a non-integer coordinate.
"""

from __future__ import annotations

import importlib
import os
import sys

from .bed import DataError, MyToolsError, UsageError

VERSION = "0.1.0"

# Pre-registered so that adding a subcommand means adding one new file, not editing
# this table. Parallel branches therefore do not collide here.
SUBCOMMANDS = {
    "sort": "mytoolslib.cmd_sort",
    "merge": "mytoolslib.cmd_merge",
    "intersect": "mytoolslib.cmd_intersect",
    "subtract": "mytoolslib.cmd_subtract",
    "closest": "mytoolslib.cmd_closest",
}

USAGE = f"""\
usage: mytools <command> [options]

mytools {VERSION} -- a small reimplementation of a subset of bedtools.
Real bedtools is the oracle; see SPEC.md.

commands:
  sort       sort by chrom (lexicographic), then start, then end
  merge      collapse overlapping and adjacent features
  intersect  report features shared between -a and -b
  subtract   remove -b regions from -a features
  closest    for each -a feature, the nearest -b feature

options:
  --version   print the version and exit
  -h, --help  print this message and exit

BED is 0-based half-open: chr1 100 200 covers bases 100..199.
"""


def parse_flags(args, *, flags=(), options=(), required=()):
    """Split bedtools-style arguments into a dict of options and the rest.

    `flags` are boolean (`-u`); `options` take a following value (`-d 10`).
    Unknown options and missing values raise UsageError (exit 2).
    """
    parsed = {name: False for name in flags}
    parsed.update({name: None for name in options})
    rest = []
    index = 0
    while index < len(args):
        arg = args[index]
        if arg in flags:
            parsed[arg] = True
            index += 1
        elif arg in options:
            if index + 1 >= len(args):
                raise UsageError(f"mytools: {arg} requires a value")
            parsed[arg] = args[index + 1]
            index += 2
        elif arg.startswith("-") and arg != "-":
            raise UsageError(f"mytools: unknown option: {arg}")
        else:
            rest.append(arg)
            index += 1
    for name in required:
        if parsed.get(name) is None:
            raise UsageError(f"mytools: {name} is required")
    return parsed, rest


def exclusive(parsed, *names):
    """Raise UsageError if more than one of `names` is set."""
    on = [name for name in names if parsed.get(name)]
    if len(on) > 1:
        raise UsageError(f"mytools: {' and '.join(on)} are mutually exclusive")


def int_option(parsed, name, default):
    value = parsed.get(name)
    if value is None:
        return default
    try:
        return int(value)
    except ValueError:
        raise UsageError(f"mytools: {name} expects an integer, got {value!r}") from None


def float_option(parsed, name, default):
    value = parsed.get(name)
    if value is None:
        return default
    try:
        return float(value)
    except ValueError:
        raise UsageError(f"mytools: {name} expects a number, got {value!r}") from None


def header_echo(enabled, out=None):
    """Return an `on_header` callable for read_bed, honouring -header."""
    stream = out if out is not None else sys.stdout
    if not enabled:
        return None

    def echo(line):
        stream.write(line + "\n")

    return echo


def _dispatch(name, args):
    module_name = SUBCOMMANDS[name]
    try:
        module = importlib.import_module(module_name)
    except ModuleNotFoundError as exc:
        # Only swallow "this subcommand does not exist yet"; a missing dependency
        # inside an existing module must still surface.
        if exc.name == module_name:
            raise UsageError(
                f"mytools: {name} is not implemented yet -- see SPEC.md and the open issues"
            ) from None
        raise
    return module.run(args)


def main(argv):
    args = argv[1:]
    if not args:
        sys.stderr.write(USAGE)
        return 2
    first = args[0]
    if first == "--version":
        print(f"mytools {VERSION}")
        return 0
    if first in ("-h", "--help"):
        sys.stdout.write(USAGE)
        return 0

    try:
        if first not in SUBCOMMANDS:
            raise UsageError(f"mytools: unknown command or option: {first}")
        return _dispatch(first, args[1:]) or 0
    except MyToolsError as exc:
        sys.stderr.write(f"{exc}\n")
        if isinstance(exc, UsageError) and first not in SUBCOMMANDS:
            sys.stderr.write(USAGE)
        return exc.exit_code
    except BrokenPipeError:
        # `mytools sort -i big.bed | head` is normal usage, not an error.
        devnull = os.open(os.devnull, os.O_WRONLY)
        os.dup2(devnull, sys.stdout.fileno())
        return 0
    except KeyboardInterrupt:
        sys.stderr.write("mytools: interrupted\n")
        return 130
    except Exception as exc:  # noqa: BLE001 -- see comment
        # SPEC.md §7: never traceback. bedtools aborts with SIGABRT on malformed
        # input and we will not. Set MYTOOLS_DEBUG=1 to get the traceback back
        # while developing.
        if os.environ.get("MYTOOLS_DEBUG"):
            raise
        sys.stderr.write(f"mytools: internal error: {type(exc).__name__}: {exc}\n")
        return 1
