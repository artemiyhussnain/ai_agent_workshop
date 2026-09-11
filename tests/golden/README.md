# Golden test cases, one file per subcommand

`tests/run_golden.sh` sources every `*.sh` in this directory. Add a subcommand by
adding a file here — `sort.sh`, `merge.sh`, `intersect.sh` — and never by editing
the harness. Several branches add cases at once, and a file each is what keeps them
from conflicting.

Each file is a fragment sourced with the harness already set up. Available:

| | |
|---|---|
| `check "<name>" -- <args>` | run `mytools` and `bedtools` with the same args, compare exit code then stdout |
| `skip "<name>" "<why>"` | record a case that cannot run here, visibly |
| `STDIN=<file> check …` | feed both commands the same stdin for one case |
| `$DATA` | the `data/` directory |
| `$SORTED_A`, `$SORTED_B`, `$SORTED_GENES` | bedtools-sorted copies (`a.bed`/`b.bed` are deliberately unsorted) |
| `$GZ_A`, `$GZ_SORTED_A` | gzipped copies, for the compression path |
| `$HEADED` | `#`/`track`/`browser` lines plus sorted data, for `-header` |
| `$EMPTY` | an empty file |
| `$tmp` | scratch directory, removed on exit |

Keep one case per flag combination, and name it after what it exercises rather
than after the flags alone.
