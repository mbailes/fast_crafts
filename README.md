# fast_crafts

Bookkeeping tools for reprocessing the FAST **CRAFTS** drift-scan survey on the OzSTAR supercomputer.

The raw data (~2 PB and growing) are copied from disk drives into scratch at
`/aphid/temp/fast_crafts/data/`, laid out as

```
data/<chunk>/<pointing>/<date>/<pointing>_arcdrift-M<beam>_<seq>_2bit.fits
e.g. data/ZD2021_1_1_2bit/Dec+3240_01_05/20220329/Dec+3240_01_05_arcdrift-M01_0001_2bit.fits
```

One **epoch** (= pointing/date) is the unit of processing: 19 beams × N files.

The aim is status you can read with `cat`: plain text files, no web servers.

## Workflow

One-time setup on OzSTAR (e.g. in `~/.bashrc`):

```
export CRAFTS_STATE_DIR=/fred/oz002/fast_crafts/state     # persistent, never purged (TBC)
git clone https://github.com/mbailes/fast_crafts.git ~/fast_crafts
```

Each batch of drives:

```
cd ~/fast_crafts && git pull

# 1. copy the drives into /aphid/temp/fast_crafts/data/ as usual

# 2. audit once the last copy has finished (only new or changed epochs are scanned)
python3 tools/crafts_audit.py /aphid/temp/fast_crafts/data 2> $CRAFTS_STATE_DIR/audit_progress.log
cat $CRAFTS_STATE_DIR/audit.txt                 # see the HEADER FIX section

# 3. fix headers: one Slurm array task per epoch directory that still needs it
slurm/submit_hdrfix.sh

# 4. watch Slurm
python3 tools/crafts_queue.py                   # or: --watch 60;  or: cat $CRAFTS_STATE_DIR/queue.txt

# 5. when the queue has drained, audit again and see what is fixed
python3 tools/crafts_audit.py /aphid/temp/fast_crafts/data 2> $CRAFTS_STATE_DIR/audit_progress.log
cat $CRAFTS_STATE_DIR/epochs_hdrfixed.txt       # epochs ready for the next stage (RFI)
```

Files in `$CRAFTS_STATE_DIR`:

| file | what |
|---|---|
| `audit.txt` | the human report: totals, per chunk, per epoch, HEADER FIX section, problems |
| `audit_epochs.tsv` | one line per epoch, incl. header state and number of stamped files |
| `epochs_hdrfixed.txt` | epochs whose every file carries the PNTFIX stamp — ready for the next stage |
| `epochs_hdrfix_todo.txt` | epochs with files never tried or that failed transiently — `submit_hdrfix.sh` submits these |
| `epochs_hdrfix_stuck.txt` | epochs whose remaining files failed for a reason a re-run won't fix (SKIP_DT, NOROOM, ...) — need a person |
| `epochs_hdrfix_nopos.txt` | epochs with no `_pos.fits` pointing file |
| `hdrfix/<chunk>__<pointing>__<date>.csv` | per-file result of every fix attempt |
| `queue.txt` | last Slurm snapshot from `crafts_queue.py` |
| `submissions.log`, `submissions/` | every array submitted, and the frozen epoch list it ran over |
| `slurm_logs/hdrfix/<job>_<task>.out` | Slurm output of each fix task |
| `audit_cache.jsonl` | audit cache (delete or use `--full` to start over) |

Header-fix states per epoch (truth = the PNTFIX stamp in each file; reasons from the per-file logs):
**DONE** every file stamped · **TODO** some files never tried / transient error · **STUCK** only permanent failures left · **NOPOS** no pointing file.

### `tools/crafts_queue.py` — Slurm monitor

Your jobs now (running / pending, with pending reasons) and how jobs ended over the last `--since` hours
(COMPLETED / FAILED / TIMEOUT / OUT_OF_MEMORY / CANCELLED), by job name, plus the progress of each array
in `submissions.log`. Prints and writes `$CRAFTS_STATE_DIR/queue.txt`; `--watch 60` refreshes.

## Tools

### `tools/crafts_audit.py`

Read-only census of the raw-data tree. It lists directories and stats files only; it never opens or modifies a file. Needs only Python 3.

```
python3 tools/crafts_audit.py /aphid/temp/fast_crafts/data -o audit.txt --tsv audit_epochs.tsv
cat audit.txt
```

Runs in two phases, with live feedback on stderr:

1. **Directory listing** (fast): lists chunks, pointings and dates only and prints a directory summary straight away (pointings, epochs and date span per chunk).
2. **Epoch census** (slow): one progress line per epoch, with running totals, files/s, an ETA and any flags (`!!`).

The report gives start/finish times, the time for each phase and the total.

Report contents:

- totals: chunks, pointings, epochs, `.fits` files, TB
- per chunk: pointings, epochs, files, size
- per epoch: files, beams present, files per beam, size, read-only count, whether `update_pointing_logs/` exists (evidence the header fix has run)
- flags: missing beams, uneven files per beam, sequence gaps/duplicates, zero-size files, unrecognised names, `.fits` at the wrong depth, empty chunks/pointings

### Persistent state and incremental re-runs

Put the report, the per-epoch TSV and the cache in a directory that is **never purged** (not `/aphid/temp`):

```
export CRAFTS_STATE_DIR=/path/to/persistent/fast_crafts_state     # e.g. in ~/.bashrc
python3 tools/crafts_audit.py /aphid/temp/fast_crafts/data 2> $CRAFTS_STATE_DIR/audit_progress.log
cat $CRAFTS_STATE_DIR/audit.txt
```

With `--state-dir` (or `$CRAFTS_STATE_DIR`) set, `audit.txt`, `audit_epochs.tsv` and `audit_cache.jsonl` go there.

Re-runs only scan what changed. The audit also reads each file's primary header (a few kB) to see whether it carries the `PNTFIX` stamp; `--no-headers` skips that. Each audited epoch is cached with its directory's modification time; raw-data epochs are never added to, and adding/removing/renaming a file always changes the directory mtime, so an unchanged epoch is skipped. Epochs that were mid-copy last time are therefore re-scanned automatically. Epochs whose header fix has started (`update_pointing_logs/` exists) but whose files are not all read-only yet are re-scanned each run and flagged `HDRFIX-INCOMPLETE`. The cache is saved every 50 epochs, so an interrupted run resumes where it stopped. Epochs that later disappear from disk (e.g. deleted after processing) stay in the cache as history.

`--full` re-scans everything; `--no-cache` neither reads nor writes the cache.

Options: `--no-size` skips the per-file `stat` (faster, but no TB or permissions); `--expect-beams N` (default 19).
Progress goes to stderr, the report to stdout (and `-o`). To keep a progress log too:
`python3 tools/crafts_audit.py /aphid/temp/fast_crafts/data -o audit.txt --tsv audit_epochs.tsv 2> audit_progress.log`

### `tools/fix_pointing_v3.py` — header (pointing) fix, with a stamp in every file

Successor to the student's `update_pointing_v2.py`. It writes the same values: RA/DEC/RA_DEG/DEC_DEG
from `<posdir>/<pointing>_<date>_pos.fits`, taking the sample nearest the file's start MJD for the
file's beam, and rejecting it if more than 1 s away. On top of v2:

- **Stamp:** every fixed file carries `PNTFIX` (tool + UTC), `PNTSRC`, `PNTMJD`, `PNTDT` and a `HISTORY` card.
  A file with `PNTFIX` is fixed; a file without it is not.
- **No whole-file rewrites:** only the primary header's own 2880-byte blocks are rewritten, in place,
  and only if the new header fits. Otherwise the file is untouched and reported `NOROOM`.
  (astropy `mode='update'` would rewrite the whole multi-GB file when a header grows.)
- Re-reads and checks the header after writing (`VERIFYFAIL` otherwise).
- **Idempotent:** stamped files are skipped (`ALREADY`); `--force` redoes them. Files fixed by v2 are
  re-fixed and stamped, and the log records how far v2's position was from v3's (`shift_arcsec`).
- **Per-file result log:** `$CRAFTS_STATE_DIR/hdrfix/<chunk>__<pointing>__<date>.csv`
  with statuses `FIXED ALREADY SKIP_DT NOROOM NOPOS NOCOL NOBEAM NOSTART VERIFYFAIL ERROR`.
- Warns if a beam's Dec is > 0.5 deg from the Dec in the pointing name.
- Files are left read-only (0444), as before.

```
python3 tools/fix_pointing_v3.py --dry-run /aphid/temp/fast_crafts/data/ZD2021_1_1_2bit/Dec+3240_01_05/20220329
python3 tools/fix_pointing_v3.py           /aphid/temp/fast_crafts/data/ZD2021_1_1_2bit/Dec+3240_01_05/20220329
```

`--posdir` (or `$CRAFTS_POSDIR`) sets where the `_pos.fits` files are; the default is the student's
`/fred/oz002/users/vmandaog/pulsar_search/FAST/crafts_2bit_miaocc01`. Needs numpy + astropy.

Batch: `slurm/fix_pointing_v3.slurm` runs one array task per epoch from a list (instructions in its header).
