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

Re-runs only scan what changed. Each audited epoch is cached with its directory's modification time; raw-data epochs are never added to, and adding/removing/renaming a file always changes the directory mtime, so an unchanged epoch is skipped. Epochs that were mid-copy last time are therefore re-scanned automatically. Epochs whose header fix has started (`update_pointing_logs/` exists) but whose files are not all read-only yet are re-scanned each run and flagged `HDRFIX-INCOMPLETE`. The cache is saved every 50 epochs, so an interrupted run resumes where it stopped. Epochs that later disappear from disk (e.g. deleted after processing) stay in the cache as history.

`--full` re-scans everything; `--no-cache` neither reads nor writes the cache.

Options: `--no-size` skips the per-file `stat` (faster, but no TB or permissions); `--expect-beams N` (default 19).
Progress goes to stderr, the report to stdout (and `-o`). To keep a progress log too:
`python3 tools/crafts_audit.py /aphid/temp/fast_crafts/data -o audit.txt --tsv audit_epochs.tsv 2> audit_progress.log`
