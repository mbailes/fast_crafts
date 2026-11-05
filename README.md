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

Options: `--no-size` skips the per-file `stat` (faster, but no TB or permissions); `--expect-beams N` (default 19).
Progress goes to stderr, the report to stdout (and `-o`). To keep a progress log too:
`python3 tools/crafts_audit.py /aphid/temp/fast_crafts/data -o audit.txt --tsv audit_epochs.tsv 2> audit_progress.log`
