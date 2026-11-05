#!/usr/bin/env python3
"""
crafts_audit.py — read-only census of the FAST CRAFTS raw-data tree.

Expected layout (from the processing scripts):

  <root>/<chunk>/<pointing>/<date>/<pointing>_<mode>-M<beam>_<seq>_<tag>.fits
  e.g. data/ZD2021_1_1_2bit/Dec+3240_01_05/20220329/Dec+3240_01_05_arcdrift-M01_0001_2bit.fits

It only lists directories and stats files (no file is opened or modified).
Anything that does not fit the layout is reported, not skipped silently.

Usage:
  crafts_audit.py [ROOT] [-o REPORT.txt] [--tsv EPOCHS.tsv] [--no-size] [--expect-beams 19]

  ROOT defaults to /aphid/temp/fast_crafts/data
"""

import argparse
import os
import re
import stat
import sys
import time
from collections import Counter, defaultdict

FITS_RE = re.compile(
    r'^(?P<ptg>.+?)_(?P<mode>[A-Za-z]+)-M(?P<beam>\d{2})_(?P<seq>\d+)(?:_(?P<tag>[^.]+))?\.fits$')

TB = 1e12


def human(nbytes):
    for unit, div in (("TB", 1e12), ("GB", 1e9), ("MB", 1e6), ("kB", 1e3)):
        if nbytes >= div:
            return f"{nbytes / div:.2f} {unit}"
    return f"{nbytes} B"


def subdirs_and_files(path, problems):
    """One scandir of path -> (sorted subdir names, list of DirEntry files)."""
    dirs, files = [], []
    try:
        with os.scandir(path) as it:
            for e in it:
                try:
                    if e.is_dir(follow_symlinks=False):
                        dirs.append(e.name)
                    else:
                        files.append(e)
                except OSError as err:
                    problems.append(f"UNREADABLE  {e.path}: {err}")
    except OSError as err:
        problems.append(f"UNREADABLE  {path}: {err}")
    return sorted(dirs), files


def compress_ranges(nums):
    """[1,2,3,7,9,10] -> '1-3,7,9-10'"""
    nums = sorted(nums)
    out, i = [], 0
    while i < len(nums):
        j = i
        while j + 1 < len(nums) and nums[j + 1] == nums[j] + 1:
            j += 1
        out.append(str(nums[i]) if i == j else f"{nums[i]}-{nums[j]}")
        i = j + 1
    return ",".join(out)


def audit_epoch(path, pointing, want_size, expect_beams, problems):
    """Census of one <pointing>/<date> directory full of .fits files."""
    dirs, files = subdirs_and_files(path, problems)
    ep = dict(nfits=0, nbytes=0, other_files=0, subdirs=dirs, beams=Counter(),
              seqs=defaultdict(list), zero=0, writable=0, badname=0,
              ptg_mismatch=0, modes=Counter(), tags=Counter(), newest=0.0)
    for e in files:
        if not e.name.endswith(".fits"):
            ep["other_files"] += 1
            continue
        ep["nfits"] += 1
        if want_size:
            try:
                st = e.stat(follow_symlinks=False)
                ep["nbytes"] += st.st_size
                ep["zero"] += st.st_size == 0
                ep["writable"] += bool(st.st_mode & stat.S_IWUSR)
                ep["newest"] = max(ep["newest"], st.st_mtime)
            except OSError as err:
                problems.append(f"UNREADABLE  {e.path}: {err}")
        m = FITS_RE.match(e.name)
        if not m:
            ep["badname"] += 1
            if ep["badname"] <= 3:
                problems.append(f"ODD NAME    {e.path}")
            continue
        beam = int(m["beam"])
        ep["beams"][beam] += 1
        ep["seqs"][beam].append(int(m["seq"]))
        ep["modes"][m["mode"]] += 1
        ep["tags"][m["tag"] or "-"] += 1
        if m["ptg"] != pointing:
            ep["ptg_mismatch"] += 1

    # --- per-epoch consistency flags ---
    flags = []
    nb = len(ep["beams"])
    if ep["nfits"] == 0:
        flags.append("EMPTY")
    elif nb != expect_beams:
        missing = sorted(set(range(1, expect_beams + 1)) - set(ep["beams"]))
        flags.append(f"BEAMS={nb}" + (f"(missing M{','.join(f'{b:02d}' for b in missing)})" if missing else ""))
    counts = set(ep["beams"].values())
    if len(counts) > 1:
        flags.append(f"UNEVEN {min(counts)}-{max(counts)}/beam")
    gaps = {}
    for b, s in ep["seqs"].items():
        full = set(range(min(s), max(s) + 1))
        missing = full - set(s)
        if missing:
            gaps[b] = missing
        if len(s) != len(set(s)):
            flags.append(f"DUPSEQ M{b:02d}")
    if gaps:
        flags.append("SEQGAPS " + " ".join(f"M{b:02d}:{compress_ranges(g)}" for b, g in sorted(gaps.items())[:4])
                     + (" ..." if len(gaps) > 4 else ""))
    if ep["zero"]:
        flags.append(f"ZEROSIZE={ep['zero']}")
    if ep["badname"]:
        flags.append(f"ODDNAMES={ep['badname']}")
    if ep["ptg_mismatch"]:
        flags.append(f"PTGNAME≠DIR={ep['ptg_mismatch']}")
    ep["hdrlog"] = "update_pointing_logs" in ep["subdirs"]
    unknown = [d for d in ep["subdirs"] if d != "update_pointing_logs"]
    if unknown:
        flags.append("SUBDIRS=" + ",".join(unknown[:4]))
    ep["flags"] = flags
    return ep


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("root", nargs="?", default="/aphid/temp/fast_crafts/data")
    ap.add_argument("-o", "--output", help="write the report here as well as to stdout")
    ap.add_argument("--tsv", help="also write one machine-readable line per epoch")
    ap.add_argument("--no-size", action="store_true", help="skip stat() of every file (faster; no TB / perms)")
    ap.add_argument("--expect-beams", type=int, default=19)
    args = ap.parse_args()

    root = os.path.abspath(args.root)
    want_size = not args.no_size
    t0 = time.time()
    problems = []
    epochs = []          # (chunk, pointing, date, ep)
    misplaced = Counter()  # fits files at the wrong depth, by level name

    chunks, f0 = subdirs_and_files(root, problems)
    misplaced["root"] += sum(e.name.endswith(".fits") for e in f0)
    for ci, chunk in enumerate(chunks, 1):
        print(f"[{ci}/{len(chunks)}] {chunk}", file=sys.stderr, flush=True)
        cpath = os.path.join(root, chunk)
        pointings, f1 = subdirs_and_files(cpath, problems)
        misplaced["chunk"] += sum(e.name.endswith(".fits") for e in f1)
        if not pointings:
            problems.append(f"EMPTY CHUNK {cpath}")
        for ptg in pointings:
            ppath = os.path.join(cpath, ptg)
            dates, f2 = subdirs_and_files(ppath, problems)
            nf2 = sum(e.name.endswith(".fits") for e in f2)
            if nf2:
                misplaced["pointing"] += nf2
                problems.append(f"FITS AT POINTING LEVEL ({nf2})  {ppath}")
            if not dates:
                problems.append(f"NO DATE DIRS {ppath}")
            for date in dates:
                ep = audit_epoch(os.path.join(ppath, date), ptg, want_size, args.expect_beams, problems)
                epochs.append((chunk, ptg, date, ep))

    # ---------------- report ----------------
    L = []
    w = L.append
    tot_files = sum(e["nfits"] for *_, e in epochs)
    tot_bytes = sum(e["nbytes"] for *_, e in epochs)
    n_ptg = len({(c, p) for c, p, _, _ in epochs})
    w(f"CRAFTS data audit   {time.strftime('%Y-%m-%d %H:%M:%S')}   root={root}")
    w(f"  chunks     {len(chunks):8d}")
    w(f"  pointings  {n_ptg:8d}")
    w(f"  epochs     {len(epochs):8d}   (pointing/date dirs)")
    w(f"  .fits      {tot_files:8d}" + (f"   {human(tot_bytes)}" if want_size else ""))
    flagged = sum(bool(e["flags"]) for *_, e in epochs)
    w(f"  epochs flagged {flagged}   misplaced .fits {sum(misplaced.values())}   problems {len(problems)}")
    w(f"  scan time {time.time() - t0:.1f} s")
    w("")

    # per chunk
    w("PER CHUNK")
    w(f"  {'chunk':<28} {'ptgs':>5} {'epochs':>6} {'fits':>9} {'size':>11}")
    bychunk = defaultdict(list)
    for c, p, d, e in epochs:
        bychunk[c].append((p, e))
    for c in chunks:
        rows = bychunk.get(c, [])
        w(f"  {c:<28} {len({p for p, _ in rows}):>5} {len(rows):>6} {sum(e['nfits'] for _, e in rows):>9}"
          f" {human(sum(e['nbytes'] for _, e in rows)) if want_size else '-':>11}")
    w("")

    # per epoch
    w("PER EPOCH  (beams = distinct M## present; f/beam = files per beam; ro = read-only files; hlog = update_pointing_logs/ exists -> header fix ran)")
    w(f"  {'chunk/pointing/date':<52} {'fits':>6} {'beams':>5} {'f/beam':>9} {'size':>10} {'ro':>6} {'hlog':>4}  flags")
    for c, p, d, e in epochs:
        cnts = e["beams"].values()
        fpb = f"{min(cnts)}-{max(cnts)}" if cnts and min(cnts) != max(cnts) else (str(min(cnts)) if cnts else "-")
        ro = f"{e['nfits'] - e['writable']}" if want_size else "-"
        w(f"  {c + '/' + p + '/' + d:<52} {e['nfits']:>6} {len(e['beams']):>5} {fpb:>9}"
          f" {human(e['nbytes']) if want_size else '-':>10} {ro:>6} {'yes' if e['hdrlog'] else '-':>4}  {' '.join(e['flags'])}")
    w("")

    # what kinds of files we saw
    modes, tags = Counter(), Counter()
    for *_, e in epochs:
        modes.update(e["modes"]); tags.update(e["tags"])
    w("FILENAME VOCABULARY")
    w("  modes: " + ", ".join(f"{k}={v}" for k, v in modes.most_common()))
    w("  tags:  " + ", ".join(f"{k}={v}" for k, v in tags.most_common()))
    other = sum(e["other_files"] for *_, e in epochs)
    if other:
        w(f"  non-.fits files inside epoch dirs: {other}")
    w("")

    if problems:
        w(f"PROBLEMS ({len(problems)})")
        for p in problems[:200]:
            w("  " + p)
        if len(problems) > 200:
            w(f"  ... {len(problems) - 200} more")

    text = "\n".join(L) + "\n"
    sys.stdout.write(text)
    if args.output:
        with open(args.output, "w") as f:
            f.write(text)
    if args.tsv:
        with open(args.tsv, "w") as f:
            f.write("chunk\tpointing\tdate\tnfits\tnbeams\tmin_per_beam\tmax_per_beam\tbytes\treadonly\thdrlog\tnewest_mtime\tflags\n")
            for c, p, d, e in epochs:
                cnts = list(e["beams"].values()) or [0]
                f.write(f"{c}\t{p}\t{d}\t{e['nfits']}\t{len(e['beams'])}\t{min(cnts)}\t{max(cnts)}\t"
                        f"{e['nbytes']}\t{e['nfits'] - e['writable'] if want_size else ''}\t{int(e['hdrlog'])}\t"
                        f"{time.strftime('%Y-%m-%dT%H:%M', time.localtime(e['newest'])) if e['newest'] else ''}\t"
                        f"{';'.join(e['flags'])}\n")


if __name__ == "__main__":
    main()
