#!/usr/bin/env python3
"""
crafts_audit.py — read-only census of the FAST CRAFTS raw-data tree.

Expected layout (from the processing scripts):

  <root>/<chunk>/<pointing>/<date>/<pointing>_<mode>-M<beam>_<seq>_<tag>.fits
  e.g. data/ZD2021_1_1_2bit/Dec+3240_01_05/20220329/Dec+3240_01_05_arcdrift-M01_0001_2bit.fits

It only lists directories and stats files (no file is opened or modified).
Anything that does not fit the layout is reported, not skipped silently.

Usage:
  crafts_audit.py [ROOT] [--state-dir DIR] [-o REPORT.txt] [--tsv EPOCHS.tsv]
                  [--no-size] [--full] [--no-cache] [--expect-beams 19]

  ROOT defaults to /aphid/temp/fast_crafts/data

Re-runs are incremental: each finished epoch is cached with its directory mtime and
skipped next time unless files were added/removed (or its header fix is mid-way).
"""

import argparse
import json
import os
import re
import stat
import sys
import time
from collections import Counter, defaultdict

FITS_RE = re.compile(
    r'^(?P<ptg>.+?)_(?P<mode>[A-Za-z]+)-M(?P<beam>\d{2})_(?P<seq>\d+)(?:_(?P<tag>[^.]+))?\.fits$')

TB = 1e12
CACHE_VERSION = 1


# ---------------------------------------------------------------------------
# Cache of finished epochs.  Raw-data epoch directories are never added to once
# copied, so an epoch whose directory mtime is unchanged need not be re-scanned.
# (Adding/removing/renaming a file always changes the directory's mtime.)
# ---------------------------------------------------------------------------

def ep_to_json(key, ep):
    d = {k: v for k, v in ep.items() if k != "seqs"}
    d["beams"] = {str(b): n for b, n in ep["beams"].items()}
    d["modes"] = dict(ep["modes"])
    d["tags"] = dict(ep["tags"])
    d["key"] = key
    d["v"] = CACHE_VERSION
    return json.dumps(d, sort_keys=True)


def ep_from_json(d):
    d = dict(d)
    d["beams"] = Counter({int(b): n for b, n in d["beams"].items()})
    d["modes"] = Counter(d["modes"])
    d["tags"] = Counter(d["tags"])
    return d


def load_cache(path):
    cache = {}
    if not path or not os.path.exists(path):
        return cache
    with open(path) as f:
        for line in f:
            try:
                d = json.loads(line)
            except ValueError:
                continue              # a torn last line from an interrupted write
            if d.get("v") == CACHE_VERSION:
                cache[d["key"]] = d
    return cache


def save_cache(path, cache):
    """Atomic rewrite: write a temp file next to it, then rename over."""
    if not path:
        return
    tmp = f"{path}.tmp.{os.getpid()}"
    with open(tmp, "w") as f:
        for d in cache.values():
            f.write(json.dumps(d, sort_keys=True) + "\n")
    os.replace(tmp, path)


def cache_still_valid(d, dir_mtime_ns, want_size):
    if d.get("dir_mtime_ns") != dir_mtime_ns:
        return False                  # files added/removed since: re-audit
    if want_size and not d.get("sized"):
        return False                  # cached by a --no-size run
    if want_size and d.get("hdrlog") and d["nfits"] - d["writable"] < d["nfits"]:
        return False                  # header fix ran but chmod 444 may have happened since
    return True


def human(nbytes):
    for unit, div in (("TB", 1e12), ("GB", 1e9), ("MB", 1e6), ("kB", 1e3)):
        if nbytes >= div:
            return f"{nbytes / div:.2f} {unit}"
    return f"{nbytes} B"


def fmt_dur(sec):
    sec = int(round(sec))
    h, rem = divmod(sec, 3600)
    m, s = divmod(rem, 60)
    return f"{h}h{m:02d}m{s:02d}s" if h else (f"{m}m{s:02d}s" if m else f"{s}s")


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
    if ep["hdrlog"] and want_size and ep["writable"]:
        flags.append(f"HDRFIX-INCOMPLETE ro={ep['nfits'] - ep['writable']}/{ep['nfits']}")
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
    ap.add_argument("--state-dir", default=os.environ.get("CRAFTS_STATE_DIR"),
                    help="persistent (never-purged) directory for the report, tsv and cache; "
                         "default $CRAFTS_STATE_DIR. With it set, -o/--tsv default to audit.txt/audit_epochs.tsv there")
    ap.add_argument("--cache", default=None,
                    help="cache of finished epochs (default: audit_cache.jsonl in the state dir, "
                         "else next to ROOT)")
    ap.add_argument("--full", action="store_true", help="ignore the cache and re-scan every epoch")
    ap.add_argument("--no-cache", action="store_true", help="neither read nor write a cache")
    args = ap.parse_args()

    root = os.path.abspath(args.root)
    want_size = not args.no_size
    problems = []
    epochs = []          # (chunk, pointing, date, ep)
    misplaced = Counter()  # fits files at the wrong depth, by level name
    t_start = time.time()
    started = time.strftime('%Y-%m-%d %H:%M:%S')

    def say(msg):
        print(f"[{fmt_dur(time.time() - t_start):>8}] {msg}", file=sys.stderr, flush=True)

    say(f"audit of {root} started {started}")
    state = args.state_dir
    if state:
        os.makedirs(state, exist_ok=True)
        args.output = args.output or os.path.join(state, "audit.txt")
        args.tsv = args.tsv or os.path.join(state, "audit_epochs.tsv")
    cache_path = None if args.no_cache else (
        args.cache or os.path.join(state or os.path.dirname(root), "audit_cache.jsonl"))
    old = load_cache(cache_path) if cache_path else {}
    cache = {} if args.full else dict(old)
    if state:
        say(f"state dir {state}  (report, tsv and cache are written here)")
    if cache_path:
        say(f"cache {cache_path}: {len(cache)} epochs loaded" + ("  (--full: ignoring them)" if args.full else ""))

    # ---- phase 1: directory skeleton only (fast: three levels of listings, no file stats) ----
    say("phase 1: listing chunk / pointing / date directories ...")
    chunks, f0 = subdirs_and_files(root, problems)
    misplaced["root"] += sum(e.name.endswith(".fits") for e in f0)
    plan = []            # (chunk, pointing, date)
    for chunk in chunks:
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
            plan += [(chunk, ptg, d) for d in dates]
    t_phase1 = time.time() - t_start

    skel = defaultdict(lambda: defaultdict(list))   # chunk -> pointing -> [dates]
    for c, p, d in plan:
        skel[c][p].append(d)
    tree = ["DIRECTORY SUMMARY  (from listings only)",
            f"  {len(chunks)} chunks, {sum(len(v) for v in skel.values())} pointings, {len(plan)} epochs",
            f"  {'chunk':<28} {'ptgs':>5} {'epochs':>6}  dates"]
    for c in chunks:
        ds = sorted(d for p in skel[c] for d in skel[c][p])
        span = f"{ds[0]} .. {ds[-1]}" if ds else "-"
        tree.append(f"  {c:<28} {len(skel[c]):>5} {len(ds):>6}  {span}")
    sys.stderr.write("\n".join(tree) + "\n")
    say(f"phase 1 done in {fmt_dur(t_phase1)}")

    # ---- phase 2: census of every epoch, with running progress ----
    t2 = time.time()
    # Decide up front which epochs need scanning (one stat per epoch dir), so the ETA is honest.
    todo, n_cached = [], 0
    for chunk, ptg, date in plan:
        key = f"{chunk}/{ptg}/{date}"
        try:
            mt = os.stat(os.path.join(root, key)).st_mtime_ns
        except OSError as err:
            problems.append(f"UNREADABLE  {os.path.join(root, key)}: {err}")
            continue
        d = cache.get(key)
        if d is not None and cache_still_valid(d, mt, want_size):
            ep = ep_from_json(d)
            problems.extend(ep.get("problems", []))
            epochs.append((chunk, ptg, date, ep))
            n_cached += 1
        else:
            todo.append((chunk, ptg, date, key, mt, d is not None))
    say(f"phase 2: {n_cached} epochs unchanged since last audit (from cache), {len(todo)} to scan"
        + ("" if want_size else " (--no-size: names only)") + " ...")
    run_files = run_bytes = 0
    last_chunk = None
    since_save = 0
    for i, (chunk, ptg, date, key, mt, was_cached) in enumerate(todo, 1):
        if chunk != last_chunk:
            say(f"chunk {chunks.index(chunk) + 1}/{len(chunks)}  {chunk}")
            last_chunk = chunk
        eprob = []
        ep = audit_epoch(os.path.join(root, key), ptg, want_size, args.expect_beams, eprob)
        problems.extend(eprob)
        ep.update(problems=eprob, dir_mtime_ns=mt, sized=want_size,
                  audited_at=time.strftime('%Y-%m-%dT%H:%M:%S'))
        epochs.append((chunk, ptg, date, ep))
        if cache_path:
            cache[key] = json.loads(ep_to_json(key, ep))
            since_save += 1
            if since_save >= 50:
                save_cache(cache_path, cache)
                since_save = 0
        run_files += ep["nfits"]
        run_bytes += ep["nbytes"]
        el = time.time() - t2
        eta = el / i * (len(todo) - i)
        say(f"  {i:>5}/{len(todo)} {ptg}/{date}  {ep['nfits']:>5} files"
            + (f" {human(ep['nbytes']):>10}" if want_size else "")
            + ("  (changed)" if was_cached else "")
            + f" | total {run_files} files" + (f" {human(run_bytes)}" if want_size else "")
            + f", {run_files / el if el else 0:.0f} files/s, ETA {fmt_dur(eta)}"
            + (f"  !! {' '.join(ep['flags'])}" if ep["flags"] else ""))
    if cache_path:
        # keep entries for epochs no longer on disk (e.g. deleted after processing) - useful history
        on_disk = {f"{c}/{p}/{d}" for c, p, d in plan}
        for k, d in old.items():
            cache.setdefault(k, d)
        gone = [k for k in cache if k not in on_disk]
        save_cache(cache_path, cache)
        say(f"cache saved: {len(cache)} epochs ({len(gone)} no longer on disk)")
    n_scanned = len(todo)
    order = {k: n for n, k in enumerate(plan)}
    epochs.sort(key=lambda t: order[(t[0], t[1], t[2])])
    t_phase2 = time.time() - t2
    t_total = time.time() - t_start
    say(f"phase 2 done in {fmt_dur(t_phase2)}; audit finished in {fmt_dur(t_total)}")

    # ---------------- report ----------------
    L = []
    w = L.append
    tot_files = sum(e["nfits"] for *_, e in epochs)
    tot_bytes = sum(e["nbytes"] for *_, e in epochs)
    n_ptg = len({(c, p) for c, p, _, _ in epochs})
    w(f"CRAFTS data audit   root={root}")
    w(f"  started {started}   finished {time.strftime('%Y-%m-%d %H:%M:%S')}   host={os.uname().nodename}")
    w(f"  chunks     {len(chunks):8d}")
    w(f"  pointings  {n_ptg:8d}")
    w(f"  epochs     {len(epochs):8d}   (pointing/date dirs)")
    w(f"  .fits      {tot_files:8d}" + (f"   {human(tot_bytes)}" if want_size else ""))
    flagged = sum(bool(e["flags"]) for *_, e in epochs)
    w(f"  epochs flagged {flagged}   misplaced .fits {sum(misplaced.values())}   problems {len(problems)}")
    w("")
    w("TIMING")
    w(f"  phase 1 (directory listing)  {fmt_dur(t_phase1):>10}")
    scanned_files = sum(e["nfits"] for *_, e in epochs if "seqs" in e)
    w(f"  phase 2 (epoch census)       {fmt_dur(t_phase2):>10}"
      + f"   {n_scanned} epochs scanned, {n_cached} from cache"
      + (f", {scanned_files / t_phase2:.0f} files/s" if t_phase2 > 0 and n_scanned else "")
      + ("" if want_size else "   (--no-size)"))
    w(f"  total                        {fmt_dur(t_total):>10}")
    w("")
    L += tree
    w("")

    # per chunk
    w("PER CHUNK")
    w(f"  {'chunk':<28} {'ptgs':>5} {'epochs':>6} {'fits':>9} {'size':>11} {'flagged':>7}")
    bychunk = defaultdict(list)
    for c, p, d, e in epochs:
        bychunk[c].append((p, e))
    for c in chunks:
        rows = bychunk.get(c, [])
        w(f"  {c:<28} {len({p for p, _ in rows}):>5} {len(rows):>6} {sum(e['nfits'] for _, e in rows):>9}"
          f" {human(sum(e['nbytes'] for _, e in rows)) if want_size else '-':>11}"
          f" {sum(bool(e['flags']) for _, e in rows):>7}")
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
