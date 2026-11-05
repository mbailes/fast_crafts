#!/usr/bin/env python3
"""
fix_pointing_v3.py — write the true beam pointing into CRAFTS PSRFITS primary headers,
and leave a stamp so every file says whether (and how) it was fixed.

Successor to viren1408/FAST_CRAFTS_Pulsar_Search preprocessing/update_pointing_v2.py.
Same positions as v2: for each file, the sample of <pointing>_<date>_pos.fits nearest to the
file's start MJD (STT_IMJD + (STT_SMJD + STT_OFFS)/86400), for the file's beam M##,
rejected if more than --max-dt seconds away.

What v3 adds
  * Stamp keywords in the primary header of every fixed file:
        PNTFIX  = 'fix_pointing_v3 <UTC>'   tool + time of the fix
        PNTSRC  = '<pos file name>'          where the position came from
        PNTMJD  = <MJD of the pos sample used>
        PNTDT   = <pos MJD - file start, s>
    plus a HISTORY card.  A file with PNTFIX is fixed; anything else is not.
  * No whole-file rewrites.  Only the primary header's own 2880-byte blocks are rewritten,
    in place, and only if the new header fits in the same number of blocks.  Otherwise the
    file is left untouched and reported NOROOM (astropy's mode='update' would silently
    rewrite the whole multi-GB file, and a killed job could leave it corrupt).
  * The header is re-read after writing and checked (VERIFYFAIL otherwise).
  * Idempotent: files already stamped are skipped (ALREADY) unless --force.
    Files fixed by v2 (RA_DEG present, no stamp) are re-fixed and stamped; the log shows
    how far the v2 position was from the v3 one (should be 0).
  * Every file gets a result line in a per-epoch CSV log in the state dir:
        <state>/hdrfix/<chunk>__<pointing>__<date>.csv
  * Files are left read-only (0444) after a successful fix, as in v2.

Usage
  fix_pointing_v3.py EPOCH_DIR [EPOCH_DIR ...]
  fix_pointing_v3.py --list epochs.txt [--task N]      # N-th line (0-based), for Slurm arrays
  options: --posdir DIR  --state-dir DIR  --max-dt 1.0  --dry-run  --force

  EPOCH_DIR = .../data/<chunk>/<pointing>/<date>   (the directory holding the .fits files)
  pos file  = <posdir>/<pointing>_<date>_pos.fits

Exit status: 0 if every file ended FIXED/ALREADY, 1 if anything else happened, 2 on usage errors.
"""

import argparse
import csv
import math
import os
import re
import stat
import sys
import time

import numpy as np
from astropy.io import fits

TOOL = "fix_pointing_v3"
BLOCK = 2880
DEFAULT_POSDIR = "/fred/oz002/users/vmandaog/pulsar_search/FAST/crafts_2bit_miaocc01"
BEAM_RE = re.compile(r"-M(\d{2})_")
DECNAME_RE = re.compile(r"Dec([+-])(\d{2})(\d{2})")

OK_STATUSES = {"FIXED", "ALREADY", "DRYRUN"}


# ---------------------------------------------------------------- raw header I/O

def read_primary_header_bytes(path):
    """Return the primary header exactly as stored: whole 2880-byte blocks up to and incl. END."""
    out = bytearray()
    with open(path, "rb") as f:
        while True:
            blk = f.read(BLOCK)
            if len(blk) < BLOCK:
                raise ValueError("file ends before primary-header END card")
            out += blk
            for i in range(0, BLOCK, 80):
                if blk[i:i + 8] == b"END     ":
                    return bytes(out)
            if len(out) > 1000 * BLOCK:
                raise ValueError("no END card in first 1000 blocks")


def write_header_in_place(path, new_bytes, old_len):
    """Overwrite the first old_len bytes of path with new_bytes (same length), fsync'd."""
    assert len(new_bytes) == old_len
    fd = os.open(path, os.O_WRONLY)
    try:
        n = os.pwrite(fd, new_bytes, 0)
        if n != old_len:
            raise OSError(f"short write {n}/{old_len}")
        os.fsync(fd)
    finally:
        os.close(fd)


# ---------------------------------------------------------------- pointing data

def read_pos_file(path):
    with fits.open(path, memmap=False) as h:
        d = h[1].data
        cols = {name: np.array(d[name]) for name in d.columns.names}
    mjd = cols["mjd"].astype(np.float64)
    order = np.argsort(mjd)
    return {k: v[order] for k, v in cols.items()}


def start_mjd(hdr):
    return hdr["STT_IMJD"] + (hdr["STT_SMJD"] + hdr.get("STT_OFFS", 0.0)) / 86400.0


def dec_from_name(pointing):
    m = DECNAME_RE.search(pointing)
    if not m:
        return None
    sign = -1 if m.group(1) == "-" else 1
    return sign * (int(m.group(2)) + int(m.group(3)) / 60.0)


def as_str(x):
    return x.decode() if isinstance(x, bytes) else str(x)


# ---------------------------------------------------------------- one file

def fix_file(path, pos, posname, max_dt, dec_nominal, dry_run, force):
    """Return a dict for the CSV log."""
    rec = dict(file=os.path.basename(path), status="", beam="", start_mjd="", pos_mjd="", dt_s="",
               ra="", dec="", ra_deg="", dec_deg="", prev_ra_deg="", prev_dec_deg="",
               shift_arcsec="", warn="", detail="")
    m = BEAM_RE.search(rec["file"])
    if not m:
        rec.update(status="NOBEAM", detail="no -M##_ in file name")
        return rec
    beam = f"M{m.group(1)}"
    rec["beam"] = beam

    raw = read_primary_header_bytes(path)
    hdr = fits.Header.fromstring(raw.decode("ascii"))

    if "PNTFIX" in hdr and not force:
        rec.update(status="ALREADY", detail=hdr["PNTFIX"],
                   ra_deg=hdr.get("RA_DEG", ""), dec_deg=hdr.get("DEC_DEG", ""))
        return rec

    try:
        t0 = start_mjd(hdr)
    except KeyError as e:
        rec.update(status="NOSTART", detail=f"missing {e}")
        return rec
    rec["start_mjd"] = f"{t0:.12f}"

    try:
        ra_deg_col, dec_deg_col = pos[f"ra_deg_{beam}"], pos[f"dec_deg_{beam}"]
        ra_col, dec_col = pos[f"ra_{beam}"], pos[f"dec_{beam}"]
    except KeyError as e:
        rec.update(status="NOCOL", detail=f"pos file has no column {e}")
        return rec

    mjd = pos["mjd"]
    i = int(np.argmin(np.abs(mjd - t0)))
    dt = (mjd[i] - t0) * 86400.0
    rec.update(pos_mjd=f"{mjd[i]:.12f}", dt_s=f"{dt:.5f}")
    if abs(dt) > max_dt:
        rec.update(status="SKIP_DT", detail=f"nearest pos sample {dt:.3f} s from file start (> {max_dt} s)")
        return rec

    ra_deg, dec_deg = float(ra_deg_col[i]), float(dec_deg_col[i])
    ra_s, dec_s = as_str(ra_col[i]).strip(), as_str(dec_col[i]).strip()
    rec.update(ra=ra_s, dec=dec_s, ra_deg=f"{ra_deg:.6f}", dec_deg=f"{dec_deg:.6f}")

    if "RA_DEG" in hdr and "DEC_DEG" in hdr:     # fixed by v2: how different is v3?
        pr, pd = float(hdr["RA_DEG"]), float(hdr["DEC_DEG"])
        rec.update(prev_ra_deg=f"{pr:.6f}", prev_dec_deg=f"{pd:.6f}")
        dra = (ra_deg - pr) * math.cos(math.radians(dec_deg))
        rec["shift_arcsec"] = f"{3600 * math.hypot(dra, dec_deg - pd):.3f}"

    warns = []
    if dec_nominal is not None and abs(dec_deg - dec_nominal) > 0.5:
        warns.append(f"DEC {dec_deg:.3f} vs name {dec_nominal:.3f}")
    rec["warn"] = "; ".join(warns)

    utc = time.strftime("%Y-%m-%dT%H:%M:%S", time.gmtime())
    hdr["RA"] = ra_s
    hdr["DEC"] = dec_s
    hdr["RA_DEG"] = (ra_deg, "[deg] beam RA at file start (from PNTSRC)")
    hdr["DEC_DEG"] = (dec_deg, "[deg] beam Dec at file start (from PNTSRC)")
    hdr["PNTFIX"] = (f"{TOOL} {utc}", "pointing fixed by")
    hdr["PNTSRC"] = (posname, "pointing source file")
    hdr["PNTMJD"] = (float(mjd[i]), "MJD of pointing sample used")
    hdr["PNTDT"] = (round(dt, 6), "[s] pointing sample minus file start")
    hdr.add_history(f"{TOOL} {utc}: RA/DEC {beam} from {posname} dt={dt:.4f}s")

    new = hdr.tostring(sep="", endcard=True, padding=True).encode("ascii")
    if len(new) != len(raw):
        rec.update(status="NOROOM",
                   detail=f"new header {len(new) // BLOCK} blocks vs {len(raw) // BLOCK}; file NOT modified")
        return rec

    if dry_run:
        rec["status"] = "DRYRUN"
        return rec

    mode = stat.S_IMODE(os.stat(path).st_mode)
    if not mode & stat.S_IWUSR:
        os.chmod(path, mode | stat.S_IWUSR)
    try:
        write_header_in_place(path, new, len(raw))
    finally:
        os.chmod(path, 0o444)

    chk = fits.Header.fromstring(read_primary_header_bytes(path).decode("ascii"))
    if chk.get("PNTFIX") != hdr["PNTFIX"] or abs(float(chk["RA_DEG"]) - ra_deg) > 1e-9:
        rec.update(status="VERIFYFAIL", detail="header re-read does not match what was written")
        return rec
    rec["status"] = "FIXED"
    return rec


# ---------------------------------------------------------------- one epoch

def fix_epoch(epoch_dir, args):
    epoch_dir = os.path.abspath(epoch_dir.rstrip("/"))
    date = os.path.basename(epoch_dir)
    pointing = os.path.basename(os.path.dirname(epoch_dir))
    chunk = os.path.basename(os.path.dirname(os.path.dirname(epoch_dir)))
    tag = f"{chunk}/{pointing}/{date}"
    posname = f"{pointing}_{date}_pos.fits"
    pospath = os.path.join(args.posdir, posname)

    files = sorted(f for f in os.listdir(epoch_dir) if f.endswith(".fits"))
    log_rows = []
    t0 = time.time()
    print(f"[{tag}] {len(files)} files, pos file {pospath}", flush=True)

    if not os.path.exists(pospath):
        log_rows = [dict(file=f, status="NOPOS", detail=pospath) for f in files]
    else:
        pos = read_pos_file(pospath)
        dec_nominal = dec_from_name(pointing)
        for k, f in enumerate(files, 1):
            try:
                rec = fix_file(os.path.join(epoch_dir, f), pos, posname, args.max_dt,
                               dec_nominal, args.dry_run, args.force)
            except Exception as e:                       # one bad file must not stop the epoch
                rec = dict(file=f, status="ERROR", detail=f"{type(e).__name__}: {e}")
            log_rows.append(rec)
            if k % 50 == 0 or k == len(files):
                print(f"[{tag}] {k}/{len(files)}  {time.time() - t0:.0f}s", flush=True)

    counts = {}
    for r in log_rows:
        counts[r["status"]] = counts.get(r["status"], 0) + 1
    summary = " ".join(f"{k}={v}" for k, v in sorted(counts.items()))
    print(f"[{tag}] done in {time.time() - t0:.0f}s: {summary}", flush=True)

    if args.state_dir and not args.dry_run:
        logdir = os.path.join(args.state_dir, "hdrfix")
        os.makedirs(logdir, exist_ok=True)
        logpath = os.path.join(logdir, f"{chunk}__{pointing}__{date}.csv")
        fields = ["run_utc", "host", "job", "epoch", "file", "status", "beam", "start_mjd", "pos_mjd",
                  "dt_s", "ra", "dec", "ra_deg", "dec_deg", "prev_ra_deg", "prev_dec_deg",
                  "shift_arcsec", "warn", "detail"]
        new = not os.path.exists(logpath)
        run = dict(run_utc=time.strftime("%Y-%m-%dT%H:%M:%S", time.gmtime()), host=os.uname().nodename,
                   job=os.environ.get("SLURM_JOB_ID", "") + (
                       f"_{os.environ['SLURM_ARRAY_TASK_ID']}" if "SLURM_ARRAY_TASK_ID" in os.environ else ""),
                   epoch=tag)
        with open(logpath, "a", newline="") as fh:
            wr = csv.DictWriter(fh, fieldnames=fields, extrasaction="ignore")
            if new:
                wr.writeheader()
            for r in log_rows:
                wr.writerow({**run, **r})
        print(f"[{tag}] log: {logpath}", flush=True)

    return all(r["status"] in OK_STATUSES for r in log_rows) and bool(log_rows)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("epochs", nargs="*", help="epoch directories (.../<chunk>/<pointing>/<date>)")
    ap.add_argument("--list", help="file with one epoch directory per line")
    ap.add_argument("--task", type=int, help="only do line TASK (0-based) of --list; for Slurm arrays")
    ap.add_argument("--posdir", default=os.environ.get("CRAFTS_POSDIR", DEFAULT_POSDIR))
    ap.add_argument("--state-dir", default=os.environ.get("CRAFTS_STATE_DIR"),
                    help="where per-epoch result logs go (<state>/hdrfix/); default $CRAFTS_STATE_DIR")
    ap.add_argument("--max-dt", type=float, default=1.0, help="max |pos sample - file start| in s (default 1)")
    ap.add_argument("--dry-run", action="store_true", help="compute and report, modify nothing, write no log")
    ap.add_argument("--force", action="store_true", help="re-fix files that already carry a PNTFIX stamp")
    args = ap.parse_args()

    epochs = list(args.epochs)
    if args.list:
        with open(args.list) as f:
            lines = [l.strip() for l in f if l.strip() and not l.startswith("#")]
        if args.task is not None:
            if not 0 <= args.task < len(lines):
                print(f"task {args.task} out of range (list has {len(lines)} lines)", file=sys.stderr)
                return 2
            lines = [lines[args.task]]
        epochs += lines
    if not epochs:
        ap.error("give epoch directories or --list")
    if not args.state_dir and not args.dry_run:
        print("WARNING: no --state-dir / $CRAFTS_STATE_DIR: files will be stamped but no result log kept",
              file=sys.stderr)

    ok = True
    for e in epochs:
        ok &= fix_epoch(e, args)
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
