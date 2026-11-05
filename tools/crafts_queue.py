#!/usr/bin/env python3
"""
crafts_queue.py — what is Slurm doing with our jobs?  Plain text, for `cat`.

  * NOW: your jobs in the queue, by job name: running / pending (with the pending reasons)
  * RECENT: how your jobs ended over the last --since hours (default 24), by job name:
    COMPLETED / FAILED / TIMEOUT / OUT_OF_MEMORY / CANCELLED / ...
  * SUBMISSIONS: the last few entries of $CRAFTS_STATE_DIR/submissions.log, with each
    array's progress (done / running / pending / failed)

Array jobs are counted per task (one task = one epoch directory).

Usage:
  crafts_queue.py                    # print, and write $CRAFTS_STATE_DIR/queue.txt
  crafts_queue.py --watch 60         # refresh every 60 s (Ctrl-C to stop)
  crafts_queue.py --since 72 --user someone
"""

import argparse
import os
import subprocess
import sys
import time
from collections import Counter, defaultdict

OK_END = {"COMPLETED"}
BAD_END = {"FAILED", "TIMEOUT", "OUT_OF_MEMORY", "NODE_FAIL", "BOOT_FAIL", "DEADLINE", "PREEMPTED", "CANCELLED"}


def run(cmd):
    try:
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
    except (OSError, subprocess.TimeoutExpired) as e:
        return None, str(e)
    if r.returncode != 0:
        return None, r.stderr.strip()
    return r.stdout, ""


def base_job(jobid):
    """'12345_7' -> '12345', '12345_[8-99%20]' -> '12345', '12345' -> '12345'"""
    return jobid.split("_")[0].split(".")[0]


def snapshot(user, since_h, state):
    L = []
    w = L.append
    w(f"CRAFTS Slurm status   {time.strftime('%Y-%m-%d %H:%M:%S')}   user={user}   host={os.uname().nodename}")
    w("")

    # ---------------- queue now ----------------
    out, err = run(["squeue", "-u", user, "-h", "-r", "-o", "%i|%j|%T|%M|%r"])
    queue = defaultdict(Counter)              # name -> state counts
    reasons = defaultdict(Counter)            # name -> pending reason counts
    by_array = defaultdict(Counter)           # base job id -> state counts
    if out is None:
        w(f"NOW: squeue failed: {err}")
    else:
        for line in out.splitlines():
            parts = line.split("|")
            if len(parts) < 5:
                continue
            jid, name, st, _, reason = parts
            queue[name][st] += 1
            by_array[base_job(jid)][st] += 1
            if st == "PENDING":
                reasons[name][reason] += 1
        tot = Counter()
        for c in queue.values():
            tot.update(c)
        w(f"NOW  running {tot['RUNNING']}   pending {tot['PENDING']}   other {sum(tot.values()) - tot['RUNNING'] - tot['PENDING']}")
        if queue:
            w(f"  {'job name':<28} {'running':>8} {'pending':>8} {'other':>6}  pending because")
            for name in sorted(queue):
                c = queue[name]
                other = sum(c.values()) - c["RUNNING"] - c["PENDING"]
                why = ", ".join(f"{r}={n}" for r, n in reasons[name].most_common(3))
                w(f"  {name:<28} {c['RUNNING']:>8} {c['PENDING']:>8} {other:>6}  {why}")
    w("")

    # ---------------- how recent jobs ended ----------------
    out, err = run(["sacct", "-u", user, "-X", "-n", "-P", "-S", f"now-{since_h}hours",
                    "-o", "JobID,JobName,State,Elapsed"])
    ended = defaultdict(Counter)
    ended_by_array = defaultdict(Counter)
    elapsed = defaultdict(list)
    if out is None:
        w(f"RECENT: sacct failed: {err}")
    else:
        for line in out.splitlines():
            parts = line.split("|")
            if len(parts) < 4:
                continue
            jid, name, st, el = parts
            st = st.split()[0]                    # 'CANCELLED by 123' -> 'CANCELLED'
            if st in ("RUNNING", "PENDING", "REQUEUED", "SUSPENDED"):
                continue
            ended[name][st] += 1
            ended_by_array[base_job(jid)][st] += 1
            if st == "COMPLETED":
                elapsed[name].append(el)
        w(f"RECENT  (ended in the last {since_h} h)")
        if not ended:
            w("  none")
        else:
            states = sorted({s for c in ended.values() for s in c}, key=lambda s: (s != "COMPLETED", s))
            w(f"  {'job name':<28} " + " ".join(f"{s[:13]:>13}" for s in states) + "   longest OK")
            for name in sorted(ended):
                c = ended[name]
                longest = max(elapsed[name]) if elapsed[name] else "-"
                w(f"  {name:<28} " + " ".join(f"{c[s]:>13}" for s in states) + f"   {longest}")
            nbad = sum(c[s] for c in ended.values() for s in BAD_END)
            if nbad:
                w(f"  !! {nbad} tasks did not complete; their .out files say why")
    w("")

    # ---------------- our submissions ----------------
    sublog = os.path.join(state, "submissions.log") if state else None
    if sublog and os.path.exists(sublog):
        with open(sublog) as f:
            subs = [l.rstrip("\n").split("\t") for l in f if l.strip()][-10:]
        w("SUBMISSIONS  (last 10 from submissions.log; per array task)")
        w(f"  {'submitted (UTC)':<17} {'job':>10} {'stage':<14} {'tasks':>6} {'done':>6} {'failed':>6} {'running':>8} {'pending':>8}")
        for s in subs:
            if len(s) < 4:
                continue
            ts, jid, stage, n = s[:4]
            e, q = ended_by_array.get(jid, Counter()), by_array.get(jid, Counter())
            failed = sum(e[x] for x in BAD_END)
            w(f"  {ts:<17} {jid:>10} {stage:<14} {n:>6} {e['COMPLETED']:>6} {failed:>6} {q['RUNNING']:>8} {q['PENDING']:>8}")
        w("  (done/failed only count tasks that ended within --since; widen it for old submissions)")
        w("")
    return "\n".join(L) + "\n"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--user", default=os.environ.get("USER"))
    ap.add_argument("--since", type=int, default=24, help="hours of job history (default 24)")
    ap.add_argument("--state-dir", default=os.environ.get("CRAFTS_STATE_DIR"))
    ap.add_argument("--watch", type=int, metavar="SEC", help="refresh every SEC seconds")
    args = ap.parse_args()
    while True:
        text = snapshot(args.user, args.since, args.state_dir)
        if args.state_dir:
            tmp = os.path.join(args.state_dir, "queue.txt.tmp")
            with open(tmp, "w") as f:
                f.write(text)
            os.replace(tmp, os.path.join(args.state_dir, "queue.txt"))
        if args.watch:
            sys.stdout.write("\033[2J\033[H")
        sys.stdout.write(text)
        sys.stdout.flush()
        if not args.watch:
            return
        time.sleep(args.watch)


if __name__ == "__main__":
    main()
