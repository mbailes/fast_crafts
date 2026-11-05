#!/bin/bash
# Submit header fixes for every epoch the last audit says still needs them.
# One Slurm array task per epoch directory, at most $MAXPAR running at once.
#
#   cd ~/fast_crafts
#   python3 tools/crafts_audit.py /aphid/temp/fast_crafts/data     # writes the lists
#   slurm/submit_hdrfix.sh                                         # submits epochs_hdrfix_todo.txt
#
# Options (environment):  LIST=other_list.txt   MAXPAR=20   FORCE=1 (submit even if hdrfix jobs are queued)
# Records every submission in $CRAFTS_STATE_DIR/submissions.log; Slurm logs go to
# $CRAFTS_STATE_DIR/slurm_logs/hdrfix/<jobid>_<task>.out
set -euo pipefail

: "${CRAFTS_STATE_DIR:?set CRAFTS_STATE_DIR to the persistent state directory}"
REPO="$(cd "$(dirname "$0")/.." && pwd)"
LIST="${LIST:-$CRAFTS_STATE_DIR/epochs_hdrfix_todo.txt}"
MAXPAR="${MAXPAR:-20}"
JOBNAME=hdrfix_v3

[ -s "$LIST" ] || { echo "nothing to submit: $LIST is empty or missing (run the audit first)"; exit 0; }

queued=$(( $(squeue -u "$USER" -h -r -n "$JOBNAME" | wc -l) ))
if [ "$queued" -gt 0 ] && [ "${FORCE:-0}" != 1 ]; then
    echo "$queued $JOBNAME tasks are already queued/running for $USER; wait for them and re-audit,"
    echo "or set FORCE=1 (stamped files are skipped, so a double submission is safe but wasteful)."
    exit 1
fi

ts=$(date -u +%Y%m%dT%H%M%SZ)
mkdir -p "$CRAFTS_STATE_DIR/submissions" "$CRAFTS_STATE_DIR/slurm_logs/hdrfix"
FROZEN="$CRAFTS_STATE_DIR/submissions/hdrfix_${ts}.txt"
cp "$LIST" "$FROZEN"                       # the array indexes into this frozen copy
N=$(( $(wc -l < "$FROZEN") ))

cd "$REPO"                                 # the .slurm finds tools/ via SLURM_SUBMIT_DIR
jobid=$(sbatch --parsable --job-name="$JOBNAME" --array="0-$((N - 1))%${MAXPAR}" \
        --output="$CRAFTS_STATE_DIR/slurm_logs/hdrfix/%A_%a.out" \
        --export=ALL,LIST="$FROZEN" slurm/fix_pointing_v3.slurm)
jobid=${jobid%%;*}

printf "%s\t%s\t%s\t%s\t%s\n" "$ts" "$jobid" "$JOBNAME" "$N" "$FROZEN" >> "$CRAFTS_STATE_DIR/submissions.log"
echo "submitted $JOBNAME job $jobid: $N epochs, max $MAXPAR at once"
echo "  list : $FROZEN"
echo "  logs : $CRAFTS_STATE_DIR/slurm_logs/hdrfix/${jobid}_*.out"
echo "  watch: python3 tools/crafts_queue.py      then re-run the audit when it has drained"
