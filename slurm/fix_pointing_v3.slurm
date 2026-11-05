#!/bin/bash
# One array task per epoch directory listed in $LIST (one path per line).
#
#   ROOT=/aphid/temp/fast_crafts/data
#   awk -F'\t' -v R=$ROOT 'NR>1{print R"/"$1"/"$2"/"$3}' $CRAFTS_STATE_DIR/audit_epochs.tsv > epochs_hdrfix.txt
#   N=$(wc -l < epochs_hdrfix.txt)
#   sbatch --array=0-$((N-1))%20 --export=ALL,LIST=$PWD/epochs_hdrfix.txt slurm/fix_pointing_v3.slurm
#
# Already-stamped files are skipped, so re-submitting the whole list is safe.
# Per-file results: $CRAFTS_STATE_DIR/hdrfix/<chunk>__<pointing>__<date>.csv
#
#SBATCH --job-name=hdrfix_v3
#SBATCH --account=oz002
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=02:00:00
#SBATCH --output=hdrfix_%A_%a.out

set -u
: "${LIST:?set LIST=/path/to/epoch_list.txt (sbatch --export=ALL,LIST=...)}"
: "${CRAFTS_STATE_DIR:?set CRAFTS_STATE_DIR to the persistent state directory}"

source /fred/oz002/users/vmandaog/vmandaog_env/bin/activate   # needs numpy + astropy

TOOLS="${SLURM_SUBMIT_DIR}/tools"          # submit from the fast_crafts checkout
echo "task ${SLURM_ARRAY_TASK_ID} on $(hostname) at $(date -u +%FT%TZ)"
python3 "${TOOLS}/fix_pointing_v3.py" --list "$LIST" --task "${SLURM_ARRAY_TASK_ID}" \
        --state-dir "$CRAFTS_STATE_DIR"
rc=$?
echo "exit ${rc} at $(date -u +%FT%TZ)"
exit $rc
