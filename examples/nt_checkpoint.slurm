#!/bin/bash
#SBATCH --job-name=clademers
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=2-00:00:00
#SBATCH --signal=B:USR1@600
#SBATCH --output=clademers-%j.log

# Edit paths, account/partition, modules, and resource requests for your HPC.
# RAM/time are placeholders. Allow enough signal lead time to drain workers
# and write a checkpoint; large target tables can take more than ten minutes.
# B: signals this shell; the trap forwards USR1 only to CladeMers, not its readers.
set -euo pipefail

CLADEMERS=/path/to/clademers
NT_PREFIX=/path/to/nt
BLAST_TAXONOMY=/path/to/blast-taxonomy
# Use shared, persistent storage and the SAME result path across submissions.
# Do not append SLURM_JOB_ID, and do not run two jobs for this path concurrently.
RESULT=/scratch/yourname/placozoa_k42

if [[ -e "$RESULT/SUCCESS" ]]; then
  echo "This result directory is already complete: $RESULT"
  exit 0
fi
if [[ -f "$RESULT/checkpoint/state.bin" ]]; then
  args=(-resume "$RESULT" -t 30 -checkpoint-seconds 300)
else
  args=(-blasttax "$BLAST_TAXONOMY" -taxid 10226 -k 42 -t 30
        -blastdb "$NT_PREFIX" -out "$RESULT" -checkpoint-seconds 300)
  # For dumps, replace -blasttax and its path with -taxnodes/-taxnames paths.
  # LCA labels are always saved. Add -rd for the summary and -wsf for per-taxon FASTAs.
  # Use -wsf-gz instead of -wsf for gzip-compressed per-taxon FASTAs.
  # Or after completion: clademers -rd -results "$RESULT" -out /path/to/report
fi

child=0
interrupted=0
forward_stop() {
  interrupted=1
  if (( child > 0 )); then
    kill -USR1 "$child" 2>/dev/null || true
  fi
}
trap forward_stop USR1 TERM INT
"$CLADEMERS" "${args[@]}" &
child=$!
if (( interrupted )); then forward_stop; fi

# A signal interrupts Bash's wait. Keep waiting until CladeMers has finished
# saving and exited, rather than ending the allocation during its checkpoint.
set +e
while true; do
  interrupted=0
  wait "$child"
  status=$?
  if (( interrupted )); then continue; fi
  break
done
if (( status == 75 )); then
  echo "CladeMers paused; inspect the log, then resubmit this same script to resume."
fi
exit "$status"
