#!/bin/bash
#SBATCH --job-name=clademers
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=2-00:00:00
#SBATCH --output=clademers-%j.log

# Edit paths, account/partition, and resource requests for your HPC.
# RAM/time above are illustrative, not measured full-NT requirements.
# Compile CladeMers on the HPC before submitting this job.
# Progress/ETA appears in this job's log; add -progress-seconds 10 for faster updates.
# Load your site's BLAST+ module here if needed.
# Compressed taxonomy is supported: use nodes.dmp.gz/names.dmp.gz below,
# with gzip on PATH (or add -gzip /path/to/gzip to the command).
set -euo pipefail

CLADEMERS=/path/to/clademers
NT_PREFIX=/path/to/nt
TAXONOMY=/path/to/taxonomy
RESULT=/scratch/yourname/placozoa_k42_${SLURM_JOB_ID}

# 30 scanner workers plus room for the producer and BLAST reader.
"$CLADEMERS" \
  -taxnodes "$TAXONOMY/nodes.dmp" \
  -taxnames "$TAXONOMY/names.dmp" \
  -taxid 10226 -k 42 -t 30 \
  -blastdb "$NT_PREFIX" -out "$RESULT"
