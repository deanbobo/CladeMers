#!/bin/bash
#SBATCH --job-name=clademers
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=2-00:00:00
#SBATCH --output=clademers-%j.log

# BLAST taxonomy mode: no nodes.dmp/names.dmp download is needed.
# Edit paths, account/partition, and resource requests for your HPC.
# RAM/time above are illustrative, not measured full-NT requirements.
# Compile CladeMers for the compute architecture before submitting this job.
# Progress/ETA appears in this job's log; add -progress-seconds 10 for faster updates.
# Load your site's BLAST+ and sqlite3 modules here if needed.
set -euo pipefail

CLADEMERS=/path/to/clademers
NT_PREFIX=/path/to/nt
# Directory containing taxonomy4blast.sqlite3, taxdb.bti, and taxdb.btd.
BLAST_TAXONOMY=/path/to/blast-taxonomy
RESULT=/scratch/yourname/placozoa_k42_${SLURM_JOB_ID}

# 30 scanner workers plus room for the producer and BLAST reader.
"$CLADEMERS" \
  -blasttax "$BLAST_TAXONOMY" \
  -taxid 10226 -k 42 -t 30 \
  -blastdb "$NT_PREFIX" -out "$RESULT"
