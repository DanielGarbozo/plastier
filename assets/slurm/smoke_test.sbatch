#!/bin/bash
#
# Smoke test: run the `test` profile end-to-end on Slurm.
#
# Submit from the repository root, passing your cluster's account/partition on
# the command line (they are deliberately not hardcoded here):
#
#     mkdir -p logs
#     sbatch --account=<your_account> --partition=<partition> assets/slurm/smoke_test.sbatch
#
# Optional environment overrides (export before sbatch):
#     PLASTIER_SCRATCH   base dir for the work dir and results   (default: /scratch/$USER)
#     PLASTIER_MODULES   `module load` line for Nextflow + Java  (default: nextflow/26.04.0 openjdk/21.0.3)
#
# The head job only orchestrates; each pipeline task is submitted as its own Slurm
# job by conf/slurm.config. Re-submit the same script to resume (-resume is set).

#SBATCH --job-name=plastier-smoke
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --output=logs/smoke_%j.log

set -euo pipefail
cd "${SLURM_SUBMIT_DIR}"

SCRATCH="${PLASTIER_SCRATCH:-/scratch/${USER}}"

# shellcheck disable=SC2086
module load ${PLASTIER_MODULES:-nextflow/26.04.0 openjdk/21.0.3}

# assets/slurm/node_scratch_bind.config binds the per-job node-local scratch dir
# (/state/partition1) into Apptainer; needed on clusters that point $TMPDIR there.
# Drop the -c line on clusters where that is not the case.
nextflow run . \
  -profile test,singularity,slurm \
  -c assets/slurm/node_scratch_bind.config \
  --slurm_account "${SLURM_JOB_ACCOUNT}" \
  --slurm_partition "${SLURM_JOB_PARTITION}" \
  -work-dir "${SCRATCH}/plastier_work_smoke" \
  --outdir "${SCRATCH}/plastier_smoke_results" \
  -resume
