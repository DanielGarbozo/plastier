#!/bin/bash
#
# Phase 3 scale-up: run the whole pipeline (stages 1-6) on the ~20 candidate S. aureus
# runs in assets/phase3/sra_ids_candidates.txt (docs/phase3_candidates.md), on Slurm.
#
# Submit from the repository root, passing your cluster's account/partition:
#
#     mkdir -p logs
#     sbatch --account=<your_account> --partition=<partition> assets/phase3/run_phase3.sbatch
#
# Optional environment overrides (export before sbatch):
#     PLASTIER_SCRATCH        base dir for the work dir and results (default: /scratch/$USER)
#     PLASTIER_MODULES        `module load` line for Nextflow + Java (default: nextflow/26.04.0 openjdk/21.0.3)
#     PLASTIER_SRA_IDS        accession list to use instead of the candidates
#     PLASTIER_MAX_TASK_TIME  longest walltime a pipeline task may request (default: 6.h, the
#                             limit of NYU's cpu_short QOS; raise it on partitions that allow more)
#
# The head job only orchestrates; every task is its own Slurm job (conf/slurm.config). If the
# head job's walltime runs out first, re-submit this same script: -resume skips finished work.

#SBATCH --job-name=plastier-phase3
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --output=logs/phase3_%j.log

set -euo pipefail
cd "${SLURM_SUBMIT_DIR}"

SCRATCH="${PLASTIER_SCRATCH:-/scratch/${USER}}"

# shellcheck disable=SC2086
module load ${PLASTIER_MODULES:-nextflow/26.04.0 openjdk/21.0.3}

# assets/slurm/node_scratch_bind.config: see its header; drop -c on clusters that do not need it.
nextflow run . \
  -profile full,singularity,slurm \
  -c assets/slurm/node_scratch_bind.config \
  --sra_ids "${PLASTIER_SRA_IDS:-assets/phase3/sra_ids_candidates.txt}" \
  --slurm_account "${SLURM_JOB_ACCOUNT}" \
  --slurm_partition "${SLURM_JOB_PARTITION}" \
  --slurm_max_time "${PLASTIER_MAX_TASK_TIME:-6.h}" \
  -work-dir "${SCRATCH}/plastier_work_phase3" \
  --outdir "${SCRATCH}/plastier_phase3_results" \
  -resume
