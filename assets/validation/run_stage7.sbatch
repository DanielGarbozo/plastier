#!/bin/bash
#
# Stage 7 validation: fetch the closed reference genomes, run stages 3-6 on them on
# Slurm, then score the tier calls against the ground truth (issue #32).
#
# Submit from the repository root, passing your cluster's account/partition:
#
#     mkdir -p logs
#     sbatch --account=<your_account> --partition=<partition> assets/validation/run_stage7.sbatch
#
# Optional environment overrides (export before sbatch):
#     PLASTIER_SCRATCH   base dir for references, work dir and results (default: /scratch/$USER)
#     PLASTIER_MODULES   `module load` line for Nextflow + Java + a python3 with pandas
#                        (default: nextflow/26.04.0 openjdk/21.0.3 anaconda3/2025.06)
#     PLASTIER_ONLY      space-separated genome accessions to run instead of all of them
#     PLASTIER_MAX_TASK_TIME  longest walltime a pipeline task may request (default: 6.h, the
#                        limit of NYU's cpu_short QOS; raise it on partitions that allow more)
#
# Re-submitting resumes (-resume): finished tasks and downloaded genomes are reused.

#SBATCH --job-name=plastier-stage7
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --output=logs/stage7_%j.log

set -euo pipefail
cd "${SLURM_SUBMIT_DIR}"

SCRATCH="${PLASTIER_SCRATCH:-/scratch/${USER}}"
REFS="${SCRATCH}/plastier_refs"
OUTDIR="${SCRATCH}/plastier_stage7_results"

# shellcheck disable=SC2086
module load ${PLASTIER_MODULES:-nextflow/26.04.0 openjdk/21.0.3 anaconda3/2025.06}

python3 assets/validation/fetch_reference_genomes.py \
  --outdir "${REFS}" \
  --samplesheet "${REFS}/assemblies.csv" \
  ${PLASTIER_ONLY:+--only ${PLASTIER_ONLY}}

# assets/slurm/node_scratch_bind.config: see its header; drop -c on clusters that do not need it.
nextflow run . \
  -profile validation,singularity,slurm \
  -c assets/slurm/node_scratch_bind.config \
  --assemblies "${REFS}/assemblies.csv" \
  --slurm_account "${SLURM_JOB_ACCOUNT}" \
  --slurm_partition "${SLURM_JOB_PARTITION}" \
  --slurm_max_time "${PLASTIER_MAX_TASK_TIME:-6.h}" \
  -work-dir "${SCRATCH}/plastier_work_stage7" \
  --outdir "${OUTDIR}" \
  -resume

mkdir -p "${OUTDIR}/stage7"
python3 bin/evaluate_metrics.py \
  --predictions "${OUTDIR}"/tier/*.tier_resolution.tsv \
  --ground-truth assets/validation/closed_genome_ground_truth.csv \
  --output "${OUTDIR}/stage7/tier_metrics.tsv" \
  --per-arg-output "${OUTDIR}/stage7/per_arg.tsv"

echo "Per-tier metrics: ${OUTDIR}/stage7/tier_metrics.tsv"
