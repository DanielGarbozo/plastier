#!/bin/bash
#
# Stage 7e sensitivity analysis (issue #35) - DIAGNOSTIC ONLY, never an accuracy figure.
# Simulates reads from closed genomes at several depths / base qualities / plasmid copy
# numbers, runs the whole pipeline on them (assembly included), and compares each gene's
# tier with its tier in the closed-genome run (bin/sensitivity_summary.py).
#
# Needs a finished assets/validation/run_stage7.sbatch run (its tier files are the
# reference, and stage7/per_arg.tsv picks the genomes). Submit from the repository root:
#
#     sbatch --account=<your_account> --partition=<partition> assets/validation/run_sensitivity.sbatch
#
# Optional environment overrides (export before sbatch):
#     PLASTIER_SCRATCH        base dir (default: /scratch/$USER); must match the stage 7 run
#     PLASTIER_MODULES        `module load` line (default: nextflow/26.04.0 openjdk/21.0.3 anaconda3/2025.06)
#     PLASTIER_SIM_GENOMES    space-separated genome accessions (default: the 5 genomes with the
#                             most ARGs truly on a plasmid, from stage7/per_arg.tsv)
#     PLASTIER_SIM_CONDITIONS d<depth>_q<qualityshift>_p<plasmidcopies> list, see simulate_reads.py
#                             (default: d10_q0_p1 d20_q0_p1 d50_q0_p1 d20_q-15_p1 d20_q0_p3)
#     PLASTIER_MAX_TASK_TIME  longest task walltime (default: 6.h, NYU cpu_short QOS limit)
#
# Re-submitting resumes (-resume).

#SBATCH --job-name=plastier-sensitivity
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --output=logs/sensitivity_%j.log

set -euo pipefail
cd "${SLURM_SUBMIT_DIR}"

SCRATCH="${PLASTIER_SCRATCH:-/scratch/${USER}}"
REFS="${SCRATCH}/plastier_refs"
REFOUT="${SCRATCH}/plastier_stage7_results"
SIMDIR="${SCRATCH}/plastier_sim_reads"
SIMOUT="${SCRATCH}/plastier_sensitivity_results"

# shellcheck disable=SC2086
module load ${PLASTIER_MODULES:-nextflow/26.04.0 openjdk/21.0.3 anaconda3/2025.06}

if [ -n "${PLASTIER_SIM_GENOMES:-}" ]; then
  GENOMES="${PLASTIER_SIM_GENOMES}"
else
  GENOMES="$(python3 - "${REFOUT}/stage7/per_arg.tsv" <<'PY'
import sys
import pandas as pd
d = pd.read_csv(sys.argv[1], sep="\t")
top = d[d.true_location == "Plasmid"].groupby("genome_id").size().sort_values(ascending=False).head(5)
if top.empty:
    sys.exit("No ARG is truly on a plasmid in per_arg.tsv - pass PLASTIER_SIM_GENOMES explicitly")
print(" ".join(top.index))
PY
)"
fi
echo "Genomes: ${GENOMES}"

python3 assets/validation/simulate_reads.py \
  --assemblies "${REFS}/assemblies.csv" --outdir "${SIMDIR}" \
  --genomes ${GENOMES} \
  --conditions ${PLASTIER_SIM_CONDITIONS:-d10_q0_p1 d20_q0_p1 d50_q0_p1 d20_q-15_p1 d20_q0_p3}

nextflow run . \
  -profile validation,singularity,slurm \
  -c assets/slurm/node_scratch_bind.config \
  --input "${SIMDIR}/sim_samplesheet.csv" \
  --slurm_account "${SLURM_JOB_ACCOUNT}" \
  --slurm_partition "${SLURM_JOB_PARTITION}" \
  --slurm_max_time "${PLASTIER_MAX_TASK_TIME:-6.h}" \
  -work-dir "${SCRATCH}/plastier_work_sensitivity" \
  --outdir "${SIMOUT}" \
  -resume

mkdir -p "${SIMOUT}/stage7"
REFERENCE=()
for g in ${GENOMES}; do REFERENCE+=("${REFOUT}/tier/${g}.tier_resolution.tsv"); done
python3 bin/sensitivity_summary.py \
  --reference "${REFERENCE[@]}" \
  --simulated "${SIMOUT}"/tier/*.tier_resolution.tsv \
  --output "${SIMOUT}/stage7/sensitivity_diagnostic.tsv" \
  --detail-output "${SIMOUT}/stage7/sensitivity_per_gene.tsv"
echo "DIAGNOSTIC ONLY. Summary: ${SIMOUT}/stage7/sensitivity_diagnostic.tsv"
