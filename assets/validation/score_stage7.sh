#!/bin/bash
#
# Score a finished stage 7 run: per-tier metrics (#32), single-tool baseline (#34) and
# PlasEval (#33). Reads only <OUTDIR>, so it can be re-run without re-running Nextflow.
#
#     assets/validation/score_stage7.sh <OUTDIR>
#
# <OUTDIR> is the --outdir of a `-profile validation ... --assemblies` run (it must hold
# tier/*.tier_resolution.tsv and mobsuite/<genome>/contig_report.txt). Results go to
# <OUTDIR>/stage7/. Needs python3 with pandas; PlasEval's own dependencies (bidict, psutil,
# networkx<3.4) are pip-installed with --user if missing.
#
# Optional environment overrides:
#     PLASEVAL_DIR   an existing PlasEval checkout (default: cloned next to <OUTDIR>)

set -euo pipefail

OUTDIR="$(cd "${1:?usage: score_stage7.sh <OUTDIR>}" && pwd)"
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
GT="${REPO}/assets/validation/closed_genome_ground_truth.csv"
PLASEVAL_DIR="${PLASEVAL_DIR:-$(dirname "${OUTDIR}")/PlasEval}"
PLASEVAL_COMMIT=8238356   # the version this scoring was tested against
OUT="${OUTDIR}/stage7"
mkdir -p "${OUT}"

shopt -s nullglob
TIER_FILES=("${OUTDIR}"/tier/*.tier_resolution.tsv)
[ ${#TIER_FILES[@]} -gt 0 ] || { echo "No tier/*.tier_resolution.tsv in ${OUTDIR}" >&2; exit 1; }
GENOMES=()
for f in "${TIER_FILES[@]}"; do GENOMES+=("$(basename "${f}" .tier_resolution.tsv)"); done

echo "== #32 per-tier precision/recall/F1 (${#GENOMES[@]} genomes)"
python3 "${REPO}/bin/evaluate_metrics.py" \
  --predictions "${TIER_FILES[@]}" --ground-truth "${GT}" \
  --output "${OUT}/tier_metrics.tsv" --per-arg-output "${OUT}/per_arg.tsv"

echo "== #34 single-tool baseline vs tier framework"
python3 "${REPO}/bin/baseline_comparison.py" \
  --predictions "${TIER_FILES[@]}" --ground-truth "${GT}" --output "${OUT}/baseline_comparison.tsv"

echo "== #33 PlasEval bin-level comparison (MOB-suite bins vs closed-genome plasmids)"
if [ ! -d "${PLASEVAL_DIR}/.git" ]; then
  git clone -q https://github.com/cchauve/PlasEval.git "${PLASEVAL_DIR}"
  git -C "${PLASEVAL_DIR}" checkout -q "${PLASEVAL_COMMIT}"
fi
python3 -c "import bidict, psutil, networkx" 2>/dev/null \
  || pip install --user -q -r "${PLASEVAL_DIR}/requirements.txt"

REPORTS=()
for g in "${GENOMES[@]}"; do
  [ -f "${OUTDIR}/mobsuite/${g}/contig_report.txt" ] \
    || { echo "Missing ${OUTDIR}/mobsuite/${g}/contig_report.txt" >&2; exit 1; }
  REPORTS+=("${OUTDIR}/mobsuite/${g}/contig_report.txt")
done
python3 "${REPO}/assets/validation/plaseval_convert.py" ground-truth --input "${GT}" --only "${GENOMES[@]}" --output "${OUT}/gt_bins.tsv"
python3 "${REPO}/assets/validation/plaseval_convert.py" predicted --input "${REPORTS[@]}" --output "${OUT}/pred_bins.tsv"
( cd "${PLASEVAL_DIR}/src" && python3 plaseval.py eval \
    --pred "${OUT}/pred_bins.tsv" --gt "${OUT}/gt_bins.tsv" \
    --out_file "${OUT}/plaseval_eval.tsv" --log_file "${OUT}/plaseval_eval.log" )

echo "Done. Results in ${OUT}/:"
ls -1 "${OUT}"
