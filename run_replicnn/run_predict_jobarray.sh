#!/usr/bin/env bash
#SBATCH --cpus-per-task=4
#SBATCH --mem=4G
#SBATCH --output=logs/predict/job_%A_%a.out
#SBATCH --error=logs/predict/job_%A_%a.err

set -euo pipefail

JOB_FILE="$1"

SECONDS=0

# Skip header
LINE=$(tail -n +2 "$JOB_FILE" | sed -n "${SLURM_ARRAY_TASK_ID}p")

# Parse task info
pred_job_id=$(echo "$LINE" | cut -f1)
model_path=$(echo "$LINE" | cut -f2)
input_file=$(echo "$LINE" | cut -f3)
out_file=$(echo "$LINE" | cut -f4)

mkdir -p "$(dirname "$out_file")"

tmp_sdf=$(mktemp)

if [[ "$model_path" == *"_delta_"* ]]; then
    left_out_chrom=$(basename "$model_path" | sed -E 's/.*_delta_(.*)/\1/')
    grep -v "^${left_out_chrom}[[:space:]]" "$input_file" > "$tmp_sdf"
else
    cp "$input_file" "$tmp_sdf"
fi

if [ ! -s "$tmp_sdf" ]; then
    echo "ERROR: temporary SDF is empty for job ${pred_job_id}"
    exit 1
fi

replicnn predict \
    -i "$tmp_sdf" \
    -m "$model_path" \
    -o "$out_file" \
    -g \
    -nl

rm "$tmp_sdf"

duration=$SECONDS
hours=$((duration/3600))
minutes=$(((duration%3600)/60))
seconds=$((duration%60))

{
  echo "=========================================="
  echo "$(date '+%Y-%m-%d %H:%M:%S') | SLURM Task ${SLURM_ARRAY_TASK_ID}/${SLURM_ARRAY_JOB_ID}"
  echo "Prediction Job ID: ${pred_job_id}"
  echo "Model: ${model_path}"
  echo "Input file: ${input_file}"
  echo "Output file: ${out_file}"
  echo "Duration: ${hours}h ${minutes}m ${seconds}s"
} >> "logs/predict/job_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.out"
