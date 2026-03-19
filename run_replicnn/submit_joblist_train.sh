#!/usr/bin/env bash
#SBATCH --job-name=replicnn_train
#SBATCH --nodelist=gpu1
#SBATCH --cpus-per-task=4
#SBATCH --mem=10G
#SBATCH --array=1-1177%11
#SBATCH --output=logs/train/job_%A_%a.out
#SBATCH --error=logs/train/job_%A_%a.err

SECONDS=0

JOB_FILE="joblist_train.tsv"
LOG_DIR="logs/train"
mkdir -p "$LOG_DIR"

# Skip the header (tail -n +2) so array indexing matches job_id order
LINE=$(tail -n +2 "$JOB_FILE" | sed -n "${SLURM_ARRAY_TASK_ID}p")

# Parse columns from joblist
job_id=$(echo "$LINE" | cut -f1)
train_org=$(echo "$LINE" | cut -f2)
train_exp=$(echo "$LINE" | cut -f3)
train_sample=$(echo "$LINE" | cut -f4)
left_out_chrom=$(echo "$LINE" | cut -f5)
chrom_file=$(echo "$LINE" | cut -f6)
model_name=$(echo "$LINE" | cut -f7)
input_file=$(echo "$LINE" | cut -f8)

# Create experiment-specific output folder
out_dir="models/${model_name}"
mkdir -p "$out_dir"

# Temporary filtered file (exclude left_out_chrom if needed)
tmp_sdf=$(mktemp)
if [ "$left_out_chrom" != "none" ]; then
    grep -v "^${left_out_chrom}[[:space:]]" "$input_file" > "$tmp_sdf"
else
    cp "$input_file" "$tmp_sdf"
fi

# Sanity check: non-empty file
if [ ! -s "$tmp_sdf" ]; then
    echo "ERROR: filtered SDF file is empty for job ${job_id} (chrom=${left_out_chrom})"
    exit 1
fi

# Run RepliCNN training (stdout/stderr already go to SLURM log files)
replicnn train \
    -i "$tmp_sdf" \
    -o "${out_dir}" \
    -g \
	-ws 201 \
	-e 300 \
	-bs 8192 \
	-v 0.1 \
	-lr 0.001

# Calculate runtime
duration=$SECONDS
hours=$((duration/3600))
minutes=$(((duration%3600)/60))
seconds=$((duration%60))

# Write metadata summary to same SLURM log file
{
  echo "=========================================="
  echo "$(date '+%Y-%m-%d %H:%M:%S') | SLURM Task ${SLURM_ARRAY_TASK_ID}/${SLURM_ARRAY_JOB_ID}"
  echo "Job ID: ${job_id}"
  echo "Model: ${model_name}"
  echo "Organism: ${train_org}"
  echo "Experiment: ${train_exp}"
  echo "Sample: ${train_sample}"
  echo "Chromosome left out: ${left_out_chrom}"
  echo "Chromsizes file: ${chrom_file}"
  echo "Input file: ${input_file}"
  echo "Filtered SDF: ${tmp_sdf}"
  echo "Duration: ${hours}h ${minutes}m ${seconds}s"
} >> "logs/train/job_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.out"
