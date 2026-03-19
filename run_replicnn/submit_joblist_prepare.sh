#!/bin/bash
#SBATCH --job-name=replicnn_prepare
#SBATCH --array=1-62%10
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G

set -euo pipefail

TASKS_FILE="joblist_prepare.tsv"

# Skip header, get the line for this array index
LINE=$(awk -v task_id=$SLURM_ARRAY_TASK_ID 'NR==task_id+1' "$TASKS_FILE")

if [[ -z "$LINE" ]]; then
  echo "ERROR: No line for task $SLURM_ARRAY_TASK_ID"
  exit 1
fi

# Parse the line into variables
IFS=$'\t' read -r job_id organism experiment_type sample_id \
    fwd_path rev_path binsize chrom_size_all path_to_tsv timing_path invert <<< "$LINE"

# Build the command
CMD="replicnn prepare \
  -fwd $fwd_path \
  -rev $rev_path \
  -bs $binsize \
  -cs $chrom_size_all \
  -o $path_to_tsv \
  "

# Add timing if available
if [[ -n "$timing_path" && "$timing_path" != "NA" ]]; then
  CMD="$CMD -t $timing_path"
fi

# Add invert flag if true/1
if [[ "$invert" == "true" || "$invert" == "1" ]]; then
  CMD="$CMD -i"
fi

echo "[$(date)] Running job $job_id ($sample_id)..."
echo "$CMD"

eval $CMD
