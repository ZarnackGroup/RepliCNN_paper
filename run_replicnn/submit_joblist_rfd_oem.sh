#!/bin/bash
#SBATCH --job-name=replicnn_rfd_oem_test
#SBATCH --array=1-62%10
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=10G
#SBATCH --exclude=jupiter3

set -euo pipefail

TASKS_FILE="../all_cross_predict/joblist_prepare.tsv"

# Get the line for this array index (skip header)
LINE=$(awk -v task_id=$SLURM_ARRAY_TASK_ID 'NR==task_id+1 {print}' "$TASKS_FILE")

if [[ -z "$LINE" ]]; then
  exit 1
fi

# Parse the line into variables
IFS=$'\t' read -r job_id organism experiment_type sample_id \
    fwd_path rev_path binsize chrom_size_all path_to_tsv timing_path invert <<< "$LINE"

# Trim possible whitespace
organism=$(echo "$organism" | xargs)
invert=$(echo "$invert" | xargs)

# Output prefix
OUTPREFIX="../../../data/oem_rfd/${sample_id}/${sample_id}"

# Set resolutions and stride based on organism
case "$organism" in
  human|mouse)
    resolutions=(50000 75000 100000 150000)
    stride=10
    ;;
  yeast)
    resolutions=(2500 5000 10000 15000)
    stride=1
    ;;
  *)
    exit 1
    ;;
esac

# Submit RFD/OEM jobs
for res in "${resolutions[@]}"; do
    for task in rfd oem; do
        CMD=(replicnn rfd_oem
			-w "$fwd_path"
			-c "$rev_path"
			-cs "$chrom_size_all"
			-o "$OUTPREFIX"
			-res $res
			-st $stride
			-t $task)

        # Append invert flag if true
		if [[ "$invert" == "true" || "$invert" == "1" ]]; then
			CMD+=(-inv)
		fi

        srun "${CMD[@]}"
    done
done
