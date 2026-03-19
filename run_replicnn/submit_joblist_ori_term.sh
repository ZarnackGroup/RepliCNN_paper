#!/bin/bash
#SBATCH --job-name=replicnn_rfd_oem
#SBATCH --array=1-62
#SBATCH --time=20:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=12G
#SBATCH --nodelist=jupiter3

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

if [[ "$organism" != "yeast" ]]; then
    echo "[$(date)] Skipping sample $sample_id (organism=$organism)"
    exit 0
fi

case "$organism" in
  human|mouse)
    ori_threshold=0.05
    ter_threshold=0.15
	window_radius=15000
	max_merge_size=15000
	n_evidence=2
	cutoff=15
	eval_resolution=75000
	resolutions=(50000 75000 100000 150000)
    stride=10
    ;;
  yeast)
    ori_threshold=0.01
    ter_threshold=0.1
	window_radius=5000
	max_merge_size=5000
	n_evidence=2
	cutoff=10
	eval_resolution=10000
	resolutions=(5000 10000 15000)
    stride=1
    ;;
  *)
    exit 1
    ;;
esac

output_prefix="../../../data/oem_rfd/${sample_id}/${sample_id}"

input=""
for res in "${resolutions[@]}"; do
    for type in rfd oem; do
        file="${output_prefix}_${type}_${stride}_${res}.bw"
        if [[ -f "$file" ]]; then
            input+=" ${file}"
        else
            echo "Warning: missing file $file" >&2
        fi
    done
done

# Trim leading space
input="${input# }"

CMD=(replicnn ori_ter
	--input ${input}
	--chromsizes ${chrom_size_all}
	--output_prefix ${output_prefix}
	--save_intermediates
	--ori-threshold ${ori_threshold}
	--ter-threshold ${ter_threshold}
	--window-radius ${window_radius}
	--max-merge-size ${max_merge_size}
	--n-evidence ${n_evidence}
	--eval_resolution ${eval_resolution}
	--cutoff ${cutoff})

echo "[$(date)] Generated CMD: ${CMD[@]}"
srun "${CMD[@]}"
