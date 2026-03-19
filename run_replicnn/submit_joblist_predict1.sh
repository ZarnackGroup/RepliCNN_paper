#!/usr/bin/env bash
set -euo pipefail

LOG_DIR="logs/predict"
mkdir -p "$LOG_DIR"

#"joblists_split/joblist_predict_1.tsv"
#"joblists_split/joblist_predict_2.tsv"
#"joblists_split/joblist_predict_5.tsv"
#"joblists_split/joblist_predict_4.tsv"

# Define your split joblists here
JOBLISTS=(
  "joblists_split/joblist_predict_3.tsv"
)

for JOB_FILE in "${JOBLISTS[@]}"; do
    # Count number of tasks (subtract 1 for header)
    NUM_TASKS=$(( $(wc -l < "$JOB_FILE") - 1 ))

    echo "Submitting predictions for $JOB_FILE with $NUM_TASKS tasks..."
	
	sleep 1

    sbatch \
        --job-name="replicnn_predict_$(basename "$JOB_FILE" .tsv)" \
		--nodelist=gpu1 \
		--gres=gpu:1 \
        --cpus-per-task=2 \
        --mem=4G \
        --array=0-1890%1 \
        --output="$LOG_DIR/job_%A_%a.out" \
        --error="$LOG_DIR/job_%A_%a.err" \
        run_predict_jobarray.sh "$JOB_FILE"
	sleep 1
done
