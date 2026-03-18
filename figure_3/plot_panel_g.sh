#!/bin/bash

SIF=~/apptainer/deeptools_3.5.6_pyhdfd78af_0.sif
BASE=PATH/scripts/04_plots_analyses/overlap_ori_term/rerun_rfd_oem_ori_ter_hct116
BIN_DIR=PATH/scripts/04_plots_analyses/iz_tss_heatmaps/new_heatmaps/score_bins

# samples: sample_name -> directory
declare -A SAMPLES=(
  [imbulrich20240401wt]=imbulrich20240401wt
  [sfbulrich20250101wt]=sfbulrich20250101wt
  [sfbulrich20250102wtrad21]=sfbulrich20250102wtrad21
)

# window sizes
WINDOWS=(50000 75000 100000 150000)

# score bins (order matters: top → bottom)
BINS=(80-100 60-80 40-60 20-40 0-20)

# region types
REGION_TYPES=(oris ters)

for SAMPLE in "${!SAMPLES[@]}"; do
  SAMPLE_DIR=${SAMPLES[$SAMPLE]}

  for REGION in "${REGION_TYPES[@]}"; do

    # build regionsFileName argument (all score bins for this region type)
    REGIONS=""
    for BIN in "${BINS[@]}"; do
      BED=${BIN_DIR}/${SAMPLE}_${REGION}_score_${BIN}.bed
      [[ -s ${BED} ]] && REGIONS="${REGIONS} ${BED}"
    done

    # skip if no bins exist (e.g. empty TERs)
    [[ -z "${REGIONS// }" ]] && continue

    for WIN in "${WINDOWS[@]}"; do
      RFD=${BASE}/${SAMPLE_DIR}/${SAMPLE}_rfd_10_${WIN}.bw
      OEM=${BASE}/${SAMPLE_DIR}/${SAMPLE}_oem_10_${WIN}.bw

      for TRACK in rfd oem; do
        if [[ $TRACK == "rfd" ]]; then
          BW=$RFD
          LABEL="RFD"
        else
          BW=$OEM
          LABEL="OEM"
        fi

        OUT=${SAMPLE}_${REGION}_${TRACK}_${WIN}_scorebins

        sbatch -c 24 --mem=48G --job-name=${OUT} --wrap "
          apptainer exec ${SIF} computeMatrix reference-point \
            --referencePoint center \
            --regionsFileName ${REGIONS} \
            --scoreFileName ${BW} \
            --beforeRegionStartLength 250000 \
            --afterRegionStartLength 250000 \
            --binSize 100 \
            --sortRegions keep \
            --sortUsing mean \
            --numberOfProcessors 24 \
            --averageTypeBins mean \
            --outFileName ${OUT}.mat.gz

          apptainer exec ${SIF} plotHeatmap \
            --heatmapWidth 5 \
            --heatmapHeight 10 \
            --regionsLabel 80-100 60-80 40-60 20-40 0-20 \
            --matrixFile ${OUT}.mat.gz \
            --outFileName ${OUT}.pdf
        "
      done
    done
  done
done
