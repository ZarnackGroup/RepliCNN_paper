#!/bin/bash

OUT_DIR=PATH/scripts/04_plots_analyses/overlap_ori_term
mkdir -p "${OUT_DIR}"

ORI_FILES=(
  PATH/data/oem_rfd/imbulrich20240401wt/imbulrich20240401wt_oris.bed
  PATH/data/oem_rfd/sfbulrich20250101wt/sfbulrich20250101wt_oris.bed
  PATH/data/oem_rfd/sfbulrich20250102wtrad21/sfbulrich20250102wtrad21_oris.bed
)

apptainer exec ~/apptainer/bedtools_2.31.1_hf5e1c6e_2.sif bedtools sort -i PATH/scripts/06_revision/okseqhmm/data/remerge/calls/all_oris_okseqhmm_replicnn.bed | \
apptainer exec ~/apptainer/bedtools_2.31.1_hf5e1c6e_2.sif bedtools merge -i - -d 50000 -c 4 -o collapse > PATH/scripts/06_revision/okseqhmm/data/remerge/calls/all_oris_okseqhmm_replicnn_merged.bed

apptainer exec ~/apptainer/bedtools_2.31.1_hf5e1c6e_2.sif bedtools sort -i PATH/scripts/06_revision/okseqhmm/data/remerge/calls/all_oris_okseqhmm.bed | \
apptainer exec ~/apptainer/bedtools_2.31.1_hf5e1c6e_2.sif bedtools merge -i - -d 50000 -c 4 -o collapse > PATH/scripts/06_revision/okseqhmm/data/remerge/calls/all_oris_okseqhmm_merged.bed

apptainer exec ~/apptainer/bedtools_2.31.1_hf5e1c6e_2.sif bedtools sort -i PATH/scripts/06_revision/okseqhmm/data/remerge/calls/all_oris_replicnn.bed | \
apptainer exec ~/apptainer/bedtools_2.31.1_hf5e1c6e_2.sif bedtools merge -i - -d 50000 -c 4 -o collapse > PATH/scripts/06_revision/okseqhmm/data/remerge/calls/all_oris_replicnn_merged.bed


sbatch -c 64 --mem=86G  --wrap " \
	apptainer exec ~/apptainer/deeptools_3.5.6_pyhdfd78af_0.sif computeMatrix reference-point \
	--referencePoint center \
	--regionsFileName \
		PATH/scripts/06_revision/okseqhmm/data/remerge/calls/oris_both.bed \
		PATH/scripts/06_revision/okseqhmm/data/remerge/calls/oris_only_replicnn.bed \
		PATH/scripts/06_revision/okseqhmm/data/remerge/calls/oris_only_okseqhmm.bed \
	--scoreFileName \
		PATH/scripts/04_plots_analyses/overlap_ori_term/rerun_rfd_oem_ori_ter_hct116/imbulrich20240401wt/imbulrich20240401wt_oem_10_75000.bw \
		PATH/scripts/04_plots_analyses/overlap_ori_term/rerun_rfd_oem_ori_ter_hct116/sfbulrich20250101wt/sfbulrich20250101wt_oem_10_75000.bw \
		PATH/scripts/04_plots_analyses/overlap_ori_term/rerun_rfd_oem_ori_ter_hct116/sfbulrich20250102wtrad21/sfbulrich20250102wtrad21_oem_10_75000.bw \
	--beforeRegionStartLength 250000 \
	--afterRegionStartLength 250000 \
	--binSize 100 \
	--sortRegions keep \
	--sortUsing mean \
	--numberOfProcessors 64 \
	--averageTypeBins mean \
	--outFileName ./heatmap.mat.gz

	apptainer exec ~/apptainer/deeptools_3.5.6_pyhdfd78af_0.sif plotHeatmap \
	--heatmapWidth 5 \
	--heatmapHeight 10 \
	--regionsLabel Both RepliCNN OkseqHMM \
	--matrixFile ./heatmap.mat.gz \
	--outFileName ./heatmap.pdf \
"
