apptainer exec PATH/tools/apptainer/deeptools_3.5.6_pyhdfd78af_0.sif \
	computeMatrix reference-point \
		--regionsFileName \
			all_called_oris.bed \
			called_not_annotated_oris.bed \
			not_called_confirmed.bed \
			not_called_likely.bed \
			not_called_dubious.bed \
		--scoreFileName \
			PATH/data/oem_rfd/GSM835651/GSM835651_rfd_1_5000.bw \
			PATH/data/oem_rfd/GSM4680468/GSM4680468_rfd_1_5000.bw \
			PATH/data/oem_rfd/imbulrich201807s16/imbulrich201807s16_rfd_1_5000.bw \
			PATH/data/oem_rfd/GSM835651/GSM835651_oem_1_5000.bw \
			PATH/data/oem_rfd/GSM4680468/GSM4680468_oem_1_5000.bw \
			PATH/data/oem_rfd/imbulrich201807s16/imbulrich201807s16_oem_1_5000.bw \
		--outFileName $(basename ${set%.bed}).mat.gz \
		--referencePoint center \
		--beforeRegionStartLength 15000 \
		--afterRegionStartLength 15000 \
		--binSize 10 \
		--sortRegions keep \
		--sortUsing mean \
		--averageTypeBins mean \
		--smartLabels \
		--numberOfProcessors 12 && \
apptainer exec PATH/tools/apptainer/deeptools_3.5.6_pyhdfd78af_0.sif \
	plotHeatmap \
		--matrixFile $(basename ${set%.bed}).mat.gz \
		--heatmapWidth 20 \
		--zMin -0.5 -0.5 -0.5 -0.5 -0.5 -0.5 \
		--zMax 0.5 0.5 0.5 0.35 0.35 0.35 \
		--outFileName $(basename ${set%.bed}).pdf
