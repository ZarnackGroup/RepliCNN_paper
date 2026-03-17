conda activate pygenometracks

awk 'OFS="\t" {
    if ($6 == "+") print $1,$2,$3,$4,$5 > "PATH/data/oem_rfd/GSM4680468/GSM4680468_oris.bed";
    else if ($6 == "-") print $1,$2,$3,$4,$5 > "PATH/data/oem_rfd/GSM4680468/GSM4680468_ters.bed";
}' PATH/data/oem_rfd/GSM4680468/GSM4680468_oris_ters.bed
awk 'OFS="\t" {
    if ($6 == "+") print $1,$2,$3,$4,$5 > "PATH/data/oem_rfd/imbulrich201807s16/imbulrich201807s16_oris.bed";
    else if ($6 == "-") print $1,$2,$3,$4,$5 > "PATH/data/oem_rfd/imbulrich201807s16/imbulrich201807s16_ters.bed";
}' PATH/data/oem_rfd/imbulrich201807s16/imbulrich201807s16_oris_ters.bed
awk 'OFS="\t" {
    if ($6 == "+") print $1,$2,$3,$4,$5 > "PATH/data/oem_rfd/GSM835651/GSM835651_oris.bed";
    else if ($6 == "-") print $1,$2,$3,$4,$5 > "PATH/data/oem_rfd/GSM835651/GSM835651_ters.bed";
}' PATH/data/oem_rfd/GSM835651/GSM835651_oris_ters.bed
awk 'OFS="\t" {
    if ($6 == "+") print $1,$2,$3,$4,$5 > "PATH/data/oem_rfd/sfbulrich20240701wtctcf/sfbulrich20240701wtctcf_oris.bed";
    else if ($6 == "-") print $1,$2,$3,$4,$5 > "PATH/data/oem_rfd/sfbulrich20240701wtctcf/sfbulrich20240701wtctcf_ters.bed";
}' PATH/data/oem_rfd/sfbulrich20240701wtctcf/sfbulrich20240701wtctcf_oris_ters.bed


awk -F'\t' '{print $1, $2, $3, $10}' OFS='\t' PATH/data/cross_prediction/models/yeast_traelseq_GSM4680468_delta_chrII/GSM4680468_pred.tsv > PATH/scripts/04_plots_analyses/yeast_summary_plot/yeast_traelseq_GSM4680468_delta_chrII_GSM4680468_pred.bg
awk -F'\t' '{print $1, $2, $3, $10}' OFS='\t' PATH/data/cross_prediction/models/yeast_gloeseq_imbulrich201807s16_delta_chrVI/imbulrich201807s16_pred.tsv > PATH/scripts/04_plots_analyses/yeast_summary_plot/yeast_gloeseq_imbulrich201807s16_delta_chrVI_imbulrich201807s16_pred.bg
awk -F'\t' '{print $1, $2, $3, $10}' OFS='\t' PATH/data/cross_prediction/models/yeast_okseq_GSM835651_delta_chrI/GSM835651_pred.tsv > PATH/scripts/04_plots_analyses/yeast_summary_plot/yeast_okseq_GSM835651_delta_chrI_GSM835651_pred.bg
awk -F'\t' '{print $1, $2, $3, $10}' OFS='\t' PATH/data/cross_prediction/models/mouse_traelseq_sfbulrich20240701wtctcf_delta_chr2/sfbulrich20240701wtctcf_pred.tsv > PATH/scripts/04_plots_analyses/yeast_summary_plot/mouse_traelseq_sfbulrich20240701wtctcf_delta_chr2_sfbulrich20240701wtctcf_pred.bg

# traelseq wt 1
region=chrII:1-813184
sample=GSM4680468
pyGenomeTracks  \
	--tracks ./example_traelseq.ini \
	--region  ${region} \
	--title "TrAEL-Seq ${sample} ${region}" \
	--height 40 \
	--width 28 \
	--outFileName ${sample}.pdf

# gloeseq
region=chrVI:1-270161
sample=imbulrich201807s16
pyGenomeTracks  \
	--tracks ./example_gloeseq.ini \
	--region  ${region} \
	--title "GLOE-Seq ${sample} ${region}" \
	--height 40 \
	--width 28 \
	--outFileName ${sample}.pdf

# okseq
region=chrI:1-230218
sample=GSM835651
pyGenomeTracks  \
	--tracks ./example_okseq.ini \
	--region  ${region} \
	--title "OK-Seq ${sample} ${region}" \
	--height 40 \
	--width 28 \
	--outFileName ${sample}.pdf
