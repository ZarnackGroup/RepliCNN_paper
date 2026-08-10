#!/bin/bash

for sample in "imb_ulrich_2024_04_01_traelseq_hct116_asy.umi.umidedup" "sfb_ulrich_2025_01_01_hct116_wt_S1.umi.umidedup" "sfb_ulrich_2025_01_02_hct116_rad21_wt_S2.umi.umidedup"; do
	apptainer exec ~/apptainer/samtools_1.23.1_ha83d96e_0.sif samtools merge \
		PATH/scripts/06_revision/okseqhmm/data/remerge/${sample}.bam \
		PATH/scripts/06_revision/okseqhmm/data/${sample}.fwd.bam \
		PATH/scripts/06_revision/okseqhmm/data/${sample}.rev.bam
	apptainer exec ~/apptainer/samtools_1.23.1_ha83d96e_0.sif samtools sort -o \
		PATH/scripts/06_revision/okseqhmm/data/remerge/${sample}_sorted.bam \
		PATH/scripts/06_revision/okseqhmm/data/remerge/${sample}.bam
	apptainer exec ~/apptainer/samtools_1.23.1_ha83d96e_0.sif samtools index PATH/scripts/06_revision/okseqhmm/data/remerge/${sample}_sorted.bam
	rm PATH/scripts/06_revision/okseqhmm/data/remerge/${sample}.bam
done

for bam in PATH/scripts/06_revision/okseqhmm/data/remerge/*_sorted.bam; do
	apptainer exec ~/apptainer/samtools_1.23.1_ha83d96e_0.sif samtools view -b -f 128 -F 16 ${bam} > tmp.fwd1.bam
	apptainer exec ~/apptainer/samtools_1.23.1_ha83d96e_0.sif samtools view -b -f 80 ${bam} > tmp.fwd2.bam
	apptainer exec ~/apptainer/samtools_1.23.1_ha83d96e_0.sif samtools merge ${bam%.bam}.fwd.bam tmp.fwd1.bam tmp.fwd2.bam
	apptainer exec ~/apptainer/samtools_1.23.1_ha83d96e_0.sif samtools index ${bam%.bam}.fwd.bam
	rm tmp.fwd1.bam tmp.fwd2.bam

	apptainer exec ~/apptainer/samtools_1.23.1_ha83d96e_0.sif samtools view -b -f 144 ${bam} > tmp.rev1.bam
	apptainer exec ~/apptainer/samtools_1.23.1_ha83d96e_0.sif samtools view -b -f 64 -F 16 ${bam} > tmp.rev2.bam
	apptainer exec ~/apptainer/samtools_1.23.1_ha83d96e_0.sif samtools merge ${bam%.bam}.rev.bam tmp.rev1.bam tmp.rev2.bam
	apptainer exec ~/apptainer/samtools_1.23.1_ha83d96e_0.sif samtools index ${bam%.bam}.rev.bam
	rm tmp.rev1.bam tmp.rev2.bam
done

for bam in PATH/scripts/06_revision/okseqhmm/data/remerge/*.fwd.bam PATH/scripts/06_revision/okseqhmm/data/remerge/*.rev.bam; do
	echo "Processing ${bam}..."
	for chromosome in chr{1..22}; do
		echo "Extracting reads for ${chromosome}..."
		apptainer exec ~/apptainer/samtools_1.23.1_ha83d96e_0.sif \
			samtools view -q 1 -f 3 -F 4 ${bam} ${chromosome} |
			awk '{print $4}' > ${bam%.bam}_${chromosome}.txt
	done
done
