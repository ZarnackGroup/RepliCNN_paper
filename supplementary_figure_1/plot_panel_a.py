import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

names = ["chromosome", "start", "end", "timing"]

hct116_timing = pd.read_csv("/PATH/data/timetables/human/human_hct116_GSE137764_hg38_rt_processed.bg", sep="\t", names=names)

human_chromosomes = ["chr"+str(i) for i in range(1, 23)]
human_samples = ["imbulrich20240401wt", "sfbulrich20250101wt", "sfbulrich20250102wtrad21"]
base_path = "/PATH/data/cross_prediction/models"

names = ["chromosome", "start", "end", "fwd", "rev", "log2", "spline", "deri", "antideri", "predicted_time"]

results_human = {
	i: {
		j: pd.concat([pd.read_csv(f"{base_path}/human_traelseq_{i}_delta_{chrom}/{j}_pred.tsv", names=names, sep="\t").query('chromosome == @chrom') for chrom in human_chromosomes])
		for j in human_samples
	}
	for i in human_samples
}

fig, ax = plt.subplots(1, 1, figsize=(5, 5))

masked = set()
for _, row in hct116_timing.query("timing==0").iterrows():
	positions = range(row.start, row.end, 10000)  # df2 bin size
	masked.update((row.chromosome, pos) for pos in positions)

sns.kdeplot(ax=ax, data=hct116_timing.query("timing!=0"), x="timing", label="GSE137764", alpha=1, color="black", linewidth=3)
for human_sample in human_samples:
	mask = [
		(chrom, start) in masked
		for chrom, start in zip(results_human[human_sample][human_sample].chromosome, results_human[human_sample][human_sample].start)
	]

	results_human[human_sample][human_sample].loc[mask, "predicted_time"] = 0

	sns.kdeplot(ax=ax, data=results_human[human_sample][human_sample].query("predicted_time!=0"), x="predicted_time", label=human_sample, alpha=0.5)

ax.set_title("HCT116 Replication Timing Distribution")
ax.legend(title="Samples", loc="upper right")

plt.tight_layout()
plt.savefig("./rt_distribution.pdf", dpi=300)
plt.show()
