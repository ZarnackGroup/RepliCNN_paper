import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

def plot_rfd_oem_bedgraph(
	chrom, 
	start, 
	end,
	df1, df2, df3, truth,
	colors=("darkred", "darkgreen", "darkblue"),
	truth_color="black",
	alpha=0.3,
	outfile=None,
	labels=None,
	max_points=100_000  # downsample if too many points
):
	x_len = end - start
	x = np.arange(start, end)
	
	dfs = [df1, df2, df3]
	n = len(dfs)
	
	if labels is None:
		labels = [f"DF{i+1}" for i in range(n)]

	# Pre-allocate truth values
	truth_values = np.zeros(x_len, dtype=float)
	sub_t = truth[(truth["chrom"] == chrom) & (truth["end"] > start) & (truth["start"] < end)]
	
	if not sub_t.empty:
		s_idx = np.maximum(sub_t["start"].values, start) - start
		e_idx = np.minimum(sub_t["end"].values, end) - start
		val = sub_t["value"].values
		for si, ei, v in zip(s_idx, e_idx, val):
			truth_values[si:ei] = v

	# Create subplots
	fig, axes = plt.subplots(nrows=n, ncols=1, sharex=True, figsize=(12, 6))
	
	for df, ax, color, label in zip(dfs, axes, colors, labels):
		values = np.zeros(x_len, dtype=float)
		sub = df[(df["chrom"] == chrom) & (df["end"] > start) & (df["start"] < end)]
		
		if not sub.empty:
			s_idx = np.maximum(sub["start"].values, start) - start
			e_idx = np.minimum(sub["end"].values, end) - start
			val = sub["value"].values
			for si, ei, v in zip(s_idx, e_idx, val):
				values[si:ei] = v

		# Downsample for plotting if needed
		if x_len > max_points:
			factor = int(np.ceil(x_len / max_points))
			values_ds = values[::factor]
			truth_ds = truth_values[::factor]
			x_ds = x[::factor]
		else:
			values_ds = values
			truth_ds = truth_values
			x_ds = x

		# Overlay ground truth
		ax.plot(x_ds, truth_ds, color=truth_color, linewidth=1, label="Ground truth")
		# Plot predicted track
		ax.plot(x_ds, values_ds, color=color, linewidth=1, label=label)
		
		ax.set_ylabel(label)
		ax.set_ylim(-1, 1)
		ax.set_xlim(start, end)
		ax.set_ylabel("")
		ax.set_xlabel("")

	fig.suptitle(f"{chrom}:{start:,}-{end:,}", fontsize=5)
	fig.supylabel("Replication timing", fontsize=5)
	fig.supxlabel(f"Genomic position on {chrom} (mb)", fontsize=5)
	plt.tight_layout(rect=[0, 0, 1, 0.95])

	if outfile:
		plt.savefig(outfile, dpi=300)
	plt.show()
	plt.close()
  
names = ["chrom","start","end","pos","neg","log2","spline","deri","antideri","value"]
true_time = pd.read_csv("PATH/data/sdfs/imbulrich20240401wt.tsv", sep="\t", names=names)
autosomes = [f"chr{i}" for i in range(1, 23)]
chromosomes = pd.read_csv("/storage/zar/shared/organisms/homo_sapiens/genomes/hg38/hg38.chrom.sizes", sep="\t", names=["chrom", "size"]).query("chrom in @autosomes")

for chromosome in chromosomes.iterrows():
	df1 = pd.read_csv(f"PATH/data/cross_prediction/models/human_traelseq_imbulrich20240401wt_delta_{chromosome[1]['chrom']}/imbulrich20240401wt_pred.tsv", sep="\t", names=names)
	df2 = pd.read_csv(f"PATH/data/cross_prediction/models/human_traelseq_sfbulrich20250101wt_delta_{chromosome[1]['chrom']}/sfbulrich20250101wt_pred.tsv", sep="\t", names=names)
	df3 = pd.read_csv(f"PATH/data/cross_prediction/models/human_traelseq_sfbulrich20250102wtrad21_delta_{chromosome[1]['chrom']}/sfbulrich20250102wtrad21_pred.tsv", sep="\t", names=names)


	plot_rfd_oem_bedgraph(
		chrom=chromosome[1]["chrom"],
		start=0,
		end=chromosome[1]["size"],
		df1=df1,
		df2=df2,
		df3=df3,
		truth=true_time,
		colors=("darkred", "darkred", "darkred"),
		truth_color="goldenrod",
		alpha=1,
		outfile=f"PATH/scripts/04_plots_analyses/examples/human_examples_{chromosome[1]['chrom']}_full.pdf",
)
