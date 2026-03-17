import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from pathlib import Path
import re
from scipy.stats import pearsonr

autosomes = [f"chr{i}" for i in range(1, 23)]

names = ["chromosome", "start", "end", "pos", "neg", "log2", "spline",  "derivative", "antiderivative", "predicted_rt"]
samples = [
	"imbulrich20240401wt",
	"sfbulrich20250101wt",
	"sfbulrich20250102wtrad21",
]
dfs = {}
for sample in samples:
	df_full = pd.DataFrame()
	for chromosome in autosomes:
		df = pd.read_csv(f"PATH/models/human_traelseq_{sample}_delta_{chromosome}/{sample}_pred.tsv", sep="\t", names=names)[["chromosome", "start", "end", "predicted_rt"]].query("chromosome == @chromosome")
		df_full = pd.concat([df_full, df], ignore_index=True).reset_index(drop=True)
	dfs[sample] = df_full.sort_values(["chromosome", "start"]).reset_index(drop=True)

true_rt = pd.read_csv("PATH/sdfs/imbulrich20240401wt.tsv", sep="\t", names=names)[["chromosome", "start", "end", "predicted_rt"]]

def significance_asterisk(p):
    if p < 0.001:
        return "***"
    elif p < 0.01:
        return "**"
    elif p < 0.05:
        return "*"
    else:
        return ""

def plot_pred_vs_truth_heatmap_4x4(dfs, true_rt, samples, outfile=None):
    all_labels = ["Truth"] + samples
    n = len(all_labels)

    matrix = np.full((n, n), np.nan)
    annot_matrix = np.full((n, n), "", dtype=object)

    truth_values = true_rt["predicted_rt"].values

    for i, label_i in enumerate(all_labels):
        vi = truth_values if label_i == "Truth" else dfs[label_i]["predicted_rt"].values
        for j, label_j in enumerate(all_labels):
            vj = truth_values if label_j == "Truth" else dfs[label_j]["predicted_rt"].values

            if i >= j:  # lower triangle and diagonal
                corr, p = pearsonr(vi, vj)
                matrix[i, j] = corr
                annot_matrix[i, j] = f"{corr:.2f}{significance_asterisk(p)}"

    plt.figure(figsize=(6,6))
    ax = sns.heatmap(
        matrix,
        xticklabels=all_labels,
        yticklabels=all_labels,
        cmap="coolwarm",
        vmin=-1, vmax=1,
        annot=annot_matrix,
        fmt="",
        square=True,
        cbar_kws={"label": "Pearson correlation"}
    )
    plt.title("Predicted RT vs Ground Truth")
    plt.tight_layout()

    if outfile:
        plt.savefig(outfile, dpi=150)
    plt.show()

samples = [
    "imbulrich20240401wt",
    "sfbulrich20250101wt",
    "sfbulrich20250102wtrad21",
]

plot_pred_vs_truth_heatmap_4x4(dfs, true_rt, samples, outfile="pred_vs_truth_heatmap.pdf")
