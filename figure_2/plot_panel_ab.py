import matplotlib as mpl
import pandas as pd
import os
from pathlib import Path
import sys
from sklearn.metrics import r2_score
from scipy.stats import pearsonr
import random

folder = Path("PATH/data/cross_prediction/models")
yeast_paths = [path for path in list(folder.iterdir()) if "yeast" in str(path)]
samples = set(["_".join(os.path.basename(sample).split("_")[1:3]) for sample in yeast_paths])
chromosomes = ["chrI", "chrII", "chrIII", "chrIV", "chrV", "chrVI", "chrVII", "chrVIII", "chrIX", "chrX", "chrXI", "chrXII", "chrXIII", "chrXIV", "chrXV", "chrXVI"]

benchmark = pd.read_csv("PATH/data/sdfs/GSM4680452.tsv",sep="\t",names=["chromosome", "start", "end", "pos", "neg", "log2", "spline", "derivative", "antiderivative", "true_time"])[["chromosome", "start", "end", "true_time"]]

# get self predictions (LOCO-CV)
data = {}
for sample in list(samples):
	if "GSM4680455" in sample: continue
	sdf = pd.DataFrame()
	for chromosome in chromosomes:
		file = f"{str(folder)}/yeast_{sample}_delta_{chromosome}/{sample.split("_")[1]}_pred.tsv"
		df = pd.read_csv(file, sep="\t", names=["chromosome", "start", "end", "pos", "neg", "log2", "spline", "derivative", "antiderivative", "pred_time"])[["chromosome", "start", "end", "pred_time"]].query("chromosome==@chromosome")
		sdf = pd.concat([sdf,df])
	sdf = sdf.assign(true_time=benchmark.true_time)
	data[sample] = sdf

scores_total = []
scores_chromosomes = []
for sample in data.keys():
	scores_total.append(tuple([
			sample,
			"Overall",
			r2_score(data[sample].true_time, data[sample].pred_time),
			pearsonr(data[sample].true_time, data[sample].pred_time)[0],
			]))
	for chromosome in chromosomes:
		scores_chromosomes.append(tuple([
			sample, 
			chromosome, 
			r2_score(data[sample].query("chromosome==@chromosome").true_time, data[sample].query("chromosome==@chromosome").pred_time),
			pearsonr(data[sample].query("chromosome==@chromosome").true_time, data[sample].query("chromosome==@chromosome").pred_time)[0],
			]))
	
scores_total = pd.DataFrame(scores_total, columns=["sample", "chromosome", "r2", "r"])
scores_chromosomes = pd.DataFrame(scores_chromosomes, columns=["sample", "chromosome", "r2", "r"])

random.seed(1)

plot_data = scores_chromosomes.copy()

family_order = {"okseq": 0, "traelseq": 1, "gloeseq": 2, "Overall": 3}

def get_family(name):
    for fam in family_order:
        if name.startswith(fam):
            return fam
    if name == "Overall":
        return "Overall"
    return "Other"

plot_data["family"] = plot_data["sample"].apply(get_family)

def trim_family(name):
    fam = get_family(name)
    if fam in ["okseq", "traelseq", "gloeseq"]:
        return name.replace(fam + "_", "", 1)
    return name

plot_data["trimmed"] = plot_data["sample"].apply(trim_family)

# compute one r per sample (use mean if multiple rows per sample)
r_per_sample = (
    plot_data
    .groupby(["family", "sample"], as_index=False)["r"]
    .mean()
)

ordered_samples = []
for fam in family_order:
    fam_samples = (
        r_per_sample[r_per_sample["family"] == fam]
        .sort_values("r", ascending=False)["sample"]
        .tolist()
    )
    ordered_samples.extend(fam_samples)

plot_data["sample"] = pd.Categorical(
    plot_data["sample"],
    categories=ordered_samples,
    ordered=True
)


mpl.rcParams.update({
    "font.family": "Arial",
    "font.size": 5,
    "axes.labelsize": 5,
    "axes.titlesize": 5,
    "xtick.labelsize": 5,
    "ytick.labelsize": 5,
    "legend.fontsize": 5,
    "pdf.fonttype": 42,   # embed TrueType fonts (critical for Affinity)
    "ps.fonttype": 42
})

plt.figure(figsize=(4.5, 1.5))

ax = sns.boxplot(
    data=plot_data,
    x="sample",
    y="r",
    hue="family",
    showfliers=False,
)

plt.tick_params(
    axis="both",
    which="major",
    labelsize=5,
    labelbottom=False,
    bottom=False,
    top=True,
    labeltop=True
)

ax.set_xticklabels(
    plot_data.set_index("sample").loc[ordered_samples, "trimmed"].unique(),
    rotation=50,
    ha="left",
    fontsize=5
)

ax.set_ylabel("")
ax.set_xlabel("")

# optional: legend tweaks for 5 pt readability
ax.legend(
    title="",
    fontsize=5,
    title_fontsize=5,
    frameon=False
)

plt.tight_layout()
plt.savefig(
    "boxplot_affinity_ready.pdf",
    format="pdf",
)
plt.show()
