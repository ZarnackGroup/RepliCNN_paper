import pandas as pd
from pathlib import Path
import bbi
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import spearmanr, pearsonr

names = ["chrom", "start", "end", "name", "score", "strand",
         "thickStart", "thickEnd", "itemRgb"]

# input files
datasets = {
    "imbulrich20240401wt_oris": "PATH/scripts/04_plots_analyses/overlap_ori_term/rerun_rfd_oem_ori_ter_hct116/ori_ter/imbulrich20240401wt_oris.bed",
    "sfbulrich20250101wt_oris": "PATH/scripts/04_plots_analyses/overlap_ori_term/rerun_rfd_oem_ori_ter_hct116/ori_ter/sfbulrich20250101wt_oris.bed",
    "sfbulrich20250102wtrad21_oris": "PATH/scripts/04_plots_analyses/overlap_ori_term/rerun_rfd_oem_ori_ter_hct116/ori_ter/sfbulrich20250102wtrad21_oris.bed",

    "imbulrich20240401wt_ters": "PATH/scripts/04_plots_analyses/overlap_ori_term/rerun_rfd_oem_ori_ter_hct116/ori_ter/imbulrich20240401wt_terms.bed",
    "sfbulrich20250101wt_ters": "PATH/scripts/04_plots_analyses/overlap_ori_term/rerun_rfd_oem_ori_ter_hct116/ori_ter/sfbulrich20250101wt_terms.bed",
    "sfbulrich20250102wtrad21_ters": "PATH/scripts/04_plots_analyses/overlap_ori_term/rerun_rfd_oem_ori_ter_hct116/ori_ter/sfbulrich20250102wtrad21_terms.bed",
}

timing = pd.read_csv("PATH/data/timetables/human/human_hct116_GSE137764_hg38_rt_processed.bg",sep="\t",names=["chromosome", "start", "end", "score"])

df = pd.read_csv('PATH/scripts/04_plots_analyses/overlap_ori_term/rerun_rfd_oem_ori_ter_hct116/ori_ter/imbulrich20240401wt_oris.bed', sep="\t", names=names)

df = df.assign(timing=df.apply(lambda row: timing[(timing["chromosome"] == row["chrom"]) & (timing["start"] <= row["end"]) & (timing["end"] >= row["start"])].score.max(), axis=1))
df = df.assign(sirt1=df.apply(lambda row: bbi.fetch("PATH/scripts/06_revision/mcm_chip/data/nextflow_out/04_reporting/igv/pSIRT1_R1.bigWig", row["chrom"], row["start"], row["end"]).sum(), axis=1))
df = df.assign(mcm2=df.apply(lambda row: bbi.fetch("PATH/scripts/06_revision/mcm_chip/data/nextflow_out/04_reporting/igv/pMCM2_R1.bigWig", row["chrom"], row["start"], row["end"]).sum(), axis=1))
df = df.assign(h4k16ac=df.apply(lambda row: bbi.fetch("PATH/scripts/06_revision/mcm_chip/data/nextflow_out/04_reporting/igv/H4K16ac_R1.bigWig", row["chrom"], row["start"], row["end"]).sum(), axis=1))

def add_mean_line(ax, x, y, order):
    means = y.groupby(x).median().reindex(order)

    ax.plot(
        range(len(order)),
        means.values,
        color="red",
        marker="o",
        linewidth=2,
        zorder=10
    )
    
def add_n_labels(ax, x, order, y=700):
    counts = x.value_counts().reindex(order, fill_value=0)

    for i, n in enumerate(counts):
        ax.text(
            i,
            y,
            f"n={n}",
            ha="center",
            va="bottom",
            fontsize=8
        )

def make_zero_quintiles(values):
    order = ["None", "Very low", "Low", "Intermediate", "High", "Very high"]
    cats = pd.Series(index=values.index, dtype="object")
    cats.loc[values == 0] = "None"
    nonzero = values != 0
    cats.loc[nonzero] = pd.qcut(
        values.loc[nonzero],
        q=5,
        labels=order[1:],
        duplicates="drop"
    ).astype(str)
    return pd.Categorical(cats, categories=order, ordered=True), order

timing_order = ["Very late", "Late", "Intermediate", "Early", "Very early"]
mcm_order = ["Very low", "Low", "Intermediate", "High", "Very high"]

fig, axes = plt.subplots(1, 4, figsize=(16, 5), sharey=True)

x1 = pd.cut(df.timing, bins=[-1,-0.5,-0.25,0.25,0.5,1], labels=timing_order, duplicates="drop")
sns.boxplot(
	x=x1,
	y=df.score,
	showfliers=False,
	ax=axes[0],
	hue=x1,
	palette=sns.color_palette("rocket")
)
add_mean_line(axes[0], x1, df.score, timing_order)
axes[0].set_xticklabels(
    axes[0].get_xticklabels(),
    rotation=-30,
    ha='left'
)
add_n_labels(axes[0], x1, timing_order)
axes[0].set_xlabel("IZ Replication timing")
axes[0].set_ylabel("RepliCNN-derived OEM score")

#plot_data = df.query("mcm2!=0")
#x3 = pd.qcut(plot_data.mcm2, q=5, labels=mcm_order, duplicates="drop")
x3, mcm_order = make_zero_quintiles(df.mcm2)
sns.boxplot(
	x=x3,
	y=df.score,
	showfliers=False,
	ax=axes[1],
	hue=x3,
	palette=sns.cubehelix_palette()
)
add_mean_line(axes[1], x3, df.score, mcm_order)
axes[1].set_xticklabels(
    axes[1].get_xticklabels(),
    rotation=-30,
    ha='left'
)
add_n_labels(axes[1], x3, mcm_order)
axes[1].set_ylabel("RepliCNN-derived OEM score")
axes[1].set_xlabel("IZ log(mcm2 ChIP-Seq signal) (Quintiles)")

#plot_data = df.query("h4k16ac!=0")
#x3 = pd.qcut(plot_data.h4k16ac, q=5, labels=mcm_order, duplicates="drop")
x3, mcm_order = make_zero_quintiles(df.h4k16ac)
sns.boxplot(
	x=x3,
	y=df.score,
	showfliers=False,
	ax=axes[2],
	hue=x3,
	palette=sns.cubehelix_palette()
)
add_mean_line(axes[2], x3, df.score, mcm_order)
axes[2].set_xticklabels(
    axes[2].get_xticklabels(),
    rotation=-30,
    ha='left'
)
add_n_labels(axes[2], x3, mcm_order)
axes[2].set_ylabel("RepliCNN-derived OEM score")
axes[2].set_xlabel("IZ log(h4k16ac ChIP-Seq signal) (Quintiles)")

#plot_data = df.query("sirt1!=0")
#x3 = pd.qcut(plot_data.sirt1, q=5, labels=mcm_order, duplicates="drop")
x3, mcm_order = make_zero_quintiles(df.sirt1)
sns.boxplot(
	x=x3,
	y=df.score,
	showfliers=False,
	ax=axes[3],
	hue=x3,
	palette=sns.cubehelix_palette()
)
add_mean_line(axes[3], x3, df.score, mcm_order)
axes[3].set_xticklabels(
    axes[3].get_xticklabels(),
    rotation=-30,
    ha='left'
)
add_n_labels(axes[3], x3, mcm_order)
axes[3].set_ylabel("RepliCNN-derived OEM score")
axes[3].set_xlabel("IZ log(sirt1 ChIP-Seq signal) (Quintiles)")
plt.tight_layout()
plt.savefig("score_timing_mcm2_boxplots.pdf", dpi=300)
plt.show()