import pandas as pd
import glob
from pathlib import Path
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

oris_okseqhmm = {
	folder.stem: pd.concat([
		pd.read_csv(file, sep="\t", names=["chromosome", "start", "end"]) for file in folder.glob("*_IZ.bed")
	]).sort_values(by=["chromosome", "start", "end"]).reset_index(drop=True)
	for folder in Path("/home/dos02bi/koenig_data/projects/rt_prediction/paper_code_for_github/scripts/06_revision/okseqhmm/data/remerge/calls").glob("*") if folder.is_dir() and (folder.stem.startswith("imb") or folder.stem.startswith("sfb"))
}

oris_replicnn = {
    "imb_ulrich_2024_04_01_traelseq_hct116_asy": pd.read_csv("/home/dos02bi/koenig_data/projects/rt_prediction/paper_code_for_github/scripts/04_plots_analyses/overlap_ori_term/rerun_rfd_oem_ori_ter_hct116/ori_ter/imbulrich20240401wt_oris.bed", sep="\t", names=["chromosome", "start", "end", "name", "score", "strand", "thickStart", "thickEnd", "itemRgb"])[["chromosome", "start", "end"]],
    "sfb_ulrich_2025_01_01_hct116_wt_S1": pd.read_csv("/home/dos02bi/koenig_data/projects/rt_prediction/paper_code_for_github/scripts/04_plots_analyses/overlap_ori_term/rerun_rfd_oem_ori_ter_hct116/ori_ter/sfbulrich20250101wt_oris.bed", sep="\t", names=["chromosome", "start", "end", "name", "score", "strand", "thickStart", "thickEnd", "itemRgb"])[["chromosome", "start", "end"]],
    "sfb_ulrich_2025_01_02_hct116_rad21_wt_S2": pd.read_csv("/home/dos02bi/koenig_data/projects/rt_prediction/paper_code_for_github/scripts/04_plots_analyses/overlap_ori_term/rerun_rfd_oem_ori_ter_hct116/ori_ter/sfbulrich20250102wtrad21_oris.bed", sep="\t", names=["chromosome", "start", "end", "name", "score", "strand", "thickStart", "thickEnd", "itemRgb"])[["chromosome", "start", "end"]],
}

for key in oris_okseqhmm.keys():

    oris_okseqhmm[key] = oris_okseqhmm[key].copy()
    oris_okseqhmm[key]["name"] = [
        f"okseqhmm_{key}_{i}" 
        for i in oris_okseqhmm[key].index
    ]

    oris_replicnn[key] = oris_replicnn[key].copy()
    oris_replicnn[key]["name"] = [
        f"replicnn_{key}_{i}" 
        for i in oris_replicnn[key].index
    ]

pd.concat(oris_replicnn.values()).reset_index(drop=True).to_csv("/home/dos02bi/koenig_data/projects/rt_prediction/paper_code_for_github/scripts/06_revision/okseqhmm/data/remerge/calls/all_oris_replicnn.bed", sep="\t", index=False, header=False)
pd.concat(oris_okseqhmm.values()).reset_index(drop=True).to_csv("/home/dos02bi/koenig_data/projects/rt_prediction/paper_code_for_github/scripts/06_revision/okseqhmm/data/remerge/calls/all_oris_okseqhmm.bed", sep="\t", index=False, header=False)
pd.concat([pd.concat(oris_okseqhmm.values()), pd.concat(oris_replicnn.values())], ignore_index=True).reset_index(drop=True).to_csv("/home/dos02bi/koenig_data/projects/rt_prediction/paper_code_for_github/scripts/06_revision/okseqhmm/data/remerge/calls/all_oris_okseqhmm_replicnn.bed", sep="\t", index=False, header=False)

for key in oris_okseqhmm.keys():
	oris_okseqhmm[key] = oris_okseqhmm[key].assign(width=oris_okseqhmm[key]["end"] - oris_okseqhmm[key]["start"])
	oris_replicnn[key] = oris_replicnn[key].assign(width=oris_replicnn[key]["end"] - oris_replicnn[key]["start"])

sns.kdeplot(oris_okseqhmm["imb_ulrich_2024_04_01_traelseq_hct116_asy"]["width"], label="OKseqHMM - TrAEL-Seq WT 1", color="C0", linestyle="solid")
sns.kdeplot(oris_replicnn["imb_ulrich_2024_04_01_traelseq_hct116_asy"]["width"], label="RepliCNN - TrAEL-Seq WT 1", color="C1", linestyle="solid")
sns.kdeplot(oris_okseqhmm["sfb_ulrich_2025_01_01_hct116_wt_S1"]["width"], label="OKseqHMM - TrAEL-Seq WT 2", color="C0", linestyle="dashed")
sns.kdeplot(oris_replicnn["sfb_ulrich_2025_01_01_hct116_wt_S1"]["width"], label="RepliCNN - TrAEL-Seq WT 2", color="C1", linestyle="dashed")
sns.kdeplot(oris_okseqhmm["sfb_ulrich_2025_01_02_hct116_rad21_wt_S2"]["width"], label="OKseqHMM - TrAEL-Seq WT 3", color="C0", linestyle="dotted")
sns.kdeplot(oris_replicnn["sfb_ulrich_2025_01_02_hct116_rad21_wt_S2"]["width"], label="RepliCNN - TrAEL-Seq WT 3", color="C1", linestyle="dotted")
plt.axvline(x=pd.concat(oris_okseqhmm.values()).width.median(), color="C0", linestyle="-", linewidth=1)
plt.text(x=pd.concat(oris_okseqhmm.values()).width.median() + 1_000, y=0.000065, s=f"OKseqHMM median: {pd.concat(oris_okseqhmm.values()).width.median():,.0f} bp", color="black", fontsize=8, rotation=0, va="bottom")
plt.text(x=pd.concat(oris_replicnn.values()).width.median() + 1_000, y=0.00007, s=f"RepliCNN median: {pd.concat(oris_replicnn.values()).width.median():,.0f} bp", color="black", fontsize=8, rotation=0, va="bottom")
plt.axvline(x=pd.concat(oris_replicnn.values()).width.median(), color="C1", linestyle="-", linewidth=1)
plt.legend(title="Method", loc="lower right")
plt.xlabel("ORI width (bp)")
plt.ylabel("Density")
plt.xlim(-10_000, 100_000)
plt.title("ORI width distribution comparison between OKseqHMM and RepliCNN")
plt.savefig("./ori_width_distribution_comparison.pdf", dpi=300)
plt.show()

dfs = []

for method_name, method in {
    "OKseqHMM": oris_okseqhmm,
    "RepliCNN": oris_replicnn
}.items():

    for sample, df in method.items():
        tmp = df.copy()
        tmp["sample"] = sample
        tmp["method"] = method_name
        dfs.append(tmp)

plot_df = pd.concat(dfs, ignore_index=True)

counts = (
    plot_df
    .groupby(["chromosome", "method", "sample"])
    .size()
    .reset_index(name="n_oris")
)

chrom_order = [f"chr{i}" for i in range(1, 23)]

fig, ax = plt.subplots(figsize=(10,4))

sns.barplot(
    data=counts,
    x="chromosome",
    y="n_oris",
    hue="method",
    order=chrom_order,
    hue_order=["OKseqHMM", "RepliCNN"],
    errorbar=None,
    ax=ax
)

sns.stripplot(
    data=counts,
    x="chromosome",
    y="n_oris",
    hue="method",
    order=chrom_order,
    hue_order=["OKseqHMM", "RepliCNN"],
    dodge=True,
    jitter=0.15,
    size=7,
    linewidth=0.5,
    edgecolor="black",
    ax=ax
)

# remove duplicate legends
handles, labels = ax.get_legend_handles_labels()
ax.legend(handles[:2], labels[:2], loc="upper right")
plt.xticks(rotation=90)
ax.set_ylabel("Number of origins")
plt.tight_layout()
plt.savefig("./ori_counts_per_chromosome.pdf", dpi=300)
plt.show()

merged_oris_all = pd.read_csv("/home/dos02bi/koenig_data/projects/rt_prediction/paper_code_for_github/scripts/06_revision/okseqhmm/data/remerge/calls/all_oris_okseqhmm_replicnn_merged.bed", sep="\t", names=["chromosome", "start", "end", "name"])
merged_oris_okseqhmm = pd.read_csv("/home/dos02bi/koenig_data/projects/rt_prediction/paper_code_for_github/scripts/06_revision/okseqhmm/data/remerge/calls/all_oris_okseqhmm_merged.bed", sep="\t", names=["chromosome", "start", "end", "name"])
merged_oris_replicnn = pd.read_csv("/home/dos02bi/koenig_data/projects/rt_prediction/paper_code_for_github/scripts/06_revision/okseqhmm/data/remerge/calls/all_oris_replicnn_merged.bed", sep="\t", names=["chromosome", "start", "end", "name"])

samples = [
	"replicnn_imb_ulrich_2024_04_01_traelseq_hct116_asy",
	"replicnn_sfb_ulrich_2025_01_01_hct116_wt_S1",
	"replicnn_sfb_ulrich_2025_01_02_hct116_rad21_wt_S2",
	"okseqhmm_imb_ulrich_2024_04_01_traelseq_hct116_asy",
	"okseqhmm_sfb_ulrich_2025_01_01_hct116_wt_S1",
	"okseqhmm_sfb_ulrich_2025_01_02_hct116_rad21_wt_S2",
]

for sample in samples:
	merged_oris_all[sample] = merged_oris_all["name"].str.contains(sample)

	if "okseqhmm" in sample:
		merged_oris_okseqhmm[sample] = merged_oris_okseqhmm["name"].str.contains(sample)
	if "replicnn" in sample:
		merged_oris_replicnn[sample] = merged_oris_replicnn["name"].str.contains(sample)
	
merged_oris_all["okseqhmm"] = merged_oris_all["name"].str.contains("okseqhmm")
merged_oris_all["replicnn"] = merged_oris_all["name"].str.contains("replicnn")

import pandas as pd
import matplotlib.pyplot as plt

rep_cols = [
    "replicnn_imb_ulrich_2024_04_01_traelseq_hct116_asy",
    "replicnn_sfb_ulrich_2025_01_01_hct116_wt_S1",
    "replicnn_sfb_ulrich_2025_01_02_hct116_rad21_wt_S2"
]

ok_cols = [
    "okseqhmm_imb_ulrich_2024_04_01_traelseq_hct116_asy",
    "okseqhmm_sfb_ulrich_2025_01_01_hct116_wt_S1",
    "okseqhmm_sfb_ulrich_2025_01_02_hct116_rad21_wt_S2"
]


# count support within each method
merged_oris_all["rep_support"] = (
    merged_oris_all[rep_cols].sum(axis=1)
)

merged_oris_all["ok_support"] = (
    merged_oris_all[ok_cols].sum(axis=1)
)


def count_support(query, column):
    return (
        merged_oris_all.query(query)[column]
        .value_counts()
        .reindex([1,2,3], fill_value=0)
    )


plot_df = pd.DataFrame({

    "RepliCNN only":
        count_support(
            "replicnn==True and okseqhmm==False",
            "rep_support"
        ),

    "Shared - RepliCNN":
        count_support(
            "replicnn==True and okseqhmm==True",
            "rep_support"
        ),

    "Shared - OKseqHMM":
        count_support(
            "replicnn==True and okseqhmm==True",
            "ok_support"
        ),

    "OKseqHMM only":
        count_support(
            "replicnn==False and okseqhmm==True",
            "ok_support"
        )
})


plot_df.index = [
    "1 sample",
    "2 samples",
    "3 samples"
]

ax = plot_df.T.plot(
    kind="bar",
    stacked=True,
    figsize=(9,5)
)

ax.set_ylabel("Number of origins")
ax.set_xlabel("")
ax.set_title(
    "Origin support across biological samples"
)

plt.xticks(rotation=45, ha="right")

ax.legend(
    title="Number of samples",
    bbox_to_anchor=(1.02,1),
    loc="upper left"
)
plt.savefig("./ori_support_across_samples.pdf", dpi=300, bbox_inches="tight")
plt.tight_layout()
plt.show()
