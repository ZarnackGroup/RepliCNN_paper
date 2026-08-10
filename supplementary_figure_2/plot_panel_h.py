import bbi
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

merged_oris_all = pd.read_csv("PATH/scripts/06_revision/okseqhmm/data/remerge/calls/all_oris_okseqhmm_replicnn_merged.bed", sep="\t", names=["chromosome", "start", "end", "name"])
timing = pd.read_csv("PATH/data/timetables/human/human_hct116_GSE137764_hg38_rt_processed.bg",sep="\t",names=["chromosome", "start", "end", "score"])

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
	
merged_oris_all["okseqhmm"] = merged_oris_all["name"].str.contains("okseqhmm")
merged_oris_all["replicnn"] = merged_oris_all["name"].str.contains("replicnn")

merged_oris_all = merged_oris_all.assign(
	replicnn_support=merged_oris_all[[sample for sample in samples if "replicnn" in sample]].sum(axis=1),
	okseqhmm_support=merged_oris_all[[sample for sample in samples if "okseqhmm" in sample]].sum(axis=1),
	timing=merged_oris_all.apply(lambda row: timing[(timing["chromosome"] == row["chromosome"]) & (timing["start"] <= row["end"]) & (timing["end"] >= row["start"])].score.max(), axis=1),
	sirt1=merged_oris_all.apply(lambda row: bbi.fetch("PATH/scripts/06_revision/mcm_chip/data/nextflow_out/04_reporting/igv/pSIRT1_R1.bigWig", row["chromosome"], int(row["end"]-row["start"]-25000), int(row["end"]-row["start"]+25000)).sum(), axis=1),
	mcm4=merged_oris_all.apply(lambda row: bbi.fetch("PATH/scripts/06_revision/mcm_chip/data/nextflow_out/04_reporting/igv/pMCM2_R1.bigWig", row["chromosome"], int(row["end"]-row["start"]-25000), int(row["end"]-row["start"]+25000)).sum(), axis=1),
	h4k16ac=merged_oris_all.apply(lambda row: bbi.fetch("PATH/scripts/06_revision/mcm_chip/data/nextflow_out/04_reporting/igv/H4K16ac_R1.bigWig", row["chromosome"], int(row["end"]-row["start"]-25000), int(row["end"]-row["start"]+25000)).sum(), axis=1),
)

plot_df = merged_oris_all.copy()

plot_df = plot_df.assign(
    width=plot_df["end"] - plot_df["start"]
)

plot_df = plot_df.assign(
    norm_sirt1=plot_df["sirt1"] / (plot_df["width"] / 1000),
    norm_mcm4=plot_df["mcm4"] / (plot_df["width"] / 1000),
    norm_h4k16ac=plot_df["h4k16ac"] / (plot_df["width"] / 1000),
)

plot_df["caller"] = np.select(
    [
        (plot_df["replicnn_support"] > 0) &
        (plot_df["okseqhmm_support"] == 0),

        (plot_df["replicnn_support"] == 0) &
        (plot_df["okseqhmm_support"] > 0),

        (plot_df["replicnn_support"] > 0) &
        (plot_df["okseqhmm_support"] > 0),
    ],
    [
        "RepliCNN only",
        "OKseqHMM only",
        "Both",
    ],
    default="None",
)

plot_df = plot_df.query("caller != 'None'").copy()

plot_long = plot_df.melt(
    id_vars=["caller"],
    value_vars=[
        "norm_mcm4",
        "norm_h4k16ac",
        "norm_sirt1",
    ],
    var_name="mark",
    value_name="signal",
)

marks = [
    "norm_mcm4",
    "norm_h4k16ac",
    "norm_sirt1",
]

mark_labels = {
    "norm_mcm4": "MCM4",
    "norm_h4k16ac": "H4K16ac",
    "norm_sirt1": "SIRT1",
}

callers = [
    "RepliCNN only",
    "OKseqHMM only",
    "Both",
]

caller_palette = {
    "RepliCNN only": "#0072B2",
    "OKseqHMM only": "#D55E00",
    "Both": "#666666",
}

fig, axes = plt.subplots(
    1,
    len(marks),
    figsize=(9, 3.5),
    sharey=False,
)

for i, mark in enumerate(marks):

    ax = axes[i]

    subset = plot_long.query(
        "mark == @mark"
    ).copy()

    subset["caller"] = pd.Categorical(
        subset["caller"],
        categories=callers,
        ordered=True,
    )

    sns.boxplot(
        data=subset,
        x="caller",
        y="signal",
        hue="caller",
        palette=caller_palette,
        order=callers,
        showfliers=False,
        width=0.6,
        ax=ax,
        legend=False,
    )

    counts = (
        subset
        .groupby("caller", observed=True)
        .size()
    )

    ymax = subset["signal"].quantile(0.99)

    for k, caller in enumerate(callers):

        if caller in counts.index:

            ax.text(
                k,
                ymax,
                f"n={counts[caller]}",
                ha="center",
                va="bottom",
                fontsize=7,
                rotation=90,
            )

    ax.set_title(mark_labels[mark])

    ax.set_xlabel("")

    ax.set_xticklabels(
        [
            "RepliCNN only",
            "OKseqHMM only",
            "Both",
        ],
        rotation=0,
    )

    if i == 0:
        ax.set_ylabel("Normalized signal")
    else:
        ax.set_ylabel("")

fig.supxlabel("ORI caller")

plt.tight_layout()

plt.savefig(
    "ori_caller_signal_boxplots.pdf",
    dpi=300,
)

plt.show()
