from sklearn.metrics import mean_absolute_error, root_mean_squared_error
from scipy.stats import pearsonr
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


def plot_rfd_oem_bedgraph(
    chrom,
    start,
    end,
    dfs_replicates,
    truth,
    colors,
    truth_color="black",
    outfile=None,
    labels=None,
    max_points=100_000,
):

    x_len = end - start
    x = np.arange(start, end)

    n_reps = len(dfs_replicates)

    if labels is None:
        labels = ["alpha", "epsilon", "watson", "crick"]

    # ---------- Ground truth ----------
    truth_values = np.full(x_len, np.nan)

    sub_t = truth[
        (truth["chrom"] == chrom)
        & (truth["end"] > start)
        & (truth["start"] < end)
    ]

    for si, ei, v in zip(
        np.maximum(sub_t["start"].values, start) - start,
        np.minimum(sub_t["end"].values, end) - start,
        sub_t["value"].values,
    ):
        truth_values[si:ei] = v


    fig, axes = plt.subplots(
        nrows=n_reps,
        figsize=(14, 3*n_reps),
        sharex=True,
        sharey=True,
    )

    if n_reps == 1:
        axes = [axes]


    for rep_idx, (dfs, ax) in enumerate(zip(dfs_replicates, axes)):

        # truth line
        if x_len > max_points:
            factor = int(np.ceil(x_len / max_points))
            x_plot = x[::factor]
            truth_plot = truth_values[::factor]
        else:
            x_plot = x
            truth_plot = truth_values


        ax.plot(
            x_plot,
            truth_plot,
            color=truth_color,
            lw=1.5,
            label="Ground truth",
        )


        # each of the four predictions
        for df, color, label in zip(dfs, colors, labels):

            values = np.full(x_len, np.nan)

            sub = df[
                (df["chrom"] == chrom)
                & (df["end"] > start)
                & (df["start"] < end)
            ]

            for si, ei, v in zip(
                np.maximum(sub["start"].values, start) - start,
                np.minimum(sub["end"].values, end) - start,
                sub["value"].values,
            ):
                values[si:ei] = v


            # metrics on full resolution
            mask = np.isfinite(truth_values) & np.isfinite(values)

            r = pearsonr(
                truth_values[mask],
                values[mask],
            )[0]

            mae = mean_absolute_error(
                truth_values[mask],
                values[mask],
            )

            rmse = root_mean_squared_error(
                truth_values[mask],
                values[mask],
            )


            # downsample only plotting
            if x_len > max_points:
                values_plot = values[::factor]
            else:
                values_plot = values


            ax.plot(
                x_plot,
                values_plot,
                color=color,
                lw=1,
                label=(
                    f"{label} "
                    f"(r={r:.2f}, MAE={mae:.2f}, RMSE={rmse:.2f})"
                ),
            )


        ax.set_ylim(-1, 1)
        ax.set_ylabel(f"rep{rep_idx+1}")
        ax.legend(
            fontsize=7,
            frameon=False,
            loc="upper right",
        )


    axes[-1].set_xlim(start, end)

    fig.suptitle(
        f"{chrom}:{start:,}-{end:,}",
        fontsize=12,
    )

    fig.supxlabel(
        f"Genomic position on {chrom}"
    )

    fig.supylabel(
        "Replication timing"
    )

    plt.tight_layout(rect=[0,0,1,0.95])


    if outfile:
        plt.savefig(
            outfile,
            dpi=300,
            bbox_inches="tight",
        )

    plt.show()
    plt.close(fig)

names = ["chrom","start","end","pos","neg","log2","spline","deri","antideri","value"]
true_time = pd.read_csv("PATH/data/sdfs/imbulrich20240401wt.tsv", sep="\t", names=names)

chromosome = "chr20"
start = 0
end = 64_444_167


alpha_rep1 = pd.read_csv(f"PATH/scripts/06_revision/pu_seq/src/models/human_puseq_alpha_rep1_delta_{chromosome}/alpha_rep1_pred.tsv", sep="\t", names=names)
epsilon_rep1 = pd.read_csv(f"PATH/scripts/06_revision/pu_seq/src/models/human_puseq_epsilon_rep1_delta_{chromosome}/epsilon_rep1_pred.tsv", sep="\t", names=names)
watson_rep1 = pd.read_csv(f"PATH/scripts/06_revision/pu_seq/src/models/human_puseq_rep1_watson_delta_{chromosome}/rep1_watson_pred.tsv", sep="\t", names=names)
crick_rep1 = pd.read_csv(f"PATH/scripts/06_revision/pu_seq/src/models/human_puseq_rep1_crick_delta_{chromosome}/rep1_crick_pred.tsv", sep="\t", names=names)

alpha_rep2 = pd.read_csv(f"PATH/scripts/06_revision/pu_seq/src/models/human_puseq_alpha_rep2_delta_{chromosome}/alpha_rep2_pred.tsv", sep="\t", names=names)
epsilon_rep2 = pd.read_csv(f"PATH/scripts/06_revision/pu_seq/src/models/human_puseq_epsilon_rep2_delta_{chromosome}/epsilon_rep2_pred.tsv", sep="\t", names=names)
watson_rep2 = pd.read_csv(f"PATH/scripts/06_revision/pu_seq/src/models/human_puseq_rep2_watson_delta_{chromosome}/rep2_watson_pred.tsv", sep="\t", names=names)
crick_rep2 = pd.read_csv(f"PATH/scripts/06_revision/pu_seq/src/models/human_puseq_rep2_crick_delta_{chromosome}/rep2_crick_pred.tsv", sep="\t", names=names)

alpha_rep3 = pd.read_csv(f"PATH/scripts/06_revision/pu_seq/src/models/human_puseq_alpha_rep3_delta_{chromosome}/alpha_rep3_pred.tsv", sep="\t", names=names)
epsilon_rep3 = pd.read_csv(f"PATH/scripts/06_revision/pu_seq/src/models/human_puseq_epsilon_rep3_delta_{chromosome}/epsilon_rep3_pred.tsv", sep="\t", names=names)
watson_rep3 = pd.read_csv(f"PATH/scripts/06_revision/pu_seq/src/models/human_puseq_rep3_watson_delta_{chromosome}/rep3_watson_pred.tsv", sep="\t", names=names)
crick_rep3 = pd.read_csv(f"PATH/scripts/06_revision/pu_seq/src/models/human_puseq_rep3_crick_delta_{chromosome}/rep3_crick_pred.tsv", sep="\t", names=names)

colors = {
    "alpha": "#0072B2",
    "epsilon": "#56B4E9",
    "watson": "#D55E00",
    "crick": "#E69F00",
}

plot_rfd_oem_bedgraph(
    chrom="chr20",
    start=0,
    end=64_444_167,
    dfs_replicates=[
        [alpha_rep1, epsilon_rep1, watson_rep1, crick_rep1],
        [alpha_rep2, epsilon_rep2, watson_rep2, crick_rep2],
        [alpha_rep3, epsilon_rep3, watson_rep3, crick_rep3],
    ],
    truth=true_time,
    colors=list(colors.values()),
    labels=list(colors.keys()),
    outfile="replication_timing_chr20.pdf",
)
