import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

def plot_rfd_oem_bedgraph(
    chrom, 
    start, 
    end,
    df1, 
    truth,
    colors=("darkred",),
    truth_color="black",
    outfile=None,
    labels=None,
    max_points=100_000
):
    x_len = end - start
    x = np.arange(start, end)

    dfs = [df1]
    n = len(dfs)

    if labels is None:
        labels = [f"DF{i+1}" for i in range(n)]

    # --- truth track ---
    truth_values = np.zeros(x_len, dtype=float)
    sub_t = truth[(truth["chrom"] == chrom) & (truth["end"] > start) & (truth["start"] < end)]

    if not sub_t.empty:
        s_idx = np.maximum(sub_t["start"].values, start) - start
        e_idx = np.minimum(sub_t["end"].values, end) - start
        val = sub_t["value"].values
        for si, ei, v in zip(s_idx, e_idx, val):
            truth_values[si:ei] = v

    # --- plotting ---
    fig, axes = plt.subplots(nrows=n, ncols=1, sharex=True, figsize=(12, 3*n))

    # ensure axes is iterable
    if n == 1:
        axes = [axes]

    for i, (df, ax, color, label) in enumerate(zip(dfs, axes, colors, labels)):

        values = np.zeros(x_len, dtype=float)
        sub = df[(df["chrom"] == chrom) & (df["end"] > start) & (df["start"] < end)]

        if not sub.empty:
            s_idx = np.maximum(sub["start"].values, start) - start
            e_idx = np.minimum(sub["end"].values, end) - start
            val = sub["value"].values
            for si, ei, v in zip(s_idx, e_idx, val):
                values[si:ei] = v

        # --- downsampling ---
        if x_len > max_points:
            factor = int(np.ceil(x_len / max_points))
            values_ds = values[::factor]
            truth_ds = truth_values[::factor]
            x_ds = x[::factor]
        else:
            values_ds = values
            truth_ds = truth_values
            x_ds = x

        # --- plotting ---
        ax.plot(x_ds, truth_ds, color=truth_color, linewidth=1, label="Ground truth")
        ax.plot(x_ds, values_ds, color=color, linewidth=1, label=label)

        ax.set_ylim(-1, 1)
        ax.set_xlim(start, end)
        ax.set_ylabel(label)
        ax.legend(loc="upper right", fontsize=6)

    fig.suptitle(f"{chrom}:{start:,}-{end:,}", fontsize=10)
    fig.supxlabel(f"Genomic position on {chrom}", fontsize=9)
    fig.supylabel("Replication timing", fontsize=9)

    plt.tight_layout(rect=[0, 0, 1, 0.95])

    if outfile:
        plt.savefig(outfile, dpi=300)

    plt.show()
    plt.close()

names = ["chrom","start","end","pos","neg","log2","spline","deri","antideri","value"]
true_time = pd.read_csv("PATH/data/sdfs/imbulrich20240401wt.tsv", sep="\t", names=names)

chromosome = "chr3"
start = 0
end = 60_000_000
train_on = "human_traelseq_sfbulrich20250101wt"
pred_on = "sfbulrich20250102wtrad21"

df1 = pd.read_csv(f"PATH/data/cross_prediction/models/{train_on}_delta_{chromosome}/{pred_on}_pred.tsv", sep="\t", names=names)

plot_rfd_oem_bedgraph(
    chrom=chromosome,
    start=start,
    end=end,
    df1=df1,
    truth=true_time,
    colors=("darkred",),
    truth_color="goldenrod",
	outfile=f"PATH/scripts/04_plots_analyses/examples/cross_examples{train_on}_{pred_on}_{chromosome}.pdf",
)

print(pearsonr(
	true_time.query("(chrom==@chromosome) & (start>=@start) & (end<=@end)").value,
	df1.query("(chrom==@chromosome) & (start>=@start) & (end<=@end)").value
)[0])

print(mean_absolute_error(
	true_time.query("(chrom==@chromosome) & (start>=@start) & (end<=@end)").value,
	df1.query("(chrom==@chromosome) & (start>=@start) & (end<=@end)").value
))
print(root_mean_squared_error(
	true_time.query("(chrom==@chromosome) & (start>=@start) & (end<=@end)").value,
	df1.query("(chrom==@chromosome) & (start>=@start) & (end<=@end)").value
))

names = ["chrom","start","end","pos","neg","log2","spline","deri","antideri","value"]
true_time = pd.read_csv("PATH/data/sdfs/imbulrich201807s16.tsv", sep="\t", names=names)

chromosome = "chrXIII"
start = 200000
end = 924431
train_on = "yeast_okseq_GSM835650"
pred_on = "imbulrich201807s16"

df1 = pd.read_csv(f"PATH/data/cross_prediction/models/{train_on}_delta_{chromosome}/{pred_on}_pred.tsv", sep="\t", names=names)

plot_rfd_oem_bedgraph(
    chrom=chromosome,
    start=start,
    end=end,
    df1=df1,
    truth=true_time,
    colors=("darkred",),
    truth_color="goldenrod",
	outfile=f"PATH/scripts/04_plots_analyses/examples/cross_examples{train_on}_{pred_on}_{chromosome}.pdf",
)

print(pearsonr(
	true_time.query("(chrom==@chromosome) & (start>=@start) & (end<=@end)").value,
	df1.query("(chrom==@chromosome) & (start>=@start) & (end<=@end)").value
)[0])

print(mean_absolute_error(
	true_time.query("(chrom==@chromosome) & (start>=@start) & (end<=@end)").value,
	df1.query("(chrom==@chromosome) & (start>=@start) & (end<=@end)").value
))
print(root_mean_squared_error(
	true_time.query("(chrom==@chromosome) & (start>=@start) & (end<=@end)").value,
	df1.query("(chrom==@chromosome) & (start>=@start) & (end<=@end)").value
))

names = ["chrom","start","end","pos","neg","log2","spline","deri","antideri","value"]
true_time = pd.read_csv("PATH/data/sdfs/imbulrich201807s16.tsv", sep="\t", names=names)

chromosome = "chrVIII"
start = 0
end = 562643
train_on = "yeast_traelseq_GSM4680460"
pred_on = "GSM4680461"

df1 = pd.read_csv(f"PATH/data/cross_prediction/models/{train_on}_delta_{chromosome}/{pred_on}_pred.tsv", sep="\t", names=names)

plot_rfd_oem_bedgraph(
    chrom=chromosome,
    start=start,
    end=end,
    df1=df1,
    truth=true_time,
    colors=("darkred",),
    truth_color="goldenrod",
	outfile=f"PATH/scripts/04_plots_analyses/examples/cross_examples{train_on}_{pred_on}_{chromosome}.pdf",
)

print(pearsonr(
	true_time.query("(chrom==@chromosome) & (start>=@start) & (end<=@end)").value,
	df1.query("(chrom==@chromosome) & (start>=@start) & (end<=@end)").value
)[0])

print(mean_absolute_error(
	true_time.query("(chrom==@chromosome) & (start>=@start) & (end<=@end)").value,
	df1.query("(chrom==@chromosome) & (start>=@start) & (end<=@end)").value
))
print(root_mean_squared_error(
	true_time.query("(chrom==@chromosome) & (start>=@start) & (end<=@end)").value,
	df1.query("(chrom==@chromosome) & (start>=@start) & (end<=@end)").value
))

names = ["chrom","start","end","pos","neg","log2","spline","deri","antideri","value"]
true_time = pd.read_csv("PATH/data/sdfs/imbulrich201807s16.tsv", sep="\t", names=names)

chromosome = "chrX"
start = 450000
end = 745751
train_on = "human_gloeseq_GSM3939127"
pred_on = "GSM4680461"

df1 = pd.read_csv(f"PATH/data/cross_prediction/models/{train_on}_delta_none/{pred_on}_pred.tsv", sep="\t", names=names)

plot_rfd_oem_bedgraph(
    chrom=chromosome,
    start=start,
    end=end,
    df1=df1,
    truth=true_time,
    colors=("darkred",),
    truth_color="goldenrod",
	outfile=f"PATH/scripts/04_plots_analyses/examples/cross_examples{train_on}_{pred_on}_{chromosome}.pdf",
)

print(pearsonr(
	true_time.query("(chrom==@chromosome) & (start>=@start) & (end<=@end)").value,
	df1.query("(chrom==@chromosome) & (start>=@start) & (end<=@end)").value
)[0])

print(mean_absolute_error(
	true_time.query("(chrom==@chromosome) & (start>=@start) & (end<=@end)").value,
	df1.query("(chrom==@chromosome) & (start>=@start) & (end<=@end)").value
))
print(root_mean_squared_error(
	true_time.query("(chrom==@chromosome) & (start>=@start) & (end<=@end)").value,
	df1.query("(chrom==@chromosome) & (start>=@start) & (end<=@end)").value
))