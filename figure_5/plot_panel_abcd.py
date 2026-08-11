import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import seaborn as sns

data = pd.read_csv("./comparisons/all_vs_all_pearsonr.tsv", sep='\t')

inclusion_list = [
	"GSM835651",
	"GSM835650",
	"GSM4680462",
	"GSM5005341",
	"GSM5005340",
	"GSM5005342",
	"imbulrich201807s16",
	"imbulrich201807s18",
	"imbulrich201807s20",
	"GSM835652",
	"GSM835653",
	"GSM4680463",
	"GSM5005319",
	"GSM5005318",
	"GSM5005335",
	"GSM5005336",
	"GSM4680453",
	"GSM5005322",
	"GSM4680452",
	"GSM4680465",
	"GSM4680464",
	"GSM4680468",
	"GSM4680469",
	"GSM5005339",
	"GSM4680461",
	"GSM4680460",
	"imbulrich20240401wt",
	"sfbulrich20250101wt",
	"sfbulrich20250102wtrad21",
	"GSM3939127",
	"GSM3939128",
	"GSM4305465",
	"GSM4305465",
]

data = data.query("train_sample in @inclusion_list").query("pred_sample in @inclusion_list")

meta = (
    data[["train_sample", "train_organism", "train_experiment"]]
    .drop_duplicates()
    .rename(columns={
        "train_sample": "sample",
        "train_organism": "organism",
        "train_experiment": "experiment"
    })
)

meta = meta.sort_values(
    by=["organism", "experiment", "sample"]
)

ordered_samples = meta["sample"].tolist()
meta = meta.set_index("sample")

heatmap_matrix = data.pivot(
    index="train_sample",
    columns="pred_sample",
    values="pearsonr"
)

heatmap_matrix = heatmap_matrix.reindex(
    index=ordered_samples,
    columns=ordered_samples
)


organisms = meta["organism"].unique()
experiments = meta["experiment"].unique()

organism_palette = dict(
    zip(organisms, sns.color_palette("Set2", len(organisms)))
)

experiment_palette = dict(
    zip(experiments, sns.color_palette("tab10", len(experiments)))
)

row_colors = pd.DataFrame({
    "Organism": meta.loc[ordered_samples, "organism"].map(organism_palette),
    "Experiment": meta.loc[ordered_samples, "experiment"].map(experiment_palette)
})

col_colors = row_colors.copy()

g = sns.clustermap(
    heatmap_matrix,
    cmap=sns.color_palette("vlag", as_cmap=True),
    center=0.5,
    vmin=0,
    vmax=1,
    row_cluster=False,
    col_cluster=False,
    row_colors=row_colors,
    col_colors=col_colors,
    linewidths=0,
    figsize=(16, 16),
    cbar_pos=(0.92, 0.3, 0.02, 0.4)  # move colorbar to right side
)

# Axis labels
g.ax_heatmap.set_xlabel("Predicted sample")
g.ax_heatmap.set_ylabel("Training sample")
g.ax_heatmap.set_title("Replication Timing Cross-Sample Prediction")

for spine in g.ax_heatmap.spines.values():
    spine.set_visible(True)
    spine.set_linewidth(1)
    spine.set_edgecolor("black")

organism_handles = [
    Patch(facecolor=color, label=org)
    for org, color in organism_palette.items()
]

experiment_handles = [
    Patch(facecolor=color, label=exp)
    for exp, color in experiment_palette.items()
]

# Place legends outside heatmap
g.ax_heatmap.legend(
    handles=organism_handles,
    title="Organism",
    bbox_to_anchor=(-0.1, 0.6),
    loc="upper left",
    frameon=False
)

g.ax_heatmap.legend(
    handles=experiment_handles,
    title="Experiment",
    bbox_to_anchor=(-0.3, 0.6),
    loc="upper left",
    frameon=False
)

plt.savefig("cross_sample_heatmap_publication_ready_small.pdf", dpi=300)
plt.show()

# self predictions
col_long = (heatmap_matrix.reset_index().melt(id_vars="train_sample", var_name="pred_sample", value_name="pearsonr"))
row_long = (heatmap_matrix.reset_index().melt(id_vars="train_sample", var_name="pred_sample", value_name="pearsonr"))

self_pred = col_long[col_long["train_sample"] == col_long["pred_sample"]]
self_pred_row = row_long[row_long["train_sample"] == row_long["pred_sample"]]

# remove them from jitter data
col_long_no_self = col_long[col_long["train_sample"] != col_long["pred_sample"]]
row_long_no_self = row_long[row_long["train_sample"] != row_long["pred_sample"]]

sample_experiment_map = meta["experiment"].to_dict()

# Add experiment annotation
col_long["experiment"] = col_long["pred_sample"].map(sample_experiment_map)
row_long["experiment"] = row_long["train_sample"].map(sample_experiment_map)
col_long_no_self["experiment"] = col_long_no_self["pred_sample"].map(sample_experiment_map)
row_long_no_self["experiment"] = row_long_no_self["train_sample"].map(sample_experiment_map)

# Palette for experiments
experiments = meta["experiment"].unique()

experiment_palette = dict(
    zip(experiments, sns.color_palette("tab10", len(experiments)))
)

plt.figure(figsize=(18, 6))

sns.boxplot(
    data=col_long,
    x="pred_sample",
    y="rmse",
    hue="experiment",
    palette=experiment_palette,
    showfliers=False,
	legend=False
)

sns.stripplot(
    data=col_long_no_self,
    x="pred_sample",
    y="rmse",
    hue="experiment",
    palette=experiment_palette,
    dodge=False,
    jitter=0.25,
    alpha=1,
    size=4,
    edgecolor="black",
    linewidth=0.5,
    legend=False
)

# self predictions
sns.scatterplot(
    data=self_pred,
    x="pred_sample",
    y="rmse",
    color="red",
    edgecolor="black",
	alpha=1,
    s=60,
    zorder=10
)

ax = plt.gca()
ax.set_yticks([0, 0.25, 0.5, 0.75, 1])
ax.set_yticklabels(["0", "0.25", "0.5", "0.75", "1"])

plt.xticks(rotation=90)
plt.xlabel("Predicted sample")
plt.ylabel("Pearson correlation")
plt.title("Distribution of correlations per Predicted Sample")
plt.tight_layout()
plt.savefig("correlation_distribution_per_predicted_sample_small.pdf", dpi=300)
plt.show()

plt.figure(figsize=(6, 18))

sns.boxplot(
    data=row_long,
    y="train_sample",
    x="rmse",
    hue="experiment",
    palette=experiment_palette,
    showfliers=False,
	legend=False
)

sns.stripplot(
    data=row_long_no_self,
    y="train_sample",
    x="rmse",
    hue="experiment",
    palette=experiment_palette,
    dodge=False,
    jitter=0.25,
    alpha=1,
    size=4,
    edgecolor="black",
    linewidth=0.5,
    legend=False
)

sns.scatterplot(
    data=self_pred_row,
    y="train_sample",
    x="rmse",
    color="red",
    edgecolor="black",
    s=60,
    zorder=10
)

ax = plt.gca()
ax.set_xticks([0, 0.25, 0.5, 0.75, 1])
ax.set_xticklabels(["0", "0.25", "0.5", "0.75", "1"])

plt.ylabel("Training sample")
plt.xlabel("Pearson correlation")
plt.title("Distribution of correlations per Training Sample")


plt.tight_layout()
plt.savefig("correlation_distribution_per_training_sample_small.pdf", dpi=300)
plt.show()

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# ---------------------------------
# 1) Remove diagonal self-prediction
# ---------------------------------

matrix_no_diag = heatmap_matrix.copy()

for sample in matrix_no_diag.index:
    if sample in matrix_no_diag.columns:
        matrix_no_diag.loc[sample, sample] = pd.NA

# ---------------------------------
# 2) Compute medians
# ---------------------------------

median_train = matrix_no_diag.median(axis=1)
median_pred = matrix_no_diag.median(axis=0)

summary = pd.DataFrame({
    "sample": matrix_no_diag.index,
    "median_training_performance": median_train.values,
    "median_prediction_performance": median_pred.values
})

# Attach metadata
summary["organism"] = meta.loc[summary["sample"], "organism"].values
summary["experiment"] = meta.loc[summary["sample"], "experiment"].values

# ---------------------------------
# 3) Define aesthetics
# ---------------------------------

# Experiment → color
experiments = summary["experiment"].unique()
experiment_palette = dict(
    zip(experiments, sns.color_palette("tab10", len(experiments)))
)

# Organism → marker
organisms = summary["organism"].unique()
marker_map = dict(zip(
    organisms,
    ["o", "s", "^", "D", "v", "P", "X"][:len(organisms)]
))

# ---------------------------------
# 4) Plot scatter
# ---------------------------------

plt.figure(figsize=(5, 5))

for org in organisms:
    for exp in experiments:

        sub = summary[
            (summary["organism"] == org) &
            (summary["experiment"] == exp)
        ]

        if sub.empty:
            continue

        plt.scatter(
            sub["median_training_performance"],
            sub["median_prediction_performance"],
            label=f"{org} | {exp}",
            color=experiment_palette[exp],
            marker=marker_map[org],
            s=100,
            edgecolors="black",
            linewidth=0.5
        )

# Reference diagonal
lims = [
    0,
    1,
]

plt.plot(lims, lims, linestyle="--", color="gray")

plt.xlim(lims)
plt.ylim(lims)

ax = plt.gca()
ax.set_aspect("equal", adjustable="box")

plt.xlabel("Median performance as training sample")
plt.ylabel("Median performance as predicted sample")
plt.title("Training vs Prediction Performance")

# Legend (outside plot)
plt.legend(loc="upper left", fontsize=8)

plt.tight_layout()
plt.savefig("median_train_vs_pred_scatter_colored_small.pdf", dpi=300)
plt.show()
