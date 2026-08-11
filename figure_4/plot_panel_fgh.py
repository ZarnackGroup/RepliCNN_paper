import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from pySankey.sankey import sankey

samples = [
	"imbulrich20240401wt",
	"sfbulrich20250101wt",
	"sfbulrich20250102wtrad21",
	"sfb_schick_2026_01_01_WT_C_Rep1",
	"sfb_schick_2026_01_02_WT_HU_Rep1",
	"sfb_schick_2026_01_03_WT_C_Rep2",
	"sfb_schick_2026_01_04_WT_HU_Rep2",
]

chromosomes = ["chr"+str(i) for i in range(1, 23)]
names = ["chromosome","start","end","pos","neg","log2","spline","derivative","antiderivative"]
true_time = \
	pd.read_csv(
		"PATH/data/sdfs/imbulrich20240401wt.tsv", 
		sep="\t", 
		names=names+["true_time_hct116"]
	)[["chromosome","start","end","true_time_hct116"]].assign(true_time_hap1=\
	pd.read_csv(
		"PATH/scripts/06_revision/sandra_traelseq/data/tsv/sfb_schick_2026_01_01_WT_C_Rep1.tsv", 
		sep="\t", 
		names=names+["true_time_hap1"]
	)[["true_time_hap1"]])

path = Path("PATH/scripts/06_revision/sandra_traelseq/src/models")

rows = []

for p in path.iterdir():

    # ignore log files
    if p.suffix == ".log":
        continue

    stem = p.stem

    # split at "_delta_"
    left, chromosome = stem.rsplit("_delta_", 1)

    parts = left.split("_")

    organism = parts[0]
    experiment = parts[1]
    sample = "_".join(parts[2:])

    rows.append({
        "organism": organism,
        "experiment": experiment,
        "sample": sample,
        "chromosome": chromosome
    })

data = (
    pd.DataFrame(rows)
      .sort_values(["organism", "experiment", "sample", "chromosome"])
      .reset_index(drop=True)
)

names = names + ["predicted_rt"]

# self predictions
for i, sample in enumerate(samples, 1):
	print(i, len(samples), sample)

	sub_data = data.query("sample == @sample")

	organism = sub_data["organism"].iat[0]
	experiment = sub_data["experiment"].iat[0]

	chromosomes = sub_data["chromosome"].unique()

	collected = []

	for chromosome in chromosomes:
		
		if chromosome == "none":
			continue

		model_dir = Path(f"{path}/{organism}_{experiment}_{sample}_delta_{chromosome}")

		pred_file = model_dir / f"{sample}_pred.tsv"

		if not pred_file.exists():
			print(f"Missing: {pred_file}")
			continue

		tmp = (
			pd.read_csv(
				pred_file,
				sep="\t",
				names=names,
			)
			.query("chromosome == @chromosome")[["chromosome", "start", "end", "predicted_rt"]]
		)

		tmp["bin"] = (
			tmp["chromosome"].astype(str)
			+ ":"
			+ tmp["start"].astype(str)
			+ "-"
			+ tmp["end"].astype(str)
		)

		pred_col = (
			f"train_{organism}_{experiment}_{sample}"
			f"_pred_{organism}_{experiment}_{sample}"
		)

		tmp = (
			tmp.rename(columns={"predicted_rt": pred_col})
			[["bin", pred_col]]
		)

		collected.append(tmp)

	if len(collected) == 0:
		print(f"No predictions found for {sample}")
		continue

	self_predictions = pd.concat(
		collected,
		ignore_index=True
	)

	if self_predictions["bin"].duplicated().any():
		raise ValueError(f"Duplicate genomic bins detected for {sample}")

	true_tmp = (
		true_time[
			["chromosome", "start", "end", "true_time_hct116", "true_time_hap1"]
		]
		.rename(columns={"chromosome": "chromosome"})
		.copy()
	)

	true_tmp["bin"] = (
		true_tmp["chromosome"]
		+ ":"
		+ true_tmp["start"].astype(str)
		+ "-"
		+ true_tmp["end"].astype(str)
	)

	true_tmp = true_tmp[
		["bin", "true_time_hct116", "true_time_hap1"]
	]

	final = pd.merge(
		self_predictions,
		true_tmp,
		on="bin",
		how="inner",
		validate="one_to_one"
	)

	final.to_csv(
		Path("./comparisons/self_predictions")
		/ f"{pred_col}.tsv",
		sep="\t",
		index=False
	)

for i, sample in enumerate(samples, 1):
	print(i, len(samples), sample)

	sub_data = data.query("sample == @sample")

	organism = sub_data["organism"].iat[0]
	experiment = sub_data["experiment"].iat[0]
	chromosomes = sub_data["chromosome"].unique()

	for pred_sample in samples:

		if pred_sample == sample:
			continue

		pred_sub = data.query("sample == @pred_sample")

		pred_organism = pred_sub["organism"].iat[0]
		pred_experiment = pred_sub["experiment"].iat[0]

		# SAME ORGANISM ONLY
		if pred_organism != organism:
			continue

		print(f"  -> Predicting {pred_sample}")

		all_chr_results = []

		for chromosome in chromosomes:

			if chromosome == "none":
				continue
				
			model_dir = Path(
				f"{path}/{organism}_{experiment}_{sample}_delta_{chromosome}"
			)

			pred_file = model_dir / f"{pred_sample}_pred.tsv"

			if not pred_file.exists():
				print(f"Missing: {pred_file}")
				continue

			tmp = pd.read_csv(
				pred_file,
				sep="\t",
				names=names,
			)

			# normalize chromosome column name
			if "chromosome" in tmp.columns:
				tmp = tmp.rename(columns={"chromosome": "chromosome"})

			tmp = (
				tmp
				.query("chromosome == @chromosome")
				[["chromosome", "start", "end", "predicted_rt"]]
			)

			# enforce correct types
			tmp["chromosome"] = tmp["chromosome"].astype(str)
			tmp["start"] = tmp["start"].astype(int)
			tmp["end"] = tmp["end"].astype(int)

			tmp["bin"] = (
				tmp["chromosome"]
				+ ":"
				+ tmp["start"].astype(str)
				+ "-"
				+ tmp["end"].astype(str)
			)

			all_chr_results.append(
				tmp[["bin", "predicted_rt"]]
			)

		if len(all_chr_results) == 0:
			continue

		final_pred = pd.concat(
			all_chr_results,
			ignore_index=True
		)

		pred_col = (
			f"train_{organism}_{experiment}_{sample}"
			f"_pred_{pred_organism}_{pred_experiment}_{pred_sample}"
		)

		final_pred = final_pred.rename(
			columns={"predicted_rt": pred_col}
		)


		# -------------------------
		# add both true RT columns
		# -------------------------

		true_tmp = (
			true_time[
				[
					"chromosome",
					"start",
					"end",
					"true_time_hct116",
					"true_time_hap1"
				]
			]
			.rename(columns={"chromosome": "chromosome"})
			.copy()
		)

		true_tmp["bin"] = (
			true_tmp["chromosome"].astype(str)
			+ ":"
			+ true_tmp["start"].astype(int).astype(str)
			+ "-"
			+ true_tmp["end"].astype(int).astype(str)
		)

		true_tmp = true_tmp[
			[
				"bin",
				"true_time_hct116",
				"true_time_hap1"
			]
		]


		final = pd.merge(
			final_pred,
			true_tmp,
			on="bin",
			how="inner",
			validate="one_to_one"
		)


		outname = (
			Path("./comparisons/same_org_predictions")
			/ f"{pred_col}.tsv"
		)

		final.to_csv(
			outname,
			sep="\t",
			index=False
		)

import pandas as pd

results = []

for sample in samples:

    filepath = (
        f"./comparisons/self_predictions/"
        f"train_human_traelseq_{sample}_pred_human_traelseq_{sample}.tsv"
    )

    data = pd.read_csv(filepath, sep="\t")

    pred_col = (
        f"train_human_traelseq_{sample}"
        f"_pred_human_traelseq_{sample}"
    )

    # extract chromosome from bin column (chr1:0-10000)
    data["chromosome"] = data["bin"].str.split(":").str[0]

    for chromosome, sub in data.groupby("chromosome"):

        results.append({
            "sample": sample,
            "chromosome": chromosome,
            "n_bins": len(sub),
            "pred_vs_hct116": sub[pred_col].corr(
                sub["true_time_hct116"]
            ),
            "pred_vs_hap1": sub[pred_col].corr(
                sub["true_time_hap1"]
            ),
            "hct116_vs_hap1": sub["true_time_hct116"].corr(
                sub["true_time_hap1"]
            )
        })


chrom_corr = pd.DataFrame(results)


import pandas as pd
from pathlib import Path

results = []

pred_path = Path("./comparisons/same_org_predictions")

for filepath in pred_path.glob("*.tsv"):

    data = pd.read_csv(filepath, sep="\t")

    # extract prediction column
    pred_cols = [
        c for c in data.columns 
        if c.startswith("train_")
    ]

    if len(pred_cols) != 1:
        print(f"Skipping {filepath.name}")
        continue

    pred_col = pred_cols[0]

    # extract training and prediction samples from column name
    # train_human_traelseq_SAMPLE_pred_human_traelseq_SAMPLE
    name = pred_col.replace("train_human_traelseq_", "")
    
    train_sample, pred_sample = name.split("_pred_human_traelseq_")

    # extract chromosome
    data["chromosome"] = (
        data["bin"]
        .str.split(":")
        .str[0]
    )

    for chromosome, sub in data.groupby("chromosome"):

        results.append({
            "train_sample": train_sample,
            "pred_sample": pred_sample,
            "chromosome": chromosome,
            "n_bins": len(sub),

            "pred_vs_hct116": sub[pred_col].corr(
                sub["true_time_hct116"]
            ),

            "pred_vs_hap1": sub[pred_col].corr(
                sub["true_time_hap1"]
            ),

            "hct116_vs_hap1": sub["true_time_hct116"].corr(
                sub["true_time_hap1"]
            )
        })


cross_chrom_corr = pd.DataFrame(results)

from sklearn.metrics import mean_absolute_error, root_mean_squared_error
from scipy.stats import pearsonr
import numpy as np
import matplotlib.pyplot as plt


def plot_rfd_oem_bedgraph(
    chrom,
    start,
    end,
    dfs_conditions,
    truth,
    colors,
    truth_color="black",
    outfile=None,
    labels=None,
    max_points=100_000,
):

    x_len = end - start
    x = np.arange(start, end)

    n_conditions = len(dfs_conditions)

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

    fig, ax = plt.subplots(
        figsize=(12, 3),
    )

    if x_len > max_points:
        factor = int(np.ceil(x_len / max_points))
        x_plot = x[::factor]
        truth_plot = truth_values[::factor]
    else:
        x_plot = x
        truth_plot = truth_values

    def get_mean_prediction(replicates):

        rep_arrays = []

        for df in replicates:

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

            rep_arrays.append(values)

        return np.nanmean(rep_arrays, axis=0)


    # ---------- calculate means ----------
    ctrl_rep1 = get_mean_prediction(dfs_conditions[0][0])
    ctrl_rep2 = get_mean_prediction(dfs_conditions[0][1])

    hu_rep1 = get_mean_prediction(dfs_conditions[1][0])
    hu_rep2 = get_mean_prediction(dfs_conditions[1][1])


    predictions = [
        (ctrl_rep1, "Control Rep1", "#0072B2", "-"),
        (ctrl_rep2, "Control Rep2", "#0072B2", "--"),
        (hu_rep1, "HU Rep1", "#D55E00", "-"),
        (hu_rep2, "HU Rep2", "#D55E00", "--"),
    ]


    # ---------- plot predictions ----------
    for values, label, color, linestyle in predictions:

        mask = np.isfinite(truth_values) & np.isfinite(values)

        if x_len > max_points:
            values_plot = values[::factor]
        else:
            values_plot = values


        ax.plot(
            x_plot,
            values_plot,
            color=color,
            linestyle=linestyle,
            lw=1.5,
            label=(
                f"{label} "
            ),
        )


    ax.set_ylim(-1, 1)
    ax.set_xlim(start, end)

    ax.set_xlabel(
        f"Genomic position on {chrom}"
    )

    ax.set_ylabel(
        "Replication timing"
    )

    ax.set_title(
        f"{chrom}:{start:,}-{end:,}"
    )

    ax.legend(
        fontsize=8,
        frameon=False,
        loc="upper right",
    )

    plt.tight_layout()

    if outfile:
        plt.savefig(
            outfile,
            dpi=300,
            bbox_inches="tight",
        )

    plt.show()
    plt.close(fig)

names = ["chrom","start","end","pos","neg","log2","spline","deri","antideri","value"]
true_time = pd.read_csv("PATH/scripts/06_revision/sandra_traelseq/data/tsv/sfb_schick_2026_01_01_WT_C_Rep1.tsv", sep="\t", names=names)

c1_HCT116_WT_1 = pd.read_csv(f"PATH/scripts/06_revision/sandra_traelseq/src/models/human_traelseq_imbulrich20240401wt_delta_{chromosome}/sfb_schick_2026_01_01_WT_C_Rep1_pred.tsv", sep="\t", names=names)
c1_HCT116_WT_2 = pd.read_csv(f"PATH/scripts/06_revision/sandra_traelseq/src/models/human_traelseq_sfbulrich20250101wt_delta_{chromosome}/sfb_schick_2026_01_01_WT_C_Rep1_pred.tsv", sep="\t", names=names)
c1_HCT116_WT_3 = pd.read_csv(f"PATH/scripts/06_revision/sandra_traelseq/src/models/human_traelseq_sfbulrich20250102wtrad21_delta_{chromosome}/sfb_schick_2026_01_01_WT_C_Rep1_pred.tsv", sep="\t", names=names)

c2_HCT116_WT_1 = pd.read_csv(f"PATH/scripts/06_revision/sandra_traelseq/src/models/human_traelseq_imbulrich20240401wt_delta_{chromosome}/sfb_schick_2026_01_03_WT_C_Rep2_pred.tsv", sep="\t", names=names)
c2_HCT116_WT_2 = pd.read_csv(f"PATH/scripts/06_revision/sandra_traelseq/src/models/human_traelseq_sfbulrich20250101wt_delta_{chromosome}/sfb_schick_2026_01_03_WT_C_Rep2_pred.tsv", sep="\t", names=names)
c2_HCT116_WT_3 = pd.read_csv(f"PATH/scripts/06_revision/sandra_traelseq/src/models/human_traelseq_sfbulrich20250102wtrad21_delta_{chromosome}/sfb_schick_2026_01_03_WT_C_Rep2_pred.tsv", sep="\t", names=names)

hu1_HCT116_WT_1 = pd.read_csv(f"PATH/scripts/06_revision/sandra_traelseq/src/models/human_traelseq_imbulrich20240401wt_delta_{chromosome}/sfb_schick_2026_01_02_WT_HU_Rep1_pred.tsv", sep="\t", names=names)
hu1_HCT116_WT_2 = pd.read_csv(f"PATH/scripts/06_revision/sandra_traelseq/src/models/human_traelseq_sfbulrich20250101wt_delta_{chromosome}/sfb_schick_2026_01_02_WT_HU_Rep1_pred.tsv", sep="\t", names=names)
hu1_HCT116_WT_3 = pd.read_csv(f"PATH/scripts/06_revision/sandra_traelseq/src/models/human_traelseq_sfbulrich20250102wtrad21_delta_{chromosome}/sfb_schick_2026_01_02_WT_HU_Rep1_pred.tsv", sep="\t", names=names)

hu2_HCT116_WT_1 = pd.read_csv(f"PATH/scripts/06_revision/sandra_traelseq/src/models/human_traelseq_imbulrich20240401wt_delta_{chromosome}/sfb_schick_2026_01_04_WT_HU_Rep2_pred.tsv", sep="\t", names=names)
hu2_HCT116_WT_2 = pd.read_csv(f"PATH/scripts/06_revision/sandra_traelseq/src/models/human_traelseq_sfbulrich20250101wt_delta_{chromosome}/sfb_schick_2026_01_04_WT_HU_Rep2_pred.tsv", sep="\t", names=names)
hu2_HCT116_WT_3 = pd.read_csv(f"PATH/scripts/06_revision/sandra_traelseq/src/models/human_traelseq_sfbulrich20250102wtrad21_delta_{chromosome}/sfb_schick_2026_01_04_WT_HU_Rep2_pred.tsv", sep="\t", names=names)

colors = {
    "HCT116_WT_1": "#0072B2",
    "HCT116_WT_2": "#1B761D",
    "HCT116_WT_3": "#720C5C",
}

dfs_replicates=[
    [c1_HCT116_WT_1, c1_HCT116_WT_2, c1_HCT116_WT_3],
    [c2_HCT116_WT_1, c2_HCT116_WT_2, c2_HCT116_WT_3],
    [hu1_HCT116_WT_1, hu1_HCT116_WT_2, hu1_HCT116_WT_3],
    [hu2_HCT116_WT_1, hu2_HCT116_WT_2, hu2_HCT116_WT_3],
]

condition_pairs = [
    (
        [c1_HCT116_WT_1, c1_HCT116_WT_2, c1_HCT116_WT_3],
        [c2_HCT116_WT_1, c2_HCT116_WT_2, c2_HCT116_WT_3],
    ),
    (
        [hu1_HCT116_WT_1, hu1_HCT116_WT_2, hu1_HCT116_WT_3],
        [hu2_HCT116_WT_1, hu2_HCT116_WT_2, hu2_HCT116_WT_3],
    ),
]


intervals =[
	("chr2", 200_000_000, 230_000_000),
	("chr8", 30_000_000, 50_000_000),
	("chr8", 70_000_000, 100_000_000),
]
for chromosome, start, end in intervals:

	plot_rfd_oem_bedgraph(
		chrom=chromosome,
		start=start,
		end=end,
		dfs_conditions=[
			(
				[c1_HCT116_WT_1, c1_HCT116_WT_2, c1_HCT116_WT_3],
				[c2_HCT116_WT_1, c2_HCT116_WT_2, c2_HCT116_WT_3],
			),
			(
				[hu1_HCT116_WT_1, hu1_HCT116_WT_2, hu1_HCT116_WT_3],
				[hu2_HCT116_WT_1, hu2_HCT116_WT_2, hu2_HCT116_WT_3],
			),
		],
		truth=true_time,
		colors=["#0072B2", "#D55E00"],
		labels=["Control", "HU"],
		outfile=f"new_rt_hap1_{chromosome}_{start}_{end}.pdf",
)

blacklist = pd.read_csv("PATH/data/timetables/human/human_hct116_GSE137764_hg38_rt_processed.bg", names=["chrom", "start", "end", "true_time"], sep="\t").query("true_time==0")

def split_bins(df, bin_size=10000):
    expanded = []

    for _, row in df.iterrows():
        starts = range(row["start"], row["end"], bin_size)

        for s in starts:
            expanded.append({
                "chrom": row["chrom"],
                "start": s,
                "end": min(s + bin_size, row["end"]),
                "true_time": row["true_time"]
            })

    return pd.DataFrame(expanded)

blacklist = split_bins(blacklist, bin_size=10000)

def remove_matching_bins(df, exclude_df, 
                         df_chrom="chrom", 
                         exclude_chrom="chrom",
                         start_col="start", 
                         end_col="end"):

    exclude_bins = set(
        zip(
            exclude_df[exclude_chrom],
            exclude_df[start_col],
            exclude_df[end_col]
        )
    )

    mask = ~df.apply(
        lambda x: (x[df_chrom], x[start_col], x[end_col]) in exclude_bins,
        axis=1
    )

    return df.loc[mask].reset_index(drop=True)

    fig = plt.figure(figsize=(6, 4))

sns.kdeplot(data=remove_matching_bins(pd.concat([c1_HCT116_WT_1, c1_HCT116_WT_2, c1_HCT116_WT_3]), blacklist), x="value", fill=False, alpha=1, linewidth=2, color="#0072B2", linestyle="solid", label="CTRL_1")
sns.kdeplot(data=remove_matching_bins(pd.concat([c2_HCT116_WT_1, c2_HCT116_WT_2, c2_HCT116_WT_3]), blacklist), x="value", fill=False, alpha=1, linewidth=2, color="#024369", linestyle="solid", label="CTRL_2")
sns.kdeplot(data=remove_matching_bins(pd.concat([hu1_HCT116_WT_1, hu1_HCT116_WT_2, hu1_HCT116_WT_3]), blacklist), x="value", fill=False, alpha=1, linewidth=2, color="#C97200", linestyle="solid", label="HU_1")
sns.kdeplot(data=remove_matching_bins(pd.concat([hu2_HCT116_WT_1, hu2_HCT116_WT_2, hu2_HCT116_WT_3]), blacklist), x="value", fill=False, alpha=1, linewidth=2, color="#714204", linestyle="solid", label="HU_2")

plt.legend(fontsize=8, frameon=False)
plt.title(f"Predicted replication timing distribution")
plt.xlabel("Replication Timing")
plt.ylabel("Density")
plt.savefig(f"rt_hap1_distribution.pdf", dpi=300, bbox_inches="tight")
plt.show()

ctrl = remove_matching_bins(pd.concat([c1_HCT116_WT_1, c1_HCT116_WT_2, c1_HCT116_WT_3, c2_HCT116_WT_1, c2_HCT116_WT_2, c2_HCT116_WT_3]), blacklist)
hu = remove_matching_bins(pd.concat([hu1_HCT116_WT_1, hu1_HCT116_WT_2, hu1_HCT116_WT_3, hu2_HCT116_WT_1, hu2_HCT116_WT_2, hu2_HCT116_WT_3]), blacklist)

ctrl = ctrl.groupby(["chrom", "start", "end"])["value"].mean().reset_index()
hu = hu.groupby(["chrom", "start", "end"])["value"].mean().reset_index()

plot_data = ctrl.assign(hu_value=hu["value"]).assign(shift=lambda x: x["hu_value"] - x["value"]).assign(time_cat=lambda x: pd.cut(x["value"], bins=[-1, -0.5, -0.25, 0.25, 0.5, 1], labels=["Very late", "Late", "Intermediate", "Early", "Very early"])).assign(hu_time_cat=lambda x: pd.cut(x["hu_value"], bins=[-1, -0.5, -0.25, 0.25, 0.5, 1], labels=["Very late", "Late", "Intermediate", "Early", "Very early"]))

sankey(
	left=plot_data["time_cat"], 
	right=plot_data["hu_time_cat"], 
	leftLabels=["Very late", "Late", "Intermediate", "Early", "Very early"],
	rightLabels=["Very late", "Late", "Intermediate", "Early", "Very early"],
	aspect=20,  
	fontsize=12)
plt.savefig(f"rt_hap1_sankey.pdf", dpi=300, bbox_inches="tight")
plt.show()