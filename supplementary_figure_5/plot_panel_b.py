import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

samples = [
	"alpha_rep1",
	"alpha_rep2",
	"alpha_rep3",
	"epsilon_rep1",
	"epsilon_rep2",
	"epsilon_rep3",
	"rep1_crick",
	"rep2_crick",
	"rep3_crick",
	"rep1_watson",
	"rep2_watson",
	"rep3_watson"
]

chromosomes = ["chr"+str(i) for i in range(1, 23)]
names = ["chromosome","start","end","pos","neg","log2","spline","derivative","antiderivative","true_time"]
true_time = pd.read_csv(
		"/home/dos02bi/koenig_data/projects/rt_prediction/paper_code_for_github/data/sdfs/imbulrich20240401wt.tsv", 
		sep="\t", 
		names=names
	)[["chromosome","start","end","true_time"]]

path = Path("/home/dos02bi/koenig_data/projects/rt_prediction/paper_code_for_github/scripts/06_revision/pu_seq/src/models/")

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

	if chromosome not in chromosomes:
		continue

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

names = ["chromosome","start","end","pos","neg","log2","spline","derivative","antiderivative","predicted_rt"]

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
			["chromosome", "start", "end", "true_time"]
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
		["bin", "true_time"]
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

import pandas as pd
from sklearn.metrics import mean_absolute_error, root_mean_squared_error

results = []

for sample in samples:

    filepath = (
        f"./comparisons/self_predictions/"
        f"train_human_puseq_{sample}_pred_human_puseq_{sample}.tsv"
    )

    data = pd.read_csv(filepath, sep="\t")

    pred_col = (
        f"train_human_puseq_{sample}"
        f"_pred_human_puseq_{sample}"
    )

    # extract chromosome from bin column (chr1:0-10000)
    data["chromosome"] = data["bin"].str.split(":").str[0]

    for chromosome, sub in data.groupby("chromosome"):

        results.append({
            "sample": sample,
            "chromosome": chromosome,
            "n_bins": len(sub),
            "pred_vs_hct116": sub[pred_col].corr(
                sub["true_time"]
            ),
			"pred_vs_hct116_mae": mean_absolute_error(sub["true_time"], sub[pred_col]),
			"pred_vs_hct116_rmse": root_mean_squared_error(sub["true_time"], sub[pred_col])
        })


chrom_corr = pd.DataFrame(results)

chrom_corr = chrom_corr.assign(
    replicate=chrom_corr["sample"].str.extract(r"(rep\d+)"),
    condition=chrom_corr["sample"].str.extract(r"^(.*?)_?rep\d+_?(.*?)$")
                               .apply(lambda x: x[0] if x[0] else x[1], axis=1)
)

fig = plt.figure(figsize=(8, 2))
sns.boxplot(
	data=chrom_corr,
	x="replicate",
	hue="condition",
	y="pred_vs_hct116",
	hue_order=["alpha", "epsilon", "watson", "crick"]
)
plt.ylim(0,1)
plt.ylabel("pearson correlation coefficient\n(predicted vs. gold-standard RT)")
plt.savefig("./pearson_correlation_boxplot.pdf", dpi=300)
plt.show()

fig = plt.figure(figsize=(8, 2))
sns.boxplot(
	data=chrom_corr,
	x="replicate",
	hue="condition",
	y="pred_vs_hct116_mae",
	hue_order=["alpha", "epsilon", "watson", "crick"]
)
plt.ylim(0,1)
plt.ylabel("mean absolute error\n(predicted vs. gold-standard RT)")
plt.savefig("./mae_boxplot.pdf", dpi=300)
plt.show()

fig = plt.figure(figsize=(8, 2))
sns.boxplot(
	data=chrom_corr,
	x="replicate",
	hue="condition",
	y="pred_vs_hct116_rmse",
	hue_order=["alpha", "epsilon", "watson", "crick"]
)
plt.ylim(0,1)
plt.ylabel("root mean squared error\n(predicted vs. gold-standard RT)")
plt.savefig("./rmse_boxplot.pdf", dpi=300)
plt.show()
