import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import upsetplot
from scipy.stats import ttest_ind
from statannotations.Annotator import Annotator

samples = [
	"sfb_schick_2026_01_01_WT_C_Rep1_oris_ters",
	"sfb_schick_2026_01_02_WT_HU_Rep1_oris_ters",
	"sfb_schick_2026_01_03_WT_C_Rep2_oris_ters",
	"sfb_schick_2026_01_04_WT_HU_Rep2_oris_ters",
]

oris = pd.read_csv("../data/ori_ter/merged_oris.bed",sep="\t",header=None,names=["chromosome","start","end","name"])

for sample in samples:
	oris[sample] = oris.name.str.contains(sample)

upset = upsetplot.UpSet(oris[samples].set_index(samples),show_counts=True,sort_by="cardinality")
upset.style_subsets(present=["sfb_schick_2026_01_01_WT_C_Rep1_oris_ters", "sfb_schick_2026_01_02_WT_HU_Rep1_oris_ters", "sfb_schick_2026_01_03_WT_C_Rep2_oris_ters", "sfb_schick_2026_01_04_WT_HU_Rep2_oris_ters"], absent=[], facecolor="blue", label="Constitutive")
upset.style_subsets(present=["sfb_schick_2026_01_01_WT_C_Rep1_oris_ters", "sfb_schick_2026_01_02_WT_HU_Rep1_oris_ters", "sfb_schick_2026_01_03_WT_C_Rep2_oris_ters", ], absent=["sfb_schick_2026_01_04_WT_HU_Rep2_oris_ters"], facecolor="blue", label="Constitutive")
upset.style_subsets(present=["sfb_schick_2026_01_01_WT_C_Rep1_oris_ters", "sfb_schick_2026_01_02_WT_HU_Rep1_oris_ters",  "sfb_schick_2026_01_04_WT_HU_Rep2_oris_ters"], absent=["sfb_schick_2026_01_03_WT_C_Rep2_oris_ters",], facecolor="blue", label="Constitutive")
upset.style_subsets(present=["sfb_schick_2026_01_01_WT_C_Rep1_oris_ters",  "sfb_schick_2026_01_03_WT_C_Rep2_oris_ters", "sfb_schick_2026_01_04_WT_HU_Rep2_oris_ters"], absent=["sfb_schick_2026_01_02_WT_HU_Rep1_oris_ters",], facecolor="blue", label="Constitutive")
upset.style_subsets(present=[ "sfb_schick_2026_01_02_WT_HU_Rep1_oris_ters", "sfb_schick_2026_01_03_WT_C_Rep2_oris_ters", "sfb_schick_2026_01_04_WT_HU_Rep2_oris_ters"], absent=["sfb_schick_2026_01_01_WT_C_Rep1_oris_ters",], facecolor="blue", label="Constitutive")
upset.style_subsets(present=["sfb_schick_2026_01_01_WT_C_Rep1_oris_ters", "sfb_schick_2026_01_03_WT_C_Rep2_oris_ters"], absent=["sfb_schick_2026_01_02_WT_HU_Rep1_oris_ters", "sfb_schick_2026_01_04_WT_HU_Rep2_oris_ters"], facecolor="darkgreen", label="WT_only")
upset.style_subsets(present=["sfb_schick_2026_01_01_WT_C_Rep1_oris_ters"], absent=["sfb_schick_2026_01_02_WT_HU_Rep1_oris_ters", "sfb_schick_2026_01_04_WT_HU_Rep2_oris_ters", "sfb_schick_2026_01_03_WT_C_Rep2_oris_ters"], facecolor="darkgreen", label="WT_only")
upset.style_subsets(present=["sfb_schick_2026_01_03_WT_C_Rep2_oris_ters"], absent=["sfb_schick_2026_01_01_WT_C_Rep1_oris_ters","sfb_schick_2026_01_02_WT_HU_Rep1_oris_ters", "sfb_schick_2026_01_04_WT_HU_Rep2_oris_ters"], facecolor="darkgreen", label="WT_only")
upset.style_subsets(present=["sfb_schick_2026_01_02_WT_HU_Rep1_oris_ters", "sfb_schick_2026_01_04_WT_HU_Rep2_oris_ters"], absent=["sfb_schick_2026_01_01_WT_C_Rep1_oris_ters", "sfb_schick_2026_01_03_WT_C_Rep2_oris_ters"], facecolor="darkred", label="HU_only")
upset.style_subsets(present=["sfb_schick_2026_01_02_WT_HU_Rep1_oris_ters"], absent=["sfb_schick_2026_01_01_WT_C_Rep1_oris_ters", "sfb_schick_2026_01_03_WT_C_Rep2_oris_ters", "sfb_schick_2026_01_04_WT_HU_Rep2_oris_ters"], facecolor="darkred", label="HU_only")
upset.style_subsets(present=["sfb_schick_2026_01_04_WT_HU_Rep2_oris_ters"], absent=["sfb_schick_2026_01_02_WT_HU_Rep1_oris_ters","sfb_schick_2026_01_01_WT_C_Rep1_oris_ters", "sfb_schick_2026_01_03_WT_C_Rep2_oris_ters"], facecolor="darkred", label="HU_only")
upset.plot()
plt.savefig("./ori_membership_upsetplot.pdf", bbox_inches="tight")
plt.show()

samples = [
	"sfb_schick_2026_01_01_WT_C_Rep1_oris_ters",
	"sfb_schick_2026_01_02_WT_HU_Rep1_oris_ters",
	"sfb_schick_2026_01_03_WT_C_Rep2_oris_ters",
	"sfb_schick_2026_01_04_WT_HU_Rep2_oris_ters",
]

wt = [
	"sfb_schick_2026_01_01_WT_C_Rep1_oris_ters",
	"sfb_schick_2026_01_03_WT_C_Rep2_oris_ters",
]

hu = [
	"sfb_schick_2026_01_02_WT_HU_Rep1_oris_ters",
	"sfb_schick_2026_01_04_WT_HU_Rep2_oris_ters",
]

df = oris.copy()

df[samples] = df[samples].astype(bool)

def classify(row):
	wt_present = row[wt].sum()
	hu_present = row[hu].sum()
	
	return_value = None

	# all four
	if wt_present == 2 and hu_present == 2:
		return_value = ["Constitutive", wt_present+hu_present]

	# three of four
	elif wt_present + hu_present == 3:
		return_value = ["Constitutive", wt_present+hu_present]

	# WT only
	elif hu_present == 0 and wt_present > 0:
		return_value = ["WT_only", wt_present+hu_present]

	# HU only
	elif wt_present == 0 and hu_present > 0:
		return_value = ["HU_only", wt_present+hu_present]

	# everything else
	else:
		return_value = ["Other", wt_present+hu_present]

	return return_value

df[["category", "support"]] = df.apply(classify, axis=1, result_type="expand")

samples = [
    "sfb_schick_2026_01_01_WT_C_Rep1_oris_ters",
    "sfb_schick_2026_01_02_WT_HU_Rep1_oris_ters",
    "sfb_schick_2026_01_03_WT_C_Rep2_oris_ters",
    "sfb_schick_2026_01_04_WT_HU_Rep2_oris_ters",
]

df = oris.copy()
df["midpoint"] = (df["start"] + df["end"]) / 2

iods = []

for sample in samples:
    present = df[df[sample]]

    for chrom, grp in present.groupby("chromosome"):
        grp = grp.sort_values("midpoint")

        distances = grp["midpoint"].diff().dropna()

        iods.extend(
            pd.DataFrame({
                "sample": sample,
                "chromosome": chrom,
                "IOD": distances
            }).to_dict("records")
        )

iods = pd.DataFrame(iods)

palette = {
    "sfb_schick_2026_01_01_WT_C_Rep1_oris_ters": "#0072B2",
    "sfb_schick_2026_01_02_WT_HU_Rep1_oris_ters": "#E69F00",
    "sfb_schick_2026_01_03_WT_C_Rep2_oris_ters": "#56B4E9",
    "sfb_schick_2026_01_04_WT_HU_Rep2_oris_ters": "#F0E442",
}

# chromosome-wise mean IOD
chrom_iod = (
    iods
    .groupby(["sample", "chromosome"], as_index=False)["IOD"]
    .mean()
)

order = [
    "sfb_schick_2026_01_01_WT_C_Rep1_oris_ters",
    "sfb_schick_2026_01_03_WT_C_Rep2_oris_ters",
    "sfb_schick_2026_01_02_WT_HU_Rep1_oris_ters",
    "sfb_schick_2026_01_04_WT_HU_Rep2_oris_ters",
]

# optional: log-transform before testing
chrom_iod["logIOD"] = np.log10(chrom_iod["IOD"])

plt.figure(figsize=(4,4))

ax = sns.boxplot(
    data=chrom_iod,
    x="sample",
    y="IOD",
    palette=palette,
	order=order,
    showfliers=False
)

sns.stripplot(
    data=chrom_iod,
    x="sample",
    y="IOD",
    color="black",
	order=order,
    alpha=0.6,
    size=4,
    ax=ax
)

pairs = [
    # WT replicates
    (
        "sfb_schick_2026_01_01_WT_C_Rep1_oris_ters",
        "sfb_schick_2026_01_03_WT_C_Rep2_oris_ters"
    ),
    # HU replicates
    (
        "sfb_schick_2026_01_02_WT_HU_Rep1_oris_ters",
        "sfb_schick_2026_01_04_WT_HU_Rep2_oris_ters"
    ),
    # Replicate 1: WT vs HU
    (
        "sfb_schick_2026_01_01_WT_C_Rep1_oris_ters",
        "sfb_schick_2026_01_02_WT_HU_Rep1_oris_ters"
    ),
    # Replicate 2: WT vs HU
    (
        "sfb_schick_2026_01_03_WT_C_Rep2_oris_ters",
        "sfb_schick_2026_01_04_WT_HU_Rep2_oris_ters"
    ),
]

annotator = Annotator(
    ax,
    pairs,
    data=chrom_iod,
    x="sample",
    y="IOD",      # perform statistics on log-transformed values
	order=order
)

annotator.configure(
    test="t-test_welch",
    text_format="star",
    loc="inside",
)

annotator.apply_and_annotate()

# restore plotting on original scale
ax.set_yscale("log")
ax.set_ylabel("Mean chromosome inter-origin distance (bp)")
ax.set_xlabel("")
ax.set_xticklabels([
    "WT\nRep1",
    "WT\nRep2",
    "HU\nRep1",
    "HU\nRep2"
])
plt.savefig("./chromosome_mean_IOD_boxplot.pdf", bbox_inches="tight")
plt.tight_layout()
plt.show()

import pandas as pd
import numpy as np

wt_samples = [
    "sfb_schick_2026_01_01_WT_C_Rep1_oris_ters",
    "sfb_schick_2026_01_03_WT_C_Rep2_oris_ters"
]

hu_samples = [
    "sfb_schick_2026_01_02_WT_HU_Rep1_oris_ters",
    "sfb_schick_2026_01_04_WT_HU_Rep2_oris_ters"
]


df = oris.copy()
df["midpoint"] = (df["start"] + df["end"]) / 2


# origins present in both HU but not WT (robust HU-specific set)
new_hu = df[
    (df[hu_samples].sum(axis=1) >= 1) &
    (df[wt_samples].sum(axis=1) == 0)
].copy()


wt_gaps = []

for chrom, grp in df[df[wt_samples].any(axis=1)].groupby("chromosome"):

    grp = grp.sort_values("midpoint")

    mids = grp["midpoint"].values

    for i in range(len(mids)-1):

        wt_gaps.append({
            "chromosome": chrom,
            "start": mids[i],
            "end": mids[i+1],
            "WT_IOD": mids[i+1]-mids[i]
        })

wt_gaps = pd.DataFrame(wt_gaps)

# Assign WT gap IOD to each HU origin
new_hu["WT_gap_IOD"] = np.nan

for idx, ori in new_hu.iterrows():

    gaps = wt_gaps[
        (wt_gaps.chromosome == ori.chromosome) &
        (wt_gaps.start < ori.midpoint) &
        (wt_gaps.end > ori.midpoint)
    ]

    if len(gaps) == 1:
        new_hu.loc[idx, "WT_gap_IOD"] = gaps.WT_IOD.values[0]


# Identify WT gaps containing HU origins
hu_gap_ids = set()

for _, ori in new_hu.dropna(subset=["WT_gap_IOD"]).iterrows():

    gaps = wt_gaps[
        (wt_gaps.chromosome == ori.chromosome) &
        (wt_gaps.start < ori.midpoint) &
        (wt_gaps.end > ori.midpoint)
    ]

    if len(gaps) == 1:
        hu_gap_ids.add(gaps.index[0])


# Build plotting dataframe
plot_df = pd.concat([
    wt_gaps[["WT_IOD"]]
        .assign(group="All WT gaps"),

    wt_gaps.loc[list(hu_gap_ids), ["WT_IOD"]]
        .assign(group="WT gaps with HU origin"),

    wt_gaps.drop(index=list(hu_gap_ids))[["WT_IOD"]]
        .assign(group="WT gaps without HU origin")
])


# log transform for statistics
plot_df["log_WT_IOD"] = np.log10(plot_df["WT_IOD"])


order = [
    "All WT gaps",
    "WT gaps with HU origin",
    "WT gaps without HU origin"
]


plt.figure(figsize=(3,4))

ax = sns.boxplot(
    data=plot_df,
    x="group",
    y="WT_IOD",
    order=order,
    showfliers=False
)


from matplotlib.ticker import FuncFormatter

def format_bp(value, pos):
    if value >= 1e6:
        return f"{value/1e6:.1f} Mb"
    elif value >= 1e3:
        return f"{value/1e3:.0f} kb"
    else:
        return f"{value:.0f} bp"

ax.yaxis.set_major_formatter(
    FuncFormatter(format_bp)
)

pairs = [
	(
		"All WT gaps",
		"WT gaps with HU origin"
	),
	(
		"All WT gaps",
		"WT gaps without HU origin"
	),
	(
		"WT gaps with HU origin",
		"WT gaps without HU origin"
	)
]

# Add statistics
annotator = Annotator(
    ax,
    pairs,
    data=plot_df,
    x="group",
    y="WT_IOD",
    order=order
)

annotator.configure(
    test="t-test_welch",
    text_format="star",
    loc="inside",
    verbose=1
)

annotator.apply_and_annotate()


ax.set_ylabel("WT inter-origin distance")
ax.set_xlabel("")

plt.xticks(
    rotation=30,
    ha="right"
)

plt.tight_layout()
plt.savefig(
    "./WT_IOD_boxplot.pdf",
    bbox_inches="tight"
)

plt.show()

rng = np.random.default_rng(42)

min_distance = 50000  # bp, adjust (e.g. 10 kb, 20 kb, 50 kb)

n_random = 1000   # per HU origin

random_positions_min_dist = []

for idx, row in hu_position.iterrows():

    gap = row["WT_gap"]

    # skip gaps that cannot accommodate the exclusion zone
    if gap <= 2 * min_distance:
        continue

    # allowed interval in absolute bp
    allowed_start = min_distance
    allowed_end = gap - min_distance

    # generate random absolute positions
    rand_abs = rng.uniform(
        allowed_start,
        allowed_end,
        n_random
    )

    # convert to relative position
    rand_rel = rand_abs / gap

    random_positions_min_dist.append(
        pd.DataFrame({
            "origin_id": idx,
            "WT_gap": gap,
            "relative_position": rand_rel
        })
    )


random_positions_min_dist = pd.concat(
    random_positions_min_dist,
    ignore_index=True
)

random_positions_min_dist.head()

obs = np.sort(
    hu_position.relative_position
)

rnd = np.sort(
    random_positions_min_dist
    .groupby("origin_id")
    .sample(1, random_state=1)
    .relative_position
)

plt.figure(figsize=(3,3))

plt.plot(
    obs,
    np.arange(len(obs))/len(obs),
    label="Observed"
)

plt.plot(
    rnd,
    np.arange(len(rnd))/len(rnd),
    label="Random + 50 kb"
)

plt.xlabel("Relative position")
plt.ylabel("ECDF")
plt.gca().set_aspect('equal', adjustable='box')
plt.savefig("./HU_origin_relative_position_ecdf.pdf", bbox_inches="tight")
plt.legend()
plt.show()