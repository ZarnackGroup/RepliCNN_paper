import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib as mpl
import matplotlib.ticker as ticker

names = ["chrom", "start", "end", "name", "score", "strand", "thickStart", "thickEnd", "itemRgb"]
threshold = 0
wt_1 = pd.read_csv("./ori_ter/imbulrich20240401wt_oris_ters.bed", sep="\t", names=names).query("score>@threshold").query("strand=='+'")
wt_2 = pd.read_csv("./ori_ter/sfbulrich20250101wt_oris_ters.bed", sep="\t", names=names).query("score>@threshold").query("strand=='+'")
wt_3 = pd.read_csv("./ori_ter/sfbulrich20250102wtrad21_oris_ters.bed", sep="\t", names=names).query("score>@threshold").query("strand=='+'")

def compute_inter_origin_distance(df):
    
    # Work on a copy
    df = df.copy()
    
    # Compute midpoint of each ORI
    df["midpoint"] = (df["start"] + df["end"]) / 2
    
    # Sort by chromosome and genomic position
    df = df.sort_values(["chrom", "midpoint"])
    
    # Compute inter-origin distance per chromosome
    df["inter_origin_distance"] = (
        df.groupby("chrom")["midpoint"]
        .diff()
    )
    
    return df

dfs = {
    "WT_1": compute_inter_origin_distance(wt_1),
    "WT_2": compute_inter_origin_distance(wt_2),
    "WT_3": compute_inter_origin_distance(wt_3),
}

# ---- Prepare long-format dataframe ----
records = []

for sample_name, df in dfs.items():
    clean = df.dropna(subset=["inter_origin_distance"])

    chrom_medians = (
        clean
        .groupby("chrom")["inter_origin_distance"]
        .median()
        .reset_index()
    )

    chrom_medians["sample"] = sample_name
    records.append(chrom_medians)

plot_df = pd.concat(records, ignore_index=True)

mpl.rcParams.update({'font.size': 10})

# ---- Plot ----
plt.figure(figsize=(8, 4))

sns.boxplot(
    data=plot_df,
    x="inter_origin_distance",
    y="sample",
    showfliers=False
)

sns.stripplot(
    data=plot_df,
    x="inter_origin_distance",
    y="sample",
    jitter=0.15,
    alpha=1,
    edgecolor="black",
    linewidth=0.5
)

plt.xlabel("Chromosome-wise median inter-IZ distance (kb)")
plt.ylabel("")
plt.xlim(0, 300_000)

plt.yticks(rotation=90)

plt.gca().xaxis.set_major_formatter(
    ticker.FuncFormatter(lambda x, pos: f"{int(x/1000)}")
)
plt.savefig("chromosome_wise_inter_origin_distance_boxplot_new_new.pdf", dpi=300)
plt.tight_layout()
plt.show()
