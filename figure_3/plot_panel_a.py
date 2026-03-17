import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from pathlib import Path
import re
from scipy.stats import pearsonr

DATA_DIR = Path("PATH/correlation_maps/averages")

files = list(DATA_DIR.glob("*.tsv"))

data = []

for f in files:
    name = f.stem
    
    # extract sample, type, resolution
    m = re.match(r"(.*)_(oem|rfd)_\d+_(\d+)", name)
    if not m:
        continue
        
    sample, signal_type, resolution = m.groups()
    
    df = pd.read_csv(f, sep="\t", header=None)
    
    values = df[5].values
    
    data.append({
        "sample": sample,
        "type": signal_type,
        "resolution": int(resolution),
        "values": values
    })

df_all = pd.DataFrame(data)

def significance_asterisk(p):
    if p < 0.001:
        return "***"
    elif p < 0.01:
        return "**"
    elif p < 0.05:
        return "*"
    else:
        return ""

def plot_combined_heatmap(df_all, resolution, outfile=None):
    df_res = df_all[df_all["resolution"] == resolution].copy()
    
    df_res["type"] = df_res["type"].str.strip().str.lower()
    df_res["sample"] = df_res["sample"].str.strip()
    
    pivot = df_res.pivot(index="sample", columns="type", values="values")
    
    samples = pivot.index.tolist()
    n_samples = len(samples)
    
    oem_matrix = np.vstack(pivot["oem"].values)
    rfd_matrix = np.vstack(pivot["rfd"].values)
    
    corr_oem = np.corrcoef(oem_matrix)
    corr_rfd = np.corrcoef(rfd_matrix)
    
    combined = np.triu(corr_oem) + np.tril(corr_rfd, -1)
    
    pvals = np.ones((n_samples, n_samples))
    for i in range(n_samples):
        for j in range(n_samples):
            if i != j:
                if i < j:
                    _, p = pearsonr(oem_matrix[i], oem_matrix[j])
                else:
                    _, p = pearsonr(rfd_matrix[i], rfd_matrix[j])
                pvals[i, j] = p

    mask = np.eye(n_samples, dtype=bool)
    
    plt.figure(figsize=(6,6))
    ax = sns.heatmap(
        combined,
        xticklabels=samples,
        yticklabels=samples,
        cmap="coolwarm",
        vmin=-1,
        vmax=1,
        square=True,
        mask=mask,
        cbar_kws={"label": "Pearson correlation"},
        annot=False
    )
    
    for i in range(n_samples):
        for j in range(n_samples):
            if not mask[i,j]:
                corr_val = combined[i,j]
                star = significance_asterisk(pvals[i,j])
                text = f"{corr_val:.2f}{star}"
                ax.text(j + 0.5, i + 0.5, text, ha="center", va="center", fontsize=10)
    
    plt.title(f"{resolution//1000} kb — Upper: OEM | Lower: RFD")
    plt.tight_layout()
    
    if outfile:
        plt.savefig(outfile, dpi=300)
    
    plt.show()
    plt.close()

for res in [50000, 75000, 100000, 150000]:
    plot_combined_heatmap(df_all, res, outfile=f"correlation_heatmap_{res}.pdf")
