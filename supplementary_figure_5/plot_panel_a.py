import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from pathlib import Path
import re

DATA_DIR = Path("PATH/scripts/06_revision/pu_seq/src/averages/files")

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

DATA_DIR = Path("PATH/scripts/04_plots_analyses/correlation_maps/averages")

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

df_all2 = pd.DataFrame(data)

df = pd.concat([df_all, df_all2], axis=0)

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import pearsonr

def significance_asterisk(p):
    if p < 0.001:
        return "***"
    elif p < 0.01:
        return "**"
    elif p < 0.05:
        return "*"
    else:
        return ""

def plot_combined_heatmap(df_all, resolution, sample_order=None, outfile=None):
    df_res = df_all[df_all["resolution"] == resolution].copy()

    df_res["type"] = df_res["type"].str.strip().str.lower()
    df_res["sample"] = df_res["sample"].str.strip()

    pivot = df_res.pivot(index="sample", columns="type", values="values")

    # Apply custom sample order if provided
    if sample_order is not None:
        # Keep only samples that actually exist
        sample_order = [s for s in sample_order if s in pivot.index]

        # Append any samples not specified in sample_order
        remaining = [s for s in pivot.index if s not in sample_order]

        pivot = pivot.reindex(sample_order + remaining)

    samples = pivot.index.tolist()
    n_samples = len(samples)

    oem_matrix = np.vstack(pivot["oem"].values)
    rfd_matrix = np.vstack(pivot["rfd"].values)

    corr_oem = np.corrcoef(oem_matrix)
    corr_rfd = np.corrcoef(rfd_matrix)

    # Upper triangle = OEM, lower triangle = RFD
    combined = np.triu(corr_oem) + np.tril(corr_rfd, -1)

    # Compute p-values
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

    plt.figure(figsize=(15, 15))
    ax = sns.heatmap(
        combined,
        xticklabels=samples,
        yticklabels=samples,
        cmap="coolwarm",
        vmin=-1,
        vmax=1,
        square=True,
        mask=mask,
        annot=False,
        linewidths=0.5,
    	linecolor="black"
    )

    # Add correlation values and significance stars
    for i in range(n_samples):
        for j in range(n_samples):
            if not mask[i, j]:
                corr_val = combined[i, j]
                star = significance_asterisk(pvals[i, j])
                ax.text(
                    j + 0.5,
                    i + 0.5,
                    f"{corr_val:.2f}",
                    ha="center",
                    va="center",
                    fontsize=10
                )

    plt.title(f"{resolution // 1000} kb — Upper: OEM | Lower: RFD")
    plt.tight_layout()

    if outfile:
        plt.savefig(outfile, dpi=300)

    plt.show()
    plt.close()

sample_order = [
    'alpha_rep1',
    'alpha_rep2',
    'alpha_rep3',
    'epsilon_rep1',
    'epsilon_rep2',
    'epsilon_rep3',
	'rep1_watson',
	'rep2_watson',
    'rep3_watson',
    'rep1_crick',
    'rep2_crick',
    'rep3_crick', 
    'imbulrich20240401wt',
    'sfbulrich20250101wt',
    'sfbulrich20250102wtrad21',
    ]

for res in [50000, 75000, 100000, 150000]:
    plot_combined_heatmap(df, res, sample_order=sample_order, outfile=f"correlation_heatmap_{res}.pdf")
