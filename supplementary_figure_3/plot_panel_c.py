import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from io import StringIO

data = """Phase	DMSO	HU	std dev dmso	std dev HU
G1	29,7	22,45	2,678307924	1,946792233
G2	10,525	13,075	0,644851404	0,95
S	59,725	64,475	2,954516317	2,745147233
"""

# Import data
df = pd.read_csv(
    StringIO(data),
    sep="\t",
    decimal=",",
)

# Order cell-cycle phases conventionally
phase_order = ["G1", "S", "G2"]
df["Phase"] = pd.Categorical(
    df["Phase"],
    categories=phase_order,
    ordered=True,
)

# Convert to long format
plot_df = pd.DataFrame({
    "Phase": list(df["Phase"]) * 2,
    "Treatment": ["DMSO"] * len(df) + ["HU"] * len(df),
    "Percentage": [
        *df["DMSO"],
        *df["HU"],
    ],
    "SD": [
        *df["std dev dmso"],
        *df["std dev HU"],
    ],
})

plot_df["Phase"] = pd.Categorical(
    plot_df["Phase"],
    categories=phase_order,
    ordered=True,
)

# Plot
fig, ax = plt.subplots(figsize=(5, 4))

sns.barplot(
    data=plot_df,
    x="Phase",
    y="Percentage",
    hue="Treatment",
    order=phase_order,
    errorbar=None,
    ax=ax,
)

# Add standard deviation error bars
for phase_idx, phase in enumerate(phase_order):

    for condition_idx, condition in enumerate(["DMSO", "HU"]):

        row = plot_df[
            (plot_df["Phase"] == phase)
            & (plot_df["Treatment"] == condition)
        ].iloc[0]

        # Seaborn places two bars approximately ±0.2 from the center
        x = phase_idx + (-0.2 if condition == "DMSO" else 0.2)

        ax.errorbar(
            x,
            row["Percentage"],
            yerr=row["SD"],
            fmt="none",
            ecolor="black",
            capsize=3,
            linewidth=1,
        )

ax.set_xlabel("Cell cycle phase")
ax.set_ylabel("Cells (%)")
ax.set_ylim(0, 100)

plt.tight_layout()
plt.savefig("cell_cycle_barplot.pdf", dpi=300)
plt.show()
