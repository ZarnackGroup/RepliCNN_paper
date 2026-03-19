import pyBigWig
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os

def plot_raw_signal(chrom, start, end, binsize,
						bigwig_file, outfile="signal_binned.png"):
	"""
	Plot binned signal from a bigWig file over a given interval.
	Uses sum of values per bin.

	Parameters
	----------
	chrom : str
		Chromosome name (e.g. 'chr1').
	start : int
		Start coordinate.
	end : int
		End coordinate.
	binsize : int
		Bin size in base pairs.
	bigwig_file : str
		Path to a bigWig file.
	outfile : str
		Output filename for the plot.
	"""
	# Define bin edges
	bins = np.arange(start, end, binsize)
	bin_centers = bins + binsize // 2

	with pyBigWig.open(bigwig_file) as bw:
		values = bw.values(chrom, start, end, numpy=True)
	values = np.nan_to_num(values, nan=0.0)

	# Truncate to multiple of binsize
	usable_len = (end - start) // binsize * binsize
	values = values[:usable_len]
	reshaped = values.reshape(-1, binsize)
	binned = reshaped.sum(axis=1)

	# Plot
	plt.figure(figsize=(8, 2))
	plt.plot(bin_centers[:len(binned)], binned, linewidth=3)
	plt.fill_between(bin_centers[:len(binned)],binned,alpha=0.8)
	plt.ylabel(f"Raw signal")
	plt.title(f"{os.path.basename(bigwig_file)}; {chrom}:{start:,}-{end:,}; n={len(binned):,}")
	plt.xlim(start,end)
	plt.ylim(0,max(binned)*1.1)
	plt.tight_layout()
	if outfile: plt.savefig(outfile, dpi=150)
	plt.show()
	plt.close()

plot_raw_signal(
	"chrII",
	275000,
	450000,
	500,
	"../../../data/bigwigs/GSM5005342.fwd.bw",
	outfile="GSM5005342.fwd.pdf",
	)

plot_raw_signal(
	"chrII",
	275000,
	450000,
	500,
	"../../../data/bigwigs/GSM5005342.rev.bw",
	outfile="GSM5005342.rev.pdf",
	)

def plot_rfd_oem(chrom, start, end, kind, bigwig_file, color, alpha, outfile="signal_unbinned.png"):
    """
    Plot unbinned signal from a bigWig file over a given interval.

    Parameters
    ----------
    chrom : str
        Chromosome name (e.g. 'chr1').
    start : int
        Start coordinate.
    end : int
        End coordinate.
    kind : str
        Label for the y-axis (e.g. "OEM", "RFD").
    bigwig_file : str
        Path to a bigWig file.
    outfile : str
        Output filename for the plot.
    """
    with pyBigWig.open(bigwig_file) as bw:
        values = bw.values(chrom, start, end, numpy=True)
    values = np.nan_to_num(values, nan=0.0)

    x = np.arange(start, end)

    # Plot
    plt.figure(figsize=(8, 2))
    plt.plot(x, values, linewidth=2, c=color, alpha=alpha)
    plt.plot(x, [0] * len(x), linewidth=1, ls="-", c="black")
    plt.fill_between(x, values, alpha=alpha-0.2, color=color)
    plt.ylabel(kind)
    plt.title(f"{os.path.basename(bigwig_file)}; {chrom}:{start:,}-{end:,}")
    plt.xlim(start, end)
    plt.ylim(min(values)*1.1, max(values)*1.1)
    plt.tight_layout()
    if outfile:
        plt.savefig(outfile, dpi=150)
    plt.show()
    plt.close()

plot_rfd_oem(
	"chrII",
	275000,
	450000,
	"OEM",
	"../../../data/oem_rfd/GSM5005342/GSM5005342_oem_1_2500.bw",
	"C1",
	0.5,
	outfile="GSM5005342_oem_1_2500.pdf",
	)

plot_rfd_oem(
	"chrII",
	275000,
	450000,
	"OEM",
	"../../../data/oem_rfd/GSM5005342/GSM5005342_oem_1_5000.bw",
	"C1",
	0.5,
	outfile="GSM5005342_oem_1_5000.pdf",
	)

plot_rfd_oem(
	"chrII",
	275000,
	450000,
	"OEM",
	"../../../data/oem_rfd/GSM5005342/GSM5005342_oem_1_10000.bw",
	"C1",
	0.75,
	outfile="GSM5005342_oem_1_10000.pdf",
	)

plot_rfd_oem(
	"chrII",
	275000,
	450000,
	"OEM",
	"../../../data/oem_rfd/GSM5005342/GSM5005342_oem_1_15000.bw",
	"C1",
	1,
	outfile="GSM5005342_oem_1_15000.pdf",
	)

plot_rfd_oem(
	"chrII",
	275000,
	450000,
	"RFD",
	"../../../data/oem_rfd/GSM5005342/GSM5005342_rfd_1_2500.bw",
	"C2",
	0.5,
	outfile="GSM5005342_rfd_1_2500.pdf",
	)

plot_rfd_oem(
	"chrII",
	275000,
	450000,
	"RFD",
	"../../../data/oem_rfd/GSM5005342/GSM5005342_rfd_1_5000.bw",
	"C2",
	0.5,
	outfile="GSM5005342_rfd_1_5000.pdf",
	)

plot_rfd_oem(
	"chrII",
	275000,
	450000,
	"RFD",
	"../../../data/oem_rfd/GSM5005342/GSM5005342_rfd_1_10000.bw",
	"C2",
	0.75,
	outfile="GSM5005342_rfd_1_10000.pdf",
	)

plot_rfd_oem(
	"chrII",
	275000,
	450000,
	"RFD",
	"../../../data/oem_rfd/GSM5005342/GSM5005342_rfd_1_15000.bw",
	"C2",
	1,
	outfile="GSM5005342_rfd_1_15000.pdf",
	)

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

def plot_rfd_oem_bedgraph(chrom, start, end,
                          df, color="C1", alpha=0.5,
                          outfile="signal_unbinned.png", label=None):
    """
    Plot unbinned signal from a single bedGraph-like DataFrame over a given interval.

    Parameters
    ----------
    chrom : str
        Chromosome name (e.g. 'chr1').
    start : int
        Start coordinate.
    end : int
        End coordinate.
    df : pandas.DataFrame
        BedGraph-like dataframe with ['chrom', 'start', 'end', 'value'].
    color : str
        Color for plotting.
    alpha : float
        Transparency for fill.
    outfile : str
        Output filename for the plot.
    label : str
        Label for y-axis.
    """
    x = np.arange(start, end)

    # Initialize signal array
    values = np.zeros(end - start, dtype=float)

    # Filter relevant rows
    sub = df[(df["chrom"] == chrom) & (df["end"] > start) & (df["start"] < end)]
    for _, row in sub.iterrows():
        s = max(start, row["start"])
        e = min(end, row["end"])
        values[s-start:e-start] = row["value"]

    # Plot
    plt.figure(figsize=(8, 2))
    plt.plot(x, values, linewidth=5, c=color)
    plt.ylabel(label if label else "signal")
    plt.title(f"{chrom}:{start:,}-{end:,}")
    plt.xlim(start, end)
    #plt.ylim(-1,1)
    plt.tight_layout()

    if outfile:
        plt.savefig(outfile, dpi=150)
    plt.show()
    plt.close()

df1 = pd.read_csv("../../../data/cross_prediction/models/yeast_traelseq_GSM5005342_delta_chrII/GSM5005342_pred.tsv",
                  sep="\t", names=["chrom","start","end","pos","neg","log2","spline","deri","antideri","value"])

df2 = pd.read_csv("../../../data/sdfs/GSM5005342.tsv",
                  sep="\t", names=["chrom","start","end","pos","neg","log2","spline","deri","antideri","value"])

plot_rfd_oem_bedgraph(
    "chrII",
    275000,
    450000,
    df1,
    color="darkred",
    alpha=0.75,
    outfile="GSM5005342_predicted.pdf",
    label="Predicted RT"
)

plot_rfd_oem_bedgraph(
    "chrII",
    275000,
    450000,
    df2,
    color="gold",
    alpha=0.75,
    outfile="GSM5005342_truth.pdf",
    label="Ground Truth RT"
)

df1 = pd.read_csv("../../../data/sdfs/GSM5005342.tsv",
                  sep="\t", names=["chrom","start","end","pos","neg","value","spline","deri","antideri","time"])
df2 = pd.read_csv("../../../data/sdfs/GSM5005342.tsv",
                  sep="\t", names=["chrom","start","end","pos","neg","log2","spline","value","antideri","time"])
df3 = pd.read_csv("../../../data/sdfs/GSM5005342.tsv",
                  sep="\t", names=["chrom","start","end","pos","neg","log2","spline","deri","value","time"])

plot_rfd_oem_bedgraph(
    "chrII",
    275000,
    450000,
    df1,
    color="darkred",
    alpha=1,
    outfile="GSM5005342_spline.pdf",
    label="Spline"
)

plot_rfd_oem_bedgraph(
    "chrII",
    275000,
    450000,
    df2,
    color="darkred",
    alpha=1,
    outfile="GSM5005342_deri.pdf",
    label="Derivative"
)

plot_rfd_oem_bedgraph(
    "chrII",
    275000,
    450000,
    df3,
    color="darkred",
    alpha=1,
    outfile="GSM5005342_antideri.pdf",
    label="Antiderivative"
)

def plot_bed_intervals(chrom, start, end,
                       df, outfile="intervals.png", label="intervals"):
    """
    Plot intervals from a BED-like DataFrame as rectangles (blocks).
    - Color comes from 'color' column (r,g,b format).
    - Alpha comes from 'score' column (scaled).
    - '+' and '-' strand intervals plotted in separate horizontal tracks.
    """
    # Filter relevant rows
    sub = df[(df["chrom"] == chrom) & (df["end"] > start) & (df["start"] < end)]

    if sub.empty:
        print("No intervals in range.")
        return

    # Normalize score to [0,1] for alpha
    scores = sub["score"].astype(float)
    if scores.max() > 0:
        alphas = scores / scores.max()
    else:
        alphas = [0.5] * len(sub)
    sub = sub.assign(alpha=alphas)

    plt.figure(figsize=(8, 2))

    for (_, row), a in zip(sub.iterrows(), sub["alpha"]):
        s = max(start, row["start"])
        e = min(end, row["end"])

        # Parse color (BED uses r,g,b)
        c = "gray"
        try:
            r, g, b = [int(x) for x in str(row["color"]).split(",")]
            c = (r/255, g/255, b/255)
        except Exception:
            pass

        # Y-position depends on strand
        if row["strand"] == "+":
            y0, y1 = 1, 2
        else:
            y0, y1 = -1, 0

        plt.fill_between([s, e], y0, y1, color=c, alpha=a)

    # Decorations
    plt.ylabel(label)
    plt.title(f"{chrom}:{start:,}-{end:,}")
    plt.xlim(start, end)
    plt.ylim(-1.5, 2.5)
    plt.yticks([1.5, -0.5], ["+", "-"])
    plt.tight_layout()

    if outfile:
        plt.savefig(outfile, dpi=150)
    plt.show()
    plt.close()

df1 = pd.read_csv(
    "../../../data/oem_rfd/GSM5005342/GSM5005342_oris_ters.bed",
    sep="\t",
    names=["chrom","start","end","name","score","strand",
           "thickStart","thickEnd","color"]
)

plot_bed_intervals(
    "chrII",
    275000,
    450000,
    df1,
    outfile="oris_ters_chrII.pdf",
    label="ORIs"
)

plot_bed_intervals(
    "chrII",
    275000,
    450000,
    df1,
    outfile="oris_ters_chrII.pdf",
    label="ORIs"
)

def plot_bed_boxes(chrom, start, end, df,
                   outfile="boxes.png", label="intervals"):
    """
    Plot intervals from a BED-like DataFrame as grey rectangles.

    Parameters
    ----------
    chrom : str
        Chromosome name (e.g. 'chr1').
    start : int
        Start coordinate.
    end : int
        End coordinate.
    df : pandas.DataFrame
        Dataframe with columns ['chrom', 'start', 'end', 'name'].
    outfile : str
        Output filename for the plot.
    label : str
        Label for y-axis.
    """
    # Filter relevant rows
    sub = df[(df["chrom"] == chrom) & (df["end"] > start) & (df["start"] < end)]

    plt.figure(figsize=(8, 1.5))

    for _, row in sub.iterrows():
        s = max(start, row["start"])
        e = min(end, row["end"])
        plt.fill_between([s, e], 0, 1, color="grey", alpha=0.8)

    plt.ylabel(label)
    plt.title(f"{chrom}:{start:,}-{end:,}")
    plt.xlim(start, end)
    plt.ylim(0-0.2, 1.2)
    plt.yticks([])
    plt.tight_layout()

    if outfile:
        plt.savefig(outfile, dpi=150)
    plt.show()
    plt.close()

df1 = pd.read_csv(
    "/home/dos02bi/koenig_data/projects/storage/data/__misc/yeast_oris/remap/oris_confirmed_saccer3.bed4",
    sep="\t",
    names=["chrom","start","end","name"]
)

plot_bed_boxes(
    "chrII",
    275000,
    450000,
    df1,
    outfile="chrII_boxes.pdf",
    label="OriDB confirmed"
)
