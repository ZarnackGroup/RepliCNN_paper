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
