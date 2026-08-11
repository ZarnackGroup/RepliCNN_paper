import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from matplotlib_venn import venn3

oris = pd.read_csv("./ori_ter/merged_oris.bed", sep="\t", header=None, names=["chrom", "start", "end", "name"])
ters = pd.read_csv("./ori_ter/merged_ters.bed", sep="\t", header=None, names=["chrom", "start", "end", "name"])

wt_1 = "imbulrich20240401wt"
wt_2 = "sfbulrich20250101wt"
wt_3 = "sfbulrich20250102wtrad21"

oris = oris.assign(wt_1= oris["name"].str.contains(wt_1)).assign(wt_2= oris["name"].str.contains(wt_2)).assign(wt_3= oris["name"].str.contains(wt_3))[["wt_1", "wt_2", "wt_3"]]
ters = ters.assign(wt_1= ters["name"].str.contains(wt_1)).assign(wt_2= ters["name"].str.contains(wt_2)).assign(wt_3= ters["name"].str.contains(wt_3))[["wt_1", "wt_2", "wt_3"]]

def plot_venn_from_boolean_df(df, cols=("wt_1", "wt_2", "wt_3"), outfile=None):
    
    A, B, C = cols
    
    # Compute set sizes
    only_A = ((df[A]) & (~df[B]) & (~df[C])).sum()
    only_B = ((~df[A]) & (df[B]) & (~df[C])).sum()
    only_C = ((~df[A]) & (~df[B]) & (df[C])).sum()
    
    A_B = ((df[A]) & (df[B]) & (~df[C])).sum()
    A_C = ((df[A]) & (~df[B]) & (df[C])).sum()
    B_C = ((~df[A]) & (df[B]) & (df[C])).sum()
    
    A_B_C = ((df[A]) & (df[B]) & (df[C])).sum()
    
    plt.figure(figsize=(6,6))
    
    venn3(
        subsets=(only_A, only_B, A_B, only_C, A_C, B_C, A_B_C),
        set_labels=cols
    )
    
    plt.title("Overlap between WT samples")
    
    if outfile:
        plt.savefig(outfile, dpi=300)
    
    plt.show()
    
plot_venn_from_boolean_df(oris, outfile="wt_ori_venn_new.pdf")

plot_venn_from_boolean_df(ters, outfile="wt_ter_venn_new.pdf")
