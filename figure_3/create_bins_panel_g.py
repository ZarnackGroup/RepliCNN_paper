import pandas as pd
from pathlib import Path

names = ["chrom", "start", "end", "name", "score", "strand",
         "thickStart", "thickEnd", "itemRgb"]

# input files
datasets = {
    "imbulrich20240401wt_oris": "PATH/scripts/04_plots_analyses/overlap_ori_term/rerun_rfd_oem_ori_ter_hct116/ori_ter/imbulrich20240401wt_oris.bed",
    "sfbulrich20250101wt_oris": "PATH/scripts/04_plots_analyses/overlap_ori_term/rerun_rfd_oem_ori_ter_hct116/ori_ter/sfbulrich20250101wt_oris.bed",
    "sfbulrich20250102wtrad21_oris": "PATH/scripts/04_plots_analyses/overlap_ori_term/rerun_rfd_oem_ori_ter_hct116/ori_ter/sfbulrich20250102wtrad21_oris.bed",

    "imbulrich20240401wt_ters": "PATH/scripts/04_plots_analyses/overlap_ori_term/rerun_rfd_oem_ori_ter_hct116/ori_ter/imbulrich20240401wt_terms.bed",
    "sfbulrich20250101wt_ters": "PATH/scripts/04_plots_analyses/overlap_ori_term/rerun_rfd_oem_ori_ter_hct116/ori_ter/sfbulrich20250101wt_terms.bed",
    "sfbulrich20250102wtrad21_ters": "PATH/scripts/04_plots_analyses/overlap_ori_term/rerun_rfd_oem_ori_ter_hct116/ori_ter/sfbulrich20250102wtrad21_terms.bed",
}

# quantile labels
labels = ["0-20", "20-40", "40-60", "60-80", "80-100"]

outdir = Path("score_bins")
outdir.mkdir(exist_ok=True)

for name, path in datasets.items():
    df = pd.read_csv(path, sep="\t", names=names)

    if name.endswith("ters"): 
        df["strand"] = "+"

    # quantile binning (robust to ties)
    df["score_bin"] = pd.qcut(
        df["score"],
        q=5,
        labels=labels,
        duplicates="drop"
    )

    # write each bin separately
    for label in labels:
        subset = df[df["score_bin"] == label]

        if subset.empty:
            continue

        outfile = outdir / f"{name}_score_{label}.bed"

        # BED output: drop helper column, no header
        subset.drop(columns="score_bin").to_csv(
            outfile,
            sep="\t",
            header=False,
            index=False
        )
