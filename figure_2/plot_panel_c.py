import upsetplot
from collections import Counter
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pybedtools
from scipy.stats import pearsonr

ids = [
    ("yeast", "traelseq", "GSM4680452"),
    ("yeast", "traelseq", "GSM4680453"),
    ("yeast", "traelseq", "GSM4680454"),
    ("yeast", "traelseq", "GSM4680457"),
    ("yeast", "traelseq", "GSM4680459"),
    ("yeast", "traelseq", "GSM4680460"),
    ("yeast", "traelseq", "GSM4680461"),
    ("yeast", "traelseq", "GSM4680462"),
    ("yeast", "traelseq", "GSM4680463"),
    ("yeast", "traelseq", "GSM4680464"),
    ("yeast", "traelseq", "GSM4680465"),
    ("yeast", "traelseq", "GSM4680468"),
    ("yeast", "traelseq", "GSM4680469"),
    ("yeast", "traelseq", "GSM5005318"),
    ("yeast", "traelseq", "GSM5005319"),
    ("yeast", "traelseq", "GSM5005322"),
    ("yeast", "traelseq", "GSM5005327"),
    ("yeast", "traelseq", "GSM5005328"),
    ("yeast", "traelseq", "GSM5005329"),
    ("yeast", "traelseq", "GSM5005334"),
    ("yeast", "traelseq", "GSM5005335"),
    ("yeast", "traelseq", "GSM5005336"),
    ("yeast", "traelseq", "GSM5005337"),
    ("yeast", "traelseq", "GSM5005338"),
    ("yeast", "traelseq", "GSM5005339"),
    ("yeast", "traelseq", "GSM5005340"),
    ("yeast", "traelseq", "GSM5005341"),
    ("yeast", "traelseq", "GSM5005342"),
    ("yeast", "traelseq", "GSM5005343"),
    ("yeast", "gloeseq", "imbulrich201807s16"),
    ("yeast", "gloeseq", "imbulrich201807s18"),
    ("yeast", "gloeseq", "imbulrich201807s20"),
    ("yeast", "okseq", "GSM835650"),
    ("yeast", "okseq", "GSM835651"),
    ("yeast", "okseq", "GSM835652"),
    ("yeast", "okseq", "GSM835653"),
]

# Read in ORIs and TERs
oris = pd.DataFrame()
ters = pd.DataFrame()

for organism, experiment, sample in ids:
    file = f"../../../data/oem_rfd/{sample}/{sample}_oris_ters.bed"
    tmp = pd.read_csv(
        file,
        sep="\t",
        names=["chrom", "start", "end", "name", "score", "strand",
               "thickStart", "thickEnd", "itemRgb"]
    )
    
    # Split ORIs and TERs, add metadata
    oris = pd.concat([
        oris,
        tmp.query("strand == '+'")
           .assign(sample=sample, organism=organism, experiment=experiment)
    ]).reset_index(drop=True)

    ters = pd.concat([
        ters,
        tmp.query("strand == '-'")
           .assign(sample=sample, organism=organism, experiment=experiment)
    ]).reset_index(drop=True)

oris_gloeseq = oris.query("experiment=='gloeseq'")[["chrom","start","end","sample"]].reset_index(drop=True)
oris_traelseq = oris.query("experiment=='traelseq'")[["chrom","start","end","sample"]].reset_index(drop=True)
oris_okseq = oris.query("experiment=='okseq'")[["chrom","start","end","sample"]].reset_index(drop=True)

oridb_confirmed = pd.read_csv("PATH/yeast_oris/remap/oris_confirmed_saccer3.bed4",sep="\t",names=["chrom","start","end","name"])
oridb_likely = pd.read_csv("PATH/yeast_oris/remap/oris_likely_saccer3.bed4",sep="\t",names=["chrom","start","end","name"])
oridb_dubious = pd.read_csv("PATH/yeast_oris/remap/oris_dubious_saccer3.bed4",sep="\t",names=["chrom","start","end","name"])

oridb_confirmed = oridb_confirmed.assign(name=oridb_confirmed["chrom"].astype(str)+":"+oridb_confirmed["start"].astype(str)+"-"+oridb_confirmed["end"].astype(str)+"_"+oridb_confirmed["name"].astype(str)+"_confirmed")
oridb_likely = oridb_likely.assign(name=oridb_likely["chrom"].astype(str)+":"+oridb_likely["start"].astype(str)+"-"+oridb_likely["end"].astype(str)+"_"+oridb_likely["name"].astype(str)+"_likely")
oridb_dubious = oridb_dubious.assign(name=oridb_dubious["chrom"].astype(str)+":"+oridb_dubious["start"].astype(str)+"-"+oridb_dubious["end"].astype(str)+"_"+oridb_dubious["name"].astype(str)+"_dubious")

assign_dist = 500

for df in [oris_gloeseq, oris_okseq, oris_traelseq]:
    for i, row in df.iterrows():
        chrom = row["chrom"]
        start = row["start"]
        end = row["end"]

        def get_overlap(db, chrom, start, end, dist):
            db_temp = db.query("chrom == @chrom")
            if db_temp.empty:
                return None

            overlap = db_temp[
                (db_temp["start"] <= end + dist) &
                (db_temp["end"] >= start - dist)
            ]

            if not overlap.empty:
                return overlap.iloc[0]["name"]
            return None

        # check in order of confidence
        if (name := get_overlap(oridb_confirmed, chrom, start, end, assign_dist)):
            kind = "confirmed"
        elif (name := get_overlap(oridb_likely, chrom, start, end, assign_dist)):
            kind = "likely"
        elif (name := get_overlap(oridb_dubious, chrom, start, end, assign_dist)):
            kind = "dubious"
        else:
            name = "unassigned"
            kind = "unassigned"

        df.at[i, "name"] = name
        df.at[i, "kind"] = kind

print(len(oris_traelseq.query("name=='unassigned'")))
print(len(oris_okseq.query("name=='unassigned'")))
print(len(oris_gloeseq.query("name=='unassigned'")))

min_reps = {
    "gloeseq": 2,
    "okseq": 3,
    "traelseq": 20
}

results = {}
for name, df in {
    "gloeseq": oris_gloeseq,
    "okseq": oris_okseq,
    "traelseq": oris_traelseq
}.items():
    tmp = (
        pybedtools.BedTool.from_dataframe(df)
        .sort()
        .merge(c=[4, 5, 6], o=["collapse", "collapse", "collapse"])
        .to_dataframe()
    )

    tmp = tmp.assign(
        n_replicates=tmp["name"].str.split(",").apply(len)
    ).query(f"n_replicates >= {min_reps[name]}")

    tmp = tmp.assign(score_list=tmp["score"].str.split(","))
    tmp = tmp.assign(
        majority_name=tmp["score_list"].apply(
            lambda vals: Counter(vals).most_common(1)[0][0]
        )
    )
    tmp = tmp.assign(strand_list=tmp["strand"].str.split(","))
    tmp = tmp.assign(
        majority_kind=tmp["strand_list"].apply(
            lambda vals: Counter(vals).most_common(1)[0][0]
        )
    )

    tmp = tmp.drop(columns=["score_list", "strand_list"])

    results[name] = tmp

annotation=pd.concat([oridb_confirmed.name,oridb_likely.name,oridb_dubious.name])

okseq_names = results["okseq"][["majority_name"]].copy()
mask = okseq_names["majority_name"] == "unassigned"
okseq_names.loc[mask, "majority_name"] = [f"unassigned_okseq_{i+1}" for i in range(mask.sum())]

traelseq_names = results["traelseq"][["majority_name"]].copy()
mask = traelseq_names["majority_name"] == "unassigned"
traelseq_names.loc[mask, "majority_name"] = [f"unassigned_traelseq_{i+1}" for i in range(mask.sum())]

gloeseq_names = results["gloeseq"][["majority_name"]].copy()
mask = gloeseq_names["majority_name"] == "unassigned"
gloeseq_names.loc[mask, "majority_name"] = [f"unassigned_gloeseq_{i+1}" for i in range(mask.sum())]

upset_sets = from_contents({
	"traelseq": traelseq_names.majority_name.drop_duplicates(keep="first"),
	"okseq": okseq_names.majority_name.drop_duplicates(keep="first"),
	"gloeseq": gloeseq_names.majority_name.drop_duplicates(keep="first"),
	"annotation":annotation
})

df = upset_sets.assign(annotation_type=np.where(
	upset_sets.id.isin(oridb_confirmed.name),
	"Confirmed",
	np.where(
		upset_sets.id.isin(oridb_likely.name),
		"Likely",
		np.where(
			upset_sets.id.isin(oridb_dubious.name),
			"Dubious",
			"Not annotated"
		)
	)
)
)

upset = UpSet(df, intersection_plot_elements=0,show_counts="{:,}")
upset.add_stacked_bars(by="annotation_type", title="Annotation", elements=10)
upset.plot()
plt.suptitle("")
plt.savefig("upset_plot.pdf")
plt.show()
