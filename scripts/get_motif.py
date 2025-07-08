from glycowork.motif.annotate import annotate_glycan, annotate_dataset, quantify_motifs
from glycowork.glycan_data.loader import df_glycan, df_species
import pandas as pd


def count_motifs_per_species(df_glycan: pd.DataFrame, condense: bool = True) -> pd.DataFrame:

    #df_sp = df_glycan["Species"].apply(lambda x: str(x).split(",")[0].strip())
    df_sp = df_species["Species"][:100]
    df_gly = df_glycan["glycan"][:100]

    df_motif = annotate_dataset(df_gly.tolist(), condense=condense)
    df_motif.index = df_sp.values
    df_sp_motif_counts = df_motif.groupby(df_motif.index).sum()

    return df_sp_motif_counts


df_sp_motif_counts = count_motifs_per_species(df_glycan, condense=True)
df_sp_motif_counts.to_csv("output/motif_counts_per_species.csv", index=True, header=True)
df_sp_motif_present= (df_sp_motif_counts > 0).astype(int)
df_sp_motif_present.to_csv("output/motif_presence_absence_per_species.csv", index=True, header=True)

print(df_sp_motif_counts.head())

