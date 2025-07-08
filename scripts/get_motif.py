import pandas as pd
from glycowork.motif.annotate import annotate_dataset
from glycowork.glycan_data.loader import df_glycan, df_species


def count_motifs_per_species(df_glycan: pd.DataFrame,
                             condense: bool = True) -> pd.DataFrame:
    df_sp = df_glycan["Species"].apply(lambda x: str(x).split(",")[0].strip())
    df_gly = df_glycan["glycan"]
    df_motif = annotate_dataset(df_gly.tolist(),
                                condense=condense)
    df_motif.index = df_sp.values
    df_sp_motif_counts = df_motif.groupby(df_motif.index
                                          ).sum()
    df_sp_motif_counts.to_csv("/output/motif_counts_per_species.csv")
    return df_sp_motif_counts

df_sp_motif_counts = count_motifs_per_species(df_glycan)
print(df_sp_motif_counts.head())

