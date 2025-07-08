from glycowork.glycan_data.loader import df_glycan, df_species
from ete3 import NCBITaxa
ncbi = NCBITaxa()


#df_sp = df_glycan["Species"][:100]
df_lineage = df_glycan["Species"][:100]
raw = df_lineage.iloc[0]
species_list = raw
species_list = [s.strip() for s in species_list]
species_list = [s.replace('_', ' ') for s in species_list]


taxids = []
for sp in species_list:
    name2taxid = ncbi.get_name_translator([sp])
    if sp in name2taxid:
        taxids.append(name2taxid[sp][0])
    else:
        print(f"Warning: Taxid not found for {sp}")


tree = ncbi.get_topology(taxids)
print(tree.get_ascii(show_internal=True))
tree.write(outfile="output/species_tree.nw")
