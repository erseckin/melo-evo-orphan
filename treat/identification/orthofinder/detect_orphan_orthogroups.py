import pandas as pd

"""
Filter Orthofinder GeneCount and Orthogroups files
to keep only Melo-specific orphan orthogroups
"""

# Input
orthofinder_genecount = pd.read_csv(
    "/path/to/Orthogroups.GeneCount.tsv",
    sep="\t",
    header=0,
    index_col=0,
    low_memory=False
)

orthofinder_orthogroups = pd.read_csv(
    "/path/to/Orthogroups.tsv",
    sep="\t",
    low_memory=False
)

# Filtering

melo_species = [
    "Megram", "Mechit", "Mehapl", "Meente",
    "Mearen", "Mejava", "Meinco", "Meluci"
]

all_species = orthofinder_genecount.columns.drop("Total")
non_melo_species = [sp for sp in all_species if sp not in melo_species]


melo_present = (orthofinder_genecount[melo_species] > 0).any(axis=1)
non_melo_absent = (orthofinder_genecount[non_melo_species] == 0).all(axis=1)

orphan_mask = melo_present & non_melo_absent

orphan_orthogroups_genecount = orthofinder_genecount[orphan_mask]
list_orphan_orthogroups = orphan_orthogroups_genecount.index.tolist()

orthofinder_orthogroups_contigs = orthofinder_orthogroups[
    orthofinder_orthogroups["Orthogroup"].isin(list_orphan_orthogroups)
]

orthofinder_orthogroups_contigs2 = orthofinder_orthogroups_contigs[
    ["Orthogroup"] + melo_species
]

# Output

orphan_orthogroups_genecount.to_csv(
    "Orphan.Orthogroups.Genecount.tsv",
    sep="\t"
)

orthofinder_orthogroups_contigs2.to_csv(
    "Orthogroups_orphans_with_id.tsv",
    sep="\t",
    index=False
)

print("Total orphan orthogroups:", orphan_orthogroups_genecount.shape[0])
print(
    "Total orphan sequences:",
    orphan_orthogroups_genecount[melo_species].sum().sum()
)

