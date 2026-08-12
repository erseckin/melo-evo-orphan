import numpy as np
import pandas as pd
import re


"""Filter Orthofinder genecount and Orthogroups files to keep only the orphan orthogroups"""

# Input files 

orthofinder_genecount = pd.read_csv('/path/to/ortholog_counts_per_species.stats.tsv',sep='\t',header=[0],index_col=[0])
orthofinder_orthogroups = pd.read_csv('/path/to/ortholog_groups.tsv',sep='\t')

# Filtering

orthofinder_genecount2 = orthofinder_genecount.drop(['Mechit.fa','Meente.fa','Megram.fa','Mehapl.fa','total'], axis=1)
orthofinder_genecount3 = orthofinder_genecount2.loc[(orthofinder_genecount2==0).all(axis=1)]
list_orphan_orthogroups = orthofinder_genecount3.index.tolist()
orphan_orthogroups_genecount = orthofinder_genecount[orthofinder_genecount.index.isin(list_orphan_orthogroups)]

print(orphan_orthogroups_genecount.shape[0])
print(orphan_orthogroups_genecount['total'].sum())

list_orphan_groups = orphan_orthogroups_genecount.index.tolist()
orthofinder_orthogroups_contigs = orthofinder_orthogroups[orthofinder_orthogroups['group_id'].isin(list_orphan_groups)]

# Output 

orphan_orthogroups_genecount.to_csv('orphan_groups_genecount.tsv', sep="\t", index=True, header=True)
orthofinder_orthogroups_contigs.to_csv('orphan_groups_contigs.tsv', sep="\t", index=False, header=True)

