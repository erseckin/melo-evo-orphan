import pandas as pd
from ete3 import NCBITaxa
ncbi = NCBITaxa()

"""Treat diamond lca results to identify genes that has an homolog outside of Meloidogyne"""

# Input

diamond_lca = pd.read_table('path/to/diamond_lca.tsv',header=None)

taxo_list = diamond_lca[1].tolist()
taxid2name = ncbi.get_taxid_translator(taxo_list)

# Filtering

diamond_lca['taxnames'] = diamond_lca[1].map(taxid2name)
diamond_lca[['taxnames']] = diamond_lca[['taxnames']].fillna('None') # None means diamond couldn't attribute any last common ancestor

diamond_lca_filtered = diamond_lca[~diamond_lca['taxnames'].str.contains("Meloidogyne|meloidogyne|None")]

# Output

diamond_lca_filtered.to_csv('diamond_lca_filtered.tsv', sep="\t", index=False, header=False) 

