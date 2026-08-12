from ete3 import Tree
import pandas as pd

"""Treat pastml results to obtain emergence node of each orthogroup"""

df = pd.read_csv("pastml_output_species.tab", sep="\t", index_col=0)  

tree = Tree("species_tree.nwk", format=1)  

for node in tree.traverse("preorder"):
    node.add_features(dist_from_root=node.get_distance(tree))

ssp_first_emergence = {}
for ssp in df.columns:
    present_nodes = df.index[df[ssp] == 1].tolist()
    if not present_nodes:
        continue
    node_dists = [(n, tree & n).dist_from_root for n in present_nodes if tree.search_nodes(name=n)]
    if node_dists:
        ssp_first_emergence[ssp] = min(node_dists, key=lambda x: x[1])[0]

pd.Series(ssp_first_emergence, name="first_emerged_at").to_csv("ssp_first_emergence.csv")

