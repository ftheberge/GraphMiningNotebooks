import h_louvain as hl
import hypernetx as hnx
import hypernetx.algorithms.hypergraph_modularity as hmod
import sys
import json

# python run_hlouvain.py edgelist_filepath dim p q seed twitch_flag

fp = open(sys.argv[1], 'r')
Lines = fp.readlines()
Edges = []
for line in Lines:
    Edges.append(set([x for x in line.strip().split(',')]))
H = hnx.Hypergraph(dict(enumerate(Edges)))


hBO = hl.hLouvainBO(H,hmod_tau="infinity",random_seed=int(sys.argv[2]))
hBO.set_params(show_bo_table=False)
result_df = hBO.hLouvain_perform_BO()
dct = hmod.part2dict(result_df["A_lstep"].tolist()[0])
with open('hlouvain_results.json', 'w') as f:
    json.dump(dct, f)