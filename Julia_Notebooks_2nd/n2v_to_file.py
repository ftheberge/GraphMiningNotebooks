import sys
import fastnode2vec as n2v
import numpy as np

# python n2v_to_file.py edgelist_filepath dim p q seed twitch_flag

if len(sys.argv) == 6:
    data = np.loadtxt(sys.argv[1])
else:
    data = np.loadtxt(sys.argv[1],skiprows=1,delimiter=",")
graph = n2v.Graph(data, directed=False, weighted=False)
nv = n2v.Node2Vec(graph, dim=int(sys.argv[2]), p=np.float64(sys.argv[3]), q=np.float64(sys.argv[4]), walk_length=80, window=5, seed=int(sys.argv[5]))
nv.train(epochs=10, verbose=False)
Y = np.array([nv.wv[i] for i in range(len(nv.wv))])
with open('_embed','w') as f:
    nid=0
    f.write(str(Y.shape[0]) + " " + str(Y.shape[1])+'\n')
    for i in range(Y.shape[0]):
        f.write(str(nid)+' ')
        for j in range(Y.shape[1]):
            f.write(str(Y[i][j])+' ')
        f.write('\n')
        nid+=1