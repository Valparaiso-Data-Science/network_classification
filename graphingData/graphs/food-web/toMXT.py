import os
import numpy as np
import networkx as nx

files = [f for f in os.listdir('.') if os.path.isfile(f)]
print(files)
for f in files:
    if f[-4:] == ".txt":
        mat = np.loadtxt(open(f, "rb"), delimiter = "\t").astype(int)
        #TL
        q2 = np.zeros((mat.shape[0], mat.shape[0]), dtype = np.int)
        #BR
        q4 = np.zeros((mat.shape[1], mat.shape[1]), dtype = np.int)
        top = np.concatenate((q2, mat), axis = 1)
        bottom = np.concatenate((mat.T, q4), axis = 1)
        adj = np.concatenate((top, bottom), axis = 0)

        G = nx.from_numpy_matrix(adj)
        nx.write_weighted_edgelist(G, "edges/" + f[:-4] + "_edges.txt")
