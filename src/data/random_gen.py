from igraph import *

er = Graph()
b = Graph()

er = er.Erdos_Renyi(n=10, m=7, directed=False, loops=False)
b = b.Barabasi(n=10, m=7)
print(er)
print(b)

#er.write_edgelist('edgelist1.csv')