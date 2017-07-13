from igraph import *
import numpy as np

er = Graph()
b = Graph()

p_list = [0.1, 0.5, 0.9]
m_list = [10, 50, 90]
#er.write_edgelist('/Volumes/Seagate Backup Plus Drive/VERUM/Network_Repository_files/edgelist1.csv')

nGraphs = range(1, 25)

#generate all the erdos_renyi and the barabasi graphs- 3 loops each generating 25 graphs
#loop 1: n=100
for n1 in nGraphs:
    n_1 = 100 + np.random.randint(-20, 20)
    p = np.random.choice(p_list)
    er = er.Erdos_Renyi(n_1, p)
    #write er to file
    m = np.random.choice(m_list)
    b = b.Barabasi(n_1, m)
    #write b to file
#loop 2 : n=1000
for n2 in nGraphs:
    n_2 = 1000 + np.random.randint(-200, 200)
    p = np.random.choice(p_list)
    er = er.Erdos_Renyi(n_2, p=p)
    #write er to file
    m = np.random.choice(m_list)
    b = b.Barabasi(n_2, m)
    # write b to file
#loop 3: n=10000
for n3 in nGraphs:
    n_3 = 10000 + np.random.randint(-2000, 2000)
    p = np.random.choice(p_list)
    er = er.Erdos_Renyi(n_3, p=p)
    #write er to file
    m = np.random.choice(m_list)
    b = b.Barabasi(n_3, m)
    # write b to file
