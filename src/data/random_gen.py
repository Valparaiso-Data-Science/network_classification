from igraph import *
import numpy as np

er = Graph()
b = Graph()

p_list = [0.01, 0.05, 0.1]
m_list = [10, 40, 60]

nGraphs = range(1, 26)
#generate all the erdos_renyi and the barabasi graphs- 3 loops each generating 25 graphs
#loop 1: n=1000
print('Lets make some graphs!')
for n1 in nGraphs:
    n_1 = 1000 + np.random.randint(-200, 200)
    p = np.random.choice(p_list)
    p_str = str(int(p * 100))
    er = er.Erdos_Renyi(n_1, p)
    #write er to file
    er.write_edgelist('er_1_'+str(n1)+ '_'+p_str+'.csv')
    m = int(np.random.choice(m_list))
    b = b.Barabasi(n_1, m)
    m_str = str(m)
    #write b to file
    b.write_edgelist('b_1_'+str(n1)+ '_'+m_str+'.csv')
    if n1 == 13:
        print('halfway done with n1!')
#loop 2 : n=10000
print('working on n2! no sweat!')
for n2 in nGraphs:
    n_2 = 10000 + np.random.randint(-2000, 2000)
    p = np.random.choice(p_list)
    p_str = str(int(p * 100))
    er = er.Erdos_Renyi(n_2, p=p)
    #write er to file
    er.write_edgelist('er_2_'+str(n2)+ '_'+p_str+'.csv')
    m = int(np.random.choice(m_list))
    b = b.Barabasi(n_2, m)
    m_str = str(m)
    # write b to file
    b.write_edgelist('b_2_'+str(n2)+'_'+m_str+'.csv')
    if n2 == 13:
        print('halfway done with n2!')
print('relax! done with n2!!')
#loop 3: n=100000
for n3 in nGraphs:
    n_3 = 100000 + np.random.randint(-20000, 20000)
    p = np.random.choice(p_list)
    p_str = str(int(p * 100))
    er = er.Erdos_Renyi(n_3, p=p)
    #write er to file
    er.write_edgelist('er_3_'+str(n3)+ '_'+p_str+'.csv')
    m = int(np.random.choice(m_list))
    b = b.Barabasi(n_3, m)
    m_str = str(m)
    # write b to file
    b.write_edgelist('b_3_'+str(n3)+'_'+m_str+'.csv')
    if n3 == 13:
        print('halfway done with n3! staying alive!')
print('done!!')