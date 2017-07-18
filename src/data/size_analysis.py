import pandas as pd

net = pd.read_csv('~/PycharmProjects/network_classification/src/data/clean_data_with_new_chem.csv', index_col=0)

net_big = net[net['Nodes']>=100000]
print('Information on Graphs with more than 100,000 Nodes:')3
print(net_big.describe())

net_med = net[net['Nodes']>=10000]
net_med = net_med[net_med['Nodes']<100000]
print('Information on Graphs with more than 10,000 Nodes and less than 100,000 Nodes:')
print(net_med.describe())

net_sm = net[net['Nodes']<10000]
print('Information on Graphs with less than 10,000 Nodes:')
print(net_sm.describe())