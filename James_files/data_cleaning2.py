import pandas as pd
import numpy as np

df = pd.read_csv( 'C:/Users/Owner/Documents/VERUM/Network stuff/network-global-stats-latest-version.txt', index_col= 0 )

df = df.drop('Chromatic Number', axis = 1)

nan_graphs = ['scc-rt-lebanon', 'scc-rt-obama', 'tech-p2p', 'soc-LiveJournal1', 'soc-sinaweibo', 'soc-student-coop', 'soc-friendster', 'soc-twitter']
for graph in nan_graphs:
    df = df[df['Graph'] != graph]



print(df.head())
#for i in df.index:
 #   df.loc[i, 'Average triangles'] = df.loc[i, 'Total triangles'] / df.loc[i, 'Nodes']