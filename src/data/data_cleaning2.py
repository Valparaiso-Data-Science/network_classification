import pandas as pd
import numpy as np

df = pd.read_csv( 'C:/Users/Owner/Documents/VERUM/Network stuff/network-global-stats-latest-version.txt')

# The attributes from the text file are read as objects, this changes them into numbers
column_list = [ 'Nodes', 'Edges', 'Density', 'Maximum degree', 'Minimum degree', 'Average degree', 'Assortativity', 'Total triangles', 'Average triangles', 'Maximum triangles', 'Avg. clustering coef.', 'Frac. closed triangles', 'Maximum k-core', 'Max. clique (lb)' ]
df[column_list] = df[column_list].apply(pd.to_numeric, errors = 'coerce')

#removes chromatic number
df = df.drop('Chromatic Number', axis = 1)

#removes graphs with NaNs that are too big to work on
nan_graphs = ['scc-rt-lebanon', 'scc-rt-obama', 'tech-p2p', 'soc-LiveJournal1', 'soc-sinaweibo', 'soc-student-coop', 'soc-friendster', 'soc-twitter']
for graph in nan_graphs:
    df = df[df.Graph != graph]

#removes the categories that are not as helpful in clustering
categories_to_exclude = [ 'DIMACS', 'DIMACS10', 'BHOSLIB', 'Temporal Reachability', 'Dynamic Networks']
for name in categories_to_exclude:
    df = df[df.Collection != name]

# removes any 0s that shouldn't be there
for i in df.index:
   df.loc[i, 'Average triangles'] = df.loc[i, 'Total triangles'] / df.loc[i, 'Nodes']

collections = np.unique( df.Collection.values )
for collection in collections:
    size = len( df[ df.Collection == collection ] )
    if size < 10:
        #print(collection, ' collection is of size ', size)
        df = df[ df.Collection != collection ]

print(df.info())

#df.to_csv( 'C:/Users/Owner/Documents/VERUM/Network stuff/git/src/data/clean_data.csv' )
#df.to_csv( 'C:/Users/Owner/Documents/VERUM/Network stuff/git/src/data/data_without_small_collections.csv' )