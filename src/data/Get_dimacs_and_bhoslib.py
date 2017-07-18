# This program retrieves the collections DIMACS, DIMACS10 and BHOSLIB
# from the raw file 'network-meta-statistics2.csv'

import pandas as pd

# Read file
df = pd.read_csv('~/Downloads/network_classification/data/raw/network-meta-statistics2.csv')

# Change column names
df.columns = (['Graph', 'Collection', 'Nodes', 'Edges', 'Density', 'Maximum degree',
       'Minimum degree', 'Average degree', 'Assortativity', 'Total triangles',
       'Average triangles', 'Maximum triangles', 'Avg. clustering coef.',
        'Frac. closed triangles', 'Maximum k-core', 'Max. clique (lb)', 'Chromatic Number'])

# Get only the dimacs and bhoslib collections
dimacs_data = df[df['Collection'] == 'DIMACS']
dimacs10_data = df[df['Collection'] == 'DIMACS10']
bhoslib_data = df[df['Collection'] == 'BHOSLIB']

df_new = pd.concat([dimacs_data, dimacs10_data, bhoslib_data])

# Get rid of Chromatic Number
del df_new['Chromatic Number']

# Drop NaN values of the frame
# Doing this drops 7 graphs
df_new = df_new.dropna()


# Move data to a csv file
df_new.to_csv('~/Downloads/network_classification/src/data/dimacs_and_bhoslib_networks.csv')