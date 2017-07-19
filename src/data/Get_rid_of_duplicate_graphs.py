# This program deletes the duplicate graphs from
# the miscellaneous collection file

import pandas as pd

# Read file
df = pd.read_csv('~/Downloads/network_classification/src/data/miscellaneous_networks.csv')
del df['Unnamed: 0']

# List of all the indeces of the graphs that we want to delete
index_list = [55,20,21,22,23,52,53,54,35]

# Delete unwanted rows
df = df.drop(df.index[index_list])

#print(len(df))
#print(df.head())
print(df[df['Graph'] == 'tech-as-skitter'])

# Copy graphs to a new csv file
df.to_csv('~/Downloads/network_classification/src/data/reduced_miscellaneous_networks.csv')