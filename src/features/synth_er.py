import pandas as pd

# This script conbines a couple csv's in order to create the full Erdos Renyi csv

infile = 'C:/Users/Owner/Documents/VERUM/Network stuff/git/src/features/synthetic_features_e1.csv'
df = pd.read_csv(infile)
print(df.head())

infile2 = 'C:/Users/Owner/Documents/VERUM/Network stuff/git/src/features/synthetic_features_e2e3.csv'
df2 = pd.read_csv(infile2)
print(df2.tail())

df_er = pd.concat([df, df2])

print(df_er.Graph.values)
#print(df.info())

df_er.to_csv('C:/Users/Owner/Documents/VERUM/Network stuff/git/src/data/synthetic_e1e2e3_complete.csv')