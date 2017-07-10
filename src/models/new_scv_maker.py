import pandas as pd
import numpy as np

infile = 'C:/Users/Owner/Documents/VERUM/Network stuff/git/src/data/data_minmaxscale.csv' # -- change for machine
df = pd.read_csv(infile, index_col=0)


minSize = 20
collections = np.unique(df.Collection.values)
for collection in collections:
    size = len(df[df.Collection == collection])
    if size < minSize:
        df = df[df.Collection != collection]

df.to_csv('C:/Users/Owner/Documents/VERUM/Network stuff/git/src/data/data_greater_than_20.csv')