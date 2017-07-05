import pandas as pd
from sklearn.preprocessing import MinMaxScaler
import numpy as np

infile = '~/PycharmProjects/network_classification/src/data/clean_data_with_new_chem.csv'
df = pd.read_csv(infile, index_col=0)

graph = list(df['Graph'])
collection = list(df['Collection'])

del df['Graph']
del df['Collection']

array = df.values

scale = MinMaxScaler()

array = scale.fit_transform(array)

#turn that array back into a dataframe!
df1 = pd.DataFrame(array)
df1.columns = df.columns

#add back the graph and collection columns!
df1['Graph']=graph
df1['Collection']=collection

#stick the last 2 columns on front
cols = list(df1.columns)
cols = cols[-2:] + cols[:-2]
df1 = df1[cols]

print(df1)

#df1.to_csv('path goes here')
