import pandas as pd
from sklearn.cluster import KMeans
import matplotlib.pyplot as plt

net = pd.read_csv('~/PycharmProjects/network_classification/src/data/data_minmaxscale.csv', index_col=0)

collection = list(net['Collection'])
graph = list(net['Graph'])

del net['Graph']
del net['Collection']
net_array = net.values

kmeans = KMeans(n_clusters=8)
labels = kmeans.fit_predict(net_array)

#***************
# TESTING TO SEE WHAT IS HAPPENING
#***************
#print(labels)
#net_ct = pd.DataFrame({'Labels':labels, 'Collection':collection})
#ct = pd.crosstab( net_ct['Collection'], net_ct['Labels'])
#print(ct)
