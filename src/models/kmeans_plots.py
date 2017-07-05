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

#************
# SEPARATING LABELS
#************

net['Label'] = labels

c0 = net[net['Label']==0]
c1 = net[net['Label']==1]
c2 = net[net['Label']==2]
c3 = net[net['Label']==3]
c4 = net[net['Label']==4]
c5 = net[net['Label']==5]
c6 = net[net['Label']==6]
c7 = net[net['Label']==7]

