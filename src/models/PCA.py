import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt

infile = '~/PycharmProjects/network_classification/src/data/data_minmaxscale.csv'
net = pd.read_csv(infile, index_col=0)

graph = list(net['Graph'])
collection = list(net['Collection'])

#i need a list of the collections as numbers-useful in coloring plot
graph_categories = []
all_categories = net['Collection'].values
all_cat_names = list(net['Collection'].values)
#creates a list of the existing categories, no repeats
for category in all_categories:
    if not category in graph_categories:
        graph_categories.append( category )
#re-assigns each collection name in all_categories to an integer, corresponding with its position in graph_categories
for i in range(len( all_categories )):
    all_categories[i] = graph_categories.index( all_categories[i] )
all_graph_names = list(net['Graph'])

del net['Graph']
del net['Collection']

samples = net.values


#intrinsic dimension calculation:
pca = PCA()
pca.fit(samples)
features = range(pca.n_components_)

#plot the feature variance to find intrinsic dimension
#plt.bar(features, pca.explained_variance_)
#plt.xticks(features)
#plt.ylabel('variance')
#plt.xlabel('PCA feature')
#plt.show()


pca_dr = PCA(n_components=2)
pca_dr.fit(samples)
transformed = pca_dr.transform(samples)

#plot data in only 2 dimensions-trying to visualize data to gain insight

xs = transformed[:,0]
ys = transformed[:,1]

plt.scatter(xs, ys, c=all_categories)
plt.show()