import pandas as pd
from sklearn.cluster import KMeans
from sklearn.preprocessing import Normalizer, Imputer
import matplotlib.pyplot as plt
from sklearn.pipeline import make_pipeline

net = pd.read_csv('~/PycharmProjects/network_classification/data/interim/data_without_dimacs_or_bhoslib.csv')

imp = Imputer(missing_values='NaN', strategy='most_frequent', axis=0)
normalize = Normalizer()

#take out cheminformatics
net = net[net['Collection']!='Cheminformatics']

del net['Graph']
del net['Collection']
del net['Chromatic Number']


net_array = net.values

ks = range(1, 20)
inertias = []

for k in ks:
    # Create a KMeans instance with k clusters: model
    model = KMeans(n_clusters=k)

    pipe = make_pipeline(imp, normalize, model)
    # Fit model to samples
    pipe.fit(net_array)

    # Append the inertia to the list of inertias
    inertias.append(model.inertia_)

# Plot ks vs inertias
plt.plot(ks, inertias, '-o')
plt.xlabel('number of clusters, k')
plt.ylabel('inertia')
plt.xticks(ks)
plt.show()