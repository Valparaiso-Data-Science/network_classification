import pandas as pd
from sklearn.cluster import KMeans
from sklearn.preprocessing import MinMaxScaler
import matplotlib.pyplot as plt

net = pd.read_csv('~/PycharmProjects/network_classification/src/data/data_without_dimacs_or_bhoslib.csv', index_col = 0)

chem = net[net['Collection'] == 'Cheminformatics']

chem_kmeans = chem.copy()

del chem_kmeans['Graph']
del chem_kmeans['Collection']
del chem_kmeans['Chromatic Number']
del chem['Chromatic Number']

chem_array = chem_kmeans.values
scale = MinMaxScaler()

#**************************
# Getting inertia graph
#**************************
def inertia_plot(data_array):
    ks = range(1, 20)
    inertias = []

    for k in ks:
        # Create a KMeans instance with k clusters: model
        model = KMeans(n_clusters=k)

        # Fit model to samples
        # Change to the data we need to check
        model.fit(data_array)

        # Append the inertia to the list of inertias
        inertias.append(model.inertia_)

    # Plot ks vs inertias
    plt.plot(ks, inertias, '-o')
    plt.xlabel('number of clusters, k')
    plt.ylabel('inertia')
    plt.xticks(ks)
    plt.title('Cheminformatics Inertia vs. k')
    plt.show()


inertia_plot(chem_array)
kmeans = KMeans(n_clusters=4)

kmeans.fit(chem_array)
labels = kmeans.predict(chem_array)

chem['Labels'] = labels

chem_0 = chem[chem['Labels']==0]
chem_0_new = chem_0.sample(frac = 0.2, random_state=42)

chem_1 = chem[chem['Labels']==1]
chem_1_new = chem_1.sample(frac = 0.2, random_state=42)

chem_2 = chem[chem['Labels']==2]
chem_2_new = chem_2.sample(frac = 0.2, random_state=42)

chem_3 = chem[chem['Labels']==3]
chem_3_new = chem_3.sample(frac = 0.2, random_state=42)

new_chem = pd.concat([chem_0_new, chem_1_new, chem_2_new, chem_3_new])
del new_chem['Labels']

#new_chem.to_csv('~/PycharmProjects/network_classification/src/data/downsampled_chem.csv')