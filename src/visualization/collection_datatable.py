import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

net = pd.read_csv('~/PycharmProjects/network_classification/src/data/clean_data_with_new_chem.csv', index_col=0)

collections = ['Brain Networks', 'Cheminformatics', 'Facebook Networks', 'Retweet Networks', 'Social Networks', 'Web Graphs']
features = ['Nodes',
           'Edges',
           'Density',
           'Maximum degree',
           'Minimum degree',
           'Average degree',
           'Assortativity',
            'Total triangles',
           'Average triangles',
           'Maximum triangles',
           'Avg. clustering coef.',
           'Frac. closed triangles',
           'Maximum k-core',
           'Max. clique (lb)'
            ]
data = pd.DataFrame(index=collections, columns=features)

mean = []
for feature in features:
    for collection in collections:
        temp = net[net['Collection']==collection] #make a smaller dataframe to work with
        mean.append(temp[feature].mean(axis=0))
    data[feature] = mean
    mean = []
print(data)