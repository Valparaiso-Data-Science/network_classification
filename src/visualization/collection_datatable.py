import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

net = pd.read_csv('~/PycharmProjects/network_classification/src/data/clean_data_with_new_chem.csv', index_col=0)

collections = [#'Biological Networks'
               'Brain Networks',
               'Cheminformatics',
               #'Collaboration Networks' 
               #'Ecology Networks' 
               'Facebook Networks',
               #'Infrastructure Networks' 
               #'Interaction Networks' 
               #'Massive Network Data'
               # 'Recommendation Networks' 
               'Retweet Networks',
               #'Scientific Computing'
                'Social Networks',
               #'Technological Networks'
               'Web Graphs']

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
data_mean = pd.DataFrame(index=collections, columns=features)

mean = []
for feature in features:
    for collection in collections:
        temp = net[net['Collection']==collection] #make a smaller dataframe to work with
        mean.append(temp[feature].mean(axis=0))
    data_mean[feature] = mean
    mean = []
print(data_mean)


data_min_max = pd.DataFrame(index=collections, columns=features)
minimum = []
maximum = []
min_max = []
for feature in features:
    for collection in collections:
        temp = net[net['Collection']==collection] #make a smaller dataframe to work with
        minimum.append(temp[feature].min(axis=0))
        maximum.append(temp[feature].max(axis=0))
    for i in np.arange(0, len(collections)):
        min_max.append(str(minimum[i]) + '-' + str(maximum[i]))
    data_min_max[feature] = min_max
    min_max = []
    minimum = []
    maximum = []
print(data_min_max)