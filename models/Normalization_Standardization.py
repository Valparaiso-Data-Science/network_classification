'''
notes about this process:
    -I used the imputer to hastily fix the NaNs (right now it is using most frequent strategy)- this will be changed once we find
    a way to calculate all the missing values.
    -Eventually I plan to use the argument parser to select which columns (statistics) to look at. Right now it drops Graph,
    Collection, and Chromatic Number.
    -I also want to use arguments to specify whether you want a cross tabulation of kmeans alg or tsne plot or both. right now it
    just does a tsne plot
    -Obviously the paths will have to be changed to run on other machines.
    --Emma
'''

import argparse
import pandas as pd
import numpy as np
from sklearn.preprocessing import Normalizer, Imputer, StandardScaler, RobustScaler
from sklearn.cluster import KMeans
from sklearn.pipeline import make_pipeline
from sklearn.manifold import TSNE
import matplotlib.pyplot as plt

def main(args):
    net = pd.read_csv('~/PycharmProjects/network_classification/data/interim/data_without_dimacs_or_bhoslib.csv')

    # getting rid of cheminformatics for fun?
    if args.minus_chem:
        net = net[net['Collection'] != 'Cheminformatics']
    if args.chem:
        net = net[net['Collection'] == 'Cheminformatics']

    #i need a list of the collections as numbers-useful in coloring tsne
    graph_categories = []
    all_categories = net['Collection'].values
    #creates a list of the existing categories, no repeats
    for category in all_categories:
        if not category in graph_categories:
            graph_categories.append( category )
    #re-assigns each collection name in all_categories to an integer, corresponding with its position in graph_categories
    for i in range(len( all_categories )):
        all_categories[i] = graph_categories.index( all_categories[i] )


    net_copy = pd.DataFrame.copy(net) #copy used for categories in the crosstab section


    del net['Graph']
    del net['Collection']
    del net['Chromatic Number']

    net_array = net.values

    imp = Imputer(missing_values='NaN', strategy='most_frequent', axis=0)

    if args.do_norm:
        print('Normalizing')
        normalize = Normalizer()
    elif args.do_ss:
        print('Standardizing')
        standard_scale = StandardScaler()
    #robust = RobustScaler()

    tsne = TSNE()

    if args.do_norm:
            pipeline = make_pipeline(imp, normalize, tsne)

    elif args.do_ss:
            pipeline = make_pipeline(imp, standard_scale, tsne)


    if args.do_tsne:
        tsne_features = pipeline.fit_transform(net_array)
        xs = tsne_features[:,0]
        ys = tsne_features[:,1]
        plt.scatter(xs, ys, c=all_categories)
        plt.show()



    #labels = pipeline.predict(net_array)

    #cross tabulation:
    #net_ct = pd.DataFrame({'Labels':labels, 'Collection':net_copy['Collection']})
    #ct = pd.crosstab(net_ct['Labels'], net_ct['Collection'])

    #print(ct)

        #ct.to_csv('~/PycharmProjects/network_classification/models/normalize_crosstab_all_data.csv')


if __name__=="__main__":

    parser=argparse.ArgumentParser(description='Command line graph inference', add_help=True)


    parser.add_argument('-n', action="store_true", dest="do_norm", default=False, help='Use normalization on data')
    parser.add_argument('-s', action="store_true", dest="do_ss", default=False, help='Use standard scaler on data')
    parser.add_argument('-ch', action='store_true', dest='minus_chem', default=False, help='Removes cheminformatics')
    parser.add_argument('-k', action='store_true', dest='do_kmeans', default=False, help='Do kmeans clustering')
    parser.add_argument('-t', action='store_true', dest='do_tsne', default=False, help='Do TSNE visualization')
    parser.add_argument('-c', action='store_true', dest='chem', default=False, help='Only cheminformatics')

    args=parser.parse_args()

    main(args)