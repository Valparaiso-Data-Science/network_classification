#notes:
    #-this program was created as a way to test how normalize vs standard scaler affected clustering (kmeans)
    #-I used the imputer to hastily fix the NaNs (right now it is using most frequent strategy)- this will be changed once we find
    #a way to calculate all the missing values.
    #-Eventually I plan to use the argument parser to select which columns (statistics) to look at. Right now it drops Graph,
    #Collection, and Chromatic Number. You can specify whether or not you want to keep cheminformatics.
    #---cheminformatics will be rescaled shortly
    #-Obviously the paths will have to be changed to run on other machines.
    #
    #--Emma


import argparse
import pandas as pd
import numpy as np
from sklearn.preprocessing import Normalizer, Imputer, StandardScaler
from sklearn.cluster import KMeans
from sklearn.pipeline import make_pipeline

import matplotlib.pyplot as plt

def main(args):
    net = pd.read_csv('~/PycharmProjects/network_classification/data/interim/data_without_dimacs_or_bhoslib.csv')

    # getting rid of cheminformatics for fun?
    if args.minus_chem:
        net = net[net['Collection'] != 'Cheminformatics']
    if args.chem:
        net = net[net['Collection'] == 'Cheminformatics']


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

    kmeans = KMeans(n_clusters=10)

    if args.do_norm:
            pipeline = make_pipeline(imp, normalize, kmeans)

    elif args.do_ss:
            pipeline = make_pipeline(imp, standard_scale, kmeans)


    labels = pipeline.predict(net_array)

    #cross tabulation:
    net_ct = pd.DataFrame({'Labels':labels, 'Collection':net_copy['Collection']})
    ct = pd.crosstab(net_ct['Labels'], net_ct['Collection'])

    print(ct)

        #ct.to_csv('~/PycharmProjects/network_classification/models/normalize_crosstab_all_data.csv')


if __name__=="__main__":

    parser=argparse.ArgumentParser(description='Command line graph inference', add_help=True)


    parser.add_argument('-n', action="store_true", dest="do_norm", default=False, help='Use normalization on data')
    parser.add_argument('-s', action="store_true", dest="do_ss", default=False, help='Use standard scaler on data')
    parser.add_argument('-ch', action='store_true', dest='minus_chem', default=False, help='Removes cheminformatics')
    parser.add_argument('-k', action='store_true', dest='do_kmeans', default=False, help='Do kmeans clustering')
    parser.add_argument('-c', action='store_true', dest='chem', default=False, help='Only cheminformatics')

    args=parser.parse_args()

    main(args)