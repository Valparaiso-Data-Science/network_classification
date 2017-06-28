#notes:
    #-this program was created as a way to test how normalize vs standard scaler affected clustering (kmeans)
    #-Eventually I plan to use the argument parser to select which columns (statistics) to look at. Right now it drops Graph,
    #Collection, and Chromatic Number. You can specify whether or not you want to keep cheminformatics.
    #---cheminformatics will be rescaled shortly
    #-Obviously the paths will have to be changed to run on other machines.
    #
    #--Emma


import argparse
import pandas as pd
import numpy as np
from sklearn.preprocessing import Normalizer, Imputer, StandardScaler, MinMaxScaler
from sklearn.cluster import KMeans
from sklearn.pipeline import make_pipeline

import matplotlib.pyplot as plt

def main(args):
    net = pd.read_csv('~/PycharmProjects/network_classification/src/data/clean_data_with_new_chem.csv', index_col = 0)
    # getting rid of cheminformatics for fun?
    if args.minus_chem:
        net = net[net['Collection'] != 'Cheminformatics']
    if args.chem:
        net = net[net['Collection'] == 'Cheminformatics']


    net_copy = pd.DataFrame.copy(net) #copy used for categories in the crosstab section


    del net['Graph']
    del net['Collection']


    net_array = net.values

    imp = Imputer(missing_values='NaN', strategy='most_frequent', axis=0)


    if args.do_ss:
        standard_scale = StandardScaler()
    elif args.do_mms:
        min_max_scale = MinMaxScaler()

    kmeans = KMeans(n_clusters=4)

    if args.do_mms:
            pipeline = make_pipeline(imp, min_max_scale, kmeans)

    elif args.do_ss:
            pipeline = make_pipeline(imp, standard_scale, kmeans)

    pipeline.fit(net_array)
    labels = pipeline.predict(net_array)

    #cross tabulation:
    net_ct = pd.DataFrame({'Labels':labels, 'Collection':net_copy['Collection']})
    ct = pd.crosstab( net_ct['Collection'], net_ct['Labels'])

    print(ct)

        #ct.to_csv('~/PycharmProjects/network_classification/models/normalize_crosstab_all_data.csv')


if __name__=="__main__":

    parser=argparse.ArgumentParser(description='Command line graph inference', add_help=True)


    parser.add_argument('-m', action="store_true", dest="do_mms", default=False, help='Use minmax scaler on data')
    parser.add_argument('-s', action="store_true", dest="do_ss", default=False, help='Use standard scaler on data')
    parser.add_argument('-ch', action='store_true', dest='minus_chem', default=False, help='Removes cheminformatics')
    parser.add_argument('-k', action='store_true', dest='do_kmeans', default=False, help='Do kmeans clustering')
    parser.add_argument('-c', action='store_true', dest='chem', default=False, help='Only cheminformatics')

    args=parser.parse_args()

    main(args)