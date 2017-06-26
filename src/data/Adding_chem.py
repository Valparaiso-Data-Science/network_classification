import pandas as pd


clean = pd.read_csv('~/PycharmProjects/network_classification/src/data/clean_data.csv', index_col = 0)
chem = pd.read_csv('~/PycharmProjects/network_classification/src/data/downsampled_chem.csv', index_col = 0)

#take out cheminformatics from data
clean = clean[clean['Collection'] != 'Cheminformatics']

#add new cheminformatics
clean_new = pd.concat([clean, chem])

#clean_new.to_csv('~/PycharmProjects/network_classification/src/data/clean_data_with_new_chem.csv')


