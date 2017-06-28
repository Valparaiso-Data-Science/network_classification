import pandas as pd
import numpy as np
from sklearn.svm import SVC
from sklearn.model_selection import train_test_split

infile = 'C:/Users/Owner/Documents/VERUM/Network stuff/git/src/data/data_minmaxscale.csv' # -- change for machine
df = pd.read_csv(infile, index_col=0)

X = df.drop(['Graph', 'Collection'], axis=1).values
y = df['Collection'].values

X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=.2)