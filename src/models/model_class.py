import pandas as pd
import numpy as np


class model:

    def __init__(self, model ):
        self.model = model()

    def test(self):
        print( self.model, '"this is a test"')


from sklearn.model_selection import train_test_split
from sklearn.svm import SVC
from sklearn.metrics import classification_report

infile = 'C:/Users/Owner/Documents/VERUM/Network stuff/git/src/data/data_minmaxscale.csv' # -- change for machine
df = pd.read_csv(infile, index_col=0)

X = df.drop(['Graph', 'Collection'], axis=1).values
y = df['Collection'].values

X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=.2)


model = SVC()  # -- change this for the model you want to use

model.fit( X_train, y_train )
y_pred = model.predict( X_test )

print(classification_report(y_test, y_pred))