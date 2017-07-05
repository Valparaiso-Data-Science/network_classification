import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split, cross_val_score
from sklearn.svm import SVC, LinearSVC
from sklearn.metrics import classification_report, confusion_matrix
from sklearn.tree import DecisionTreeClassifier
from sklearn.naive_bayes import GaussianNB
from sklearn.linear_model import LogisticRegression

infile = 'C:/Users/Owner/Documents/VERUM/Network stuff/git/src/data/data_minmaxscale.csv' # -- change for machine
df = pd.read_csv(infile, index_col=0)


# -- gets rid of small collections -- in the future we will import a file that already has this, so we will delete this block later
collections = np.unique( df.Collection.values )
dropoff = 0 #only looks at collections with more than dropoff graphs
for collection in collections:
    size = len( df[ df.Collection == collection ] )
    if size < dropoff:
        df = df[ df.Collection != collection ]

drop_list = ['Graph', 'Collection',
            'Nodes',
             #'Edges', 'Density',
             'Maximum degree', #'Minimum degree', 'Average degree', 'Assortativity',
             'Total triangles',
             'Average triangles', 'Maximum triangles', #'Avg. clustering coef.', 'Frac. closed triangles',
             'Maximum k-core', 'Max. clique (lb)'
            ]
X = df.drop(drop_list, axis=1).values
y = df['Collection'].values
split = .25 # -- percent of data to hold for test set
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size = split)


model = DecisionTreeClassifier()  # -- change this for the model you want to use

def modelFitTest(model, X, y, split=.25, cv=4):

    cvscores = cross_val_score(model, X, y, cv = cv)
    print(cv, '-fold cross validation')
    print('cv scores: ', cvscores)
    print('cv average: ', np.mean(cvscores))

    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=split)
    model.fit( X_train, y_train )
    y_pred = model.predict( X_test )

    print('Test size: ', split)
    print(classification_report(y_test, y_pred))
    print('Confusion matrix')
    print(confusion_matrix(y_test, y_pred))
    print('score of prediction: ', model.score(X_test, y_test))