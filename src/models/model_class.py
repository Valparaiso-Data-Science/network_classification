import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split, cross_val_score, StratifiedKFold
from sklearn.metrics import classification_report, confusion_matrix


infile = 'C:/Users/Owner/Documents/VERUM/Network stuff/git/src/data/data_minmaxscale.csv' # -- change for machine
df = pd.read_csv(infile, index_col=0)

def modelFitTest(model, df, minSize=20, dropList=['Graph', 'Collection'], split=.25, cv=4, feat_comp = False):


    collections = np.unique(df.Collection.values)
    for collection in collections:
        size = len(df[df.Collection == collection])
        if size < minSize:
            df = df[df.Collection != collection]
    print('Using collections of size >', minSize)

    X = df.drop(dropList, axis=1).values
    y = df['Collection'].values
    print('Excluding categories: ', dropList)

    iterator = StratifiedKFold(n_splits = cv, shuffle = True,# random_state = 42
                                )
    cvscores = cross_val_score(model, X, y, cv = iterator)
    print(cv, '-fold stratified cross validation')

    if feat_comp == False:
        print('cv scores: ', cvscores)
        print('cv average: ', np.mean(cvscores))

        X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=split, stratify=y)
        model.fit( X_train, y_train )
        y_pred = model.predict( X_test )

        print('Test size: ', split)
        print(classification_report(y_test, y_pred))
        print('Confusion matrix')
        print(confusion_matrix(y_test, y_pred))
        print('score of prediction: ', model.score(X_test, y_test))

    else:

        #feature comparison -- compares the current cv average from input to the cv average from using all features
        Xold = df.drop(['Graph', 'Collection'], axis = 1).values
        oldCVscores = cross_val_score(model, Xold, y, cv = iterator)
        print('New cv average - Old cv average = ', np.mean(cvscores) - np.mean(oldCVscores))
        return np.mean(cvscores) - np.mean(oldCVscores)