import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split, cross_val_score, StratifiedKFold, cross_val_predict
from sklearn.metrics import classification_report, confusion_matrix


infile = 'C:/Users/Owner/Documents/VERUM/Network stuff/git/src/data/data_minmaxscale.csv' # -- change for machine
df = pd.read_csv(infile, index_col=0)


class DataFrame():

    # constructor that assigns a pandas dataframe as internal data
    def __init__(self, df):
        self.df = df


    def modelFitTest(self, model, minSize=20, dropList=['Graph', 'Collection'], cv=5, feat_comp=False, prnt=True):

        dfNew = self.df.copy()
        collections = np.unique(dfNew.Collection.values)
        for collection in collections:
            size = len(dfNew[dfNew.Collection == collection])
            if size < minSize:
                dfNew = dfNew[dfNew.Collection != collection]


        X = dfNew.drop(dropList, axis=1).values
        y = dfNew['Collection'].values


        iterator = StratifiedKFold(n_splits = cv, shuffle = True,# random_state = 42
                                        )
        cvscores = cross_val_score(model, X, y, cv = iterator)


        if prnt == True:
            print('Using collections of size >', minSize)
            print('Excluding categories: ', dropList)
            print(cv, '-fold stratified cross validation')

        if feat_comp == False:

            #X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=split, stratify=y)
            #model.fit( X_train, y_train )
            #y_pred = model.predict( X_test )

            cv_pred = cross_val_predict(model, X, y, cv = iterator)

            # To include prints of the outputs
            if prnt == True:
                print('cv scores: ', cvscores)
                print('cv average: ', np.mean(cvscores))
                #print('Test size: ', split)
                print(classification_report(y, cv_pred))
                print('Confusion matrix')
                print(confusion_matrix(y, cv_pred))
                #print('score of prediction: ', model.score(X_test, y_test))


        else:
            #feature comparison -- compares the current cv average from input to the cv average from using all features
            Xold = dfNew.drop(['Graph', 'Collection'], axis = 1).values
            oldCVscores = cross_val_score(model, Xold, y, cv = iterator)

            if prnt == True:
                print('New cv average - Old cv average = ', np.mean(cvscores) - np.mean(oldCVscores))
                print( 'New cv average: ', np.mean(cvscores))

            return np.mean(cvscores) - np.mean(oldCVscores)

