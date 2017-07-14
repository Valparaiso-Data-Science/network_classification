import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split, cross_val_score, StratifiedKFold, cross_val_predict
from sklearn.metrics import classification_report, confusion_matrix


infile = 'C:/Users/Owner/Documents/VERUM/Network stuff/git/src/data/data_minmaxscale.csv' # -- change for machine
df = pd.read_csv(infile, index_col=0)


class ModelTester():

    # constructor that assigns a pandas dataframe as internal data
    def __init__(self, df):
        self.df = df

    # returns self.df when printed
    def __str__(self):
        return str(self.df)

# Returns classification report and confusion matrix from cv over entire dataframe. If feat_comp == True, returns the difference in cv score from df with all features included
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

            return cvscores, np.mean(cvscores), classification_report(y, cv_pred), confusion_matrix(y, cv_pred)

        else:
            #feature comparison -- compares the current cv average from input to the cv average from using all features
            Xold = dfNew.drop(['Graph', 'Collection'], axis = 1).values
            oldCVscores = cross_val_score(model, Xold, y, cv = iterator)

            if prnt == True:
                print('New cv average - Old cv average = ', np.mean(cvscores) - np.mean(oldCVscores))
                print( 'New cv average: ', np.mean(cvscores))

            return np.mean(cvscores) - np.mean(oldCVscores)


# Returns a dataframe with the names of mislabeled graphs (from one full cv run of the model), the actual collections of these graphs, and the predicted collections
    def get_mislabel_analysis(self, model, minSize=20, dropList=['Graph', 'Collection'], cv=5, prnt=True):

        dfNew = self.df.copy()
        collections = np.unique(dfNew.Collection.values)
        for collection in collections:
            size = len(dfNew[dfNew.Collection == collection])
            if size < minSize:
                dfNew = dfNew[dfNew.Collection != collection]

        X = dfNew.drop(dropList, axis=1).values
        y = dfNew['Collection'].values
        names = dfNew['Graph'].values

        iterator = StratifiedKFold(n_splits=cv, shuffle=True,  # random_state = 42
                                   )
        cv_pred = cross_val_predict(model, X, y, cv=iterator)
        cv_results_dict = {'Name': names, 'Actual': y, 'Predicted': cv_pred}

        column_order = ['Name', 'Actual', 'Predicted']
        results = pd.DataFrame(cv_results_dict, columns=column_order, copy=True)

        if prnt == True:
            print('Using collections of size >', minSize)
            print('Excluding features: ', dropList)
            print(classification_report(y, cv_pred))
            print(confusion_matrix(y, cv_pred))


        return results[results.Actual != results.Predicted]

    # returns a new dataframe consisting of the graphs that were mislabeled 'percent_return' - percent of the time after 'repeat' tests of the model
    # includes in new dataframe how many times each graph was mislabeled, and in which category
    def get_mislabeled_graphs(self, model, repeat=100, percent_return=0.5, minSize=20, dropList=['Graph', 'Collection'],
                              cv=5):

        dfNew = self.df.copy()
        allNameDict = {'CountTot': np.zeros(len(dfNew.Graph.values)),
                       'CountWeb': np.zeros(len(dfNew.Graph.values)),
                       'CountFb': np.zeros(len(dfNew.Graph.values)),
                       'CountSoc': np.zeros(len(dfNew.Graph.values)),
                       'CountRt': np.zeros(len(dfNew.Graph.values)),
                       'CountChem': np.zeros(len(dfNew.Graph.values)),
                       'CountBn': np.zeros(len(dfNew.Graph.values)), }
        order = ['CountTot', 'CountBn', 'CountChem', 'CountFb', 'CountRt', 'CountSoc', 'CountWeb']
        all_names = pd.DataFrame(allNameDict, columns=order, index=dfNew.Graph.values, copy=True)

        for i in range(repeat):
            analysis = self.get_mislabel_analysis(model, minSize, dropList, cv, prnt=False)
            for graph in analysis.Name.values:
                all_names.loc[graph, 'CountTot'] += 1

                if analysis[analysis.Name == graph].iloc[0, 2] == 'Brain Networks':
                    all_names.loc[graph, 'CountBn'] += 1
                elif analysis[analysis.Name == graph].iloc[0, 2] == 'Cheminformatics':
                    all_names.loc[graph, 'CountChem'] += 1
                elif analysis[analysis.Name == graph].iloc[0, 2] == 'Facebook Networks':
                    all_names.loc[graph, 'CountFb'] += 1
                elif analysis[analysis.Name == graph].iloc[0, 2] == 'Retweet Networks':
                    all_names.loc[graph, 'CountRt'] += 1
                elif analysis[analysis.Name == graph].iloc[0, 2] == 'Social Networks':
                    all_names.loc[graph, 'CountSoc'] += 1
                elif analysis[analysis.Name == graph].iloc[0, 2] == 'Web Graphs':
                    all_names.loc[graph, 'CountWeb'] += 1

        percentOver = all_names[all_names['CountTot'] >= percent_return * repeat]
        return percentOver