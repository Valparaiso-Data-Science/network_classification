import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split, cross_val_score, StratifiedKFold, cross_val_predict
from sklearn.metrics import classification_report, confusion_matrix
from sklearn.preprocessing import MinMaxScaler


#infile = 'C:/Users/Owner/Documents/VERUM/Network stuff/git/src/data/data_minmaxscale.csv' # -- change for machine
#df = pd.read_csv(infile, index_col=0)


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

        #scales the data
        graphs = list(dfNew['Graph'])
        categories = list(dfNew['Collection'])
        del dfNew['Graph']
        del dfNew['Collection']

        array = dfNew.values
        scale = MinMaxScaler()
        scaledArray = scale.fit_transform(array)
#create the new df
        scaleDF = pd.DataFrame(scaledArray)
        scaleDF.columns = dfNew.columns
        scaleDF['Graph'] = graphs
        scaleDF['Collection'] = categories
        cols = list(scaleDF.columns)
        cols = cols[-2:] + cols[:-2]
        scaleDF = scaleDF[cols]



        collections = np.unique(scaleDF.Collection.values)
        for collection in collections:
            size = len(scaleDF[scaleDF.Collection == collection])
            if size < minSize:
                scaleDF = scaleDF[scaleDF.Collection != collection]


        X = scaleDF.drop(dropList, axis=1).values
        y = scaleDF['Collection'].values


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
            Xold = scaleDF.drop(['Graph', 'Collection'], axis = 1).values
            oldCVscores = cross_val_score(model, Xold, y, cv = iterator)

            if prnt == True:
                print('New cv average - Old cv average = ', np.mean(cvscores) - np.mean(oldCVscores))
                print( 'New cv average: ', np.mean(cvscores))

            return np.mean(cvscores) - np.mean(oldCVscores)


# Returns a dataframe with the names of mislabeled graphs (from one full cv run of the model), the actual collections of these graphs, and the predicted collections
    def get_mislabel_analysis(self, model, minSize=20, dropList=['Graph', 'Collection'], cv=5, prnt=True):

        dfNew = self.df.copy()

        # scales the data
        graphs = list(dfNew['Graph'])
        categories = list(dfNew['Collection'])
        del dfNew['Graph']
        del dfNew['Collection']

        array = dfNew.values
        scale = MinMaxScaler()
        scaledArray = scale.fit_transform(array)
        # create the new df
        scaleDF = pd.DataFrame(scaledArray)
        scaleDF.columns = dfNew.columns
        scaleDF['Graph'] = graphs
        scaleDF['Collection'] = categories
        cols = list(scaleDF.columns)
        cols = cols[-2:] + cols[:-2]
        scaleDF = scaleDF[cols]

        collections = np.unique(scaleDF.Collection.values)
        for collection in collections:
            size = len(scaleDF[scaleDF.Collection == collection])
            if size < minSize:
                scaleDF = scaleDF[scaleDF.Collection != collection]

        X = scaleDF.drop(dropList, axis=1).values
        y = scaleDF['Collection'].values
        names = scaleDF['Graph'].values

        iterator = StratifiedKFold(n_splits=cv, shuffle=True,  # random_state = 42
                                   )
        cv_pred = cross_val_predict(model, X, y, cv=iterator)
        cv_results_dict = {'Name': names, 'Actual': y, 'Predicted': cv_pred}

        column_order = ['Name', 'Actual', 'Predicted']
        results = pd.DataFrame(cv_results_dict, columns=column_order, copy=True)

        if prnt == True:
            print('Using collections of size >', minSize)
            print('Excluding features: ', dropList)




        return results[results.Actual != results.Predicted]




    # returns a new dataframe consisting of the graphs that were mislabeled 'percent_return' - percent of the time after 'repeat' tests of the model
    # includes in new dataframe how many times each graph was mislabeled, and in which category
    # externalData is used to run this code testing other data, not doing cv (like with misc graphs)
    def get_mislabeled_graphs(self, model, externalData=False, repeat=100, percent_return=0.5, minSize=20, dropList=['Graph', 'Collection'],
                              cv=5):


        dfNew = self.df.copy()

#scales the data
        graphs = list(dfNew['Graph'])
        categories = list(dfNew['Collection'])
        del dfNew['Graph']
        del dfNew['Collection']

        array = dfNew.values
        scale = MinMaxScaler()
        scaledArray = scale.fit_transform(array)
#create the new df
        scaleDF = pd.DataFrame(scaledArray)
        scaleDF.columns = dfNew.columns
        scaleDF['Graph'] = graphs
        scaleDF['Collection'] = categories
        cols = list(scaleDF.columns)
        cols = cols[-2:] + cols[:-2]
        scaleDF = scaleDF[cols]

# This removes the collections that are too small
        collections = np.unique(scaleDF.Collection.values)
        for collection in collections:
            size = len(scaleDF[scaleDF.Collection == collection])
            if size < minSize:
                scaleDF = scaleDF[scaleDF.Collection != collection]

# Creates a dictionary that will be used to create a new dataframe. It is different if working with externalData
        collectionNames = np.unique(scaleDF.Collection.values)
        if type(externalData) == bool:
            allNameDict = {'TotalCount': np.zeros(len(scaleDF.Graph.values))}
            for i in range( len(collectionNames) ):
                allNameDict[collectionNames[i]] = np.zeros(len(scaleDF.Graph.values))
            order = ['TotalCount'] + list(collectionNames)
# for working with externalData
        # this one is specifically for misc graphs
        elif 'Collection Hypothesis' in externalData.columns:
            allNameDict = {'Hypothesis' : externalData['Collection Hypothesis'].values, 'TotalCount': np.zeros(len(externalData.Graph.values))}
            for i in range(len(collectionNames)):
                allNameDict[collectionNames[i]] = np.zeros(len(externalData.Graph.values))
            order = ['Hypothesis', 'TotalCount'] + list(collectionNames)
        else:
            allNameDict = {'TotalCount': np.zeros(len(externalData.Graph.values))}
            for i in range(len(collectionNames)):
                allNameDict[collectionNames[i]] = np.zeros(len(externalData.Graph.values))
            order = ['TotalCount'] + list(collectionNames)


# Creates DF with organzied column order, and index depending on if we are running the method on self.df data or external data (misc graphs)
        if type(externalData) == bool:
            all_names = pd.DataFrame(allNameDict, columns=order, index=scaleDF.Graph.values, copy=True)

   # For external data
        else:
            all_names = pd.DataFrame(allNameDict, columns=order, index=externalData.Graph.values, copy=True)

     #This has the method run on self.df data, with cv
        if type(externalData) == bool:
          # This trains and test the model 'repeat' number of times, keeping track of how each graph was mislabeled each time
            for i in range(repeat):
                analysis = self.get_mislabel_analysis(model, minSize, dropList, cv, prnt=False)

    # Keeps the count
                for graph in analysis.Name.values:
                    all_names.loc[graph, 'TotalCount'] += 1

                    mislabel = analysis[analysis.Name == graph].iloc[0, 2]
                    all_names.loc[graph, mislabel] += 1

    # This is for externalData == a dataframe, so the method runs testing some external data (aka misc graphs)
        else:
            # This trains the model and tests on the externalData 'repeat' number of times, keeping track of how each graph was mislabeled each time
            for i in range(repeat):
                analysis = self.train_predict(model, externalData, dropList, minSize)
                mislabeled = analysis[ analysis.Hypothesis != analysis.Predicted ]

                # Keeps the count
                for graph in mislabeled.Name.values:
                    all_names.loc[graph, 'TotalCount'] += 1

                    mislabel = mislabeled[mislabeled.Name == graph].iloc[0, 2]
                    all_names.loc[graph, mislabel] += 1

# This returns a Df with only the graphs that were labeled 'percent_return' -percent of the time
        percentOver = all_names[all_names['TotalCount'] >= percent_return * repeat]
        return percentOver

# returns new df with collections in collectionList combined into one new collection: newName
    def combine_collections(self, collectionList, newName):
        dfNew = self.df.copy()
        collectionDFs = []
        for label in collectionList:
            dfSpec = dfNew[ dfNew.Collection == label ]
            collectionDFs.append( dfSpec )
            for i in dfSpec.index:
                dfNew.loc[i, 'Collection'] = newName

#for the misc df
            if 'Collection Hypothesis' in dfNew.columns:
                dfSpec = dfNew[ dfNew['Collection Hypothesis'] == label ]
                for i in dfSpec.index:
                    dfNew.loc[i, 'Collection Hypothesis'] = newName

        return dfNew

# Trains a model on full set of data, and then predicts on new data 'testDF'
    def train_predict(self, model, testDF, dropList=['Graph', 'Collection'], minSize=20):
        dfNew = self.df.copy()
        graphs = list(dfNew['Graph'])
        categories = list(dfNew['Collection'])
        del dfNew['Graph']
        del dfNew['Collection']

        array = dfNew.values
        scale = MinMaxScaler()
        scaledArray = scale.fit_transform(array)
        # create the new df
        scaleDF = pd.DataFrame(scaledArray)
        scaleDF.columns = dfNew.columns
        scaleDF['Graph'] = graphs
        scaleDF['Collection'] = categories
        cols = list(scaleDF.columns)
        cols = cols[-2:] + cols[:-2]
        scaleDF = scaleDF[cols]

        # This removes the collections that are too small
        collections = np.unique(scaleDF.Collection.values)
        for collection in collections:
            size = len(scaleDF[scaleDF.Collection == collection])
            if size < minSize:
                scaleDF = scaleDF[scaleDF.Collection != collection]

# Creates the training data
        X =scaleDF.drop(dropList, axis=1).values
        y = scaleDF['Collection'].values

# Sets up the testing data
        testNames = testDF.Graph.values
        hypothesis = testDF['Collection Hypothesis'].values
        if 'Collection Hypothesis' in testDF.columns:
            testFeatures = testDF.drop(['Collection Hypothesis'] + dropList, axis=1)
        else:
            testFeatures = testDF.drop(dropList, axis=1)

        model.fit(X, y)
        scaledTestFeatures = scale.transform(testFeatures)
        pred = model.predict(scaledTestFeatures)

        resultsDict = {'Name' : testNames, 'Hypothesis' : hypothesis, 'Predicted' : pred}
        colOrder = ['Name', 'Hypothesis', 'Predicted']
        results = pd.DataFrame(resultsDict, columns=colOrder)
        return results