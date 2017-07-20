import pandas as pd
import numpy as np
from sklearn.tree import DecisionTreeClassifier
from sklearn.naive_bayes import GaussianNB
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import GridSearchCV, StratifiedKFold, LeaveOneOut
from sklearn.preprocessing import MinMaxScaler

from git.src.models.model_class import ModelTester

infile = 'C:/Users/Owner/Documents/VERUM/Network stuff/git/src/data/clean_data_with_new_chem.csv' # -- change for machine
df = pd.read_csv(infile, index_col=0)



tester = ModelTester(df)
model = RandomForestClassifier(n_estimators=30)
print(tester.train_test(model))







# Everything below is for finding the best RandomForest parameter, commented out for convenience
# best param changed each time, seemed to be most often around 30

#
#dfNew = df.copy()
## scales the data and removes collections< 20
#graphs = list(dfNew['Graph'])
#categories = list(dfNew['Collection'])
#del dfNew['Graph']
#del dfNew['Collection']
#
#array = dfNew.values
#scale = MinMaxScaler()
#scaledArray = scale.fit_transform(array)
## create the new df
#scaleDF = pd.DataFrame(scaledArray)
#scaleDF.columns = dfNew.columns
#scaleDF['Graph'] = graphs
#scaleDF['Collection'] = categories
#cols = list(scaleDF.columns)
#cols = cols[-2:] + cols[:-2]
#scaleDF = scaleDF[cols]
#
#collections = np.unique(scaleDF.Collection.values)
#for collection in collections:
#    size = len(scaleDF[scaleDF.Collection == collection])
#    if size < 20:
#        scaleDF = scaleDF[scaleDF.Collection != collection]
#
#X = scaleDF.drop(['Graph','Collection'], axis=1).values
#y = scaleDF['Collection'].values


#numberTrees = np.arange(8,70)
#gridParam = { 'n_estimators' : numberTrees }
#
#model = RandomForestClassifier()
#iterator = StratifiedKFold(n_splits = 5, shuffle = True)
#RandForestcv = GridSearchCV(model, gridParam, cv=5)
#
#RandForestcv.fit(X,y)
#bestParam = RandForestcv.best_params_['n_estimators']
#print(bestParam)
#print(RandForestcv.best_score_)
#print(RandForestcv.get_params())
#print(RandForestcv.grid_scores_)