import pandas as pd
import numpy as np
from sklearn.tree import DecisionTreeClassifier
from sklearn.naive_bayes import GaussianNB
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.svm import SVC, LinearSVC
from sklearn.model_selection import GridSearchCV, StratifiedKFold, LeaveOneOut
from sklearn.preprocessing import MinMaxScaler
import sys
import time

from git.src.models.model_class import ModelTester

infile = 'C:/Users/Owner/Documents/VERUM/Network stuff/git/src/data/clean_data_with_new_chem.csv' # -- change for machine
df = pd.read_csv(infile, index_col=0)

remove = ['Graph', 'Collection',
          # 'Nodes',
          # 'Edges',
          # 'Density',
           'Maximum degree',
           'Minimum degree',
          # 'Average degree',
           'Assortativity',
            'Total triangles',
          # 'Average triangles',
          # 'Maximum triangles',
          # 'Avg. clustering coef.',
          # 'Frac. closed triangles',
          # 'Maximum k-core',
          # 'Max. clique (lb)'
            ]
tester = ModelTester(df)
forestModel = RandomForestClassifier(30)
NBModel = GaussianNB()


combined = tester.combine_collections(['Brain Networks', 'Biological Networks'], 'Brain-Bio')
combTesterBnBio = ModelTester(combined)

combined = tester.combine_collections(['Web Graphs', 'Technological Networks'], 'Web-Tech')
combTesterWbTch = ModelTester(combined)

combined = combTesterBnBio.combine_collections(['Web Graphs', 'Technological Networks'], 'Web-Tech')
combTesterBoth = ModelTester(combined)



scoresArr = np.zeros(3)
for i in range(100):
    webTech = combTesterWbTch.modelFitTest(forestModel, prnt=False, LOO=True)
    brainBio = combTesterBnBio.modelFitTest(forestModel, prnt=False, LOO=True)
    both = combTesterBoth.modelFitTest(forestModel, prnt=False, LOO=True)

    currentScores = [ webTech,
                      brainBio,
                      both
                      ]
    arr = np.array(currentScores)
    scoresArr += arr
    print( 'Iteration ', i )


scoresAv = scoresArr / 100
print('WebTech: ', scoresAv[0])
print('BrainBio: ', scoresAv[1])
print('Both: ', scoresAv[2])

sys.exit('Im done')
start = time.time()

scoresArr = np.zeros(12)
for i in range(100):
    RF = tester.modelFitTest(forestModel, prnt=False)
    #RFLOO = tester.modelFitTest(forestModel, LOO=True, prnt=False)
    NB = tester.modelFitTest(NBModel, dropList=remove, prnt=False)
    #NBLOO = tester.modelFitTest(NBModel, LOO=True, dropList=remove, prnt=False)
    DT = tester.modelFitTest(DecisionTreeClassifier(), prnt=False)
    #DTLOO = tester.modelFitTest(DecisionTreeClassifier(), LOO=True, prnt=False)
    Lin = tester.modelFitTest(LinearSVC(), prnt=False)
    #LinLOO = tester.modelFitTest(LinearSVC(), LOO=True, prnt=False)
    svc = tester.modelFitTest(SVC(), prnt=False)
    #svcLOO = tester.modelFitTest(SVC(), LOO=True, prnt=False)
    Log = tester.modelFitTest(LogisticRegression(), prnt=False)
    #LogLOO = tester.modelFitTest(LogisticRegression(), LOO=True, prnt=False)

    currentScores = [RF,# RFLOO,
                 NB,# NBLOO,
                 DT,# DTLOO,
                 Lin,# LinLOO,
                 svc,# svcLOO,
                 Log,# LogLOO
                 ]
    arr = np.array(currentScores)
    scoresArr += arr
scoresAv = scoresArr / 100

print(scoresAv)

end = time.time()
print(end-start)

print('RF: ', scoresAv[0])
print('RFLOO: ', scoresAv[1])
print('NB: ', scoresAv[2])
print('NBLOO: ', scoresAv[3])
print('DT: ', scoresAv[4])
print('DTLOO: ', scoresAv[5])
print('Lin: ', scoresAv[6])
print('LinLOO: ', scoresAv[7])
print('svc: ', scoresAv[8])
print('svcLOO: ', scoresAv[9])
print('Log: ', scoresAv[10])
print('LogLOO: ', scoresAv[11])



















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