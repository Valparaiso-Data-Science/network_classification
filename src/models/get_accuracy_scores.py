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
           'Edges',
          # 'Density',
           'Maximum degree',
           'Minimum degree',
          # 'Average degree',
          # 'Assortativity',
            'Total triangles',
           'Average triangles',
           'Maximum triangles',
          # 'Avg. clustering coef.',
          # 'Frac. closed triangles',
          # 'Maximum k-core',
          # 'Max. clique (lb)'
            ]
tester = ModelTester(df)
forestModel = RandomForestClassifier(30)
NBModel = GaussianNB()



#**************************************************************
# For getting accuracy of 6 models
#**************************************************************
start = time.time()
# Currently changed to find std. dev.
scoresArr = np.array([[0]*12])
for i in range(100):
    print('iteration ', i)
    RF = tester.modelFitTest(forestModel, prnt=False)
    RFLOO = tester.modelFitTest(forestModel, LOO=True, prnt=False)
    NB = tester.modelFitTest(NBModel, dropList=remove, prnt=False)
    NBLOO = tester.modelFitTest(NBModel, LOO=True, dropList=remove, prnt=False)
    DT = tester.modelFitTest(DecisionTreeClassifier(), prnt=False)
    DTLOO = tester.modelFitTest(DecisionTreeClassifier(), LOO=True, prnt=False)
    Lin = tester.modelFitTest(LinearSVC(), prnt=False)
    LinLOO = tester.modelFitTest(LinearSVC(), LOO=True, prnt=False)
    svc = tester.modelFitTest(SVC(), prnt=False)
    svcLOO = tester.modelFitTest(SVC(), LOO=True, prnt=False)
    Log = tester.modelFitTest(LogisticRegression(), prnt=False)
    LogLOO = tester.modelFitTest(LogisticRegression(), LOO=True, prnt=False)

    currentScores = [RF, RFLOO,
                 NB, NBLOO,
                 DT, DTLOO,
                 Lin, LinLOO,
                 svc, svcLOO,
                 Log, LogLOO
                 ]
    arr = np.array([currentScores])
    scoresArr = np.concatenate((scoresArr, arr), axis=0)
    #print(scoresArr)
    #scoresArr += currentScores

    if i == 49:
        mid = time.time()
        print('half way time: ', mid-start)
#scoresAv = scoresArr / 100

#print(scoresAv)
scoresArr = np.delete(scoresArr, 0, 0)
scoresStdDev = np.std(scoresArr, axis=0)

print('RF: ', scoresStdDev[0])
print('RFLOO: ', scoresStdDev[1])
print('NB: ', scoresStdDev[2])
print('NBLOO: ', scoresStdDev[3])
print('DT: ', scoresStdDev[4])
print('DTLOO: ', scoresStdDev[5])
print('Lin: ', scoresStdDev[6])
print('LinLOO: ', scoresStdDev[7])
print('svc: ', scoresStdDev[8])
print('svcLOO: ', scoresStdDev[9])
print('Log: ', scoresStdDev[10])
print('LogLOO: ', scoresStdDev[11])

#print('RF: ', scoresAv[0])
#print('RFLOO: ', scoresAv[1])
#print('NB: ', scoresAv[2])
#print('NBLOO: ', scoresAv[3])
#print('DT: ', scoresAv[4])
#print('DTLOO: ', scoresAv[5])
#print('Lin: ', scoresAv[6])
#print('LinLOO: ', scoresAv[7])
#print('svc: ', scoresAv[8])
#print('svcLOO: ', scoresAv[9])
#print('Log: ', scoresAv[10])
#print('LogLOO: ', scoresAv[11])


end = time.time()
print('total time: ', end-start)


#*******************************************************************
# For LOO scores using different sized collections, different models
#*******************************************************************
start = time.time()
allScores = np.zeros(2)
for i in range(100):
    print('iteration ', i)
    #currentRF0 = tester.modelFitTest(forestModel, LOO=True, prnt=False, minSize=0 )
    #currentRF1 = tester.modelFitTest(forestModel, LOO=True, prnt=False, minSize=10)
    #currentRF2 = tester.modelFitTest(forestModel, LOO=True, prnt=False, minSize=20)
    currentGNB0 = tester.modelFitTest(NBModel, dropList=remove, LOO=True, prnt=False, minSize=0 )
    #currentGNB1 = tester.modelFitTest(NBModel, dropList=remove, LOO=True, prnt=False, minSize=10)
   # currentGNB2 = tester.modelFitTest(NBModel, dropList=remove, LOO=True, prnt=False, minSize=20)
    #currentDT0 = tester.modelFitTest(DecisionTreeClassifier(), LOO=True, prnt=False, minSize=0 )
    currentDT1 = tester.modelFitTest(DecisionTreeClassifier(), LOO=True, prnt=False, minSize=10)
    #currentDT2 = tester.modelFitTest(DecisionTreeClassifier(), LOO=True, prnt=False, minSize=20)
    #currentLin0 = tester.modelFitTest(LinearSVC(), LOO=True, prnt=False, minSize=0 )
    #currentLin1 = tester.modelFitTest(LinearSVC(), LOO=True, prnt=False, minSize=10)
    #currentLin2 = tester.modelFitTest(LinearSVC(), LOO=True, prnt=False, minSize=20)

    ar = np.array([#currentRF0,
          #currentRF1 ,
          #currentRF2 ,
          currentGNB0,
          #currentGNB1,
          #currentGNB2,
          #currentDT0 ,
          currentDT1 ,
          #currentDT2 ,
          #currentLin0,
          #currentLin1,
          #currentLin2
                     ])

    if i == 49:
        mid = time.time()
        print('halfway time: ', mid-start)

    allScores += ar
avScores = allScores / 100

print('RF')
#print('>0 -- ', avScores[0])
#print('>10 -- ', avScores[1])
#print('>20 -- ', avScores[2])
print('GNB')
print('>0 -- ', avScores[0])
#print('>10 -- ', avScores[4])
#print('>20 -- ', avScores[5])
print('DT')
#print('>0 -- ', avScores[6])
print('>10 -- ', avScores[1])
#print('>20 -- ', avScores[8])
print('LinSVC')
#print('>0 -- ', avScores[9])
#print('>10 -- ', avScores[10])
#print('>20 -- ', avScores[11])

end = time.time()
print('total time: ', end-start)
sys.exit()
#*************************************************************
# For getting f1 score among collections for different models
#*************************************************************
start = time.time()

scoresRF = np.zeros(6)
scoresGNB = np.zeros(6)
scoresDT = np.zeros(6)
scoresLin = np.zeros(6)
for i in range(100):
    print('iteration ', i)

   # currentRF = tester.modelFitTest(forestModel, f1Score=True, LOO=True, prnt=False)
   # currentGNB = tester.modelFitTest(NBModel, dropList=remove, LOO=True, f1Score=True, prnt=False)
   # currentDT = tester.modelFitTest(DecisionTreeClassifier(), LOO=True, f1Score=True, prnt=False)
    currentLin = tester.modelFitTest(LinearSVC(), LOO=True, f1Score=True, prnt=False)

   # scoresRF += np.array(currentRF)
   # scoresGNB += np.array(currentGNB)
   # scoresDT += np.array(currentDT)
    scoresLin += np.array(currentLin)


averageRF = scoresRF/100
averageGNB = scoresGNB/100
averageDT = scoresDT/100
averageLin = scoresLin/100

#print('Random Forest')
#print('Brain: ', averageRF[0])
#print('Chem: ', averageRF[1])
#print('Facebook: ', averageRF[2])
#print('Retweet: ', averageRF[3])
#print('Social: ', averageRF[4])
#print('Web: ', averageRF[5])
#print('Gaussian Naive Bayes')
print('Brain: ', averageGNB[0])
print('Chem: ', averageGNB[1])
print('Facebook: ', averageGNB[2])
print('Retweet: ', averageGNB[3])
print('Social: ', averageGNB[4])
print('Web: ', averageGNB[5])

#print('Decision Tree')
#print('Brain: ', averageDT[0])
#print('Chem: ', averageDT[1])
#print('Facebook: ', averageDT[2])
#print('Retweet: ', averageDT[3])
#print('Social: ', averageDT[4])
#print('Web: ', averageDT[5])

print('Linear SVC')
print('Brain: ', averageLin[0])
print('Chem: ', averageLin[1])
print('Facebook: ', averageLin[2])
print('Retweet: ', averageLin[3])
print('Social: ', averageLin[4])
print('Web: ', averageLin[5])

end = time.time()
print('time: ', end-start)

sys.exit("I'm done")


#***********************************************
# For getting accuracy of combined collections
#***********************************************

combined = tester.combine_collections(['Brain Networks', 'Biological Networks'], 'Brain-Bio')
combTesterBnBio = ModelTester(combined)

combined = tester.combine_collections(['Web Graphs', 'Technological Networks'], 'Web-Tech')
combTesterWbTch = ModelTester(combined)

combined = combTesterBnBio.combine_collections(['Web Graphs', 'Technological Networks'], 'Web-Tech')
combTesterBoth = ModelTester(combined)

scoresArr = np.zeros(3)
for i in range(100):
    webTech = combTesterWbTch.modelFitTest(NBModel, dropList=remove, prnt=False, LOO=True)
    brainBio = combTesterBnBio.modelFitTest(NBModel, dropList=remove, prnt=False, LOO=True)
    both = combTesterBoth.modelFitTest(NBModel, dropList=remove, prnt=False, LOO=True)

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