''' This script was mostly used to test out the model class, and get mislabel analyses. There is not much organization, and things were constantly being added and changed'''

import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split, cross_val_score, StratifiedKFold
from sklearn.svm import SVC, LinearSVC
from sklearn.metrics import classification_report, confusion_matrix
from sklearn.tree import DecisionTreeClassifier
from sklearn.naive_bayes import GaussianNB
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
import sys

#from git.src.models.model_class import modelFitTest
from git.src.models.model_class import ModelTester

infile = 'C:/Users/Owner/Documents/VERUM/Network stuff/git/src/data/clean_data_with_new_chem.csv' # -- change for machine
df = pd.read_csv(infile, index_col=0)

#remove Interaction Networks
df = df[ df.Collection != 'Interaction Networks' ]
df = df[ df.Collection != 'Collaboration Networks' ]

#collections to remove from model
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

forestModel = RandomForestClassifier(n_estimators=30)
GNBmodel = GaussianNB()

tester = ModelTester(df)

#tester.modelFitTest(forestModel, LOO=True)
tester.modelFitTest(GNBmodel, LOO=True, minSize=10)
tester.modelFitTest(GNBmodel, LOO=True, minSize=0)

#print(tester.modelFitTest(forestModel, f1Score=True, prnt=False))

#analysisNB = tester.get_mislabel_analysis(GNBmodel, dropList=remove, LOO=True, proba=True)
#print(analysisNB)
#mislabelRF = tester.get_mislabeled_graphs(forestModel, LOO=True)
#mislabelRF.to_csv('C:/Users/Owner/Documents/VERUM/Network stuff/git/src/models/RandomForest_mislabeled_LOOcv.csv')

#mislabelGNB = tester.get_mislabeled_graphs(GNBmodel, dropList=remove, LOO=True)
#mislabelGNB.to_csv('C:/Users/Owner/Documents/VERUM/Network stuff/git/src/models/GaussianNB_mislabeled_LOOcv.csv')


#analysisNB.to_csv('C:/Users/Owner/Documents/VERUM/Network stuff/git/src/models/GNB_vs_RF_mislabel_comparison.csv')
sys.exit("I'm done'")

combined = tester.combine_collections(['Brain Networks', 'Biological Networks'], 'Brain-Bio')
combTester = ModelTester(combined)

combined = combTester.combine_collections(['Web Graphs', 'Technological Networks'], 'Web-Tech')
combTester = ModelTester(combined)

#combTester.modelFitTest(forestModel, minSize=16, LOO=True)

#combined = combTester.combine_collections(['Ecology Networks', 'Scientific Computing'], 'Sci-Eco')
#combTester = ModelTester(combined)

#combTester.modelFitTest(forestModel, minSize=16, LOO=True)
#tester.modelFitTest(forestModel, minSize=16, LOO=True)













# Start to experiment with misc graphs
miscFile = 'C:/Users/Owner/Documents/VERUM/Network stuff/git/src/data/reduced_miscellaneous_networks.csv'
misc = pd.read_csv(miscFile, index_col=0)
miscObject = ModelTester(misc)

renamedMisc = miscObject.combine_collections(['Web Graphs', 'Technological Networks'], 'Web-Tech')
miscObject = ModelTester(renamedMisc)
renamedMisc = miscObject.combine_collections(['Brain Networks', 'Biological Networks'], 'Brain-Bio')
#miscObject = ModelTester(renamedMisc)
#renamedMisc = miscObject.combine_collections(['Ecology Networks', 'Scientific Computing'], 'Sci-Eco')

#miscComb = combTester.train_predict(forestModel, renamedMisc)




miscScore = []
for i in range(100):
    miscTest = combTester.train_predict(forestModel, renamedMisc)
    correct = len(miscTest[miscTest.Hypothesis == miscTest.Predicted].Name.values)
    score = correct / 32
    miscScore.append(score)

print(miscTest[miscTest.Hypothesis == miscTest.Predicted])
print(miscScore)
print('average score: ', np.mean(miscScore))

#miscTest = combTester.train_predict(forestModel, renamedMisc, minSize=16)
#print(miscTest)
#print(miscTest.info())
#print(miscTest[ miscTest.Hypothesis == miscTest.Predicted])
#miscAnalysis = combTester.get_mislabeled_graphs(forestModel, externalData=renamedMisc, minSize=16)
#print(miscAnalysis)


infile_er = 'C:/Users/Owner/Documents/VERUM/Network stuff/git/src/data/synthetic_e1e2e3_complete.csv'
df_er = pd.read_csv(infile_er, index_col=0)

dfSyn = pd.concat([df, df_er])
testerSyn = ModelTester(dfSyn)
tester.modelFitTest(forestModel)
tester.modelFitTest(GNBmodel, dropList=remove)
testerSyn.modelFitTest(forestModel, LOO=False)
testerSyn.modelFitTest(GNBmodel, dropList=remove)
#print(tester.get_mislabeled_graphs(GNBmodel, dropList=remove))
#print(tester.get_mislabeled_graphs(forestModel))




