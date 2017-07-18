import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split, cross_val_score, StratifiedKFold
from sklearn.svm import SVC, LinearSVC
from sklearn.metrics import classification_report, confusion_matrix
from sklearn.tree import DecisionTreeClassifier
from sklearn.naive_bayes import GaussianNB
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier

#from git.src.models.model_class import modelFitTest
from git.src.models.model_class import ModelTester

infile = 'C:/Users/Owner/Documents/VERUM/Network stuff/git/src/data/clean_data_with_new_chem.csv' # -- change for machine
df = pd.read_csv(infile, index_col=0)

#collections to remove from model
remove = ['Graph', 'Collection',
          # 'Nodes',
          # 'Edges',
          # 'Density',
          # 'Maximum degree',
          # 'Minimum degree',
          # 'Average degree',
          # 'Assortativity',
          #  'Total triangles',
          # 'Average triangles',
          # 'Maximum triangles',
          # 'Avg. clustering coef.',
          # 'Frac. closed triangles',
          # 'Maximum k-core',
          # 'Max. clique (lb)'
            ]

tester = ModelTester(df)
#print(tester.get_mislabel_analysis(RandomForestClassifier(), dropList=remove))


combined = tester.combine_collections(['Web Graphs', 'Technological Networks'], 'Web-Tech')
combTester = ModelTester(combined)
combined = combTester.combine_collections(['Brain Networks', 'Biological Networks'], 'Brain-Bio')
combTester = ModelTester(combined)
#combined = combTester.combine_collections([#'Collaboration Networks',
#                                           'Interaction Networks',
#                                           'Recommendation Networks'
#                                           ], 'Collab_Inter_Rec')
#combTester = ModelTester(combined)
combined = combTester.combine_collections(['Ecology Networks', 'Scientific Computing'], 'Sci-Eco')
combTester = ModelTester(combined)

combMislabel = combTester.get_mislabel_analysis(RandomForestClassifier(), minSize=16)
print(combMislabel
      [ combMislabel.Actual == 'Sci-Eco']
     )
combTester.modelFitTest(RandomForestClassifier(), minSize=16)
tester.modelFitTest(RandomForestClassifier(), minSize=0)

miscFile = 'C:/Users/Owner/Documents/VERUM/Network stuff/git/src/data/miscellaneous_networks.csv'
misc = pd.read_csv(miscFile, index_col=0)
miscObject = ModelTester(misc)

renamedMisc = miscObject.combine_collections(['Web Graphs', 'Technological Networks'], 'Web-Tech')
miscObject = ModelTester(renamedMisc)
renamedMisc = miscObject.combine_collections(['Brain Networks', 'Biological Networks'], 'Brain-Bio')
miscObject = ModelTester(renamedMisc)
renamedMisc = miscObject.combine_collections(['Ecology Networks', 'Scientific Computing'], 'Sci-Eco')



miscTest = combTester.train_predict(RandomForestClassifier(), renamedMisc, minSize=0)
print(miscTest)
print(miscTest[ miscTest.Hypothesis == miscTest.Predicted])






#print('Decision Tree')
#modelFitTest(DecisionTreeClassifier(), df, cv=5, dropList=remove, feat_comp=False )
##print('Decision Tree')
##modelFitTest(DecisionTreeClassifier(), df, cv=5, split=.3, minSize=10, dropList=remove, feat_comp=False)
##print('Decision Tree')
##modelFitTest(DecisionTreeClassifier(), df, cv=5, split=.3, minSize=20, dropList=remove, feat_comp=False)

#print('SVC')
#modelFitTest(SVC(), df, cv=5, split=.3, dropList=remove, feat_comp=False)
#print('SVC')
#modelFitTest(SVC(), df, cv=5, split=.3, minSize=10, dropList=remove, feat_comp=False)
#print('SVC')
#modelFitTest(SVC(), df, cv=5, split=.3, minSize=20, dropList=remove, feat_comp=False)

#print('Linear SVC')
#modelFitTest(LinearSVC(), df, cv=5, dropList=remove, feat_comp=False )
#print('Linear SVC')
#modelFitTest(LinearSVC(), df,cv=5, split=.3, minSize=10, dropList=remove, feat_comp=False)
#print('Linear SVC')
#modelFitTest(LinearSVC(), df,cv=5, split=.3, minSize=20, dropList=remove, feat_comp=False)

#print('Gaussian Naive Bayes')
#modelFitTest(GaussianNB(), df, cv=5, dropList=remove, feat_comp=False )
#print('Gaussian Naive Bayes')
#modelFitTest(GaussianNB(), df, cv=5, split=.3, minSize=10, dropList=remove, feat_comp=False)
#print('Gaussian Naive Bayes')
#modelFitTest(GaussianNB(), df, cv=5, split=.3, minSize=20, dropList=remove, feat_comp=False)

#print('Logistic Regression')
#modelFitTest(LogisticRegression(), df, cv=5, split=.3, dropList=remove, feat_comp=False)
#print('Logistic Regression')
#modelFitTest(LogisticRegression(), df, cv=5, split=.3, minSize=10, dropList=remove, feat_comp=False)
#print('Logistic Regression')
#modelFitTest(LogisticRegression(), df, cv=5, split=.3, minSize=20, dropList=remove, feat_comp=False)

#print('RandomForestClassifier')
#modelFitTest(RandomForestClassifier(), df, cv=5, dropList=remove, feat_comp=False )
#print('RandomForestClassifier')
#modelFitTest(RandomForestClassifier(), df, cv=5, split=.3, minSize=10, dropList=remove, feat_comp=False)
#print('RandomForestClassifier')
#modelFitTest(RandomForestClassifier(), df, cv=5, split=.3, minSize=20, dropList=remove, feat_comp=False)


