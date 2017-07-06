import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split, cross_val_score, StratifiedKFold
from sklearn.svm import SVC, LinearSVC
from sklearn.metrics import classification_report, confusion_matrix
from sklearn.tree import DecisionTreeClassifier
from sklearn.naive_bayes import GaussianNB
from sklearn.linear_model import LogisticRegression

from git.src.models.model_class import modelFitTest

infile = 'C:/Users/Owner/Documents/VERUM/Network stuff/git/src/data/data_minmaxscale.csv' # -- change for machine
df = pd.read_csv(infile, index_col=0)

#collections to remove from model
remove = ['Graph', 'Collection',
            #'Nodes',
           'Edges',
            #'Density', 'Maximum degree', 'Minimum degree', 'Average degree',
             #'Assortativity',
             #'Total triangles', 'Average triangles', 'Maximum triangles', 'Avg. clustering coef.',
             #'Frac. closed triangles', 'Maximum k-core',
             #'Max. clique (lb)'
            ]


print('Decision Tree')
modelFitTest(DecisionTreeClassifier(), df, cv=5, split=.3, dropList=remove, feat_comp=True)
print('Decision Tree')
modelFitTest(DecisionTreeClassifier(), df, cv=5, split=.3, minSize=10, dropList=remove, feat_comp=True)
print('Decision Tree')
modelFitTest(DecisionTreeClassifier(), df, cv=5, split=.3, minSize=20, dropList=remove, feat_comp=True)

print('SVC')
modelFitTest(SVC(), df, cv=5, split=.3, dropList=remove, feat_comp=True)
print('SVC')
modelFitTest(SVC(), df, cv=5, split=.3, minSize=10, dropList=remove, feat_comp=True)
print('SVC')
modelFitTest(SVC(), df, cv=5, split=.3, minSize=20, dropList=remove, feat_comp=True)

print('Linear SVC')
modelFitTest(LinearSVC(), df,cv=5, split=.3, dropList=remove, feat_comp=True)
print('Linear SVC')
modelFitTest(LinearSVC(), df,cv=5, split=.3, minSize=10, dropList=remove, feat_comp=True)
print('Linear SVC')
modelFitTest(LinearSVC(), df,cv=5, split=.3, minSize=20, dropList=remove, feat_comp=True)

print('Gaussian Naive Bayes')
modelFitTest(GaussianNB(), df, cv=5, split=.3, dropList=remove, feat_comp=True)
print('Gaussian Naive Bayes')
modelFitTest(GaussianNB(), df, cv=5, split=.3, minSize=10, dropList=remove, feat_comp=True)
print('Gaussian Naive Bayes')
modelFitTest(GaussianNB(), df, cv=5, split=.3, minSize=20, dropList=remove, feat_comp=True)

print('Logistic Regression')
modelFitTest(LogisticRegression(), df, cv=5, split=.3, dropList=remove, feat_comp=True)
print('Logistic Regression')
modelFitTest(LogisticRegression(), df, cv=5, split=.3, minSize=10, dropList=remove, feat_comp=True)
print('Logistic Regression')
modelFitTest(LogisticRegression(), df, cv=5, split=.3, minSize=20, dropList=remove, feat_comp=True)