import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split, cross_val_score, cross_val_predict, StratifiedKFold
from sklearn.svm import SVC, LinearSVC
from sklearn.metrics import classification_report, confusion_matrix
from sklearn.tree import DecisionTreeClassifier
from sklearn.naive_bayes import GaussianNB
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier

infile = 'C:/Users/Owner/Documents/VERUM/Network stuff/git/src/data/data_minmaxscale.csv' # -- change for machine
df = pd.read_csv(infile, index_col=0)

techGraphs = df[ df.Collection == 'Technological Networks' ]
webGraphs = df[ df.Collection == 'Web Graphs' ]
web_tech = pd.concat([webGraphs, techGraphs])
dfNew = df.copy()
for i in web_tech.index:
    dfNew.loc[i, 'Collection'] = 'Web-Tech'


from git.src.models.model_class import modelFitTest

#modelFitTest(RandomForestClassifier(), dfNew, cv=5)

from git.src.models.classification_analysis import get_mislabel_analysis, get_mislabeled_graphs

get_mislabel_analysis(RandomForestClassifier(), dfNew)
print(get_mislabel_analysis(RandomForestClassifier(), dfNew))