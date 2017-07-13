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

# Returns a dataframe showing all mislabeled graphs, their actual collection, and which collection they were classified as
def get_mislabel_analysis(model, df, minSize=20, dropList=['Graph', 'Collection'], cv=5, prnt=True):


    collections = np.unique(df.Collection.values)
    for collection in collections:
        size = len(df[df.Collection == collection])
        if size < minSize:
            df = df[df.Collection != collection]

    X = df.drop(dropList, axis=1).values
    y = df['Collection'].values
    names = df['Graph'].values

    iterator = StratifiedKFold(n_splits=cv, shuffle=True,  # random_state = 42
                               )
    cv_pred = cross_val_predict(model, X, y, cv=iterator)
    cv_results_dict = {'Name' : names, 'Actual' : y, 'Predicted' : cv_pred}

    column_order = ['Name', 'Actual', 'Predicted']
    results = pd.DataFrame(cv_results_dict, columns = column_order, copy=True)

    if prnt == True:
        print('Using collections of size >', minSize)
        print('Excluding categories: ', dropList)
        print(classification_report(y, cv_pred))
        print(confusion_matrix(y, cv_pred))

    return results[ results.Actual != results.Predicted ]



# returns a dataframe consisting of the graphs that were mislabeled 'percent_return' - percent of the time after 'repeat' tests of the model
def get_mislabeled_graphs(model, df, repeat=100, percent_return=0.5, minSize=20, dropList=['Graph', 'Collection'], cv=5):
    allNameDict = {'CountTot' : np.zeros( len(df.Graph.values) ),
                   'CountWeb' : np.zeros( len(df.Graph.values) ),
                   'CountFb': np.zeros( len(df.Graph.values) ),
                   'CountSoc': np.zeros( len(df.Graph.values) ),
                   'CountRt': np.zeros( len(df.Graph.values) ),
                   'CountChem': np.zeros( len(df.Graph.values) ),
                   'CountBn': np.zeros( len(df.Graph.values) ),}
    order = ['CountTot', 'CountBn', 'CountChem', 'CountFb', 'CountRt', 'CountSoc', 'CountWeb']
    all_names = pd.DataFrame(allNameDict, columns=order, index=df.Graph.values, copy=True)

    for i in range( repeat ):
        analysis = get_mislabel_analysis(model, df, minSize, dropList, cv, prnt=False)
        for graph in analysis.Name.values:
            all_names.loc[graph, 'CountTot'] += 1

            if analysis[analysis.Name == graph].iloc[0, 2] == 'Brain Networks':
                all_names.loc[graph, 'CountBn'] += 1
            elif analysis[ analysis.Name == graph ].iloc[0,2] == 'Cheminformatics':
                all_names.loc[ graph, 'CountChem' ] += 1
            elif analysis[ analysis.Name == graph ].iloc[0,2] == 'Facebook Networks':
                all_names.loc[ graph, 'CountFb' ] += 1
            elif analysis[ analysis.Name == graph ].iloc[0,2] == 'Retweet Networks':
                all_names.loc[ graph, 'CountRt' ] += 1
            elif analysis[ analysis.Name == graph ].iloc[0,2] == 'Social Networks':
                all_names.loc[ graph, 'CountSoc' ] += 1
            elif analysis[ analysis.Name == graph ].iloc[0,2] == 'Web Graphs':
                all_names.loc[ graph, 'CountWeb' ] += 1


    percentOver = all_names[ all_names['CountTot'] >= percent_return * repeat ]
    return percentOver


Dtree50 = get_mislabeled_graphs(RandomForestClassifier(), df )
print( Dtree50, #DtreeAnalysis
 )