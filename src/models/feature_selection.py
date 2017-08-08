# This code uses feature selection algorithms with differente models
# to obtain the most important features of our data
# Algorithms used are stability selection, rfe and rfecv

# Import relevant packages
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.feature_selection import RFE, RFECV
from sklearn.linear_model import LogisticRegression, Lasso, LinearRegression, RandomizedLasso, RandomizedLogisticRegression
from sklearn.svm import LinearSVC
from sklearn.ensemble import RandomForestClassifier
from sklearn import tree
from sklearn.model_selection import StratifiedKFold
import numpy as np


# Read file
df = pd.read_csv('~/Downloads/network_classification/src/data/data_greater_than_20.csv', index_col='Unnamed: 0')

# Delete categorical column
del df['Graph']


# Change the collection names to numbers
df['Collection'] = df['Collection'].astype('category')
df['Collection'] = df['Collection'].cat.codes

# Get all column names
col_names = df.keys()

# Create array of the values and define X and Y
df_array = df.values
X = df_array[:, 1:15]
Y = df_array[:, 0] #target variable (collection)

# Number of important features we want to extract
# Using 1 will give us the ranking of all of the features
num_of_features = 1



#**************************************
# RFE feature selection function
#**************************************

# Define function that runs RFE in different models
def rfe_function(model, model_name, nu
    m_of_features, X, Y):
    rfe = RFE(model, num_of_features)
    fit = rfe.fit(X,Y)
    print("\n" + model_name)
    print("Num Features: " + str(fit.n_features_))
    print("Selected Features: " + str(fit.support_))
    print("Feature Ranking: " + str(fit.ranking_))
    # Print the most important feature
    print("\nMost Important Feature in " + model_name + ": ")
    for i, val in enumerate(fit.ranking_):
        if val == True:
            print(col_names[i+1])


# Calling function to run RFE with different models
rfe_function(LogisticRegression(), "Logistic Regression", num_of_features, X, Y)
rfe_function(Lasso(alpha = 0.1), "Lasso Regression", num_of_features, X, Y)
rfe_function(LinearSVC(), "Linear SVC", num_of_features, X, Y)
rfe_function(LinearRegression(), "Linear Regression", num_of_features, X, Y)
rfe_function(tree.DecisionTreeClassifier(random_state=42), "Decision Tree", num_of_features, X, Y)
rfe_function(RandomForestClassifier(random_state=42), "Random Forest", num_of_features, X, Y)



# *******************************************
# Stability Selection Method
# *******************************************

def stability_select(model, model_name, X, Y)
    model.fit(X,Y)
    print(model_name)
    print("Features sorted by their score:")
    print(sorted(zip(map(lambda x: round(x, 4), model_1.scores_), col_names), reverse=True))


stability_select(RandomizedLasso(alpha=0.01), "\nRandomized Lasso Regression", X, Y)
stability_select(RandomizedLogisticRegression(), "\nRandomized Logistic Regression", X, Y)

#**************************************
# RFECV feature selection
#**************************************
list = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
df_features = pd.DataFrame(index=pd.Series(range(14)))
iterations = 100
total_num_features = 14

def rfecv(model, name):
    for j in range(iterations):
        rfecv = RFECV(estimator=model, step=1, cv=StratifiedKFold(n_splits=5, shuffle=True),
                      scoring='accuracy')
        rfecv.fit(X, Y)
        # Uncomment this to see the ranking of the features in rfecv
        #print("Optimal number of features : %d" % rfecv.n_features_)
        #print("Selected Features: " + str(rfecv.support_))
        #print("Feature Ranking: " + str(rfecv.ranking_))

        # Get grid scores
        g_scores = rfecv.grid_scores_
        indices = np.argsort(g_scores)[::-1]
        print('Printing RFECV results:')

        for f in range(X.shape[1]):
            # print("%d. Number of features: %d, Grid_Score: %f" % (f + 1, indices[f] + 1, g_scores[indices[f]]))
            for i in range(total_num_features):
                if indices[f]  == i :
                    list[i] = list[i] + g_scores[indices[f]]
        #print('List: ' + str(list))

    df_features['scores'] = list
    print(df_features.columns)


    for m in range(total_num_features):
        list[m] = list[m]/iterations

    print('New List: ' + str(list))



    # Plot number of features VS. cross-validation scores
    plt.figure()
    plt.title(name, {'size': '22'})
    plt.xlabel("Number of features selected", {'size': '18'})
    plt.ylabel("Cross validation score \n(nb of correct classifications)", {'size': '18'})
    plt.plot(range(1, len(list) + 1), list)
    plt.show()

    return list



# Calling RFECV function for all three models
rfecv(LinearSVC(), "RFECV - Linear SVC")
rfecv(tree.DecisionTreeClassifier(), "RFECV - Decision Tree")
rfecv(RandomForestClassifier(), "RFECV - Random Forest")


# Used to save rfecv data
# Change file or path name as desired
list = pd.DataFrame(list)
#list.to_csv('~/Downloads/network_classification/src/data/rfecv_random_tree_data.csv')