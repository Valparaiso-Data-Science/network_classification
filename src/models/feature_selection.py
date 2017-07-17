# This code uses feature selection algorithms with differente models
# to obtain the most important features of our data


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
Y = df_array[:, 0]

# Number of important features we want to extract
# Using 1 will give us the ranking of all of the features
num_of_features = 1



#**************************************
# RFE feature selection function
#**************************************

# Define function that runs RFE in different models
def rfe_function(model, model_name, num_of_features, X, Y):
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

#  Create Randomized Lasso Regression Model with RFE
model_1 = RandomizedLasso(alpha=0.01)
model_1.fit(X, Y)
print("\nRandomized Lasso Regression")
print("Features sorted by their score:")
print(sorted(zip(map(lambda x: round(x, 4), model_1.scores_), col_names), reverse=True))


#  Create Randomized Logistic Regression Model with RFE
model_2 = RandomizedLogisticRegression()
model_2.fit(X, Y)
print("\nRandomized Logistic Regression")
print("Features sorted by their score:")
print(sorted(zip(map(lambda x: round(x, 4), model_2.scores_), col_names), reverse=True))



#**************************************
# RFECV feature selection
#**************************************

def rfecv(model, name):
    rfecv = RFECV(estimator=model, step=1, cv=StratifiedKFold(n_splits=5, shuffle=True, random_state=42),
                  scoring='accuracy')
    rfecv.fit(X, Y)
    print("Optimal number of features : %d" % rfecv.n_features_)
    # print("Selected Features: " + str(rfecv.support_))
    print("Feature Ranking: " + str(rfecv.ranking_))

    # Get grid scores
    g_scores = rfecv.grid_scores_
    indices = np.argsort(g_scores)[::-1]
    print('Printing RFECV results:')
    for f in range(X.shape[1]):
        print("%d. Number of features: %d, Grid_Score: %f" % (f + 1, indices[f] + 1, g_scores[indices[f]]))

    # Plot number of features VS. cross-validation scores
    plt.figure()
    plt.title(name, {'size': '16'})
    plt.xlabel("Number of features selected", {'size': '11'})
    plt.ylabel("Cross validation score (nb of correct classifications)", {'size': '11'})
    plt.plot(range(1, len(rfecv.grid_scores_) + 1), rfecv.grid_scores_)
    plt.show()


# Calling RFECV function for all three models
rfecv(model_3, "RFECV - Linear SVC")
rfecv(model_5, "RFECV - Decision Tree")
rfecv(model_6, "RFECV - Random Forest")