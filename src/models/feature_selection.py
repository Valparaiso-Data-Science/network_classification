import pandas as pd
import matplotlib.pyplot as plt
from sklearn.feature_selection import SelectFromModel, SelectKBest, RFE, RFECV
from sklearn.linear_model import LogisticRegression, Lasso, LinearRegression, RandomizedLasso, RandomizedLogisticRegression
from sklearn.svm import LinearSVC
from sklearn.ensemble import RandomForestClassifier
from sklearn import tree
from sklearn.model_selection import StratifiedKFold
import numpy as np

# Read file
df = pd.read_csv('~/Downloads/network_classification/src/data/data_minmaxscale.csv', index_col='Unnamed: 0')

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
num_of_features = 1

# RFE method
# Create various models to compare
#*******************************************
# Create Logistic Regression Model with RFE
#*******************************************
model_1 = LogisticRegression()
rfe_1 = RFE(model_1, num_of_features)
fit_1 = rfe_1.fit(X, Y)
print("\nLogistic Regression")
print("Num Features: " + str(fit_1.n_features_))
print("Selected Features: " + str(fit_1.support_))
print("Feature Ranking: " + str(fit_1.ranking_))

#*******************************************
# Create Lasso Regression Model with RFE
#*******************************************
model_2 = Lasso(alpha = 0.1)
rfe_2 = RFE(model_2, num_of_features)
fit_2 = rfe_2.fit(X, Y)
print("\nLasso Regression")
print("Num Features: " + str(fit_2.n_features_))
print("Selected Features: " + str(fit_2.support_))
print("Feature Ranking: " + str(fit_2.ranking_))

#*******************************************
# Create SVM Model with RFE
#*******************************************
#use svc linear
model_3 = LinearSVC()
rfe_3 = RFE(model_3, num_of_features)
fit_3 = rfe_3.fit(X, Y)
print("\nLinear SVC")
print("Num Features: " + str(fit_3.n_features_))
print("Selected Features: " + str(fit_3.support_))
print("Feature Ranking: " + str(fit_3.ranking_))

#*******************************************
# Create Linear Regression Model with RFE
#*******************************************
model_4 = LinearRegression()
rfe_4 = RFE(model_4, num_of_features)
fit_4 = rfe_4.fit(X, Y)
print("\nLinear Regression")
print("Num Features: " + str(fit_4.n_features_))
print("Selected Features: " + str(fit_4.support_))
print("Feature Ranking: " + str(fit_4.ranking_))

#*******************************************
# Create Decision Tree Model with RFE
#*******************************************
model_5 = tree.DecisionTreeClassifier(random_state=42)
rfe_5 = RFE(model_5, num_of_features)
fit_5 = rfe_5.fit(X, Y)
print("\nDecision Tree")
print("Num Features: " + str(fit_5.n_features_))
print("Selected Features: " + str(fit_5.support_))
print("Feature Ranking: " + str(fit_5.ranking_))

#*******************************************
# Create Random Forest Model with RFE
#*******************************************
model_6 = RandomForestClassifier(random_state=42)
rfe_6 = RFE(model_6, num_of_features)
fit_6 = rfe_6.fit(X, Y)
print("\nRandom Forest")
print("Num Features: " + str(fit_5.n_features_))
print("Selected Features: " + str(fit_5.support_))
print("Feature Ranking: " + str(fit_5.ranking_))


# Print all the keys
model_list = ["Logistic Regression", "Lasso Regression", "Linear SVC", "Linear Regression",
              "Decision Tree", "Random Forest"]
fit_list = [fit_1.ranking_, fit_2.ranking_, fit_3.ranking_, fit_4.ranking_, fit_5.ranking_, fit_6.ranking_,
            ]
for i, fit in enumerate(fit_list):
    print("\nMost Important Features in " + model_list[i] + ": ")
    for i, val in enumerate(fit):
        if val == True:
            print(col_names[i+1])


# Stability Selection Method
# *******************************************
#  Create Randomized Lasso Regression Model with RFE
# *******************************************
model_7 = RandomizedLasso(alpha=0.01)
model_7.fit(X, Y)
print("\nRandomized Lasso Regression")
print("Features sorted by their score:")
print(sorted(zip(map(lambda x: round(x, 4), model_7.scores_), col_names), reverse=True))

# *******************************************
#  Create Randomized Logistic Regression Model with RFE
# *******************************************
model_8 = RandomizedLogisticRegression()
model_8.fit(X, Y)
print("\nRandomized Logistic Regression")
print("Features sorted by their score:")
print(sorted(zip(map(lambda x: round(x, 4), model_8.scores_), col_names), reverse=True))


#**************************************
# RFECV feature selection
#**************************************
# Can change the model and cv
rfecv = RFECV(estimator = model_6, step = 1, cv = StratifiedKFold(n_splits=5, shuffle = True, random_state=42),
              scoring = 'accuracy')
rfecv.fit(X, Y)
print("Optimal number of features : %d" % rfecv.n_features_)
#print("Selected Features: " + str(rfecv.support_))
print("Feature Ranking: " + str(rfecv.ranking_))

# Get grid scores
g_scores = rfecv.grid_scores_
indices = np.argsort(g_scores)[::-1]
print('Printing RFECV results:')
for f in range(X.shape[1]):
    print("%d. Number of features: %d, Grid_Score: %f" % (f + 1, indices[f]+1, g_scores[indices[f]]))

# Plot number of features VS. cross-validation scores
plt.figure()
plt.title("RFECV - Random Forest", {'size':'16'})
plt.xlabel("Number of features selected", {'size':'11'})
plt.ylabel("Cross validation score (nb of correct classifications)", {'size':'11'})
plt.plot(range(1, len(rfecv.grid_scores_) + 1), rfecv.grid_scores_)
plt.show()
