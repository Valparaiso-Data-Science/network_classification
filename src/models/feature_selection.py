import pandas as pd
import numpy as np
from sklearn.feature_selection import SelectFromModel, SelectKBest, RFE
from sklearn.linear_model import LogisticRegression, Lasso
from sklearn.svm import LinearSVC

# Read file
df = pd.read_csv('~/Downloads/network_classification/src/data/clean_data_with_new_chem.csv', index_col='Unnamed: 0')

# Delete categorical column
del df['Graph']

# Change the collection names to numbers
df['Collection'] = df['Collection'].astype('category')
df['Collection'] = df['Collection'].cat.codes

# Get all column names
col_names = df.keys()

# Create array of the values and define X and Y
df_array = df.values
X = df_array[:, 0:14]
Y = df_array[:, 0]

# Number of important features we want to extract
num_of_features = 6

# Create various models to compare
#*******************************************
# Create Logistic Regression Model with RFE
#*******************************************
model_1 = LogisticRegression()
rfe_1 = RFE(model_1, num_of_features)
fit_1 = rfe_1.fit(X, Y)
print("Logistic Regression")
print("Num Features: " + str(fit_1.n_features_))
print("Selected Features: " + str(fit_1.support_))
print("Feature Ranking: " + str(fit_1.ranking_))

#*******************************************
# Create Lasso Regression Model with RFE
#*******************************************
model_2 = Lasso(alpha = 0.1)
rfe_2 = RFE(model_2, num_of_features)
fit_2 = rfe_2.fit(X, Y)
print("Lasso Regression")
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
print("Linear SVC")
print("Num Features: " + str(fit_3.n_features_))
print("Selected Features: " + str(fit_3.support_))
print("Feature Ranking: " + str(fit_3.ranking_))



# Print all the keys
fit_list = [fit_1.ranking_, fit_2.ranking_, fit_3.ranking_]
for i, fit in enumerate(fit_list):
    print("Fit " + str(i+1))
    for i, val in enumerate(fit):
        if val == True:
            print(col_names[i])



# Temporarily here; probably will delete this chunk of code later
#******************
# Select K Best
#******************
#test = SelectKBest(score_func=chi2, k=10)
#fit = test.fit(X, Y)
# summarize scores
#np.set_printoptions(precision=3)
#print(fit.scores_)
#features = fit.transform(X)
# summarize selected features
#print(features[0:5,:])