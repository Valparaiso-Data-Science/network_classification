import pandas as pd
import numpy as np
from sklearn.feature_selection import SelectFromModel, SelectKBest, RFE
from sklearn.linear_model import LogisticRegression, Lasso
from sklearn.svm import SVC
import matplotlib.pyplot as plt

# Read file
df = pd.read_csv('~/Downloads/network_classification/src/data/clean_data_with_new_chem.csv', index_col='Unnamed: 0')

# Delete categorical column
del df['Graph']

# Change the collection names to numbers
df['Collection'] = df['Collection'].astype('category')
df['Collection'] = df['Collection'].cat.codes


# Create array of the values and define X and Y
df_array = df.values
X = df_array[:, 0:14]
Y = df_array[:, 0]

# Number of important features we want to extract
num_of_features = 6

#*******************************************
# Create Logistic Regression Model with RFE
#*******************************************
model1 = LogisticRegression()
rfe1 = RFE(model1, num_of_features)
fit1 = rfe1.fit(X, Y)
print("Num Features: " + str(fit1.n_features_))
print("Selected Features: " + str(fit1.support_))
print("Feature Ranking: " + str(fit1.ranking_))

#*******************************************
# Create Lasso Regression Model with RFE
#*******************************************
model2 = Lasso(alpha = 0.1)
rfe2 = RFE(model2, num_of_features)
fit2 = rfe2.fit(X, Y)
print("Num Features: " + str(fit2.n_features_))
print("Selected Features: " + str(fit2.support_))
print("Feature Ranking: " + str(fit2.ranking_))

# Print all the keys
print(df.keys())

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