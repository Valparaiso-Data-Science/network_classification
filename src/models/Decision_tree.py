import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import randint
from sklearn.tree import DecisionTreeClassifier
from sklearn.model_selection import RandomizedSearchCV, train_test_split


df = pd.read_csv('C:/Users/Owner/Documents/VERUM/Network stuff/git/src/data/clean_data.csv', index_col=0)
X = df.drop(['Graph', 'Collection'], axis=1).values
y = df['Collection'].values
# Setup the parameters and distributions to sample from: param_dist
param_dist = {"max_depth": [3, None],
              "max_features": randint(1, 9),
              "min_samples_leaf": randint(1, 9),
              "criterion": ["gini", "entropy"]}

# Instantiate a Decision Tree classifier: tree
tree = DecisionTreeClassifier()

# Instantiate the RandomizedSearchCV object: tree_cv
tree_cv = RandomizedSearchCV(tree, param_dist, cv= 5)

# Fit it to the data
tree_cv.fit(X, y)

# Print the tuned parameters and score
print("Tuned Decision Tree Parameters: {}".format(tree_cv.best_params_))
print("Best score is {}".format(tree_cv.best_score_))

tree2 = DecisionTreeClassifier( max_features=7, min_samples_leaf=2)
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=.2)
tree2.fit(X_train, y_train)

print( tree2.decision_path(X_train))