import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import randint
from sklearn.tree import DecisionTreeClassifier, export_graphviz
from sklearn.model_selection import RandomizedSearchCV, train_test_split


#df = pd.read_csv('C:/Users/Owner/Documents/VERUM/Network stuff/git/src/data/clean_data.csv', index_col=0)
df = pd.read_csv('C:/Users/Owner/Documents/VERUM/Network stuff/git/src/data/clean_data_with_new_chem.csv', index_col=0) #Data with down-sampled Cheminfo
X = df.drop(['Graph', 'Collection'], axis=1).values
y = df['Collection'].values
# Setup the parameters and distributions to sample from: param_dist
param_dist = {"max_depth": [3, None],
              "max_features": randint(1, 9),
              "min_samples_leaf": randint(1, 9),
              "criterion": ["gini", "entropy"]}


X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=.2)

# Instantiate a Decision Tree classifier: tree
tree = DecisionTreeClassifier()

# Instantiate the RandomizedSearchCV object: tree_cv
tree_cv = RandomizedSearchCV(tree, param_dist, cv= 5)

# Fit it to the data
tree_cv.fit(X_train, y_train)

# Print the tuned parameters and score
print("Tuned Decision Tree Parameters: {}".format(tree_cv.best_params_))
print("Best score is {}".format(tree_cv.best_score_))


tree2 = DecisionTreeClassifier( max_features=7, min_samples_leaf=5)

tree2.fit(X_train, y_train)

#print( tree2.decision_path(X_train))
#tree2.decision_path(X_train)

#export_graphviz(tree2, out_file='full_cleandata_tree.dot')