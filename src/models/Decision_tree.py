import numpy as np
import pandas as pd
from scipy.stats import randint
from sklearn.tree import DecisionTreeClassifier, export_graphviz
from sklearn.model_selection import RandomizedSearchCV, GridSearchCV, train_test_split, cross_val_score

infile = 'C:/Users/Owner/Documents/VERUM/Network stuff/git/src/data/data_minmaxscale.csv' # -- change for machine
df = pd.read_csv(infile, index_col=0)

# This adjusts the dataframe so that we only look at four distinct collections, by removing ones we don't want
names = ['Biological Networks','Collaboration Networks',
         'Ecology Networks', 'Infrastructure Networks',
         'Interaction Networks', 'Massive Network Data',
         'Recommendation Networks', 'Scientific Computing',
         'Social Networks', 'Technological Networks', 'Web Graphs'
         #'Facebook Networks', 'Brain Networks', 'Cheminformatics', 'Retweet Networks'
 ]
#for name in names:
   # df = df[df['Collection'] != name ]

#To be removed later, once there is a file that already has this
collections = np.unique( df.Collection.values )
for collection in collections:
    size = len( df[ df.Collection == collection ] )
    if size < 10:
        df = df[ df.Collection != collection ]

X = df.drop(['Graph', 'Collection'], axis=1).values
y = df['Collection'].values



# Setup the parameters and distributions to sample from: param_dist
# This one sets up the parameters for the Randomized Search
param_dist_rand = {"max_depth": [20, None],
              "max_features": randint(1, 9),
              "min_samples_leaf": randint(1, 9),
              "criterion": ["gini", "entropy"]}

# This sets up the parameters for the Grid Search
param_dist_grid = {"max_depth": np.arange(2, 25),
              "max_features": np.arange(1,14),
              "min_samples_leaf": np.arange(1, 20),
              "criterion": ["gini", "entropy"]}


X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=.2)

# Instantiate a Decision Tree classifier: used for the Random/Grid Search
tree = DecisionTreeClassifier()

#  ---- The following lines are commented out as they are only used to determine the best parameters
# Instantiate the RandomizedSearchCV object: tree_cv
#tree_cv = RandomizedSearchCV(tree, param_dist_rand, cv= 5)
#tree_cv = GridSearchCV(tree, param_dist_grid, cv= 5)

# Fit it to the data
#tree_cv.fit(X_train, y_train)

# Print the tuned parameters and score
#print("Tuned Decision Tree Parameters: {}".format(tree_cv.best_params_))
#print("Best score is {}".format(tree_cv.best_score_))

# This classifier is set up with the best parameters found by the grid search
tree2 = DecisionTreeClassifier( criterion='entropy', max_depth=24, max_features=8, min_samples_leaf=1)

cv_scores = cross_val_score(tree2, X_train, y_train, cv = 5)
print('average cv score from training set: ', np.mean( cv_scores ))

tree2.fit(X_train, y_train)
tree2.predict(X_test)
test_score = tree2.score(X_test, y_test)
print('test_score: ', test_score)


#print( tree2.decision_path(X_train))
#tree2.decision_path(X_train)

#export_graphviz(tree2, out_file='full_cleandata_tree.dot')