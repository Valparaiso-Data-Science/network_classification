import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.linear_model import Lasso
import matplotlib.pyplot as plt
from sklearn.model_selection import GridSearchCV

infile = '~/PycharmProjects/network_classification/src/data/data_minmaxscale.csv' # --- change this for specific machine
df = pd.read_csv(infile, index_col=0)

X = df.drop(['Graph', 'Collection'], axis=1).values
y1 = df['Collection'].values

y = df.Collection.astype('category').cat.codes.values

lasso = Lasso(alpha=0.05)

lasso_coef = lasso.fit(X, y).coef_

name = df.columns.values
index = [0, 1]

names = np.delete(name, index)

plt.plot(range(len(names)), lasso_coef)
plt.xticks(range(len(names)), names, rotation=90)
plt.ylabel = ('Coefficients')
plt.show()



#doing grid search cv to tune hyperparameter alpha
#found that alpha=0.05 is our best bet
#param_grid = {'alpha': np.arange(0.01, 3, 0.005)}
#lasso_cv = GridSearchCV(lasso, param_grid, cv=5)
#lasso_cv.fit(X, y)
#print("Tuned Alpha Parameter: {}".format(lasso_cv.best_params_))
#print("Best score is {}".format(lasso_cv.best_score_))
