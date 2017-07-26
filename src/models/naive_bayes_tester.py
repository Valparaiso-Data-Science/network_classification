import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split, cross_val_score, StratifiedKFold
from sklearn.metrics import classification_report, confusion_matrix
from sklearn.naive_bayes import GaussianNB


from git.src.models.model_class import ModelTester

infile = 'C:/Users/Owner/Documents/VERUM/Network stuff/git/src/data/clean_data_with_new_chem.csv' # -- change for machine
df = pd.read_csv(infile, index_col=0)

scores = [] # sets up list for all cv score averages
scoresDiff = [] # sets up list for all differences of cv score averages
iterations = 100
toDrop = ['Graph', 'Collection',
          # 'Nodes',
          # 'Edges',
          # 'Density',
           'Maximum degree',
           'Minimum degree',
          # 'Average degree',
           'Assortativity',
           'Total triangles',
          # 'Average triangles',
           'Maximum triangles',
          # 'Avg. clustering coef.',
          # 'Frac. closed triangles',
          # 'Maximum k-core',
          # 'Max. clique (lb)'
            ]
tester = ModelTester(df)

for i in range(iterations):
    score = tester.modelFitTest(GaussianNB(), dropList=toDrop, cv=5, prnt=False)
    scores.append(score)

print('average score: ', np.mean(scores))
#print('average score difference: ', np.mean(scoresDiff))
