import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split, cross_val_score, StratifiedKFold
from sklearn.metrics import classification_report, confusion_matrix
from sklearn.naive_bayes import GaussianNB
import sys


from git.src.models.model_class import ModelTester

infile = 'C:/Users/Owner/Documents/VERUM/Network stuff/git/src/data/clean_data_with_new_chem.csv' # -- change for machine
df = pd.read_csv(infile, index_col=0)

scores = [] # sets up list for all cv score averages
scoresDiff = [] # sets up list for all differences of cv score averages

toDrop = ['Graph', 'Collection',
          # 'Nodes',
           'Edges',
          # 'Density',
           'Maximum degree',
           'Minimum degree',
          # 'Average degree',
          # 'Assortativity',
           'Total triangles',
           'Average triangles',
           'Maximum triangles',
          # 'Avg. clustering coef.',
          # 'Frac. closed triangles',
          # 'Maximum k-core',
          # 'Max. clique (lb)'
            ]
tester = ModelTester(df)
features = ['Nodes', 'Edges', 'Density', 'Maximum degree', 'Minimum degree', 'Average degree', 'Assortativity', 'Total triangles', 'Average triangles', 'Maximum triangles', 'Avg. clustering coef.', 'Frac. closed triangles', 'Maximum k-core', 'Max. clique (lb)'
            ]
def rfe(model, iterations=100, initialFeatures=['Nodes']):
    sFeatures = set(features)
    sInitial = set(initialFeatures)
    diff = sFeatures.difference(sInitial)

# gets an initial bestScore, using only features in initialFeatures
    scores = []
    for i in range(iterations):
        score = tester.modelFitTest(model, dropList=['Graph', 'Collection'] + list( diff ), prnt=False)
        scores.append(score)
    bestScore = np.mean(scores)

# checks if adding a feature improves score
    bestFeature = None
    for feature in list( diff ):
        newDiff = diff.difference({feature})

        scores =[]
        for i in range(iterations):
            score = tester.modelFitTest(model, dropList=['Graph', 'Collection'] + list(newDiff), prnt=False)
            scores.append(score)
        newScore = np.mean(scores)

     # Every time it finds a better score than using the initialFeatures, it redefines bestScore
        if newScore > bestScore:
            print('found better score adding ', feature, 'score: ', newScore)
            bestScore = newScore
            bestFeature = feature

    if bestFeature == None:
        # checks if taking away one feature improves score
        worstFeature = None
        for feature in initialFeatures:
            newInitial = sInitial.difference({feature})
            newDiff = sFeatures.difference(newInitial)
            scores = []
            for i in range(iterations):
                score = tester.modelFitTest(model, dropList=['Graph', 'Collection'] + list(newDiff), prnt=False)
                scores.append(score)
            newScore = np.mean(scores)

# Every time it finds a better score than using the initialFeatures, it redefines bestScore
            if newScore > bestScore:
                print('found better score removing ', feature, 'score: ', newScore)
                bestScore = newScore
                worstFeature = feature

# if taking away one feature didn't improve the score at all, function finally finishes
        if worstFeature == None:
            print('best score: ', bestScore)
            print('features used: ', sFeatures.difference(diff))

# if at least one feature was found to improve score when removed, the function removes that feature and repeats
        else:
            print('removed ', worstFeature)
            nextInitial = sInitial.difference({worstFeature})
            bestScore = rfe(model, initialFeatures=list(nextInitial))

# if at least one feature was found to improve the score when added, function adds that feature and repeats
    else:
        print('added ', bestFeature)
        bestScore = rfe(model, initialFeatures=initialFeatures + [bestFeature])

    return bestScore

#rfe(GaussianNB(), initialFeatures=[])

#result = rfe(GaussianNB(), initialFeatures=['Density', 'Avg. clustering coef.', 'Frac. closed triangles', 'Maximum k-core', 'Max. clique (lb)'])
result = rfe(GaussianNB(), initialFeatures= list( set(features).difference( set(toDrop) ) ))
#result = rfe(GaussianNB(), initialFeatures=features)
#print(result)

sys.exit('done')

for i in range(100):
    print('iteration ', i)
    score = tester.modelFitTest(GaussianNB(), dropList=toDrop, LOO=True, prnt=False)
    scores.append(score)

print('average score: ', np.mean(scores))
##print('average score difference: ', np.mean(scoresDiff))
