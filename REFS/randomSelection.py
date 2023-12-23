#This script generates a list of random indexes to be searched in the "features_0.csv" file and store the features selected in a new file.
#The intention of doing this is because we need to validate that the features selected by REFS work comparong results with the ones obtained
#in a random way, after select this random features, the data associated with these features is also gathered.
#The output are three csv files: 1) features randomly selected, 2) the data associated to these features, and 3) the indexes randomly selected

import csv
import pandas as pd
import random
import numpy as np

features_List = pd.read_csv('features_0.csv', header=None) #Reading the features list
data_Inf = pd.read_csv('data_0_Nochim.csv', header=None)

features_Size = len(features_List) #to establish the maximun number allowed in the random
num_Features = 16 #number of features wanted
indexes = [] #list of indexes
features_Selected = [] #to store the features
data_Selected = [] #to store the data associated to the feature
i = 0

while len(indexes) <= num_Features-1:
    band = 0
    index_Tmp = random.randint (0,features_Size-1)
    if i == 0:
        indexes.append(index_Tmp)
    if index_Tmp in indexes:
        band = 1
    if band != 1:
        indexes.append(index_Tmp)
    i = i+1

for j in range (0, len(indexes)):
    index_Tmp = indexes[j]
    ft = features_List.values[index_Tmp-1]
    dt = data_Inf.values[:,index_Tmp-1]
    features_Selected.append(ft)
    data_Selected.append(dt)
    
indexes_Final = pd.DataFrame(indexes)
feature_Final = pd.DataFrame(features_Selected)
data_Final = pd.DataFrame(data_Selected)
dataTransposed = data_Final.T

indexes_Final.to_csv('indexes.csv', header=None, index=None)
feature_Final.to_csv('features.csv', header=None, index=None)
dataTransposed.to_csv('data.csv', header=None, index=None)

print("done")