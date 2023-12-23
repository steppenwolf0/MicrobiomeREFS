import copy
import datetime
import logging
import numpy as np
import os
import sys

# used for normalization
from sklearn.preprocessing import StandardScaler

# this is an incredibly useful function
from pandas import read_csv
import pandas as pd
#import feature selection with SelectKBest
from sklearn.feature_selection import SelectKBest, f_regression, f_classif, chi2, mutual_info_classif, mutual_info_regression

#load data
def loadDataset() :
	dfData = pd.read_csv(r".\data\data.csv", header=None, sep=',', dtype=float).values
	dfLabels = pd.read_csv(r".\data\labels.csv", header=None)
	biomarkers = pd.read_csv(r".\data\features_0.csv", header=None)
	dfData=np.transpose(dfData)
	return dfData, dfLabels.values.ravel(), biomarkers.values.ravel() # to have it in the format that the classifiers like

def featureSelection() :
	#declare selector with 4 features using F-score
	selector=SelectKBest(f_classif, k=9)
	#load Dataset
	X, y, biomarkerNames = loadDataset()
	#Normalize Data
	scaler = StandardScaler()
	X = scaler.fit_transform(X)
	#Calculate Scores
	X_new = selector.fit_transform(X, y)
	#Get positions of Best Scores
	selected=selector.get_support(indices=True)
	#Print ANOVA F-Values
	print("ANOVA F-value")
	print(selector.scores_[selected])
	#Print P-values
	print("p values")
	print(selector.pvalues_[selected])
	#Print Resulting Features
	print("features names")
	print(biomarkerNames[selected])
	print("features index")
	#Print Features Index
	print(selected)
	# create folder
	folderName ="./results/"
	if not os.path.exists(folderName) : os.makedirs(folderName)
	#Print reduce Dataset
	pd.DataFrame(X_new).to_csv(folderName+"data_"+str(0)+".csv", header=None, index =None)
	pd.DataFrame(biomarkerNames[selected]).to_csv(folderName+"features_"+str(0)+".csv", header=None, index =None)
	pd.DataFrame(y).to_csv(folderName+"labels.csv", header=None, index =None)
	return 

if __name__ == "__main__" :
	sys.exit( featureSelection() )