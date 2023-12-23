# Script that makes use of more advanced feature selection techniques
# by Alberto Tonda, 2017

import copy
import datetime
import graphviz
import logging
import numpy as np
import os
import sys
import pandas as pd 

from sklearn.ensemble import AdaBoostClassifier
from sklearn.ensemble import BaggingClassifier
from sklearn.ensemble import ExtraTreesClassifier
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.ensemble import RandomForestClassifier

from sklearn.linear_model import ElasticNet
from sklearn.linear_model import ElasticNetCV
from sklearn.linear_model import Lasso
from sklearn.linear_model import LassoCV
from sklearn.linear_model import LogisticRegression
from sklearn.linear_model import LogisticRegressionCV
from sklearn.linear_model import PassiveAggressiveClassifier
from sklearn.linear_model import RidgeClassifier
from sklearn.linear_model import RidgeClassifierCV
from sklearn.linear_model import SGDClassifier

from sklearn.multiclass import OneVsOneClassifier 
from sklearn.multiclass import OneVsRestClassifier
from sklearn.multiclass import OutputCodeClassifier

from sklearn.naive_bayes import BernoulliNB
from sklearn.naive_bayes import GaussianNB
from sklearn.naive_bayes import MultinomialNB 

from sklearn.neighbors import KNeighborsClassifier
from sklearn.neighbors import NearestCentroid
from sklearn.neighbors import RadiusNeighborsClassifier

from sklearn.svm import SVC

from sklearn.tree import DecisionTreeClassifier
from sklearn.tree import ExtraTreeClassifier
from sklearn.neural_network import MLPClassifier
# used for normalization
from sklearn.preprocessing import  Normalizer
from sklearn.preprocessing import StandardScaler
from sklearn.preprocessing import MinMaxScaler
# used for cross-validation
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_curve, auc
from sklearn.metrics import matthews_corrcoef
from scipy import interp
import matplotlib.pyplot as plt
# this is an incredibly useful function
from pandas import read_csv
from sklearn.utils import resample

def loadDataset() :
	index=0
	# data used for the predictions
	dfData = read_csv("./best/data_"+str(index)+".csv", header=None, sep=',')
	dfLabels = read_csv("./best/labels.csv", header=None)
		
	return dfData.values, dfLabels.values.ravel() # to have it in the format that the classifiers like


def runFeatureReduce() :
	
	# a few hard-coded values
	numberOfFolds = 10
	
	# list of classifiers, selected on the basis of our previous paper "
	classifierList = [
		
			# [GradientBoostingClassifier(n_estimators=300), "GradientBoostingClassifier(n_estimators=300)"],
			# [RandomForestClassifier(n_estimators=300), "RandomForestClassifier(n_estimators=300)"],
			# [LogisticRegression(), "LogisticRegression"],
			# [PassiveAggressiveClassifier(),"PassiveAggressiveClassifier"],
			# [SGDClassifier(), "SGDClassifier"],
			# [SVC(kernel='linear'), "SVC(linear)"],
			# [RidgeClassifier(), "RidgeClassifier"],
			# [BaggingClassifier(n_estimators=300), "BaggingClassifier(n_estimators=300)"],
			
			[AdaBoostClassifier(n_estimators=300), "AdaBoostClassifier(n_estimators=300)"],
            [ExtraTreesClassifier(n_estimators=300), "ExtraTreesClassifier(n_estimators=300)"],
            [KNeighborsClassifier(), "KNeighborsClassifier"],
            [LassoCV(), "LassoCV"],
			[MLPClassifier(), "MLPClassifier"]
			]
	
	# this is just a hack to check a few things
	#classifierList = [
	#		[RandomForestClassifier(), "RandomForestClassifier"]
	#		]

	print("Loading dataset...")
	X, y = loadDataset()
	
	print(len(X))
	print(len(X[0]))
	print(len(y))

	
	labels=np.max(y)+1
	# prepare folds
	skf = StratifiedKFold(n_splits=numberOfFolds, shuffle=True)
	indexes = [ (training, test) for training, test in skf.split(X, y) ]

	
	# iterate over all classifiers
	classifierIndex = 0
	
	for originalClassifier, classifierName in classifierList :
		
		print("\nClassifier " + classifierName)
		classifierPerformance = []		
		fig1, ax1 = plt.subplots()
		tprs = []
		aucs = []
		mccs = []
		mean_fpr = np.linspace(0, 1, 100)
		i = 0
		j = 0
		yTest=[]
		yNew=[]
	
		for train_index, test_index in indexes :
			
			X_train, X_test = X[train_index], X[test_index]
			y_train, y_test = y[train_index], y[test_index]
			
			Majority=1
			Minority=0
			
		
			# let's normalize, anyway
			# MinMaxScaler StandardScaler Normalizer
			scaler = StandardScaler()
			X_train = scaler.fit_transform(X_train)
			X_test = scaler.transform(X_test)

			if False:	
				# Separate majority and minority classes
				df_majority = X_train[y_train==Majority]
				df_minority = X_train[y_train==Minority]
				
				print("df_majority "+str(len(df_majority)))
				print("df_minority "+str(len(df_minority)))
				
					# Upsample minority class
				df_minority_upsampled = resample(df_minority, 
												 replace=True,     # sample with replacement
												 n_samples=len(df_majority),    # to match majority class
												 random_state=123) # reproducible results
				print("df_minority_upsampled "+str(len(df_minority_upsampled)))
				# Combine majority class with upsampled minority class
				df_upsampled = np.concatenate([df_majority, df_minority_upsampled])
				print('df_upsampled '+str(len(df_upsampled)))
				print('df_upsampled2D  '+str(len(df_upsampled[0])))
				
				arr=[]
				arr = [Minority for i in range(len(df_majority))] 
				print('arr '+str(len(arr)))
				
				arr2=[]
				arr2 = [Majority for i in range(len(df_majority))] 
				print('arr '+str(len(arr2)))
				
				y_train = np.concatenate([arr2, arr])
				print('y_train '+str(len(y_train)))
				X_train=df_upsampled
		
			
			classifier = copy.deepcopy(originalClassifier)
			classifier.fit(X_train, y_train)
			scoreTraining = classifier.score(X_train, y_train)
			scoreTest = classifier.score(X_test, y_test)
			
			y_new = classifier.predict(X_test)
            
			
			for i in range(0, len(y_new)):
				y_new[i]=round(y_new[i])
                
			
			
			fpr, tpr, thresholds = roc_curve(y_test, y_new)
			tprs.append(interp(mean_fpr, fpr, tpr))
			tprs[-1][0] = 0.0
			roc_auc = auc(fpr, tpr)
			aucs.append(roc_auc)
			mcc = matthews_corrcoef(y_test, y_new)
			print("mcc %0.4f"%(mcc))
			mccs.append(mcc)
			
			ax1.plot(fpr, tpr, lw=1, alpha=0.3,
					 label='ROC fold %d (AUC = %0.2f)' % (j, roc_auc))
			j += 1
			i += 1
			
				
			
			
			print("\ttraining: %.4f, test: %.4f" % (scoreTraining, scoreTest))
			classifierPerformance.append( scoreTest )
		classifierIndex+=1
		line ="%s \t %.4f \t %.4f \n" % (classifierName, np.mean(classifierPerformance), np.std(classifierPerformance))
		
		print(line)
	
		ax1.plot([0, 1], [0, 1], linestyle='--', lw=2, color='r',
			 label='Chance', alpha=.8)
		
		mean_mcc = np.mean(mccs)
		print(mccs)
		std_mcc = np.std(mccs)
		print("mean_mcc %0.4f std_mcc %0.4f"%(mean_mcc,std_mcc))
		mean_tpr = np.mean(tprs, axis=0)
		mean_tpr[-1] = 1.0
		mean_auc = auc(mean_fpr, mean_tpr)
		std_auc = np.std(aucs)
		ax1.plot(mean_fpr, mean_tpr, color='b',
				 label=r'Mean ROC (AUC = %0.2f $\pm$ %0.2f)' % (mean_auc, std_auc),
				 lw=2, alpha=.8)

		std_tpr = np.std(tprs, axis=0)
		tprs_upper = np.minimum(mean_tpr + std_tpr, 1)
		tprs_lower = np.maximum(mean_tpr - std_tpr, 0)
		ax1.fill_between(mean_fpr, tprs_lower, tprs_upper, color='grey', alpha=.2,
						 label=r'$\pm$ 1 std. dev.')

		ax1.axis(xmin=-0.05,xmax=1.05)
		ax1.axis(ymin=-0.05,ymax=1.05)

		ax1.set_xlabel('False Positive Rate')
		ax1.set_ylabel('True Positive Rate')
		ax1.set_title('ROC %s'% ( classifierName))
		ax1.legend(loc="lower right")
		#plt.show()
		plt.savefig("%s.pdf"% ( classifierName))
		plt.savefig("%s.png"% ( classifierName), dpi=300)
	return

if __name__ == "__main__" :
	sys.exit( runFeatureReduce() )