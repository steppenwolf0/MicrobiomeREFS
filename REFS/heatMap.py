import seaborn as sn
import pandas as pd
import matplotlib.pyplot as plt
from pandas import read_csv
from sklearn.preprocessing import StandardScaler,MinMaxScaler
import sys
from matplotlib.colors import LinearSegmentedColormap
directory="best"
import numpy as np
from numpy import array
def loadDataset() :
	index=0
	# data used for the predictions
	dfData = read_csv("./"+directory+"/data_"+str(index)+".csv", header=None, sep=',')
	# data used for the predictions

	dfLabels = read_csv("./"+directory+"/labels_0.csv", header=None)
	dfFeats = read_csv("./"+directory+"/features_"+str(index)+".csv", header=None)
	dfClasses = read_csv("./"+directory+"/classes.csv", header=None)		
	return dfData.values, dfLabels.values.ravel(), dfFeats.values.ravel(), dfClasses.values.ravel() # to have it in the format that the classifiers like

def test():

	X,y,feats,classes=loadDataset()
	
	indexY=np.argsort(y, axis=0) 
	print(indexY)
	X=X[indexY]
	y=y[indexY]
	
	pd.DataFrame(y).to_csv("./"+directory+"/y.csv", header=None, index =None)
	
	scaler = MinMaxScaler()
	X = scaler.fit_transform(X)
	dfData=X
	
	ylabels = ["" for i in range(len(y))]
	if y[0]==0:
		ylabels[0]=str(classes[0])
	for  i in range (1, len(y)):
		if y[i]==1 and y[i-1]==0:
			ylabels[i]=str(classes[1])
		if y[i]==2 and y[i-1]==1:
			ylabels[i]=str(classes[2])
	
	xlabels = feats
	f, ax = plt.subplots(figsize=(20,5))
	#ax = sn.heatmap(dfData, annot=True, fmt="d", cmap=sn.color_palette("Blues"),  cbar=False)
	#palette = sn.color_palette("GnBu_d",30)
	palette = sn.color_palette("Blues",30)
	#palette.reverse()
	sn.set(font_scale = 1.0)
	# Define colors

	ax = sn.heatmap(dfData, annot=False, annot_kws={"size": 10}, fmt=".4f", cmap=palette,    
	#cbar_kws={'label': 'Accuracy 'orientation': 'horizontal'},
	cbar_kws={'label': '', "shrink": 0.50 }, 
	cbar=True, xticklabels=xlabels, yticklabels=ylabels)
	ax.tick_params(axis='both', which='major', labelsize=10)
	ax.set_xticklabels(ax.get_xticklabels(), rotation=90, horizontalalignment='center')
	ax.set_yticklabels(ax.get_yticklabels(), rotation=0, horizontalalignment='right')
	labels_Final = ax.get_xticklabels()
	pd.DataFrame(labels_Final).to_csv("./"+directory+"/axisX.csv", header=None, index =None)
	

	
	#x.set_xlabel('Cancer Type', fontsize=12)
	#ax.set_ylabel("Classifier", fontsize=12)
	plt.show()
	#plt.savefig("test")		
	print(np.max(X))
	return	

if __name__ == "__main__" :
	sys.exit( test() )

