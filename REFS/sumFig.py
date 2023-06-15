# -*- coding: utf-8 -*-
"""
Created on Thu Jan  7 12:15:27 2021

@author: alber
"""
# convert float to percentage string
def convert_to_percentage(f) :
    if f != "" :
        print("Converting \"%.4f\"..." % f)
        percentage = f * 100.0
        return "%.1f%%" % percentage
    else :
        return ""

# script to create a figure
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.cm
#df = pd.read_csv("Results-Variant-UK.csv")
#df = pd.read_csv("Results-Variant-AU.csv")

file2 = open('./best/sumA.csv', 'w')

file2.write("features, run0, run1, run2, run3, run4, run5, run6, run7, run8, run9\n")

filepath = './best/sum.csv'
with open(filepath) as file1:
   line = file1.readline()
   cnt = 1
   while line:
       #print("Line {}: {}".format(cnt, line.strip()))
       file2.write(line)
       line = file1.readline()
       cnt += 1
file1.close()
file2.close()



df = pd.read_csv("./best/sumA.csv")

x = df['features'].values

features=pd.read_csv("./best/features_0.csv", header=None)
print("len features:"+str(len(features.values)))

maxValue=len(features.values)


runs = [r for r in list(df) if r != 'features']

#We declare 15 because numbers 0-5 are almost white.
n_lines=15
c = np.arange(1, n_lines + 1)

norm = mpl.colors.Normalize(vmin=c.min(), vmax=c.max())
cmap = mpl.cm.ScalarMappable(norm=norm, cmap=mpl.cm.Blues)
cmap.set_array([])

import matplotlib.pyplot as plt
fig = plt.figure(figsize=(8,6))
ax = fig.add_subplot(111)
ax.set_xscale('log')

#starts in 5 because numbers 0-5 are almost white.
i=5
for r in runs :
    y = df[r].values
    ax.plot(x, y, label=r, c=cmap.to_rgba(i))
    i=i+1

plt.xticks(x, [str(a) for a in list(x)], rotation=90, fontsize=6)
y_locs = ax.get_yticks()
print(y_locs)
plt.yticks(y_locs, [convert_to_percentage(p) for p in y_locs])

ax.axvline(linewidth=2, color='r', x=maxValue)
ax.grid(linestyle='--')    
ax.legend(loc='best')
ax.set_xlabel("Number of features (log scale)")
ax.set_ylabel("Ensemble accuracy")
ax.set_title("Accuracy vs number of features in REFS runs")
plt.savefig("sumFig.pdf")
plt.savefig("sumFig.png", dpi=300)
