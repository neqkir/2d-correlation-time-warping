##########################################################################################
##########################################################################################
####
####             Correlation Optimized Warping testing
#### 
##########################################################################################
 
import cow
import cow2d
import see2D

import matplotlib.pyplot as plt
 
import math
import numpy as np


def visualizeSignals(mat):
 
    nrows = len(mat)
    ncols = len(mat[0])
    
    fig = plt.figure('2d profiles')
    ax = fig.gca(projection='3d')
    
    xr = np.arange(nrows)
    yr = np.arange(ncols)
    
    for i in range(nrows):
        x = xr[i]
        z = np.array(mat[i])
        ax.plot(yr, z, zs = x, zdir='z')
    plt.show()
    
#############################################################################
## test with synthetic data (built from 1-D chromatogram)
#############################################################################

with open("A1.txt","r") as infile:
    a1 = [float(line) for line in infile]
with open("A2.txt","r") as infile:
    a2 = [float(line) for line in infile]

a1 = a1[9700:9900] 
a2 = a2[9700:9900] 
            
 
Y = []
Y.append(a1)
Y.append(a2)
cow2 = cow2d.COW2D(Y, 16, 16, 6, 6, 1, 1)
##sample1=[cow2.sample[0][j] for j in range(len(cow2.sample[0]))]
##sample2=[cow2.sample[1][j] for j in range(len(cow2.sample[0]))]
##plt.plot(sample1, 'b')
##plt.plot(sample2, 'r')
##plt.show()

##cow2.align()

##target1=[cow2.target[0][j] for j in range(len(cow2.target[0]))]
##target2=[cow2.target[1][j] for j in range(len(cow2.target[0]))]
##plt.plot(target1, 'b')
##plt.plot(target2, 'r')
##plt.show()

##visualizeSignals(cow2.target)
##visualizeSignals(cow2.sample)

Y = []
with open("A1.txt","r") as infile:
    a1 = [float(line) for line in infile]
with open("A2.txt","r") as infile:
    a2 = [float(line) for line in infile]
a1 = a1[9700:9900] 
a2 = a2[9700:9900]
for i in range(10):
    Y.append(a1)
    Y.append(a2)
cow2 = cow2d.COW2D(Y, 16, 16, 6, 6, 1, 1)
cow2.align()
# before alignment
visualizeSignals(cow2.sample)
# after alignment
visualizeSignals(cow2.target)

#######################################################################################
##  We show limitations of 2D comprehensive cow 
