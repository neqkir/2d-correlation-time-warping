###################################################################################################
###################################################################################################
####
####             Correlation Optimized Warping 2D
#### 
###################################################################################################

## A Comprehensive Two-Dimensional Retention Time Alignment Algorithm To Enhance Chemometric Analysis
## of Comprehensive Two-Dimensional Separation Data,
## Karisa M. Pierce, Lianna F. Wood, Bob W. Wright, Robert E. Synovec, 2005

## Algorithm is adapted from this article
##        - replacing algo for 1D alignment with COW algorithm
##        - building target matrix replacing target matrix current row with sample matrix rows averaged with
##            window_size target matrix rows before and after current row



import math, copy, operator, random, os, sys
import numpy as np
sys.path.append(os.getcwd())
import cow

class COW2D:

    n1, n2 = 0, 0  # nb of frames, dimension 1 and dimension 2
    sample, target = [], [] # sample matrix, reference/target matrix ; sample matrix is 2D chromatogram to be arranged
    s1, s2 = 0, 0  # slack parameters, maximum shift for a segment length, dimension 1 and dimension 2
    L1, L2 = 0, 0 # window parameters, dimension 1 and dimension 2

    def __init__(self, sampleV = None, nbFram1 = None, nbFram2 = None, sl1 = None, sl2 = None, win1 = None, win2 = None):
        """
        . slack: maximum length for decrease or increase in segment lenght
        . X: sample
        . Y: reference
        . n: number of frames in sample and reference
        """
        if (nbFram1 < 4 or nbFram2 < 4):
            print 'WARNING: Number of frames should be at least 4'
            return None

        if(sampleV == None):
            return None
        
        self.s1, self.s2 = sl1, sl2
        self.sample = sampleV
        self.n1, self.n2 = nbFram1, nbFram2
        self.L1, self.L2 = win1, win2

        self.target = [[self.sample[i][j] for j in range(len(self.sample[0]))] for i in range(len(self.sample))]  # at first, target is sample

     
    def align(self):
        """
        Align 2d sample chromatogram with successive averages in target chromatogram
        """
        # align along secondary chromatograms (dimension 2)
        for i in range(len(self.sample)): # while there are some rows
            self.updateAvgTargetDim2(i)
            aligned = self.alignSampleToTargetDim2(i)
            if(aligned != False):
                self.updateTargetAfterAlignDim2(i, aligned)
            

    def alignDim2(self):
        """
        Align 2d sample chromatogram with successive averages in target chromatogram
        """
        # align along secondary chromatograms (dimension 2)
        for i in range(len(self.sample)): # while there are some rows
            self.updateAvgTargetDim2(i)
            aligned = self.alignSampleToTargetDim2(i)
            self.updateTargetAfterAlignDim2(i, aligned)
        

    def alignDim1(self):
        """
        Align 2d sample chromatogram with successive averages in target chromatogram
        """
        # align along dimension 1
        self.target = [list(x) for x  in zip(*self.target)]
        self.sample = [list(x) for x  in zip(*self.sample)]

        self.alignDim2()
        
        self.target = [list(x) for x  in zip(*self.target)]
        self.sample = [list(x) for x  in zip(*self.sample)]
            
    def updateAvgTargetDim2(self, r):
        """
        Average row r in target matrix with self.L2 rows before
        and self.L2 rows after r.
        """
        M, N = len(self.target[0]), len(self.target)
        avg = [self.target[r][j] for j in range(M)]

        i = r-1
        while(i >= 0 and r-i <= self.L2):
            for j in range(M):
                avg[j] = avg[j] + self.target[i][j] 
            i = i - 1
        imin = i + 1 # lowest row taken into account
  
        i = r + 1
        while(i < N and i-r <= self.L2):
            for j in range(M):
                avg[j] = avg[j] + self.target[i][j] 
            i = i + 1
        imax = i - 1 # highest row taken into account

        avg = [avg[i]/(float)(imax-imin+1) for i in range(len(avg))]        

        for j in range(M):
            self.target[r][j] = avg[j]  # replace row r in target with average


            
    def updateAvgTargetDim1(self, c):
        """
        Average col c in target matrix with self.L1 rows before
        and self.L1 rows after c.
        """
        self.target = [list(x) for x  in zip(*self.target)]
        self.updateAvgTargetDim2(c)
        self.target = [list(x) for x  in zip(*self.target)]
      


    def alignSampleToTargetDim2(self,r):
        """
        Align 1D vector at row r in sample to vector at row r in target
        """
        sampleVect = [self.sample[r][j] for j in range(len(self.target[0]))]
        targetVect = [self.target[r][j] for j in range(len(self.target[0]))]
        
        align = cow.COW(sampleVect, targetVect, self.n2, self.s2)
        aligned = align.warp_sample_to_target()
        
        return aligned
        


    def alignSampleToTargetDim1(self,c):
        """
        Align 1D vector at col c in sample to vector at col c in target
        """
        self.target = [list(x) for x  in zip(*self.target)]
        aligned = self.alignSampleToTargetDim2(c)
        self.target = [list(x) for x  in zip(*self.target)]
        
        return aligned



    def updateTargetAfterAlignDim2(self, r, aligned):
        """
        Update target after aligning 1D chromatograms
        """
        for j in range(len(self.target[0])):
            self.target[r][j] = aligned[j] 



    def updateTargetAfterAlignDim1(self, c, aligned):
        """
        Update target after aligning 1D chromatograms
        """
        for i in range(len(self.target)):
            self.target[i][c] = aligned[i]
        
        
def bisect_left(a, x, lo=0, hi=None):
    """
    Return the index where to insert item x in list a, assuming a is sorted.
    Search is bound to lo:hi. Default for hi is len(a). Returned index is used
    directly for insertion
    """
    if (hi == None):
        hi = len(a)
    while lo < hi:
        mid = (lo+hi)/2
        if a[mid] < x:
            lo = mid+1
        else:
            hi = mid
    return lo



