###################################################################################################
###################################################################################################
####
####             Correlation Optimized Warping
#### 
###################################################################################################

## Based on "Correlation optimized warping and dynamic
## time warping as preprocessing methods
## for chromatographic data"
## GiorgioTomasi*,Frans van den Bergand Claus Andersson

import math, copy, operator, random 
import numpy as np

class COW:

    n = 0  #nb of frames
    X, Y = [], [] # sample vector, reference vector
    s = 0  # slack parameter, maximum shift for a segment length
    bounds = []
  

    def __init__(self, sampleV = None , refV = None, nbFrames = None, slack = None):
        """
        . slack: maximum length for decrease or increase in segment lenght
        . X: sample
        . Y: reference
        . n: number of frames in sample and reference
        """
        if (sampleV == None or refV == None or nbFrames == None or slack == None): 
            return None
        if (nbFrames < 4):
            print 'WARNING: Number of frames should be at least 4'
            return None
        self.s = slack
        self.X = sampleV
        self.Y = refV
        self.n = nbFrames
        if(self.X == None or self.Y == None):
            return None
        self.bounds = self.setBounds() # 2x(self.n+1) contains boundaries for segments in X and Y
  


    def correlation(self, frameX, frameY):
        """
        Returns correlation between frameX and frameY if segments
        are same length
        """
        pX, pY = len(frameX), len(frameY)
        
        if(pX != pY):
            return 0
        
        correl,normX,normY = 0,0,0
        
        # substract average
        frameX = frameX - np.mean(frameX)
        frameY = frameY - np.mean(frameY)
        
        for i in range(pX):
            correl += frameX[i] * frameY[i]
            normX  += frameX[i]**2
            normY  += frameY[i]**2 
         
     
        if(math.sqrt(normX) == 0 and math.sqrt(normY) == 0):
            return 1 # two constants
        
        prod = math.sqrt(normX) * math.sqrt(normY)
        
        if(prod < 10e-5):
            prod = 10e-5

        if(prod == 0):
            return 0
        
        return correl /(float)(prod)



    def resample(self, segX, N, beforeX = None, afterX = None):
        """
        Returns segX resampled with N points. If len(segX)== N, return segX.
            . segX: segment from sample vector
            . N: number of points in reference vector
            . beforeX: value in sample vector just before segX, None if first segment in X
            . afterX: value in sample vector just after segX, None if last segment in X
        """
        pX = len(segX)   #points in segX and segY
        sampledSegX = [] # segment X after sampling
         
        if(pX == N):
            return segX # no resampling is needed

        Xabs = map(lambda x: (x+1)/(float)(pX+1), range(pX))  

        ## fill in sampledSegX vector
        
        if(pX > N):  # subsample segX
           for i in range(N):
               currentXabs = (i+1)/(float)(N+1) # abscissa in resampled segX
               idx = bisect_left(Xabs,currentXabs)
               if (idx > 0):
                   slope = (segX[idx] - segX[idx-1])/(float)(Xabs[idx] - Xabs[idx-1])
                   sampledSegX.append(segX[idx-1] + slope * (currentXabs - Xabs[idx-1])) # value for resampled segX at currentXabs
               else:  # insertion failed
                   return []  
                
        elif(pX < N):  # first complete segX with before and after points, then subsample segX
            
            Xabs = [0] + Xabs + [1]
            
            if(beforeX == None): # no interpolation with previous value is possible
                beforeX = segX[0] # repeat first value in segment
            if(afterX == None): # no interpolation with next value is possible
                afterX = segX[-1] # repeat last value in segment
            
            segX = [beforeX] + segX + [afterX]
            
            # print segX 
            for i in range(N):
                currentXabs = (i+1)/(float)(N+1) # abscissa in resampled segX
                idx = bisect_left(Xabs,currentXabs)
                # print "abscissa for insertion" + str(currentXabs) + ", " + str(idx)
                if (idx > 0):
                    slope = (segX[idx] - segX[idx-1])/(float)(Xabs[idx] - Xabs[idx-1])
                    sampledSegX.append(segX[idx-1] + slope * (currentXabs - Xabs[idx-1])) # value for resampled segX at currentXabs
                else: # insertion failed
                    return []

        return sampledSegX



 
    def setBounds(self):
        """
        Initialize segment boundaries.
        """
        nX, nY = len(self.X), len(self.Y)
        bounds = []
        nSegX  = nX/self.n
        nSegY  = nY/self.n
        bX, bY = [], [] 
        # example
        # [..][...][....] corresponds to boundaries [022559] that we abbreviate in [0259]
        for i in range((self.n)):
            bX.append(i*nSegX)
            bY.append(i*nSegY)
        bX.append(len(self.X)) # last segment contains also remainder in nX/self.n
        bY.append(len(self.Y))
        bounds.append(bX)
        bounds.append(bY)
        
        return bounds
      

    def forward_search(self):
        """
        Fill in matrix with maximal values or cumulative correlation.
        Fill in matrix with optimal paths / optimal shifts.
        """
        bX, bY, pX = self.bounds[0], self.bounds[1], len(self.X)
        if(self.n != (len(bX)-1) or self.n != (len(bY)-1)): # test number of bounds
            return False
        
        F, U = [], []  # matrix with cumulative correlation, matrix with optimal path

        # init matrix
        for i in range(self.n):
            tp1, tp2 = [], []
            for j in range(pX):
                tp1.append(float(0))
                tp2.append([])
            F.append([tp1[i] for i in range(pX)])
            U.append([tp2[i] for i in range(pX)])

        bIdx = range(self.n+1)[::-1] # indices for bounds

        ##############################################################################
        # consider each segment right to left:

        # first segment 
        segY = self.Y[bY[bIdx[1]]:bY[bIdx[0]]]  
        jm, jM = self.possible_positions_bounds(self.n, pX/self.n)
        
        for j in range(jm, jM+1): # loop over possible left bound position for first segment
            resampledSegX = self.resample(self.X[j:], len(segY), self.X[j-1], None) # resample segX to size of segY   
            F[-1][j] = self.correlation(resampledSegX, segY) # update  with correlation btw resampled segX and segY
            U[-1][j] = [j] # remember position
        # next segments
        for i in bIdx[1:-2]: # loop over segments
            
            segY = self.Y[bY[i-1]:bY[i]]   

            jm, jM = self.possible_positions_bounds(i, pX/self.n) # possible left bound position  
  
            for j in range(jm, jM+1): # loop over possible left bound position for each segment
                # compute maximum correlation and corresponding path for current position
                correl_comp = self.max_correl_over_combinations(F, U, i, j, segY)
                if(correl_comp == False):
                    continue
                else:
                    max_corr, max_path = correl_comp

                F[i-1][j] = max_corr
                U[i-1][j] = max_path  # remember position  
        # last segment 
        correl_comp = self.max_correl_over_combinations(F, U, 1, 0, self.Y[:bY[bIdx[-2]]]) # maximum correlation and corresponding path for last position, (0,0)
        
        if(correl_comp != False):
            max_corr, max_path = correl_comp

        F[0][0] = max_corr
        U[0][0] = max_path 
             
        #
        #################################################################################
        
        return F, U

    

    def max_correl_over_combinations(self, F, U, segIdX, position, segY):
        """
        For given position in matrix F with cumative correlations,
        explore possible combinations of a previous segment left boundary
        and a current segment length, expressed as a value for the slack parameter.
        
        Return maximum cumulative correlation for this position over combinations.
        Returns path corresponding to this maximum, a path is concatenation of past path
        and slack parameter value leading to maximum.
        
        . segIdx: index for current segment, in 1:self.n ; segIdx is matrix row index
        . position: position in sample vector for current boundary position
        . F: self.n x len(self.X) matrix with cumulative correlation
        . segX: X segment  
        . segY: Y segment 
        """
        if (segIdX == self.n):
            print 'error: index for previous array exceeds number of segments'
            return

        segL = len(self.X)/self.n # length for majority of segments in X
        s = self.s 

        max_corr = float(0)
        max_corr_path = []
        
        for j in range(segL - s, segL + s + 1): # loop over possible segment lengths
            
            prevj = position + j # column for previous left boundary position
            # print 'prevj ' + str(prevj) + '    segIdx   ' + str(segIdX) + '   lenF  '  + str(len(F)) + ' ' + str(len(F[0])) 
            if((prevj < len(F[0])) and F[segIdX][prevj] > float(0) ): # examine if position was possible at previous iteration

                # print 'seg X ' + str(self.X[position:prevj]) + 'before ' + str(self.X[position-1]) + 'after ' + str(self.X[prevj]) + 'X '  + str(self.X)

                if(position > 0 ): # all segment but last, r to l
                    resampledSegX = self.resample(self.X[position:prevj], len(segY), self.X[position-1], self.X[prevj])
                else: # last segment
                    resampledSegX = self.resample(self.X[position:prevj], len(segY), self.X[0], self.X[prevj])
                    
                # print 'resampledSegX ' + str(resampledSegX) + 'segIdx '  + str(segIdX) + ' prevj ' + str(prevj) + 'F[segIdX][prevj]  ' +  str(F[segIdX][prevj]) + ' corr(resampledSegX,segY)' + str(self.correlation(resampledSegX,segY))
                
                if(F[segIdX][prevj] + self.correlation(resampledSegX,segY) > max_corr):         
                    max_corr = F[segIdX][prevj] + self.correlation(resampledSegX,segY) # update max_corr if current combination is improving cumulative correlation
                    max_corr_path = U[segIdX][prevj] +  [position] # update path accordingly
                    
            else:
                # print 'WARNING: no reachable previous position '
                continue
            
        return max_corr, max_corr_path



    def possible_positions_bounds(self, i, m):
        """
        For current row in F and U matrices, returns maximum and minimum.
        Note that jm and jM are column indices starting with 0
        . i: current row index in matrix
        . m: length for most segments in sample signal
        . slack: maximum increase or decrease in segment length
        """
        pX, s = len(self.X), self.s
        
        # lowest col index to inspect at row i
        jm = max([(i - 1)*(m - s),(pX - (pX % self.n) - (self.n - i + 1)*(m + s))])   
        jm = max(0, jm) # should not exceed matrix dimension
         
        # largest col index to inspect at row i
        jM = min([(i - 1)*(m + s),(pX - (pX % self.n) - (self.n - i + 1)*(m - s))])  
        jM = min(jM, (len(self.X)-1)) # should not exceed matrix dimension

        return jm, jM
        
        
    def score(self, sX, sY, eX, eY):
        """
        Score X[sX:eX] w.r.t. Y[sY, eY] and return score.
        . sX, sY = segment start 
        . eX, eY = segment end 
        """
        segX, segY = self.X[sX:eX], self.Y[sY:eY]
        pX, pY = len(segX), len(segY)
        
        if (pX != pY): # resample if needed
            
            if(eX == pX):
                afterX = None# last segment in X 
            else:
                afterX = self.X[eX]
            if(sX==0):
                beforeX = None # first segment in X
            else:
                beforeX = self.X[sX-1]
            
            segX = self.resample(segX, pY, beforeX, afterX) # resample
            
        return self.correlation(segX, segY) # compute correlation
                

    def warp_sample_to_target(self):
        
        fwdsrch = self.forward_search()
        sY = len(self.Y)/self.n # length for majority of segments in Y
        pLastY = len(self.Y)%self.n + sY

        if(fwdsrch != False):
            path = fwdsrch[1][0][0]
            path = path[::-1]
        else:
            return
 
        if(len(path) < 3):
            print 'WARNING: optimal path do not have required number of elements, there is some problem.'
            return False
        
        warpedX = []
        
        #first segment left  
        segX = self.X[:path[1]]
        #print 'segX ' + str(segX) + 'before None' + 'after  '  + str(self.X[path[1]]) + 'resampledSegX ' + str(resampledSegX)
        resampledSegX = self.resample(segX, sY, None, self.X[path[1]])
        warpedX = warpedX + resampledSegX
        
        #next segments
        for i in range(len(path)-2):
            segX = self.X[path[i+1]:path[i+2]]
            resampledSegX = self.resample(segX, sY, self.X[(path[i+1]-1)], self.X[path[i+2]])
            #print 'segX ' + str(segX) + 'before ' + str( self.X[(path[i+1]-1)]) + 'after ' + str(self.X[path[i+2]]) + 'resampledSegX ' + str(resampledSegX)
            warpedX = warpedX + resampledSegX
            
        #last segment right
        segX = self.X[path[-1]:] 
        resampledSegX = self.resample(segX, pLastY, self.X[(path[-1]-1)], None)
        warpedX = warpedX + resampledSegX
        # print 'segX ' + str(segX) + 'before ' + str(self.X[(path[-1]-1)]) + 'after None' + 'resampledSegX ' + str(resampledSegX) + str(warpedX) +  str(len(warpedX))
        
        return warpedX


    
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
