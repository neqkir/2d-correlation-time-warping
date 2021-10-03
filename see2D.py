##########################################################################################
##########################################################################################
####
####            See 2D chromatograms   
#### 
##########################################################################################

import math
import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import axes3d
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
 

class GC2D:

    gc2d = []
    retT = []
    matGC2D = []
    matRetT = []
    mod = 0 # modulation period
    N = 0 # number of secondary chromatograms
    n = 0 # number of data points in 1-D raw GCxGC data
    
    def __init__(self, pathAbundance, pathRetTimes, modulation):
        
        with open(pathAbundance,"r") as infile:
            self.gc2d = [float(line) for line in infile]
            
        self.mod = modulation
        self.n = len(self.gc2d)
        self.N = len(self.gc2d)/self.mod
        
        with open(pathRetTimes,"r") as infile:
            self.retT = [float(line) for line in infile]
 
             
    def modulate(self, imin, imax):
        """
        Build modulated 2-D chromatogram.
        Indices between imin and imax are taken into account
        """
        times = self.retT[imin:imax]
        gc2d = self.gc2d[imin:imax]

        if(len(times)!=len(gc2d)):
            return

        nrows = (int)((times[(len(times)-1)] - times[0])/self.mod) + 1

        self.matGC2D = []
        self.matRetT = []
        
        for i in range(nrows):
            self.matGC2D.append([])
            self.matRetT.append([])

        for i in range(len(times)-1):
            t = times[i] - times[0]
            x = (int)(t/self.mod) # row in 2-D GC2D matrix
            (self.matRetT[x]).append(t%self.mod) # to this row, append remaining ret time
            (self.matGC2D[x]).append(gc2d[i]) # to this row, append corresponding abundance

          
    def chromatogramIndex(self, retT):
        """
        Return chromatogram to which the point at 1D retention time retT belongs
        """

        return retT/self.mod

    def show_2D(self, mini, maxi, peaks):
        """
        Rough peak detect 2D map show
        Returns retention times for peaks (first and second dimension)
        . peaks: contains indices for peaks in gc2d[mini:maxi]
        """
        times = self.retT[mini:maxi]
        # match each peak to its retention time in the first dimension
        peakT1 = [((times[i]/self.mod)*self.mod) for i in peaks]
        peakT2 = []
        for i in peaks:  
            peakT2.append(times[i]%self.mod)
        return peakT1, peakT2

    def peakDetect(self):
        return self.peakDetectSeg(0,self.n-1)

    def roughPeakDetectSeg(self, imin, imax):
        """
        Returns list of roughly identified peaks in signal between indices imin and imax.
        """
        GCGC = self.gc2d[imin:imax]
        print GCGC
        size = len(GCGC)
        i = 1
        peaks = []
        
        while(i < size-2):

            if(GCGC[i-1] < GCGC[i] and GCGC[i+1] < GCGC[i] ): # apex
                jmax, j, jminleft = i, i-2, i-1
    
                #localNoiseMin = self.estimNoise() ## TODO

                # look for min left
                while((j>0) and (GCGC[j] < GCGC[jminleft]+ 10**(-5))): # + thshld * localNoiseMin)):
                    if(GCGC[j] < GCGC[jminleft]):
                        jminleft = j
                        #localNoiseMin = self.estimNoise()
                    j = j-1

                j,jminright = i+2, i+1
                #localNoiseMin = self.estimNoise()

                # look for min right
                while((j < size-1) and (GCGC[j] < GCGC[jminright]+ 10**(-5))): #+ thshld * localNoiseMin)):
                    if(GCGC[j] < GCGC[jminright]):
                        jminright = j
                        #localNoiseMin = self.estimNoise()
                    j = j+1
         
                # valley - valley
                valleysDiff = abs(GCGC[jminleft] - GCGC[jminright])
                # minimum max - valley difference
                maxDiff =  min(abs(GCGC[jmax] - GCGC[jminright]), abs(GCGC[jmax] - GCGC[jminleft]))

                if(valleysDiff < 2.0 * maxDiff):
                    #print 'valleyDiff ' + str(valleysDiff) +'maxDiff ' + str(maxDiff)
                    #print 'jminleft ' + str(jminleft) +'jminright ' + str(jminright)
                    peaks.append(jmax)
                    
                i = jminright + 1
 
            else:
                i = i+1
                
        return peaks


    def peakDetectSeg(self, imin, imax):
        """
        Return refined list of peaks in signal between imin and imax.
        
        """
        gc2d, times = self.gc2d[imin:imax], self.retT[imin:imax]  # segments
        
        rPeaks = self.roughPeakDetectSeg(imin, imax) # first list of potential peaks (indices)
        i = 1
        peaks = []
        while(i < len(rPeaks) - 2):
        
            if(gc2d[rPeaks[i-1]] < gc2d[rPeaks[i]] and gc2d[rPeaks[i+1]] < gc2d[rPeaks[i]] ): # apex
                peaks.append(rPeaks[i])
            i = i + 1         
        return peaks

    
    def visualizeSignals(self, imin, imax):
        """
        3d plot signal between imin and imax
        """
        # don't like incomplete secondary chromatogram
        # max_ = bisect_left(self.retT[imin:imax],(self.retT[imax]/self.mod)*self.mod-10e-5)+imin
        times = self.retT[imin:imax]
        nrows = (int)((times[(len(times)-1)] - times[0])/self.mod) + 1
        
        self.modulate(imin, imax)
        
        fig = plt.figure('2d profiles')
        ax = fig.gca(projection='3d')
        for i in range(nrows-1):
            x = self.matRetT[i][0] + self.mod * i
            y = np.array(self.matRetT[i])
            z = np.array(self.matGC2D[i])
            ax.plot(y, z, zs = x, zdir='z')
        plt.show()



    def visualize(self, imin, imax, typ_ = "wireframe"):
        """
        3d plot signal between imin and imax
        . typ_: type of plot, "wireframce", "surface"
        """
        # don't like incomplete secondary chromatogram
        times = self.retT[imin:imax]
        nrows = (int)((times[(len(times)-1)] - times[0])/self.mod) + 1
        
        self.modulate(imin, imax)
        
        fig = plt.figure('3d view')
        ax = fig.gca(projection='3d')
        
        x = []
        for i in range(nrows):
            x.append(self.matRetT[i][0] + self.mod * i)

 
        y = []
        for i in range(len(self.matRetT[0])):
            y.append(self.matRetT[0][i])
        y = y[:-1]
 
        X,Y = np.meshgrid(x,y)
 
        z = [tuple(self.matGC2D[i]) for i in range(len(self.matGC2D))]
        
        zzip = zip(*z)
 
        if(typ_ == "wireframe"):
            ax.plot_wireframe(X,Y,zzip)
            plt.show()
        elif(typ_ == "contour"):
            cset = ax.contour(X, Y, zzip, zdir='z', offset=0)
            plt.show()
        elif(typ_ == "surf_contours"):
            surf = ax.plot_surface(X, Y, zzip, rstride=1, cstride=1, alpha=0.3)
            cset = ax.contour(X, Y, zzip, zdir='z', offset=-40)
            cset = ax.contour(X, Y, zzip, zdir='x', offset=-40)
            cset = ax.contour(X, Y, zzip, zdir='y', offset=-40)
            plt.show()




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
