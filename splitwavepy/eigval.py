"""
The eigenvalue method of Silver and Chan (1991)
"""

import core as c
from . import plotting as p
import numpy as np
import matplotlib.pyplot as plt
from scipy import signal, stats

class Measurement:
    
    def __init__(self,data=None,degs=None,lags=None,lam1=None,lam2=None,window=None):
        
        self.method = 'Eigenvalue'
        
        if data is None:           
            # generate synthetic
            self.data, self.degs, self.lags, self.lam1, self.lam2, self.window = _synthM()      
        elif data is not None and (degs is None or lags is None or lam1 is None or lam2 is None or window is None):                                              
            # generate measurement using defaults
            self.data,self.degs, self.lags, self.lam1, self.lam2, self.window = _grideigval(data)           
        else:           
            # everything provided
            self.data = data
            self.degs = degs
            self.lags = lags
            self.lam1 = lam1
            self.lam2 = lam2
            self.window = window
            
        # ensure data is a Pair for convenience functions
        if not isinstance(self.data,c.Pair):
            self.data = c.Pair(self.data)
        
        # get some measurement attributes
        # uses ratio lam1/lam2 to find optimal fast and lag parameters
        maxloc = c.max_idx(self.lam1/self.lam2)
        self.fast = self.degs[maxloc]
        self.lag  = self.lags[maxloc]
        
        # get some useful stuff
        self.data_corr = c.unsplit(self.data.data,self.fast,self.lag)
        self.srcpol = c.pca(self.data_corr)
        self.srcpoldata = c.rotate(self.data.data,-self.srcpol)
        self.srcpoldata_corr = c.rotate(self.data_corr,-self.srcpol)
        
        # signal to noise ratio estimates
        # self.snr = c.snr(c.window(self.srcpoldata_corr,self.window))
        # self.snrRH = c.snrRH(c.window(self.srcpoldata_corr,self.window))
        self.snr = np.max(self.lam1/self.lam2) 

        # number degrees of freedom
        self.ndf = ndf(c.window(self.srcpoldata_corr[1,:],self.window))
        # value of lam2 at 95% confidence contour
        self.lam2_95 = Ftest(self.lam2,self.ndf,alpha=0.05)

        # convert traces to Pair class for convenience
        self.data_corr = c.Pair(self.data_corr)
        self.srcpoldata = c.Pair(self.srcpoldata)
        self.srcpoldata_corr = c.Pair(self.srcpoldata_corr)
        

    def plot(self,vals=None,cmap='magma',lam2_95=True):
        """
        plot the measurement.
        by default plots lam1/lam2 with the lambda2 95% confidence interval overlaid
        """
        if vals is None:
            vals = self.lam1 / self.lam2
            
        plt.contourf(self.lags,self.degs,vals,cmap=cmap)
        
        if lam2_95 is True:
            plt.contour(self.lags,self.degs,self.lam2,levels=[self.lam2_95])
        
        plt.show()

    # def save():
    #     """
    #     Save Measurement for future referral
    #     """
        
    
def eigvalcov(data):
    """
    return sorted eigenvalues of covariance matrix
    lambda1 first, lambda2 second
    """
    return np.sort(np.linalg.eigvals(np.cov(data)))
    
    
def _grideigval(data, maxlag=None, window=None, stepang=None, steplag=None):

    # set some defaults
    if maxlag is None:
        maxlag = int(data[0].size / 10)
        maxlag = maxlag if maxlag%2==0 else maxlag + 1
    if steplag is None:
        steplag = 2 * int(np.max([1,maxlag/80]))
    if stepang is None:
        stepang = 2
    if window is None:
        # by default whatevers smaller,
        # half trace length or 10 * max shift
        # ensure window is odd length
        window = int(np.min([data.shape[1] * 0.5,maxlag * 10]))
        window = window if window%2==1 else window + 1

    degs, lags = np.meshgrid(np.arange(0,180,stepang),
                             np.arange(0,maxlag,steplag).astype(int))

    shape = degs.shape
    lam1 = np.zeros(shape)
    lam2 = np.zeros(shape)
    for ii in np.arange(shape[1]):
        temp = c.rotate(data,degs[0,ii])
        for jj in np.arange(shape[0]):
            # remove splitting so use inverse operator (negative lag)
            temp2 = c.lag(temp,-lags[jj,ii])
            temp3 = c.window(temp2,window)
            lam2[jj,ii], lam1[jj,ii] = eigvalcov(temp3)
            
    return data,degs,lags,lam1,lam2,window

def grideigval(data, maxlag=None, window=None, stepang=None, steplag=None):
    data,degs,lags,lam1,lam2,window = _grideigval(data)
    return Measurement(data=data,degs=degs,lags=lags,lam1=lam1,lam2=lam2,window=window)
    

def ndf(y,taper=False,detrend=True):
    """
    Estimates number of degrees of freedom using noise trace y.
    Uses the improvement found by Walsh et al (2013).
    By default will detrend data to ensure zero mean
    and will taper edges using a Tukey filter affecting amplitudes of data at edges (extreme 5%)
    """

    if taper is True:
        y = y * signal.tukey(y.size,0.05)
        
    if detrend is True:
        # ensure no trend on the noise trace
        y = signal.detrend(y)

  
    Y = np.fft.fft(y)
    amp = np.absolute(Y)
    
    # estimate E2 and E4 following Walsh et al (2013)
    a = np.ones(Y.size)
    a[0] = a[-1] = 0.5
    E2 = np.sum( a * amp**2)
    E4 = (np.sum( (4 * a**2 / 3) * amp**4))
    
    ndf = 2 * ( 2 * E2**2 / E4 - 1 )
    
    return ndf
    
def Ftest(lam2,ndf,alpha=0.05):
    """
    returns lambda2 value at 100(1-alpha)% confidence interval
    by default alpha=0.05 yielding the value of lambda2 at 95% confidence interval
    """
    lam2min = lam2.min()
    k = 2 # two parameters, phi and dt.
    R = ((lam2 - lam2min)/k) /  (lam2min/(ndf-k))
    F = stats.f.ppf(1-alpha,k,ndf)
    lam2alpha = lam2min * ( 1 + (k/(ndf-k)) * F)
    return lam2alpha
    
def _synthM(deg=25,lag=10):
    P = c.Pair()
    P.split(deg,lag)
    return _grideigvalcov(P.data)
    

     
