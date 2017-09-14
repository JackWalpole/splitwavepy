"""
The eigenvalue method of Silver and Chan (1991)
Uses Pair to do high level work
"""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

# from . import Pair
from . import eigval

class EigenM:
    
    """
    Silver and Chan (1991) eigenvalue method measurement.
    """
    
    def __init__(self,*args,**kwargs):
        
        if 
        Pair = Pair(*args,**kwargs)
        # if data is None:
        #     # generate synthetic Pair
        #     Pair = pair.Pair(**kwargs)
        #     #
        #     self.data,
        #     self.degs,
        #     self.lags, self.lam1, self.lam2, self.window = _synthM()
        # elif data is not None and (degs is None or
        #                            lags is None or
        #                            lam1 is None or
        #                            lam2 is None or
        #                            window is None):
        #     # generate measurement using defaults
        #     self.data,self.degs, self.lags, self.lam1, self.lam2, self.window = _grideigval(data)
        # else:
        #     # everything provided from low level (except delta?)
        #     for key,arg in kwargs.items():
        #         self.key = arg
            
            # self.data = data
            # self.degs = degs
            # self.lags = lags
            # self.lam1 = lam1
            # self.lam2 = lam2
            # self.window = window
        #
        # # ensure data is a Pair for convenience functions
        # if not isinstance(self.data,pair.Pair):
        #     self.data = c.Pair(self.data)
        
        # get some measurement attributes
        # uses ratio lam1/lam2 to find optimal fast and lag parameters
        maxloc = c.max_idx(self.lam1/self.lam2)
        self.fast = self.degs[maxloc]
        self.lag  = self.lags[maxloc]
        # generate "squashed" profiles
        self.fastprofile = np.sum(self.lam1/self.lam2, axis=0)
        self.lagprofile = np.sum(self.lam1/self.lam2, axis=1)
        # generate redefined "NI" value
        self.ni = ni(self)
        
        # get some useful stuff
        self.data_corr = c.unsplit(self.data.data,self.fast,self.lag)
        self.srcpol = c.pca(self.data_corr)
        self.srcpoldata = c.rotate(self.data.data,-self.srcpol)
        self.srcpoldata_corr = c.rotate(self.data_corr,-self.srcpol)
        
        # signal to noise ratio estimates
        # self.snr = c.snr(c.window(self.srcpoldata_corr,self.window))
        # self.snrRH = c.snrRH(c.window(self.srcpoldata_corr,self.window))
        # self.snr = np.max(self.lam1/self.lam2)
        ### if total energy = signal + noise = lam1 + lam2
        ### lam1 = signal + 1/2 noise
        ### lam2 = 1/2 noise
        ### then signal / noise = 
        self.snr = np.max((self.lam1-self.lam2)/(2*self.lam2))

        # number degrees of freedom
        self.ndf = ndf(c.window(self.srcpoldata_corr[1,:],self.window))
        # value of lam2 at 95% confidence contour
        self.lam2_95 = ftest(self.lam2,self.ndf,alpha=0.05)

        # convert traces to Pair class for convenience
        self.data_corr = c.Pair(self.data_corr)
        self.srcpoldata = c.Pair(self.srcpoldata)
        self.srcpoldata_corr = c.Pair(self.srcpoldata_corr)
        

    def plot(self,vals=None,cmap='viridis',lam2_95=False,polar=True):
        """
        plot the measurement.
        by default plots lam1/lam2 with the lambda2 95% confidence interval overlaid
        """
        
        
        
        if vals is None:
            vals = self.lam1 / self.lam2
        
        if polar is True:
            rads = np.deg2rad(np.column_stack((self.degs,self.degs+180,self.degs[:,0]+360)))
            lags = np.column_stack((self.lags,self.lags,self.lags[:,0]))
            vals = np.column_stack((vals,vals,vals[:,0]))
            fig, ax = plt.subplots(subplot_kw=dict(projection='polar'))
            ax.contourf(rads,lags,vals,50,cmap=cmap)
            ax.set_theta_direction(-1)
            ax.set_theta_offset(np.pi/2.0)
            if lam2_95 is True:
                lam2 = np.column_stack((self.lam2,self.lam2,self.lam2[:,0]))
                plt.contour(rads,lags,lam2,levels=[self.lam2_95])
        else:
            plt.contourf(self.lags,self.degs,vals,50,cmap=cmap)        
            if lam2_95 is True:
                plt.contour(self.lags,self.degs,self.lam2,levels=[self.lam2_95])
            

            
        
        plt.show()

    # def save():
    #     """
    #     Save Measurement for future referral
    #     """
    
# def _synthM(deg=25,lag=10):
#     P = c.Pair()
#     P.split(deg,lag)
#     return eigval.grideigval(P.data)

def ni(M):
    """
    measure of self-similarity in measurements at 90 degree shift in fast direction
    """
    halfway = int(M.degs.shape[1]/2)
    diff = M.fastprofile - np.roll(M.fastprofile,halfway)
    mult = M.fastprofile * np.roll(M.fastprofile,halfway)
    sumdiffsq = np.sum(diff**2)
    summult = np.sum(mult)
    return sumdiffsq/summult