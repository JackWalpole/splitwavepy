import numpy as np

class Data:
    
    def __init__(self,*args,**kwargs):
        
        self.x = np.arange(10)
            
        
# class Grid:
#
#     def __init__(self):
#
#         self.lags = np.arange(30)
#         self.degs = np.arange(-90,90,2)
        
    # def gridsearch():

class Pair(Data):
    
    def __init__(self, *args, **kwargs):
        
        # Derive from Data
        Data.__init__(self, *args, **kwargs)
        
        # Derive from Grid
        # Grid. __init__(self, *args, **kwargs)
        
    # def eigenM(self, *args, **kwargs):
    #
    #     self.EigenM = EigenM(self)
    
    def eigen(self):
        
        self.eigm = {}
        self.eigm['']
        self.eigm['lam1'] = self.x * 2    
        
    # class EigenM(Grid):
    #
    #     def __init__(self, **kwargs):
    #
    #         Measure.__init__(self)
    #         # self.lam1 = self.x * 2
        

