class Animal:
    
    # def __init__(self):     
        
    @property
    def name(self):
        return self.__name
    
    @name.setter
    def name(self, name, surname):
        if not isinstance(name, str):
            raise Exception('nah bruv')
        self.__name = name + surname
            
    # @property
    # def delta(self):
    #     return self.__delta
    #
    # @delta.setter
    # def delta(self, delta):
    #     if delta <= 0: raise ValueError('delta must be positive')
    #     self.__delta = float(delta)
    #
    # @delta.getter
    # def