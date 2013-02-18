'''
Created on 15 Feb 2013

@author: rjs10
'''

class AeroelasticOps:
    """@brief Options specific to aeroelastic calculations
    @param ElasticAxis Same as Theodorsen's a variable in aerofoil section."""
    
    def __init__(self, ElasticAxis, gForce=9.81, AirDensity = 1.0):
        self.ElasticAxis = ElasticAxis
        self.gForce = gForce
        self.AirDensity = AirDensity

if __name__ == '__main__':
    pass