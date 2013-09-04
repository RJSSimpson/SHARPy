'''
Created on 15 Feb 2013

@author: rjs10
'''

class AeroelasticOps:
    """@brief Options specific to aeroelastic calculations.
    @param ElasticAxis Same as Theodorsen's a variable in aerofoil section.
    @param InertialAxis Same as Theodorsen's a variable in aerofoil section.
    @param Tight Attempt tightly-coupled aeroelastic simulation (untested)
    @param ImpStart Start nonlinear aeroelastic simulation with impulsive
                    start in aero solver. TODO: move to aero inputs.
    """
    
    def __init__(self, ElasticAxis, InertialAxis, 
                  AirDensity = 1.0,
                  Tight = False, ImpStart = False):
        
        assert (-1.0 <= ElasticAxis and ElasticAxis <= 1.0),\
                'Elastic axis prescribed outside section'
        assert (-1.0 <= InertialAxis and InertialAxis <= 1.0),\
                'Inertial axis prescribed outside section'
        assert (AirDensity >= 0.0),'Negative air density requested'
        
        self.ElasticAxis = ElasticAxis
        self.InertialAxis = InertialAxis
        self.AirDensity = AirDensity
        self.Tight = Tight
        self.ImpStart = ImpStart

if __name__ == '__main__':
    AeroelasticOps(0.2,0.4,1.20)