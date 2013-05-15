'''PyAero.UVLM.Utils.DerivedTypes
@brief      Classes containing data for UVLM solver.
@author     Rob Simpson
@contact    r.simpson11@imperial.ac.uk 
@version    0.0
@date       30/01/2013
@pre        None
@warning    None
'''

import ctypes as ct
import numpy as np

class VMopts:
    """@brief options for UVLM
    @param M Chordwise panels.
    @param N Spanwise panels.
    @param ImageMethod Use image method across y-z plane if True.
    @param Mstar Number of wake panels.
    @param Steady Flag for specifying Steady(True) or Unsteady(False) solution.
    @param KJMeth If True use Joukowski force calculation method, if false use
           generalised Katz and Plotkin.
    @param NewAIC Calculate a new AIC matrix, required if the aerodynamic 
           surface has deformed significantly.
    @param DelTime Desired time step.
    @param Rollup Flag to specify wake roll-up due to induced velocities.
    @param NumCores Number of cores for parallel computation using openMP."""
    
    def __init__(self, M, N, ImageMethod = False, Mstar = 1, Steady = True,\
                 KJMeth = True, NewAIC = True, DelTime = 0.0, Rollup = False,\
                 NumCores = 1):
        self.M = ct.c_uint(M)
        self.N = ct.c_uint(N)
        self.ImageMethod = ct.c_bool(ImageMethod)
        self.Mstar = ct.c_int(Mstar)
        self.Steady = ct.c_bool(Steady)
        self.KJMeth = ct.c_bool(KJMeth)
        self.NewAIC = ct.c_bool(NewAIC)
        self.DelTime = ct.c_double(DelTime)
        self.Rollup = ct.c_bool(Rollup)
        self.NumCores = ct.c_uint(NumCores)
        

class VMinput:
    """@brief input data for UVLM.
    @param c chord.
    @param b span.
    @param area wing area.
    @param U_infty Free-stream flow magnitude.
    @param alpha free-stream direction (AoA).
    @param theta root twist angle.
    """
    
    def __init__(self, c, b, U_mag, alpha, theta, ZetaDotTest=0.0, \
                 WakeLength = 50.0):
        self.c = c
        self.b = b
        self.area = c*b
        self.U_infty = np.zeros((3))
        self.U_infty[:] = [0.0, -U_mag*np.cos(alpha), U_mag*np.sin(alpha)]
        self.alpha = alpha
        self.theta = theta
        self.ZetaDotTest = ZetaDotTest
        self.WakeLength = WakeLength
    

if __name__ == '__main__':
    c = 1
    b = 8
    U_mag = 25.0
    alpha = 2*np.pi/180.0
    VMINPUT = VMinput(c, b, U_mag, alpha)
    print(VMINPUT.U_infty)