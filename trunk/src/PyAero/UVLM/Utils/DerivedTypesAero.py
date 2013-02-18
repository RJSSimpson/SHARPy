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
    @param Mstar Number of wake panels."""
    
    def __init__(self, M, N, ImageMethod = False, Mstar = 1, Steady = True,\
                 KJMeth = True):
        self.M = ct.c_uint(M)
        self.N = ct.c_uint(N)
        self.ImageMethod = ct.c_bool(ImageMethod)
        self.Mstar = ct.c_int(Mstar)
        self.Steady = ct.c_bool(Steady)
        self.KJMeth = ct.c_bool(KJMeth)
        

class VMinput:
    """@brief input data for UVLM.
    @param c chord.
    @param b span.
    @param area wing area.
    @param U_infty Free-stream flow magnitude.
    @param alpha free-stream direction (AoA).
    @param theta root twist angle.
    """
    
    def __init__(self, c, b, U_mag, alpha, theta,ZetaDotTest=0.0):
        self.c = c
        self.b = b
        self.area = c*b
        self.U_infty = np.zeros((3))
        self.U_infty[:] = [0.0, -U_mag*np.cos(alpha), U_mag*np.sin(alpha)]
        self.alpha = alpha
        self.theta = theta
        self.ZetaDotTest = ZetaDotTest
    

if __name__ == '__main__':
    c = 1
    b = 8
    U_mag = 25.0
    alpha = 2*np.pi/180.0
    VMINPUT = VMinput(c, b, U_mag, alpha)
    print(VMINPUT.U_infty)