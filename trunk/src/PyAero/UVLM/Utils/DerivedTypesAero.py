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

class VMopts:
    """@brief data for UVLM
    @param M Chordwise panels.
    @param N Spanwise panels."""
    
    def __init__(self, M, N, ImageMethod = False):
        self.M = ct.c_uint(M)
        self.N = ct.c_uint(N)
        self.ImageMethod = ct.c_bool(ImageMethod)
    

if __name__ == '__main__':
    pass