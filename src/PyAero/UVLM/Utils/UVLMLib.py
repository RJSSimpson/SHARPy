'''@package PyAero.UVLM.Utils.Vorticity
@brief      Low level functions for evaluating the biot-savart law.
@author     Rob Simpson
@contact    r.simpson11@imperial.ac.uk 
@version    0.0
@date       29/01/2013
@pre        None
@warning    None
'''

import ctypes as ct
import SharPySettings as Settings
import numpy as np

UVLMPath = Settings.UVLMLibDir + Settings.UVLMLibName
UVLMLib = ct.cdll.LoadLibrary(UVLMPath)

cpp_vorticity_biotsegment = UVLMLib.cpp_wrap_vorticity_biotsegment
cpp_vorticity_biotsegment_map = UVLMLib.cpp_wrap_vorticity_biotsegment_map
c_vorticity_biotsegment = UVLMLib.c_wrap_vorticity_biotsegment
cpp_test_biotsegment = UVLMLib.cpp_wrap_test_biotsegment
c_test_biotsegment = UVLMLib.c_wrap_test_biotsegment
cpp_solver_vlm = UVLMLib.cpp_wrap_solver_vlm
cpp_AIC = UVLMLib.cpp_wrap_AIC

"""ctypes does not check whether the correct number OR type of input arguments
are passed to each of these functions - great care must be taken to ensure the
number of arguments and the argument types are correct.
TODO: developer provided argtypes which must be a sequence of C data types.
http://docs.python.org/2/library/ctypes.html"""

cpp_vorticity_biotsegment.restype = None
cpp_vorticity_biotsegment_map.restype = None
c_vorticity_biotsegment.restype = None
cpp_test_biotsegment.restype = None
c_test_biotsegment.restype = None
cpp_solver_vlm.restype = None
cpp_AIC.restype = None

def AreEqual(arg1,arg2):
    """@brief Returns true if argumants are equal."""
    if arg1 == arg2:
        return True
    else:
        return False

def Cpp_Vorticity_BiotSegment(xP,x1,x2,Gamma,Uind):
    """@brief Wrapper for c_vorticity_biotsegment."""
    
    cpp_vorticity_biotsegment(xP.ctypes.data_as(ct.POINTER(ct.c_double)), \
                            x1.ctypes.data_as(ct.POINTER(ct.c_double)), \
                            x2.ctypes.data_as(ct.POINTER(ct.c_double)), \
                            ct.byref(ct.c_double(Gamma)), \
                            Uind.ctypes.data_as(ct.POINTER(ct.c_double)) )
    
    
def Cpp_Vorticity_BiotSegment_Map(xP,x1,x2,Gamma,Uind):
    """@brief Wrapper for c_vorticity_biotsegment."""
    
    cpp_vorticity_biotsegment_map(xP.ctypes.data_as(ct.POINTER(ct.c_double)), \
                            x1.ctypes.data_as(ct.POINTER(ct.c_double)), \
                            x2.ctypes.data_as(ct.POINTER(ct.c_double)), \
                            ct.byref(ct.c_double(Gamma)), \
                            Uind.ctypes.data_as(ct.POINTER(ct.c_double)) )
    
    
def C_Vorticity_BiotSegment(xP,x1,x2,Gamma,Uind):
    """@brief Wrapper for c_vorticity_biotsegment."""
    
    c_vorticity_biotsegment(xP.ctypes.data_as(ct.POINTER(ct.c_double)), \
                            x1.ctypes.data_as(ct.POINTER(ct.c_double)), \
                            x2.ctypes.data_as(ct.POINTER(ct.c_double)), \
                            ct.byref(ct.c_double(Gamma)), \
                            Uind.ctypes.data_as(ct.POINTER(ct.c_double)) )
    
    
def Cpp_test_biotsegment(NumTests):
    """@brief Wrapper for cpp_test_biotsegment."""
    
    cpp_test_biotsegment(ct.byref(ct.c_int(NumTests)))
    
    
def C_test_biotsegment(NumTests):
    """@brief Wrapper for cpp_test_biotsegment."""
    
    c_test_biotsegment(ct.byref(ct.c_int(NumTests)))


def Cpp_Solver_VLM(Zeta, ZetaDot, Uext, ZetaStar, VMOPTS, Forces, \
                   Gamma, GammaStar, AIC = None, BIC = None):
    """@details wrapper for cpp_solver_vlm."""
    
    "If memory for AIC and BIC is not allocated, allocate here."
    if AIC == None and BIC == None:
        AIC = np.zeros((VMOPTS.M.value*VMOPTS.N.value,\
                        VMOPTS.M.value*VMOPTS.N.value), \
                        ct.c_double,'C')
        BIC = np.zeros((VMOPTS.M.value*VMOPTS.N.value,\
                        VMOPTS.M.value*VMOPTS.N.value), \
                        ct.c_double,'C')
    
    cpp_solver_vlm(Zeta.ctypes.data_as(ct.POINTER(ct.c_double)),\
                   ZetaDot.ctypes.data_as(ct.POINTER(ct.c_double)),\
                   Uext.ctypes.data_as(ct.POINTER(ct.c_double)),\
                   ZetaStar.ctypes.data_as(ct.POINTER(ct.c_double)),\
                   ct.byref(VMOPTS.M),\
                   ct.byref(VMOPTS.N),\
                   ct.byref(VMOPTS.ImageMethod), \
                   ct.byref(VMOPTS.Mstar),\
                   ct.byref(VMOPTS.Steady),\
                   ct.byref(VMOPTS.KJMeth),\
                   ct.byref(VMOPTS.NewAIC),\
                   ct.byref(VMOPTS.DelTime),\
                   ct.byref(VMOPTS.Rollup),\
                   ct.byref(VMOPTS.NumCores),\
                   Forces.ctypes.data_as(ct.POINTER(ct.c_double)), \
                   Gamma.ctypes.data_as(ct.POINTER(ct.c_double)), \
                   GammaStar.ctypes.data_as(ct.POINTER(ct.c_double)), \
                   AIC.ctypes.data_as(ct.POINTER(ct.c_double)), \
                   BIC.ctypes.data_as(ct.POINTER(ct.c_double)))

def Colloc(zeta_G_panel):
    """@brief Collocation from grid position.
        
    @param zeta_G_panel A (2,2,3) array of (i,j) panel corner points.
    @return colloc panel collocation point.
    """
    
    if zeta_G_panel.shape[0] > 2 or zeta_G_panel.shape[1] > 2:
        raise TypeError("Extent of corner point array greater than 2x2(x3).")
    elif zeta_G_panel.shape[0] < 2 or zeta_G_panel.shape[1] < 2:
        raise TypeError("Extent of corner point array less than 2x2(x3).")
    
    return ( 0.25 * (zeta_G_panel[0,0,:] +
                     zeta_G_panel[0,1,:] +
                     zeta_G_panel[1,1,:] +
                     zeta_G_panel[1,0,:]  ) )
    
def Cpp_AIC(zetaSrc,mSrc,nSrc,zetaTgt,mTgt,nTgt,AIC):
    """@details Wrapper for cpp_wrap_AIC.
    @param zetaSrc Source grid points.
    @param mSrc Chordwise panels.
    @param nSrc Spanwise.
    @param zetaTgt Target grid points.
    @param mSrc Chordwise panels.
    @param nSrc Spanwise.
    @return AIC influence coefficient matrix."""
    cpp_AIC(zetaSrc.ctypes.data_as(ct.POINTER(ct.c_double)),
            ct.byref(ct.c_uint(mSrc)),
            ct.byref(ct.c_uint(nSrc)),
            zetaTgt.ctypes.data_as(ct.POINTER(ct.c_double)),
            ct.byref(ct.c_uint(mTgt)),
            ct.byref(ct.c_uint(nTgt)),
            AIC.ctypes.data_as(ct.POINTER(ct.c_double)))
    
        

if __name__ == '__main__':
    xP = np.array([0.0,0.0,-1.0], ct.c_double, order='C')
    x1 = np.array([-0.5,0.0,0.0], ct.c_double, order='C')
    x2 = np.array([0.5,0.0,0.0], ct.c_double, order='C')
    Gamma = 1.0
    Uind = np.array([0.0,0.0,0.0], ct.c_double, order='C')
    Cpp_Vorticity_BiotSegment_Map(xP, x1, x2, Gamma, Uind)
    print(Uind)