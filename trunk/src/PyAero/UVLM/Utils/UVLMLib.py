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



if __name__ == '__main__':
    xP = np.array([0.0,0.0,-1.0], ct.c_double, order='C')
    x1 = np.array([-0.5,0.0,0.0], ct.c_double, order='C')
    x2 = np.array([0.5,0.0,0.0], ct.c_double, order='C')
    Gamma = 1.0
    Uind = np.array([0.0,0.0,0.0], ct.c_double, order='C')
    Cpp_Vorticity_BiotSegment_Map(xP, x1, x2, Gamma, Uind)
    print(Uind)