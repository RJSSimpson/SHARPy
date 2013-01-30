'''@package PyAero.UVLM.Utils.Vorticity
@brief      Low level functions for evaluating the biot-savart law.
@author     Rob Simpson
@contact    r.simpson11@imperial.ac.uk 
@version    0.0
@date       29/01/2013
@pre        None
@warning    None
'''

import numpy as np
import ctypes as ct
from math import pow
from numpy.linalg import norm
import UVLMLib
from Misc import Timer
import matplotlib.pyplot as plt

def BiotSegment(xP,x1,x2,Gamma,Uind):
    """@brief Velocity induced by segment x1->x2 at point xP."""
    
    if Gamma == 0.0:
        return np.array([0.0, 0.0, 0.0])
    
    
    "Define position vectors"
    r0 = x2 - x1
    r1 = xP - x1
    r2 = xP - x2
    rx = np.cross(r1,r2)
    
    
    "Check if within core radius"
    if norm(r1) <= 1.0e-5 or norm(r2) <= 10e-5 or norm(rx) <= 1.0e-5:     
        return np.array([0.0, 0.0, 0.0])
    
    
    "Calculate 'K' (Katz and Plotkin)"
    r0d1 = np.dot(r0,r1)
    r0d2 = np.dot(r0,r2)
    K = (Gamma/(4.0*np.pi*pow(norm(rx),2.0))) \
        * (r0d1/norm(r1) - r0d2/norm(r2))
    
    
    "overwrite induced velocity"
    Uind[0] = K*rx[0]
    Uind[1] = K*rx[1]
    Uind[2] = K*rx[2]
    
    
if __name__ == '__main__':
    "test speed of BiotSavart"
    
    NumPanels = [10, 40, 160]
    
    x1 = np.array([-0.5, 0.0, 0.0])
    x2 = np.array([0.5, 0.0, 0.0])
    xP = np.array([0.0, 0.0, -1.0])
    Gamma = 1.0
    Uind = np.array([0.0, 0.0, 0.0])
    
    "Initialise timing result array"
    Times = np.zeros((3,6))
              
    Counter = 0
    for Test in range(len(NumPanels)):
        
        TotalCalls = int((NumPanels[Test]*4.0)**2.0)
        #test Py
        with Timer() as t_Py:
            for i in range(TotalCalls):
                BiotSegment(xP,x1,x2,Gamma,Uind)
                
        print('Request took %.03f sec.' % t_Py.interval)
        
        #test C++ from Py
        with Timer() as t_Cpp:
            for i in range(TotalCalls):
                UVLMLib.Cpp_Vorticity_BiotSegment(xP,x1,x2,Gamma,Uind)
    
        print('Request took %.03f sec.' % t_Cpp.interval)
        
        """#test C++ using Eigen::Map: NOTE: slightly worse actually!
        with Timer() as t_CppMap:
            for i in range(TotalCalls):
                UVLMLib.Cpp_Vorticity_BiotSegment_Map(xP,x1,x2,Gamma,Uind)
    
        print('Request took %.03f sec.' % t_CppMap.interval)"""
        
        #test with C from Py
        with Timer() as t_C:
            for i in range(TotalCalls):
                UVLMLib.C_Vorticity_BiotSegment(xP,x1,x2,Gamma,Uind)
                
        print('Request took %.03f sec.' % t_C.interval)
                
        """#test with C direct call: NOTE: this makes no difference
        with Timer() as t_Cdirect:
            for i in range(TotalCalls):
                UVLMLib.c_vorticity_biotsegment(xP.ctypes.data_as(ct.POINTER(ct.c_double)), \
                                x1.ctypes.data_as(ct.POINTER(ct.c_double)), \
                                x2.ctypes.data_as(ct.POINTER(ct.c_double)), \
                                ct.byref(ct.c_double(Gamma)), \
                                Uind.ctypes.data_as(ct.POINTER(ct.c_double)) )
                
        print('Request took %.03f sec.' % t_Cdirect.interval)"""
        
        #test with Cpp loop in cpp
        with Timer() as t_CppCloop:
            UVLMLib.Cpp_test_biotsegment(TotalCalls)
            
        print('Request took %.03f sec.' % t_CppCloop.interval)
        
        #test with C loop in cpp
        with Timer() as t_CCloop:
            UVLMLib.C_test_biotsegment(TotalCalls)
            
        print('Request took %.03f sec.' % t_CCloop.interval)
        
        Times[Test,:] = [NumPanels[Test], t_Py.interval, t_Cpp.interval, \
                         t_C.interval, t_CppCloop.interval, t_CCloop.interval]
    #END test loop
    print(Times)
    plt.figure(1)
    plt.plot(Times[:,0],Times[:,1],'ko-')
    plt.plot(Times[:,0],Times[:,2],'ks-')
    plt.plot(Times[:,0],Times[:,3],'k^-')
    plt.plot(Times[:,0],Times[:,4],'rs-')
    plt.plot(Times[:,0],Times[:,5],'r^-')
    plt.xlabel('Number of Panels')
    plt.ylabel('Time, s.')
    plt.legend(('Python',\
               'Python calling C++ library',\
               'Python calling C library',\
               'C++ code',\
               'C code'),'upper left')
    plt.show()