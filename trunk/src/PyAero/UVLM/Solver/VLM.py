'''@package PyAero.UVLM.VLM
@brief      quasi-steady VLM solution with fixed-wake.
@author     Rob Simpson
@contact    r.simpson11@imperial.ac.uk 
@version    0.0
@date       30/01/2013
@pre        None
@warning    None
'''

import UVLMLib
import numpy as np
from DerivedTypesAero import VMopts
import ctypes as ct
from XbeamLib import Psi2TransMat

def Run_Cpp_Solver_VLM():
    "Create inputs"
    M = 1
    N = 1
    VMOPTS = VMopts(M,N,True)
    
    "define wing geom"
    c = 1
    semi_b = 8
    area = c*semi_b
    
    "Define AoA using CRV"
    AoA = 10.0
    Psi = np.zeros((3))
    Psi[0] = AoA*np.pi/180.0
    "get transformation matrix"
    R = Psi2TransMat(Psi)
    
    
    "define grid nodes"
    DeltaC = c/M
    DeltaS = semi_b/N
    Zeta = np.zeros((VMOPTS.M.value+1,VMOPTS.N.value+1,3),ct.c_double,'C')
    for i in range(VMOPTS.M.value+1):
        for j in range(VMOPTS.N.value+1):
            Zeta[i][j][:] = np.dot(R,[j*DeltaS,-i*DeltaC,0.0]) #Aero grid defined differently
        #END for j
    #END for i   
    
    "define free stream velocity"
    Uinfinity = np.array([0.0, -100.0, 0.0])
    
    "define wake grid"
    ZetaStar = np.zeros((2,VMOPTS.N.value+1,3),ct.c_double,'C');
    ZetaStar[0,:] = Zeta[VMOPTS.M.value,:]
    ZetaStar[1,:] = ZetaStar[0,:]
    for j in range(VMOPTS.N.value+1):
        ZetaStar[1,j] += 10000 * (Uinfinity/np.linalg.norm(Uinfinity)) * \
                               -(Zeta[VMOPTS.M.value][0][1] - Zeta[0][0][1])
    #END for j
    
    "define zeros for surface motion"
    ZetaDot = np.zeros((VMOPTS.M.value+1,VMOPTS.N.value+1,3),ct.c_double,'C')
    
    "define external velocities (gusts)"
    Uext = np.zeros((VMOPTS.M.value+1,VMOPTS.N.value+1,3),ct.c_double,'C')
    for i in range(VMOPTS.M.value+1):
        for j in range(VMOPTS.N.value+1):
            Uext[i,j,:] = Uinfinity
        #END for j
    #END for i
    
    "define arrays for outputs"
    Forces = np.zeros((VMOPTS.M.value+1,VMOPTS.N.value+1,3),ct.c_double,'C')
    Gamma = np.zeros((VMOPTS.M.value,VMOPTS.N.value),ct.c_double,'C')
    GammaStar = np.zeros((1,VMOPTS.N.value),ct.c_double,'C')
    
    UVLMLib.Cpp_Solver_VLM(Zeta, ZetaDot, Uext, ZetaStar, VMOPTS, \
                           Forces, Gamma, GammaStar)
    
    Coeff = np.zeros((3))
    for i in range(VMOPTS.M.value+1):
        for j in range(VMOPTS.N.value+1):
            Coeff[:] += Forces[i,j,:]/(0.5*np.linalg.norm(Uinfinity)**2.0*area)
        #end for j
    #end for i
    
    print(Coeff)

if __name__ == '__main__':
    Run_Cpp_Solver_VLM()
    