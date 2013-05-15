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
from DerivedTypesAero import VMopts, VMinput
import ctypes as ct
from XbeamLib import Psi2TransMat
import DerivedTypes
import PostProcess
import SharPySettings as Settings
import os

Settings.OutputDir = os.getcwd() + '/'
Settings.OutputFileRoot = ''

def InitSteadyGrid(VMOPTS,VMINPUT):
    """@brief Initialise steady grid and zero grid velocities."""
    
    "Define theta using CRV, used for wing twist"
    "use EA here"
    Psi = np.zeros((3))
    Psi[0] = VMINPUT.theta
    "get transformation matrix"
    R = Psi2TransMat(Psi)
    
    "define grid nodes"
    DeltaC = VMINPUT.c/VMOPTS.M.value
    DeltaS = VMINPUT.b/VMOPTS.N.value
    Zeta = np.zeros((VMOPTS.M.value+1,VMOPTS.N.value+1,3),ct.c_double,'C')
    for i in range(VMOPTS.M.value+1):
        for j in range(VMOPTS.N.value+1):
            Zeta[i][j][:] = np.dot(R,[j*DeltaS,-i*DeltaC,0.0])
        #END for j
    #END for i
    
    "define zeros for surface motion"
    ZetaDot = np.zeros((VMOPTS.M.value+1,VMOPTS.N.value+1,3),ct.c_double,'C')
    
    "apply y-dir motion for tests"
    if VMINPUT.ZetaDotTest != 0.0:
        for i in range(VMOPTS.M.value+1):         
            for j in range(VMOPTS.N.value+1):
                ZetaDot[i,j,1] = VMINPUT.ZetaDotTest
            # END for j
        #END for j
    #END if
    
    "define gamma"
    Gamma = np.zeros((VMOPTS.M.value,VMOPTS.N.value),ct.c_double,'C')
    
    "define arrays for outputs"
    Forces = np.zeros((VMOPTS.M.value+1,VMOPTS.N.value+1,3),ct.c_double,'C')
    
    return Zeta, ZetaDot, Gamma, Forces


def InitSteadyWake(VMOPTS,VMINPUT,Zeta,VelA_G = None):
    """@brief Initialise steady wake.
    @param VMOPTS options.
    @param VMINPUT inputs.
    @param Zeta Surface grid.
    @param VelA_G Velocity of surface FoR in G-frame."""
    
    #TODO: Wake length is a function of inputs delta_tprime & Mstar

    # Create empty array for wake grid."
    ZetaStar = np.zeros((VMOPTS.Mstar.value+1,VMOPTS.N.value+1,3),
                        ct.c_double,'C');
    
    # Set first row equal to trailing-edge."
    ZetaStar[0,:] = Zeta[VMOPTS.M.value,:]
    
    # Calculate incremental distance vector from trailing-edge to far-field."
    if VelA_G == None:
        DeltaX = (VMINPUT.WakeLength * 
                 (VMINPUT.U_infty / np.linalg.norm(VMINPUT.U_infty)))
    elif VelA_G != None:
        DeltaX = (VMINPUT.WakeLength * 
                 ((VMINPUT.U_infty - VelA_G) / 
                 np.linalg.norm(VMINPUT.U_infty - VelA_G)))
            
    # Divide total distance by number of wake panels."
    DeltaX = DeltaX/float(VMOPTS.Mstar.value)
    
    # Populate wake grid according to increment from trailing-edge.
    for i in range(1, VMOPTS.Mstar.value+1):         
        for j in range(VMOPTS.N.value+1):
            ZetaStar[i,j] = Zeta[VMOPTS.M.value, j] + i*DeltaX
        #END for j
    #END for i
    
    # Declare empty array for wake gamma."
    GammaStar = np.zeros((VMOPTS.Mstar.value,VMOPTS.N.value),ct.c_double,'C')
    
    return ZetaStar, GammaStar


def InitSteadyExternalVels(VMOPTS,VMINPUT):
    "@brief Initialse external velocities (free-stream only)"
    Uext = np.zeros((VMOPTS.M.value+1,VMOPTS.N.value+1,3),ct.c_double,'C')
    for i in range(VMOPTS.M.value+1):
        for j in range(VMOPTS.N.value+1):
            Uext[i,j,:] = VMINPUT.U_infty[:]
        #END for j
    #END for i
    
    return Uext


def Run_Cpp_Solver_VLM(VMOPTS,VMINPUT):
    
    "init grid"
    Zeta, ZetaDot, Gamma, Forces = InitSteadyGrid(VMOPTS,VMINPUT)
    
    "init wake"
    ZetaStar, GammaStar = InitSteadyWake(VMOPTS,VMINPUT, Zeta)
    
    "init external velocities"  
    Uext = InitSteadyExternalVels(VMOPTS,VMINPUT)
    
    
    "Solve"
    UVLMLib.Cpp_Solver_VLM(Zeta, ZetaDot, Uext, ZetaStar, VMOPTS, \
                           Forces, Gamma, GammaStar)
    
    
    "Print tecplot file to check wake and grid etc"
    Variables=['X', 'Y', 'Z','Gamma']
    
    "write header"
    
    Filename = Settings.OutputDir + Settings.OutputFileRoot + 'AeroGrid.dat'
    
    FileObject = PostProcess.WriteAeroTecHeader(Filename,\
                                                'Default',\
                                                Variables)
    
    "write surface zone data"
    PostProcess.WriteAeroTecZone(FileObject, 'Surface', Zeta,\
                                -1, 0,\
                                0.0, Variables, False, Gamma)
    
    "write wake data"
    PostProcess.WriteAeroTecZone(FileObject, 'Wake', ZetaStar,\
                                -1, 0,\
                                0.0, Variables, False, GammaStar)
    
    "close file"
    PostProcess.CloseAeroTecFile(FileObject)
    
    "post process to get coefficients"
    return PostProcess.GetCoeffs(VMOPTS, Forces, VMINPUT)


if __name__ == '__main__':
    "Create inputs"
    M = 1
    N = 1
    Mstar = 25
    VMOPTS = VMopts(M,N,False,Mstar,True,True)
    
    "define wing geom"
    c = 1
    b = 1000 #semi-span
    
    "free stream conditions"
    U_mag = 25.0
    alpha = 1.0*np.pi/180.0
    theta = 0.0*np.pi/180.0
    VMINPUT = VMinput(c, b, U_mag, alpha, theta)
    
    "run solver"
    Coeffs = Run_Cpp_Solver_VLM(VMOPTS,VMINPUT)
    
    print(Coeffs)
    