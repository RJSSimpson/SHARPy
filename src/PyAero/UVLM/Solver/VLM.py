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
import PostProcess
import SharPySettings as Settings
import PyBeam.Utils.DerivedTypes as DerivedTypes
from PyBeam.Utils import BeamInit
from PyFSI.Beam2UVLM import CoincidentGrid
from PyFSI.Beam2UVLM import InitSection
import PyAero.UVLM.Utils.DerivedTypesAero as DerivedTypesAero

def InitSteadyGrid(VMOPTS,VMINPUT):
    """@brief Initialise steady grid and zero grid velocities."""
    
    # Define theta using CRV, used for wing twist
    # use EA here"
    Psi = np.zeros((3))
    Psi[0] = VMINPUT.theta
    # Get rotation matrix
    R = Psi2TransMat(Psi)
    
    # Define grid nodes
    DeltaC = VMINPUT.c/VMOPTS.M.value
    DeltaS = VMINPUT.b/VMOPTS.N.value
    Zeta = np.zeros((VMOPTS.M.value+1,VMOPTS.N.value+1,3),ct.c_double,'C')
    for i in range(VMOPTS.M.value+1):
        for j in range(VMOPTS.N.value+1):
            Zeta[i][j][:] = np.dot(R,[j*DeltaS,-i*DeltaC,0.0])
        #END for j
    #END for i
    
    # Define zeros for surface motion
    ZetaDot = np.zeros((VMOPTS.M.value+1,VMOPTS.N.value+1,3),ct.c_double,'C')
    
    # Apply y-dir motion for tests
    if VMINPUT.ZetaDotTest != 0.0:
        for i in range(VMOPTS.M.value+1):         
            for j in range(VMOPTS.N.value+1):
                ZetaDot[i,j,1] = VMINPUT.ZetaDotTest
            # END for j
        #END for j
    #END if
    
    # Define gamma.
    Gamma = np.zeros((VMOPTS.M.value,VMOPTS.N.value),ct.c_double,'C')
    
    # Define arrays for outputs.
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

class MyStruct():
    """@brief Default empty class for use as a general purpose struct."""
    def __init__(self):
        """Init with no attributes."""
        


def Run_Cpp_Solver_VLM(VMOPTS, VMINPUT, VMUNST = None, AELOPTS = None):
    
    # init grid
    if VMUNST == None:
        # Set up variables assuming no unsteady motion
        VMUNST = MyStruct()
        VMUNST.PsiA_G = np.array([0.0, 0.0,0.0])
        VMUNST.VelA_G = np.array([0.0, 0.0, 0.0])
        VMUNST.OmegaA_G = np.array([0.0, 0.0, 0.0])
        VMUNST.OriginA_G = np.array([0.0, 0.0, 0.0])
        
    if AELOPTS == None:
        # Set up defaults for grid generation.
        AELOPTS = MyStruct()
        AELOPTS.ElasticAxis = 0.0   
     
    # Initialise DerivedTypes for PyBeam initialization. 
    XBOPTS = DerivedTypes.Xbopts()
    XBINPUT = DerivedTypes.Xbinput(2, VMOPTS.N.value, VMINPUT.b)  
    # Create a dummy stiffness matrix to avoid errors.
    XBINPUT.BeamStiffness = np.eye(6, dtype = ct.c_double)
    # Use PyBeam to define reference beam (elastic axis).
    XBINPUT, XBOPTS, NumNodes_tot, XBELEM, PosIni, PsiIni,\
            XBNODE, NumDof = BeamInit.Static(XBINPUT, XBOPTS)
    # Delete unnecessary variables.
    del XBELEM, XBNODE, NumDof  
    # Copy reference beam to current (deformed) variables.
    PosDefor = PosIni.copy(order = 'F')
    PsiDefor = PsiIni.copy(order = 'F')
    # Declare empty array for beam DoF rates.
    PosDotDef = np.zeros_like(PosDefor, ct.c_double, 'F')
    PsiDotDef = np.zeros_like(PsiDefor, ct.c_double, 'F')
       
    # Check the specified inputs from the PyAero have been properly applied.
    assert NumNodes_tot.value - 1 == VMOPTS.N.value, "Initialisation wrong"
    
    # Initialise section coordinates.
    Section = InitSection(VMOPTS,VMINPUT,AELOPTS.ElasticAxis)
    
    # Initialise origin and orientation of surface velocities in a-frame
    CGa = Psi2TransMat(VMUNST.PsiA_G)
    VelA_A = np.dot(CGa.T,VMUNST.VelA_G)
    OmegaA_A = np.dot(CGa.T,VMUNST.OmegaA_G)
    
    # Declare empty array for aerodynamic grid and velocities.
    Zeta = np.zeros((VMOPTS.M.value + 1,VMOPTS.N.value + 1,3), ct.c_double, 'C')
    ZetaDot = np.zeros((VMOPTS.M.value + 1, VMOPTS.N.value + 1, 3),
                       ct.c_double, 'C')
    
    # Initialise aerodynamic grid and velocities.
    CoincidentGrid(PosDefor, PsiDefor, Section,
                   VelA_A, OmegaA_A, PosDotDef, PsiDotDef,
                   XBINPUT, Zeta, ZetaDot, VMUNST.OriginA_G, VMUNST.PsiA_G,
                   VMINPUT.ctrlSurf)
    
    # init wake
    ZetaStar, GammaStar = InitSteadyWake(VMOPTS,VMINPUT, Zeta)
    
    # init external velocities  
    Uext = InitSteadyExternalVels(VMOPTS,VMINPUT)
    
    # Solver variables.
    Gamma = np.zeros((VMOPTS.M.value,VMOPTS.N.value),ct.c_double,'C')
    Forces = np.zeros((VMOPTS.M.value+1,VMOPTS.N.value+1,3),ct.c_double,'C')
    
    # Solve
    UVLMLib.Cpp_Solver_VLM(Zeta, ZetaDot, Uext, ZetaStar, VMOPTS,
                           Forces, Gamma, GammaStar)
    
    
    # Print tecplot file to check wake and grid etc
    Variables=['X', 'Y', 'Z','Gamma']
    
    # write header
    
    Filename = Settings.OutputDir + Settings.OutputFileRoot + 'AeroGrid.dat'
    
    FileObject = PostProcess.WriteAeroTecHeader(Filename,
                                                'Default',
                                                Variables)
    
    # write surface zone data
    PostProcess.WriteAeroTecZone(FileObject, 'Surface', Zeta,
                                -1, 0,
                                0.0, Variables, False, Gamma)
    
    # write wake data
    PostProcess.WriteAeroTecZone(FileObject, 'Wake', ZetaStar,
                                -1, 0,
                                0.0, Variables, False, GammaStar)
    
    # close file
    PostProcess.CloseAeroTecFile(FileObject)
    
    if Settings.WriteUVLMdebug == True:
        debugFile = Settings.OutputDir + 'gamma_1degBeta.dat'
        np.savetxt(debugFile, Gamma.flatten('C'))
    
    # post process to get coefficients
    return PostProcess.GetCoeffs(VMOPTS, Forces, VMINPUT, VMUNST.VelA_G)


if __name__ == '__main__':
    # Create inputs
    from Goland import *
    
    # Re-define control Surface
    ctrlSurf = DerivedTypesAero.ControlSurf(iMin = M - M/4,
                                            iMax = M,
                                            jMin = N - N/4,
                                            jMax = N,
                                            typeMotion = 'cos',
                                            betaBar = 1.0*np.pi/180.0,
                                            omega = 0.0)
    # Inputs.
    WakeLength = 100.0
    VMINPUT = DerivedTypesAero.VMinput(c = c,
                                   b = XBINPUT.BeamLength,
                                   U_mag = Umag,
                                   alpha = 0.0*np.pi/180.0,
                                   theta = 0.0,
                                   WakeLength = WakeLength,
                                   ctrlSurf = ctrlSurf)
    
    # Redefine VMOPTS
    VMOPTS.Mstar = ct.c_int(1)
    VMOPTS.Steady = ct.c_bool(True)
    VMOPTS.ImageMethod = ct.c_bool(True)
    # run solver
    Coeffs = Run_Cpp_Solver_VLM(VMOPTS,VMINPUT)
    
    print(Coeffs)
    