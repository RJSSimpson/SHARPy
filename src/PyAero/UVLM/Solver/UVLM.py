'''
@package    PyAero.UVLM.UVLM
@brief      UVLM solution.
@author     Rob Simpson
@contact    r.simpson11@imperial.ac.uk 
@version    0.0
@date       18/02/2013
@pre        None
@warning    None
'''

import UVLMLib
import numpy as np
from DerivedTypesAero import VMopts, VMinput, ControlSurf
import ctypes as ct
from XbeamLib import Psi2TransMat
import DerivedTypes
import PostProcess
import SharPySettings as Settings
import os
from VLM import InitSteadyExternalVels, InitSteadyWake
import BeamInit
from PyFSI.Beam2UVLM import InitSection, CoincidentGrid
from PyCoupled.Utils.DerivedTypesAeroelastic import AeroelasticOps

class VMUnsteadyInput:
    """@brief Contains data for unsteady run of UVLM.
    @param WakeLength Length of wake in chordlengths.
    @param DelS non-dim timestep s = omega*c/U.
    @param NumChordLengths Number of chord lengths to travel 
    in prescribed simulation.
    @param VelA_G Velocity of reference frame.
    @param OmegaA_G Initial angular vel of reference frame.
    @param OriginA_G Origin of reference frame in G-frame.
    @param PsiA_G Orientation of reference frame in G-frame."""
    
    def __init__(self, VMOPTS, VMINPUT, WakeLength,\
                 DelS, NumChordLengths,\
                 VelA_G, OmegaA_G,\
                 OriginA_G = np.zeros((3),ct.c_double),\
                 PsiA_G = np.zeros((3),ct.c_double)):
        
        self.WakeLength = WakeLength
        self.NumChordLengths = NumChordLengths
        if np.linalg.norm(VelA_G) != 0.0:
            self.VelA_G = VelA_G/np.linalg.norm(VelA_G)
        else:
            self.VelA_G = np.array([0.0,0.0,0.0])
        self.OmegaA_G = OmegaA_G
        self.OriginA_G = OriginA_G
        self.PsiA_G = PsiA_G
        
        # Physical time step information."
        self.NumSteps = NumChordLengths/DelS 
        self.FinalTime = NumChordLengths*VMINPUT.c / \
                    np.linalg.norm(VMINPUT.U_infty-VelA_G)
        self.DelTime = self.FinalTime/self.NumSteps
        
        # Check Mstar is high enough to discretize wake.
        if self.DelTime != self.WakeLength*VMINPUT.c/VMOPTS.Mstar.value:
            print("DelTime requested is ",self.DelTime,\
                  "\nDelWakePanel is ",self.WakeLength*VMINPUT.c/VMOPTS.Mstar.value,\
                  "\nChanging Mstar to ", self.WakeLength*VMINPUT.c/self.DelTime)
            VMOPTS.Mstar.value = int(self.WakeLength*VMINPUT.c/self.DelTime)
            
        # Set DelTime for VMOPTS
        VMOPTS.DelTime = ct.c_double(self.DelTime)
        

def Run_Cpp_Solver_UVLM(VMOPTS,VMINPUT,VMUNST,AELOPTS):
    """@brief UVLM solver with prescribed inputs."""
    
    
    # Initialise DerivedTypes for PyBeam initialisation. 
    XBOPTS = DerivedTypes.Xbopts()
    XBINPUT = DerivedTypes.Xbinput(2, VMOPTS.N.value, VMINPUT.b)  
    # Create a dummy stiffness matrix to avoid errors.
    XBINPUT.BeamStiffness = np.eye(6,6,0,ct.c_double)
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
    VelA_A = np.dot(CGa,VMUNST.VelA_G)
    OmegaA_A = np.dot(CGa,VMUNST.OmegaA_G)
    
    # Declare empty array for aerodynamic grid and velocities.
    Zeta = np.zeros((VMOPTS.M.value + 1,VMOPTS.N.value + 1,3), ct.c_double, 'C')
    ZetaDot = np.zeros((VMOPTS.M.value + 1, VMOPTS.N.value + 1, 3),
                       ct.c_double, 'C')
    
    # Initialise aerodynamic grid and velocities.
    CoincidentGrid(PosDefor, PsiDefor, Section,
                   VelA_A, OmegaA_A, PosDotDef, PsiDotDef,
                   XBINPUT, Zeta, ZetaDot, VMUNST.OriginA_G, VMUNST.PsiA_G,
                   VMINPUT.ctrlSurf)
    
    # Initialise wake for unsteady solution.
    ZetaStar, GammaStar = InitSteadyWake(VMOPTS, VMINPUT, Zeta,
                                         VMUNST.VelA_G)
    
    # Initialise external velocities.  
    Uext = InitSteadyExternalVels(VMOPTS,VMINPUT)
    
    # Declare empty solver variables.
    Forces = np.zeros_like(Zeta, ct.c_double, 'C')
    Gamma = np.zeros((VMOPTS.M.value,VMOPTS.N.value),ct.c_double,'C')
    AIC = np.zeros((VMOPTS.M.value*VMOPTS.N.value,
                   VMOPTS.M.value*VMOPTS.N.value),
                   ct.c_double,'C')
    BIC = np.zeros((VMOPTS.M.value*VMOPTS.N.value,
                   VMOPTS.M.value*VMOPTS.N.value),
                   ct.c_double,'C')
    
    # Open tecplot file object."
    Variables=['X', 'Y', 'Z','Gamma']
    Filename = Settings.OutputDir + Settings.OutputFileRoot + 'AeroGrid.dat'
    FileObject = PostProcess.WriteAeroTecHeader(Filename,
                                                'Default',
                                                Variables)
    
    # Initialise vector of time steps.
    Time = np.arange(0.0,VMUNST.FinalTime+VMUNST.DelTime,VMUNST.DelTime)
    # Create Array for storing time and coefficient data.
    CoeffHistory = np.zeros((len(Time),4))
    
    # Loop through time steps.
    for iTimeStep in range(len(Time)):
        
        # Set forces array to zero (+= operator used in C++ library).
        Forces[:,:,:] = 0.0
        
        if iTimeStep > 0:
            
            # Update geometry.
            VMUNST.OriginA_G[:] += VMUNST.VelA_G[:]*VMUNST.DelTime
            # TODO: update OmegaA_A in pitching problem.
            
            # Update control surface defintion.
            if VMINPUT.ctrlSurf != None:
                VMINPUT.ctrlSurf.update(Time[iTimeStep])
            
            # Update aerodynamic surface.
            CoincidentGrid(PosDefor, PsiDefor, Section,\
                       VelA_A, OmegaA_A, PosDotDef, PsiDotDef,
                       XBINPUT, Zeta, ZetaDot, VMUNST.OriginA_G, VMUNST.PsiA_G,
                       VMINPUT.ctrlSurf)
            
            # Convect wake downstream.        
            ZetaStar = np.roll(ZetaStar,1,axis = 0)
            GammaStar = np.roll(GammaStar,1,axis = 0)
            # Overwrite 1st row with with new trailing-edge position.
            ZetaStar[0,:] = Zeta[VMOPTS.M.value,:]
            # Overwrite Gamma with TE value from previous time step.
            GammaStar[0,:] = Gamma[VMOPTS.M.value-1,:]
            
            # set NewAIC to false.
            if VMINPUT.ctrlSurf == False:
                VMOPTS.NewAIC = ct.c_bool(False)
            
        # END if iTimeStep > 1
        
        
        # Calculate forces on aerodynamic grid.
        UVLMLib.Cpp_Solver_VLM(Zeta, ZetaDot, Uext, ZetaStar, VMOPTS,
                               Forces, Gamma, GammaStar, AIC, BIC)
         
        #print(PostProcess.GetCoeffs(VMOPTS, Forces, VMINPUT, VMUNST.VelA_G))
        
        CoeffHistory[iTimeStep,0] = Time[iTimeStep]
        CoeffHistory[iTimeStep,1:] = PostProcess.GetCoeffs(VMOPTS, Forces,
                                                           VMINPUT,
                                                           VMUNST.VelA_G)
        
        # Write aerodynamic surface data as tecplot zone data.
        PostProcess.WriteAeroTecZone(FileObject, 'Surface', Zeta,\
                                     iTimeStep, len(Time),\
                                     Time[iTimeStep], Variables, False, Gamma)
        
        # Write wake data as tecplot zone data.
        PostProcess.WriteAeroTecZone(FileObject, 'Wake', ZetaStar,\
                                    iTimeStep, len(Time),\
                                    Time[iTimeStep], \
                                    Variables, False, GammaStar)
        
        # Rollup due to external velocities.
        ZetaStar[:,:] += VMINPUT.U_infty * VMOPTS.DelTime.value
        
    # END for iTimeStep
    
    # Close tecplot file object.
    PostProcess.CloseAeroTecFile(FileObject)
    
    return CoeffHistory
        

if __name__ == '__main__':
    # Define options for aero solver.
    M = 4
    N = 10
    Mstar = 80
    VMOPTS = VMopts(M,N,True,Mstar,False,False,True,0.0,True,4)
    
    # Define wing geometry parameters.
    c = 1
    b = 4 #semi-span
    
    # Define Control Surface
    ctrlSurf = ControlSurf(3, 4, 6, 9, 'sin', 1*np.pi/180.0, 2*np.pi)
    
    # Define free stream conditions.
    U_mag = 1.0
    alpha = 0.0*np.pi/180.0
    theta = 10.0*np.pi/180.0
    
    VMINPUT = VMinput(c, b, U_mag, alpha, theta, 0.0, 15.0,  ctrlSurf)
    
    # Define unsteady solver parameters.
    WakeLength = 10.0
    DeltaS = 1.0/4.0
    NumChordLengths = 10.0
    VelA_A = np.array([0, 0, 0])
    OmegaA_A = np.array([0, 0, 0])
    VMUNST = VMUnsteadyInput(VMOPTS,VMINPUT,\
                             WakeLength,\
                             DeltaS,\
                             NumChordLengths,\
                             VelA_A, OmegaA_A)
    
    # Define 'aeroelastic' options.
    AELOPTS = AeroelasticOps(ElasticAxis = -0.5,
                             gForce =  0.0,
                             AirDensity = 1.0)
    
    # Run C++ solver
    Coeffs = Run_Cpp_Solver_UVLM(VMOPTS,VMINPUT,VMUNST,AELOPTS)
    
    print(Coeffs)
    