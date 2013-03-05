'''@package PyAero.UVLM.UVLM
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
from DerivedTypesAero import VMopts, VMinput
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

Settings.OutputDir = os.getcwd() + '/'
Settings.OutputFileRoot = ''

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
        
        "physical timestep info"
        self.NumSteps = NumChordLengths/DelS 
        self.FinalTime = NumChordLengths*VMINPUT.c / \
                    np.linalg.norm(VMINPUT.U_infty-VelA_G)
        self.DelTime = self.FinalTime/self.NumSteps
        
        "check Mstar is big enough for simulation"
        if self.DelTime != self.WakeLength*VMINPUT.c/VMOPTS.Mstar.value:
            print("DelTime requested is ",self.DelTime,\
                  "\nDelWakePanel is ",self.WakeLength*VMINPUT.c/VMOPTS.Mstar.value,\
                  "\nChanging Mstar to ", self.WakeLength*VMINPUT.c/self.DelTime)
            VMOPTS.Mstar.value = int(self.WakeLength*VMINPUT.c/self.DelTime)
            
        "Set DelTime for VMOPTS"
        VMOPTS.DelTime = ct.c_double(self.DelTime)
        

def Run_Cpp_Solver_UVLM(VMOPTS,VMINPUT,VMUNST,AELOPTS):
    """@brief UVLM solver with prescribed inputs."""
    
    "define reference line as in PyBeam"
    XBOPTS = DerivedTypes.Xbopts()
    
    XBINPUT = DerivedTypes.Xbinput(2, VMOPTS.N.value, VMINPUT.b)
    "dummy stiffness matrix"
    XBINPUT.BeamStiffness = np.eye(6,6,0,ct.c_double)
    
    XBINPUT, XBOPTS, NumNodes_tot, XBELEM, PosIni, PsiIni,\
            XBNODE, NumDof = BeamInit.Static(XBINPUT, XBOPTS)
            
    PosDefor = PosIni.copy(order = 'F')
    PsiDefor = PsiIni.copy(order = 'F')
    PosDotDef = np.zeros_like(PosDefor, ct.c_double, 'F')
    PsiDotDef = np.zeros_like(PsiDefor, ct.c_double, 'F')   
    
    assert NumNodes_tot.value - 1 == VMOPTS.N.value, 'initialisatiion wrong'
    
    "init grid based section + definition of reference line"
    Section = InitSection(VMOPTS,VMINPUT,AELOPTS.ElasticAxis)
    
    Zeta = np.zeros((VMOPTS.M.value + 1,VMOPTS.N.value + 1,3),ct.c_double,'C')
    ZetaDot = np.zeros((VMOPTS.M.value + 1,VMOPTS.N.value + 1,3),ct.c_double,'C')
    
    
    "Initialise origin and orientation of VEls in a-frame"
    CGa = Psi2TransMat(VMUNST.PsiA_G)
    
    VelA_A = np.dot(CGa,VMUNST.VelA_G)
    OmegaA_A = np.dot(CGa,VMUNST.OmegaA_G)
    
    "get aerogrid and velocities"
    CoincidentGrid(PosDefor, PsiDefor, Section,\
                   VelA_A, OmegaA_A, PosDotDef, PsiDotDef,
                   XBINPUT, Zeta, ZetaDot, VMUNST.OriginA_G, VMUNST.PsiA_G)
    
#    print("Section=\n",Section,"\n")
#    print("Zeta=\n",Zeta,"\n")
#    print("ZetaDot=\n",ZetaDot,"\n")
    
    "TODO: Wake length is a function of inputs delta_tprime & Mstar"
    "init wake for unsteady solution"
    ZetaStar, GammaStar = InitSteadyWake(VMOPTS, VMINPUT, Zeta,\
                                         VMUNST.VelA_G)
    
    "init external velocities"  
    Uext = InitSteadyExternalVels(VMOPTS,VMINPUT)
    
    
    "solver variables"
    Forces = np.zeros_like(Zeta, ct.c_double, 'C')
    Gamma = np.zeros((VMOPTS.M.value,VMOPTS.N.value),ct.c_double,'C')
    AIC = np.zeros((VMOPTS.M.value*VMOPTS.N.value,VMOPTS.M.value*VMOPTS.N.value), \
                   ct.c_double,'C')
    BIC = np.zeros((VMOPTS.M.value*VMOPTS.N.value,VMOPTS.M.value*VMOPTS.N.value), \
                   ct.c_double,'C')

    Time = np.arange(0.0,VMUNST.FinalTime+VMUNST.DelTime,VMUNST.DelTime)
    
    
    "Print tecplot file to check wake and grid etc"
    Variables=['X', 'Y', 'Z','Gamma']
    Filename = Settings.OutputDir + Settings.OutputFileRoot + 'AeroGrid.dat'
    FileObject = PostProcess.WriteAeroTecHeader(Filename,\
                                                'Default',\
                                                Variables)
    
    "initial geometry and wake"
    ZetaStar, GammaStar = InitSteadyWake(VMOPTS, VMINPUT, Zeta, VMUNST.VelA_G)
    
    "Init time-history"
    CoeffHistory = np.zeros((len(Time),4))
    
    "loop through timesteps"
    for iTimeStep in range(len(Time)):
        
        "zero solver variables"
        Forces[:,:,:] = 0.0
        
        if iTimeStep > 0:
            "update geometry"
            VMUNST.OriginA_G[:] += VMUNST.VelA_G[:]*VMUNST.DelTime
            
            "TODO: update OmegaA_A in pitching problem."
            
            "surface"
            CoincidentGrid(PosDefor, PsiDefor, Section,\
                       VelA_A, OmegaA_A, PosDotDef, PsiDotDef,
                       XBINPUT, Zeta, ZetaDot, VMUNST.OriginA_G, VMUNST.PsiA_G)
            
            "update wake geom"        
            "'roll' data"
            ZetaStar = np.roll(ZetaStar,1,axis = 0)
            GammaStar = np.roll(GammaStar,1,axis = 0)
            "overwrite grid points with TE"
            ZetaStar[0,:] = Zeta[VMOPTS.M.value,:]
            "overwrite Gamma with TE value from previous timestep"
            GammaStar[0,:] = Gamma[VMOPTS.M.value-1,:]
            
            
            "set NewAIC to false"
            VMOPTS.NewAIC = ct.c_bool(False)
            
        # END if iTimeStep > 1
        
        
        "Solve (TODO: Rollup)"
        UVLMLib.Cpp_Solver_VLM(Zeta, ZetaDot, Uext, ZetaStar, VMOPTS, \
                               Forces, Gamma, GammaStar, AIC, BIC)
        
        
        #print(PostProcess.GetCoeffs(VMOPTS, Forces, VMINPUT, VMUNST.VelA_G))
        
        CoeffHistory[iTimeStep,0] = Time[iTimeStep]
        CoeffHistory[iTimeStep,1:] = PostProcess.GetCoeffs(VMOPTS, Forces, \
                                                 VMINPUT, VMUNST.VelA_G)
        
        "write surface zone data"
        PostProcess.WriteAeroTecZone(FileObject, 'Surface', Zeta,\
                                     iTimeStep, len(Time),\
                                     Time[iTimeStep], Variables, False, Gamma)
        
        "write wake data"
        PostProcess.WriteAeroTecZone(FileObject, 'Wake', ZetaStar,\
                                    iTimeStep, len(Time),\
                                    Time[iTimeStep], \
                                    Variables, False, GammaStar)
        
        "Rollup due to external velocities"
        ZetaStar[:,:] += VMINPUT.U_infty*VMOPTS.DelTime.value
        
    # END for iTimeStep
    
    "close file"
    PostProcess.CloseAeroTecFile(FileObject)
    
    "return time history"
    return CoeffHistory

if __name__ == '__main__':
    "Create inputs"
    M = 4
    N = 13
    Mstar = 80
    VMOPTS = VMopts(M,N,True,Mstar,False,False,True,0.0,True,4)
    
    "define wing geom"
    c = 1
    b = 4 #semi-span
    
    "free stream conditions"
    U_mag = 1.0
    alpha = 0.0*np.pi/180.0
    theta = 10.0*np.pi/180.0
    VMINPUT = VMinput(c, b, U_mag, alpha, theta)
    
    "unsteady stuff"
    WakeLength = 10.0
    DeltaS = 1.0/16.0
    NumChordLengths = 10.0
    VelA_A = np.array([0, 0, 0])
    OmegaA_A = np.array([0, 0, 0])
    VMUNST = VMUnsteadyInput(VMOPTS,VMINPUT,\
                             WakeLength,\
                             DeltaS,\
                             NumChordLengths,\
                             VelA_A, OmegaA_A)
    
    "Aeroelastic options"
    AELOPTS = AeroelasticOps(ElasticAxis = -0.5, gForce =  0.0, AirDensity = 1.0)
    
    "run solver"
    Coeffs = Run_Cpp_Solver_UVLM(VMOPTS,VMINPUT,VMUNST,AELOPTS)
    
    print(Coeffs)
    