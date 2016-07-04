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
from VLM import InitSteadyExternalVels, InitSteadyWake
import BeamInit
from PyFSI.Beam2UVLM import InitSection, CoincidentGrid
from PyCoupled.Utils.DerivedTypesAeroelastic import AeroelasticOps
from PyAero.UVLM.Utils.DerivedTypesAero import VMUnsteadyInput
from scipy.io import savemat
import getpass

def Run_Cpp_Solver_UVLM(VMOPTS,VMINPUT,VMUNST,AELOPTS,vOmegaHist=None):
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
    if vOmegaHist == None:
        VelA_A = np.dot(CGa.T,VMUNST.VelA_G)
        OmegaA_A = np.dot(CGa.T,VMUNST.OmegaA_G)
    else:
        VelA_A = np.dot(CGa.T,vOmegaHist[0,1:4])
        OmegaA_A = np.dot(CGa.T,vOmegaHist[0,4:7])
    
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
    if vOmegaHist == None:
        velForWake = VMUNST.VelA_G
    else:
        velForWake = vOmegaHist[0,1:4]
        
    ZetaStar, GammaStar = InitSteadyWake(VMOPTS, VMINPUT, Zeta, velForWake)
    
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
            if vOmegaHist == None:
                VMUNST.OriginA_G[:] += VMUNST.VelA_G[:]*VMUNST.DelTime
                # TODO: update OmegaA_A, PsiA_A in pitching problem.
            else:
                VMUNST.OriginA_G[:] += vOmegaHist[iTimeStep,1:4]*VMUNST.DelTime
                VelA_A = vOmegaHist[iTimeStep,1:4]
                # TODO: update OmegaA_A, PsiA_G in pitching problem.
            
            # Update control surface defintion.
            if VMINPUT.ctrlSurf != None:
                VMINPUT.ctrlSurf.update(Time[iTimeStep])
            
            # Update aerodynamic surface.
            CoincidentGrid(PosDefor,
                           PsiDefor,
                           Section,
                           VelA_A,
                           OmegaA_A,
                           PosDotDef,
                           PsiDotDef,
                           XBINPUT,
                           Zeta,
                           ZetaDot,
                           VMUNST.OriginA_G,
                           VMUNST.PsiA_G,
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
        if iTimeStep == 0:
            sav = VMOPTS.NewAIC
            VMOPTS.NewAIC = ct.c_bool(1)
        
        UVLMLib.Cpp_Solver_VLM(Zeta, ZetaDot, Uext, ZetaStar, VMOPTS,
                                   Forces, Gamma, GammaStar, AIC, BIC)
        
        if iTimeStep == 0:
            VMOPTS.NewAIC = sav
            del sav
         
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
    
    # write geometry and circ dist. to file for matlab
    if iTimeStep == len(Time)-1 and True:
        savemat('/home/rjs10/Desktop/uvlmOut.mat',
                {'zeta':Zeta.flatten('C'),
                 'zetaW':ZetaStar.flatten('C'),
                 'gamma':Gamma.flatten('C'),
                 'gammaW':GammaStar.flatten('C')},
                False,
                oned_as='column'
               )
    
    return CoeffHistory
        

if __name__ == '__main__':
    Settings.OutputDir = '/home/' + getpass.getuser() + '/Documents/MATLAB/newUVLM/nonZeroAerofoil/UVLM/'
    writeToMat = True
    
    # Define options for aero solver.
    M = 20
    N = 1
    Mstar = 10*M
    U_mag = 1.0
    alpha = 0.0*np.pi/180.0
    theta = 1.0*np.pi/180.0
    imageMeth=True
    
    # Define wing geometry parameters.
    c = 1
    b = 2000 #semi-span
    
    VMOPTS = VMopts(M,N,imageMeth,Mstar,False,True,False,c/(M*U_mag),False,4)
    
    # Define Control Surface
    ctrlSurf = None #ControlSurf(3, 4, 6, 9, 'sin', 1*np.pi/180.0, 2*np.pi)
    
    # Define free stream conditions.
    
    VMINPUT = VMinput(c, b, U_mag, alpha, theta, Mstar/M,  ctrlSurf)
    
    # Define unsteady solver parameters.
    WakeLength = Mstar/M
    DeltaS = 1/M # define wrt to #CHORD
    NumChordLengths = 200.0
    VelA_A = np.array([0, 0, 0])
    OmegaA_A = np.array([0, 0, 0])
    VMUNST = VMUnsteadyInput(VMOPTS,VMINPUT,\
                             WakeLength,\
                             DeltaS,\
                             NumChordLengths,\
                             VelA_A, OmegaA_A)
    
    # Generate history of axis sytem motions (in inertial frame)
    hBar=0.01
    k=1.0
    omegaY=2*U_mag*k/c
    Time = np.arange(0.0,VMUNST.FinalTime+VMUNST.DelTime,VMUNST.DelTime)
    vOmegaHist = np.zeros((len(Time),7))
    vOmegaHist[:,0] = Time
    vOmegaHist[:,1] = 0.0 # along A-frame x-axis
    vOmegaHist[:,2] = omegaY*hBar*np.cos(omegaY*Time) # y-axis (fore/aft)
    vOmegaHist[:,3] = 0.0 # z-axis (up/down)
    vOmegaHist[:,4] = 0.0 # about A-frame x-axis
    vOmegaHist[:,5] = 0.0 # about y
    vOmegaHist[:,6] = 0.0 # about z
    
    # Define 'aeroelastic' options.
    AELOPTS = AeroelasticOps(ElasticAxis = -0.5)
    
    # Run C++ solver
    Coeffs = Run_Cpp_Solver_UVLM(VMOPTS,VMINPUT,VMUNST,AELOPTS,vOmegaHist)
    
    print(Coeffs[-50:,:])
    
    if writeToMat == True:
        fileName = Settings.OutputDir + 'UVLMrectAR' + str(b/c) + '_m' + str(M) + 'mW' + str(Mstar) + 'n' + str(N) + 'delS' + str(2*DeltaS) + '_alpha' + str(theta)[:7]
        if imageMeth != False:
            fileName += '_half'
        savemat(fileName,
                {'vOmegaHist':vOmegaHist,
                 'Coeffs':Coeffs},
                True)
    