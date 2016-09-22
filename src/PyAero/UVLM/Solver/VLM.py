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
from PyAero.UVLM.Utils.Linear import genSSuvlm, nln2linStates, runLinearAero
import getpass
from scipy.io.matlab.mio import savemat, loadmat
import matplotlib.pyplot as plt

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
    if VelA_G is None:
        DeltaX = (VMINPUT.WakeLength * 
                 (VMINPUT.U_infty / np.linalg.norm(VMINPUT.U_infty)))
    elif VelA_G is not None:
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
    Coeffs = PostProcess.GetCoeffs(VMOPTS, Forces, VMINPUT, VMUNST.VelA_G)
    return Coeffs, Zeta, ZetaStar, Gamma, GammaStar, Forces, Uext


if __name__ == '__main__':
    
    Settings.OutputDir = '/home/' + getpass.getuser() + '/Documents/MATLAB/newUVLM/nonZeroAerofoil/'
    writeToMat = False
    runLinear = True
    
    # Inputs.
    m=20
    n=1
    Umag = 1.0
    alpha = 1.0*np.pi/180.0
    chord = 1.0
    span=2.e3
    WakeLength = 10.0
    imageMeth = True
    VMINPUT = DerivedTypesAero.VMinput(chord ,
                                   b = span,
                                   U_mag = Umag,
                                   alpha = alpha,
                                   theta = 0.0,
                                   WakeLength = WakeLength,
                                   ctrlSurf = None)
    
    
    VMOPTS = VMopts(m,
                    n,
                    imageMeth,
                    Mstar = 1,
                    Steady = True,
                    KJMeth = True)
    # run solver
    Coeffs, Zeta, ZetaStar, Gamma, GammaStar, foo, Uext = Run_Cpp_Solver_VLM(VMOPTS,VMINPUT)[0:7]
    del foo
    
    # unsteady params
    mW=10*m
    delS=2/m
    
    # transform states/inputs 
    gam, gamW, gamPri, zeta, zetaW, zetaPri, nu, beam2aero = nln2linStates(Zeta, ZetaStar, Gamma, GammaStar, Uext, m, n, mW, chord)
    
    # generate linear model
    E,F,G,C,D = genSSuvlm(gam,gamW,gamPri,zeta,zetaW,zetaPri,nu,m,n,mW,delS,imageMeth)
    
    # matrices for aerofoil DoFs
    e=0.25+zeta[0]-0.25/m
    f=0.75+zeta[0]-0.25/m
    
    # convert inputs from general kinematics to aerofil DoFs with heave
    T = np.zeros((9*(m+1)*(n+1),6))
    rot=np.zeros((3,3))
    rot[0,0]=np.cos(alpha)
    rot[0,2]=-np.sin(alpha)
    rot[1,1]=1.0
    rot[2,0]=np.sin(alpha)
    rot[2,2]=np.cos(alpha)
    for i in range(m+1):
        for j in range(n+1):
            q=i*(n+1)+j
            # alpha, alphaPrime
            T[3*(m+1)*(n+1)+3*q+2,0] = -(zeta[3*q]+0.25/m-e)
            T[3*q+2,1] = -(zeta[3*q]+0.25/m-e)
            # plunge velocity
            T[3*q:3*q+3,2]=np.dot(rot,np.array([0, 0, -1]))
            # in-plane velocity (+ nose direction)
            T[3*q:3*q+3,5] = np.dot(rot,np.array([-1, 0, 0]))
            # beta, betaPrime
            if zeta[3*q]+0.25/m > f:
                T[3*(m+1)*(n+1)+3*q+2,3] = -(zeta[3*q]+0.25/m-f)
                T[3*q+2,4] = -(zeta[3*q]+0.25/m-f)
                
                
    G_s = np.dot(G,T)
    D_s = np.dot(D,T)
    
    # get coefficients as output
    T_coeff = np.zeros((3,3*(m+1)*(n+1)))
    T_coeff[0,0::3] = 1.0 #drag
    T_coeff[1,2::3] = 1.0 #lift
    for i in range(m+1):
        for j in range(n+1):
            q=i*(n+1)+j
            # moment = r*L, +ve nose-up, at quarter chord
            T_coeff[2,3*q+2] = -(zeta[3*q]-e)
    
    C_coeff = np.dot(T_coeff,C)
    D_coeff = np.dot(T_coeff,D)
    D_s_coeff = np.dot(T_coeff,D_s)
    
    # spanwise lift distribution as output
    T_span=np.zeros((n+1,3*(m+1)*(n+1)))
    for jj in range(n+1):
        T_span[jj,3*jj+2::3*(n+1)] = 1.0

    if writeToMat == True:
        fileName = Settings.OutputDir + 'rectWingAR' + str(span/chord) + '_m' + str(m) + 'mW' + str(mW) + 'n' + str(n) + 'delS' + str(delS) + '_alpha' + str(alpha)
        if e != 0.25:
            fileName += 'e'+str(e)
        if f != 0.75:
            fileName += 'f'+str(f)
        if imageMeth != False:
            fileName += 'half'
        savemat(fileName,
                {'E':E, 'F':F, 'G':G, 'C':C, 'D':D, 'm':m, 'mW':mW, 'delS':delS,
                 'G_s':G_s, 'D_s':D_s,
                 'C_coeff':C_coeff, 'D_coeff':D_coeff, 'D_s_coeff':D_s_coeff,
                 'T_coeff':T_coeff, 'T_span':T_span,
                 'AR':span/chord, 'm':m, 'mW':mW, 'n':n, 'zeta':zeta},
                True)
        
    if runLinear == True:
        nT = 4001 # number of time steps
        u = np.zeros((nT,G_s.shape[1])) # inputs
        k = 1.0 # reduced frequency
        hBar = 0.01 # 1% of chord
        tau = np.linspace(0.0, nT*delS, nT)
        surge=hBar*np.sin(k*tau)
        u[:,5]=k*hBar*np.cos(k*tau) # in-plane vibrations
        tOut, yOut = runLinearAero(E, F, G_s, C_coeff, D_s_coeff, delS, nT, u)[0:2]
        # load UVLM data
        dataDir = '/home/rjs10/Documents/MATLAB/newUVLM/nonZeroAerofoil/UVLM/'
        data = loadmat(dataDir + 'UVLMrectAR2000.0_m20mW200n1delS0.1_alpha0.01745_half.mat')
        # Plot
        plt.plot(tOut,surge,'b--',
                 tOut, 2*yOut[:,1]/(2*np.pi*np.pi/180)/2000.0,'r-',
                 data['Coeffs'][:,0]*2,(data['Coeffs'][:,3]-Coeffs[2])/(2*np.pi*np.pi/180),'k--')
        plt.xlabel(r'$\tau$')
        plt.ylabel(r'$C_l / 2\pi\alpha_0$')
        plt.title('Lift perturbations due to surging motion')
        plt.grid(True)
        plt.ylim((-0.05,0.05))
        #plt.savefig("test.png")
        plt.show()
    # end if
    