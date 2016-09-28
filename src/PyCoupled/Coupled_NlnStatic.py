'''@package PyCoupled.Coupled_NlnStatic
@brief      NonlinearStatic Beam + UVLM.
@author     Rob Simpson
@contact    r.simpson11@imperial.ac.uk 
@version    0.0
@date       11/02/2013
@pre        None
@warning    None
'''

import sys
import SharPySettings as Settings
import DerivedTypes
import BeamIO
import BeamLib
import BeamInit
import numpy as np
import ctypes as ct
from PyFSI.Beam2UVLM import InitSection
from PyFSI.Beam2UVLM import CoincidentGrid
from PyAero.UVLM.Utils import UVLMLib
from PyAero.UVLM.Utils.Linear import nln2linStates
from PyFSI.Beam2UVLM import CoincidentGridForce
from PyAero.UVLM.Utils import DerivedTypesAero
from PyAero.UVLM.Utils import PostProcess
from PyAero.UVLM.Solver.VLM import InitSteadyExternalVels, InitSteadyWake
from PyCoupled.Utils.DerivedTypesAeroelastic import AeroelasticOps
from PyBeam.Utils.XbeamLib import AddGravityLoads, Skew, Psi2TransMat, Tangential
from DerivedTypesAero import ControlSurf
from PyBeam.Utils.Linear import genSSbeam
from scipy.io.matlab.mio import savemat
from PyAero.UVLM.Utils.Linear import genSSuvlm
from Misc import iNode2iElem
from getpass import getuser

def Solve_Py(XBINPUT,XBOPTS,VMOPTS,VMINPUT,AELAOPTS):
    """Nonlinear static solver using Python to solve aeroelastic
    equation. Assembly of structural matrices is carried out with 
    Fortran subroutines. Aerodynamics solved using PyAero.UVLM."""
    
    assert XBOPTS.Solution.value == 112, ('NonlinearStatic requested' +
                                              ' with wrong solution code')
    
    # Initialize beam.
    XBINPUT, XBOPTS, NumNodes_tot, XBELEM, PosIni, PsiIni, XBNODE, NumDof \
                = BeamInit.Static(XBINPUT,XBOPTS)
    # Set initial conditions as undef config.
    PosDefor = PosIni.copy(order='F')
    PsiDefor = PsiIni.copy(order='F')
    if XBOPTS.PrintInfo.value==True:
        sys.stdout.write('Solve nonlinear static case in Python ... \n')
    # Initialise structural eqn tensors.
    KglobalFull = np.zeros((NumDof.value,NumDof.value),
                           ct.c_double, 'F'); ks = ct.c_int()
    FglobalFull = np.zeros((NumDof.value,NumDof.value),
                           ct.c_double, 'F'); fs = ct.c_int()
    DeltaS  = np.zeros(NumDof.value, ct.c_double, 'F')
    Qglobal = np.zeros(NumDof.value, ct.c_double, 'F')
    x       = np.zeros(NumDof.value, ct.c_double, 'F')
    dxdt    = np.zeros(NumDof.value, ct.c_double, 'F')
    # Beam Load Step tensors
    Force = np.zeros((XBINPUT.NumNodesTot,6),ct.c_double,'F')
    iForceStep     = np.zeros((NumNodes_tot.value,6), ct.c_double, 'F')
    iForceStep_Dof = np.zeros(NumDof.value, ct.c_double, 'F')
    
    # Initialze Aero.    
    Section = InitSection(VMOPTS,VMINPUT,AELAOPTS.ElasticAxis)
    # Declare memory for Aero grid and velocities.
    Zeta = np.zeros((Section.shape[0],PosDefor.shape[0],3),ct.c_double,'C')
    ZetaDot = np.zeros((Section.shape[0],PosDefor.shape[0],3),ct.c_double,'C')
    # Additional Aero solver variables.
    AeroForces = np.zeros((VMOPTS.M.value+1,VMOPTS.N.value+1,3),ct.c_double,'C')
    Gamma = np.zeros((VMOPTS.M.value,VMOPTS.N.value),ct.c_double,'C')
    # Init external velocities.
    Uext = InitSteadyExternalVels(VMOPTS,VMINPUT)
    # Create zero triads for motion of reference frame.
    VelA_A = np.zeros((3))
    OmegaA_A = np.zeros((3))
    # Create zero vectors for structural vars not used in static analysis.
    PosDotDef = np.zeros_like(PosDefor, ct.c_double, 'F')
    PsiDotDef = np.zeros_like(PsiDefor, ct.c_double, 'F')      
    
    # Define tecplot stuff.
    FileName = Settings.OutputDir + Settings.OutputFileRoot + 'AeroGrid.dat'
    Variables = ['X', 'Y', 'Z','Gamma']        
    FileObject = PostProcess.WriteAeroTecHeader(FileName, 
                                                'Default',
                                                Variables)
    
    # Start Load Loop.
    for iLoadStep in range(XBOPTS.MaxIterations.value):
        # Reset convergence parameters and loads.
        Iter = 0
        ResLog10 = 0.0
        Force[:,:] = 0.0
        AeroForces[:,:,:] = 0.0
        oldPos = PosDefor.copy(order='F')
        oldPsi = PsiDefor.copy(order='F')    
                
        # Calculate aero loads.
        if hasattr(XBINPUT, 'ForcedVel'):
            CoincidentGrid(PosDefor,
                           PsiDefor,
                           Section,
                           XBINPUT.ForcedVel[0,:3],
                           XBINPUT.ForcedVel[0,3:],
                           PosDotDef,
                           PsiDotDef,
                           XBINPUT,
                           Zeta,
                           ZetaDot,
                           ctrlSurf = VMINPUT.ctrlSurf)
        else:
            CoincidentGrid(PosDefor, PsiDefor, Section, VelA_A, 
                           OmegaA_A, PosDotDef, PsiDotDef, XBINPUT,
                           Zeta, ZetaDot,
                           ctrlSurf = VMINPUT.ctrlSurf)
        
        if hasattr(XBINPUT,'ForcedVel'):
            ZetaStar, GammaStar = InitSteadyWake(VMOPTS,VMINPUT,Zeta,XBINPUT.ForcedVel[0,:3])
        else:
            ZetaStar, GammaStar = InitSteadyWake(VMOPTS,VMINPUT,Zeta)
            
        # Define AICs here for debgugging - Rob 16/08/2016
        AIC = np.zeros((VMOPTS.M.value*VMOPTS.N.value,
                        VMOPTS.M.value*VMOPTS.N.value),
                        ct.c_double,'C')
        BIC = np.zeros((VMOPTS.M.value*VMOPTS.N.value,
                        VMOPTS.M.value*VMOPTS.N.value),
                        ct.c_double,'C')
        
        # Solve for AeroForces.
        UVLMLib.Cpp_Solver_VLM(Zeta, ZetaDot, Uext, ZetaStar, VMOPTS, 
                               AeroForces, Gamma, GammaStar, AIC, BIC)
        
        AeroForces[:,:,:] = AELAOPTS.AirDensity*AeroForces[:,:,:]
        
        # Write solution to tecplot file.
        PostProcess.WriteUVLMtoTec(FileObject,
                                   Zeta,
                                   ZetaStar,
                                   Gamma,
                                   GammaStar,
                                   iLoadStep,
                                   XBOPTS.NumLoadSteps.value,
                                   iLoadStep*1.0,
                                   Text = True)
        
        # Map AeroForces to beam.
        CoincidentGridForce(XBINPUT, PsiDefor, Section, AeroForces,
                            Force)
        
        # Add gravity loads.
        AddGravityLoads(Force,XBINPUT,XBELEM,AELAOPTS,
                        PsiDefor,VMINPUT.c)
        
        # Apply factor corresponding to force step.
        if iLoadStep < XBOPTS.NumLoadSteps.value:
            iForceStep = (Force + XBINPUT.ForceStatic)*float( (iLoadStep+1) ) / \
                                                XBOPTS.NumLoadSteps.value
        else:
            # continue at full loading until equilibrium 
            iForceStep = Force + XBINPUT.ForceStatic
        
        if XBOPTS.PrintInfo.value == True:
            sys.stdout.write('  iLoad: %-10d\n' %(iLoadStep+1))
            sys.stdout.write('   SubIter DeltaF     DeltaX     ResLog10\n')
            
        # Start Newton Iteration.
        while( (ResLog10 > np.log10(XBOPTS.MinDelta.value)) 
             & (Iter < XBOPTS.MaxIterations.value) ):
            
            Iter += 1
            if XBOPTS.PrintInfo.value == True:
                sys.stdout.write('   %-7d ' %(Iter))
                
            # Set structural eqn tensors to zero
            KglobalFull[:,:] = 0.0; ks = ct.c_int()
            FglobalFull[:,:] = 0.0; fs = ct.c_int()
            Qglobal[:] = 0.0
            
            # Assemble matrices for static problem
            BeamLib.Cbeam3_Asbly_Static(XBINPUT,
                                        NumNodes_tot,
                                        XBELEM,
                                        XBNODE,
                                        PosIni,
                                        PsiIni,
                                        PosDefor,
                                        PsiDefor,
                                        iForceStep,
                                        NumDof,
                                        ks,
                                        KglobalFull,
                                        fs,
                                        FglobalFull,
                                        Qglobal,
                                        XBOPTS)
            
            # Get state vector from current deformation.
            PosDot = np.zeros((NumNodes_tot.value,3), ct.c_double, 'F')
            PsiDot = np.zeros((XBINPUT.NumElems,Settings.MaxElNod,3),
                              ct.c_double, 'F')
            
            BeamLib.Cbeam_Solv_Disp2State(NumNodes_tot,
                                          NumDof,
                                          XBINPUT,
                                          XBNODE,
                                          PosDefor,
                                          PsiDefor,
                                          PosDot,
                                          PsiDot,
                                          x,
                                          dxdt)
            
            
            # Get forces on unconstrained nodes.
            BeamLib.f_fem_m2v(ct.byref(NumNodes_tot),
                              ct.byref(ct.c_int(6)),
                              iForceStep.ctypes.data_as(ct.POINTER(ct.c_double)),
                              ct.byref(NumDof),
                              iForceStep_Dof.ctypes.data_as(ct.POINTER(ct.c_double)),
                              XBNODE.Vdof.ctypes.data_as(ct.POINTER(ct.c_int)) )
            
            # Calculate \Delta RHS.
            Qglobal = Qglobal - np.dot(FglobalFull, iForceStep_Dof)
            
            
            # Calculate \Delta State Vector. 
            DeltaS = - np.dot(np.linalg.inv(KglobalFull), Qglobal)
                
            if XBOPTS.PrintInfo.value == True:
                sys.stdout.write('%-10.4e %-10.4e ' 
                                 % (max(abs(Qglobal)),max(abs(DeltaS))))
            
            # Update Solution.
            BeamLib.Cbeam3_Solv_Update_Static(XBINPUT,
                                              NumNodes_tot,
                                              XBELEM,
                                              XBNODE,
                                              NumDof,
                                              DeltaS,
                                              PosIni,
                                              PsiIni,
                                              PosDefor,
                                              PsiDefor)
            
            # Record residual at first iteration.
            if(Iter == 1):
                Res0_Qglobal = max(abs(Qglobal)) + 1.e-16
                Res0_DeltaX  = max(abs(DeltaS)) + 1.e-16

            # Update residual and compute log10
            Res_Qglobal = max(abs(Qglobal)) + 1.e-16
            Res_DeltaX  = max(abs(DeltaS)) + 1.e-16
            ResLog10 = max([np.log10(Res_Qglobal/Res0_Qglobal), 
                            np.log10(Res_DeltaX/Res0_DeltaX)]  )
            
            if XBOPTS.PrintInfo.value == True:
                sys.stdout.write('%8.4f\n' %(ResLog10))
            
            # Stop the solution.                
            if(ResLog10 > 10.):
                sys.stderr.write(' STOP\n')
                sys.stderr.write(' The max residual is %e\n' %(ResLog10))
                exit(1)
            elif Res_DeltaX < 1.e-14:
                break
        # END Newton iteration
        
        # After incremental loading continue until equilibrium reached
        if iLoadStep >= XBOPTS.NumLoadSteps.value:
            Pos_error = PosDefor - oldPos
            Psi_error = PsiDefor - oldPsi
            if( (np.linalg.norm(Pos_error)<=XBOPTS.MinDelta) & \
                (np.linalg.norm(Psi_error)<=XBOPTS.MinDelta) ):
                break
    # END Load step loop


    if XBOPTS.PrintInfo.value==True:
        sys.stdout.write(' ... done\n')    
     
    # Write deformed configuration to file.
    ofile = Settings.OutputDir + Settings.OutputFileRoot + '_SOL112_def.dat'
    if XBOPTS.PrintInfo.value==True:
        sys.stdout.write('Writing file %s ... ' %(ofile))
    fp = open(ofile,'w')
    fp.write('TITLE="Non-linear static solution: deformed geometry"\n')
    fp.write('VARIABLES="iElem" "iNode" "Px" "Py" "Pz" "Rx" "Ry" "Rz"\n')
    fp.close()
    if XBOPTS.PrintInfo.value==True:
        sys.stdout.write('done\n')
    WriteMode = 'a'
    
    BeamIO.OutputElems(XBINPUT.NumElems, NumNodes_tot.value, XBELEM, 
                       PosDefor, PsiDefor, ofile, WriteMode)
    
    # Print deformed configuration.
    if XBOPTS.PrintInfo.value==True:
        sys.stdout.write('--------------------------------------\n')
        sys.stdout.write('NONLINEAR STATIC SOLUTION\n')
        sys.stdout.write('%10s %10s %10s\n' %('X','Y','Z'))
        for inodi in range(NumNodes_tot.value):
            sys.stdout.write(' ')
            for inodj in range(3):
                sys.stdout.write('%12.5e' %(PosDefor[inodi,inodj]))
            sys.stdout.write('\n')
        sys.stdout.write('--------------------------------------\n')
        
    # Close Tecplot ascii FileObject.
    PostProcess.CloseAeroTecFile(FileObject)
    
    if True:
        # print and save deformed wing total force coefficients (may include gravity
        # and other applied static loads). Coefficients in wind-axes.
        cF = np.zeros((3))
        for i in range(VMOPTS.M.value+1):
            for j in range(VMOPTS.N.value+1):
                cF += AeroForces[i,j,:] 
        
        Calpha=Psi2TransMat(np.array([VMINPUT.alpha, 0.0, 0.0]))
        cF=np.dot(Calpha,cF)
        cF=cF/(0.5*AELAOPTS.AirDensity*np.linalg.norm(VMINPUT.U_infty)**2.0*VMINPUT.c*XBINPUT.BeamLength)
        print("Reference condition total force coefficients: {}".format(cF))
        fileName = Settings.OutputDir + Settings.OutputFileRoot + 'refCf'
        savemat(fileName,{'cF':cF})
        
    # Return solution
    return PosDefor, PsiDefor, Zeta, ZetaStar, Gamma, GammaStar, iForceStep, NumNodes_tot, NumDof, XBELEM, XBNODE, PosIni, PsiIni, Uext

if __name__ == '__main__':
    ## solve nlnstatic problem.
    Settings.OutputDir='/home/' + getuser() + '/Documents/MATLAB/Patil_HALE/incrementalTests/infiniteAlpha1SurgingWing/'
    
    # beam options.
    XBOPTS = DerivedTypes.Xbopts(FollowerForce = ct.c_bool(False),
                                 MaxIterations = ct.c_int(50),
                                 PrintInfo = ct.c_bool(True),
                                 NumLoadSteps = ct.c_int(10),
                                 Solution = ct.c_int(112),
                                 MinDelta = ct.c_double(1e-5))
     
    # beam inputs.
    XBINPUT = DerivedTypes.Xbinput(3,5)
    XBINPUT.BeamLength = 2000.0
    XBINPUT.BeamStiffness[0,0] = 1.0e+09
    XBINPUT.BeamStiffness[1,1] = 1.0e+09
    XBINPUT.BeamStiffness[2,2] = 1.0e+09
    XBINPUT.BeamStiffness[3,3] = 1.0e+04
    XBINPUT.BeamStiffness[4,4] = 2.0e+04
    XBINPUT.BeamStiffness[5,5] = 5.0e+06
      
    XBINPUT.BeamStiffness[:,:] = XBINPUT.BeamStiffness[:,:]*1.0e12
      
    XBINPUT.BeamMass[0,0] = 0.75
    XBINPUT.BeamMass[1,1] = 0.75
    XBINPUT.BeamMass[2,2] = 0.75
    XBINPUT.BeamMass[3,3] = 0.1
    XBINPUT.BeamMass[4,4] = 0.001
    XBINPUT.BeamMass[5,5] = 0.001
    
    # Goland
#     XBINPUT = DerivedTypes.Xbinput(3,10)
#     XBINPUT.BeamLength = 6.096
#     XBINPUT.BeamStiffness[0,0] = 1.0e+09
#     XBINPUT.BeamStiffness[1,1] = 1.0e+09
#     XBINPUT.BeamStiffness[2,2] = 1.0e+09
#     XBINPUT.BeamStiffness[3,3] = 0.99e+06
#     XBINPUT.BeamStiffness[4,4] = 9.77e+06
#     XBINPUT.BeamStiffness[5,5] = 1.0e+09
#     XBINPUT.BeamStiffness[:,:] = 1.0*XBINPUT.BeamStiffness[:,:]
#     XBINPUT.BeamMass[0,0] = 35.71
#     XBINPUT.BeamMass[1,1] = 35.71
#     XBINPUT.BeamMass[2,2] = 35.71
#     XBINPUT.BeamMass[3,3] = 8.64
#     XBINPUT.BeamMass[4,4] = 0.001
#     XBINPUT.BeamMass[5,5] = 0.001
#     # Off diagonal terms (in Theodorsen sectional coordinates)
#     ElasticAxis = -0.34
#     InertialAxis = -7.0/50.0
#     x_alpha = InertialAxis - ElasticAxis
#     # pitch-plunge coupling term (b-frame coordinates)
#     c = 1.8288
#     cgLoc = 0.5*c*np.array([0.0, x_alpha, 0.0])
#     cgSkew = Skew(cgLoc)
#     XBINPUT.BeamMass[:3,3:] = XBINPUT.BeamMass[0,0] * cgSkew
#     XBINPUT.BeamMass[3:,:3] = XBINPUT.BeamMass[:3,3:].T
#     XBINPUT.BeamMass[4,4] = 0.001 + XBINPUT.BeamMass[0,0]*pow(cgLoc[2],2.0)
#     XBINPUT.BeamMass[5,5] = 0.001 + XBINPUT.BeamMass[0,0]*pow(cgLoc[1],2.0)
    
    
    # Gravity
    XBINPUT.g = 0.0 
     
    # aero options. 
    M = 4
    N = XBINPUT.NumNodesTot - 1
    VMOPTS = DerivedTypesAero.VMopts(M = M,
                                     N = N,
                                     ImageMethod = True,
                                     Steady = True,
                                     KJMeth = True)
     
    # aero inputs.
    U_mag = 25
    ctrlSurf = None
#     U_mag=170
    chord = 1.0 #m
#     chord = c
    wakeLength = 10 # in chords
    alphaVec = (1,)#range(1); #deg
    for alpha in alphaVec:
        VMINPUT = DerivedTypesAero.VMinput(c = chord,
                                           b = XBINPUT.BeamLength,
                                           U_mag = U_mag*1.0,
                                           alpha = alpha*np.pi/180.0,
                                           theta = 0.0,
                                           ctrlSurf = ctrlSurf,
                                           WakeLength = wakeLength)
         
        # aeroelastic opts.
        # Density due to US standard atmosphere at 20km.
        EApos = 0.0 # Elastic axis position in Theodorsens' coords
        IApos = 0.0
#         EApos=ElasticAxis
#         IApos=InertialAxis
#         rho=1.02
        AELAOPTS = AeroelasticOps(EApos,IApos,0.08891)
#         AELAOPTS = AeroelasticOps(EApos,IApos,rho)
        
        Settings.OutputFileRoot='M'+str(M)+'N'+str(N)+'_V'+str(U_mag)+'_alpha'+str(alpha)
        
        # solve static problem
        PosDefor, PsiDefor, Zeta, ZetaStar, Gamma, GammaStar, iForceStep, NumNodes_tot, NumDof, XBELEM, XBNODE, PosIni, PsiIni, Uext = Solve_Py(XBINPUT,XBOPTS,VMOPTS,VMINPUT,AELAOPTS)
        
        # save reference gamma
        if True:
            fileName=Settings.OutputDir+Settings.OutputFileRoot+'_Gamma0'
            savemat(fileName,
                {'Gamma0':Gamma},
                True
            )
        
        ### linearized beam model
        A,B,C,kMat,phiSort = genSSbeam(XBINPUT,
                             NumNodes_tot,
                             NumDof,
                             XBELEM,
                             XBNODE,
                             PosIni,
                             PsiIni,
                             PosDefor,
                             PsiDefor,
                             XBOPTS,
                             chord,
                             U_mag,
                             modal=10)
        
        ### linearized aerodynamics
        # Unsteady aero parameters
        mW=wakeLength*M #10 chord-lengths
        delS=2.0/M
        
        # transform states/inputs 
        gam, gamW, gamPri, zeta, zetaW, zetaPri, nu, beam2aero = nln2linStates(Zeta, ZetaStar, Gamma, GammaStar, Uext, M, N, mW, chord, U_mag)
        
        # generate model
        Ea,Fa,Ga,Ca,Da = genSSuvlm(gam,gamW,gamPri,zeta,zetaW,zetaPri,nu,M,N,mW,delS,imageMeth=True)
        
        ### linearized interface
        xiZeta = np.zeros((3*(M+1)*(N+1),6*(N+1)))
        #transformation from beam axes to aero axes
        e=EApos/2.0+0.5 # EA position
        for ii in range(M+1):
            for jj in range(N+1):
                q=(N+1)*ii+jj;
                jElem,jjElem=iNode2iElem(jj, N+1, XBINPUT.NumNodesElem)
                Cj=Psi2TransMat(PsiDefor[jElem,jjElem,:])
                Tang=Tangential(PsiDefor[jElem,jjElem,:])
                xiZeta[3*q:3*q+3,6*jj:6*jj+3]=beam2aero
                xiZeta[3*q:3*q+3,6*jj+3:6*jj+6]=-np.dot(beam2aero,
                                                 np.dot(Cj,
                                                 np.dot(Skew([0, chord*(e-float(ii)/float(M)), 0]),
                                                        Tang)))
        # structural DoFs to aero inputs
        xiUa = np.zeros((9*(M+1)*(N+1),12*(N+1)))
        xiUa[:3*(M+1)*(N+1),:6*(N+1)] = xiZeta/chord
        xiUa[3*(M+1)*(N+1):6*(M+1)*(N+1),6*(N+1):12*(N+1)] = xiZeta/chord
        
        if True:
            fileName=Settings.OutputDir+Settings.OutputFileRoot+'_SSbeam'
            savemat(fileName,
                    {'A':A,'B':B,'C':C,'kMat':kMat,'phiSort':phiSort},
                    True)
            fileName=Settings.OutputDir+Settings.OutputFileRoot+'_SSaero'
            savemat(fileName,
                    {'Ea':Ea,'Fa':Fa,'Ga':Ga,'Ca':Ca,'Da':Da,'delS':delS,
                     'm':M, 'mW':mW, 'n':N, 'zeta':Zeta},
                    True)
            fileName=Settings.OutputDir+Settings.OutputFileRoot+'_SScoupl'
            savemat(fileName,
                    {'xiZeta':xiZeta,'xiUa':xiUa},
                    True)
            