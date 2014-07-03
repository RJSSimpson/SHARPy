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
from PyFSI.Beam2UVLM import CoincidentGridForce
from PyAero.UVLM.Utils import DerivedTypesAero
from PyAero.UVLM.Utils import PostProcess
from PyAero.UVLM.Solver.VLM import InitSteadyExternalVels
from PyAero.UVLM.Solver.VLM import InitSteadyWake
from PyCoupled.Utils.DerivedTypesAeroelastic import AeroelasticOps
from PyBeam.Utils.XbeamLib import AddGravityLoads
from DerivedTypesAero import ControlSurf

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
#         for iLoadStep in range(XBOPTS.NumLoadSteps.value):
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
            
            # Solve for AeroForces.
            UVLMLib.Cpp_Solver_VLM(Zeta, ZetaDot, Uext, ZetaStar, VMOPTS, 
                                   AeroForces, Gamma, GammaStar)
            
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
            
        # Return solution
        return PosDefor, PsiDefor, Zeta, ZetaStar, Gamma, GammaStar, iForceStep

if __name__ == '__main__':
#      solve nlnstatic problem.
#     
#      beam options.
#     XBOPTS = DerivedTypes.Xbopts(FollowerForce = ct.c_bool(False),
#                                  MaxIterations = ct.c_int(50),
#                                  PrintInfo = ct.c_bool(True),
#                                  NumLoadSteps = ct.c_int(25),
#                                  Solution = ct.c_int(112),
#                                  MinDelta = ct.c_double(1e-4))
#     
#      beam inputs.
#     XBINPUT = DerivedTypes.Xbinput(3,20)
#     XBINPUT.BeamLength = 16.0
#     XBINPUT.BeamStiffness[0,0] = 1.0e+09
#     XBINPUT.BeamStiffness[1,1] = 1.0e+09
#     XBINPUT.BeamStiffness[2,2] = 1.0e+09
#     XBINPUT.BeamStiffness[3,3] = 1.0e+04
#     XBINPUT.BeamStiffness[4,4] = 2.0e+04
#     XBINPUT.BeamStiffness[5,5] = 5.0e+06
#     
#     XBINPUT.BeamStiffness[:,:] = 1.0*XBINPUT.BeamStiffness[:,:]
#     
#     XBINPUT.BeamMass[0,0] = 0.75
#     XBINPUT.BeamMass[1,1] = 0.75
#     XBINPUT.BeamMass[2,2] = 0.75
#     XBINPUT.BeamMass[3,3] = 0.1
#     XBINPUT.BeamMass[4,4] = 0.001
#     XBINPUT.BeamMass[5,5] = 0.001
#     
#     
#      aero options. 
#     M = 4
#     N = XBINPUT.NumNodesTot - 1
#     VMOPTS = DerivedTypesAero.VMopts(M = M,
#                                      N = N,
#                                      ImageMethod = True,
#                                      Steady = True,
#                                      KJMeth = True)
#     
#      aero inputs.
#     ctrlSurf = DerivedTypesAero.ControlSurf(iMin = M - M/4,
#                                             iMax = M,
#                                             jMin = N - N/4,
#                                             jMax = N,
#                                             typeMotion = 'cos',
#                                             betaBar = 5.0*np.pi/180.0,
#                                             omega = 0.0)
#     VMINPUT = DerivedTypesAero.VMinput(c = 1.0,
#                                        b = XBINPUT.BeamLength,
#                                        U_mag = 25.0,
#                                        alpha = 0.0*np.pi/180.0,
#                                        theta = 0.0,
#                                        ctrlSurf = ctrlSurf)
#     
#      aeroelastic opts.
#      Density due to US standard atmosphere at 20km.
#     AELAOPTS = AeroelasticOps(0.0,0.0,0.08891)
#     
#     Solve_Py(XBINPUT,XBOPTS,VMOPTS,VMINPUT,AELAOPTS)
    
    
    
    ### Goland Static with Henrik ###
    
    # Beam options.
    XBOPTS = DerivedTypes.Xbopts(FollowerForce = ct.c_bool(False),
                                 MaxIterations = ct.c_int(50),
                                 PrintInfo = ct.c_bool(True),
                                 NumLoadSteps = ct.c_int(5),
                                 Solution = ct.c_int(112),
                                 MinDelta = ct.c_double(1e-4))
    # beam inputs.
    XBINPUT = DerivedTypes.Xbinput(2,28)
    XBINPUT.BeamLength = 6.096
    XBINPUT.BeamStiffness[0,0] = 1.0e+09
    XBINPUT.BeamStiffness[1,1] = 1.0e+09
    XBINPUT.BeamStiffness[2,2] = 1.0e+09
    XBINPUT.BeamStiffness[3,3] = 0.99e+06
    XBINPUT.BeamStiffness[4,4] = 9.77e+06
    XBINPUT.BeamStiffness[5,5] = 1.0e+09
    XBINPUT.BeamStiffness[:,:] = 1.0*XBINPUT.BeamStiffness[:,:]
    XBINPUT.BeamMass[0,0] = 35.71
    XBINPUT.BeamMass[1,1] = 35.71
    XBINPUT.BeamMass[2,2] = 35.71
    XBINPUT.BeamMass[3,3] = 8.64
    XBINPUT.BeamMass[4,4] = 0.001
    XBINPUT.BeamMass[5,5] = 0.001
    # Off diagonal terms (in Theodorsen sectional coordinates)
    ElasticAxis = -0.34
    InertialAxis = -7.0/50.0
    x_alpha = InertialAxis - ElasticAxis
    # pitch-plunge coupling term (b-frame coordinates)
    c = 1.8288
    mOff = x_alpha*(c/2)*XBINPUT.BeamMass[0,0]
    XBINPUT.BeamMass[2,3] = -mOff
    XBINPUT.BeamMass[0,5] = mOff
    XBINPUT.BeamMass[3:,:3] = XBINPUT.BeamMass[:3,3:].T
    # Get suggested panelling.
    Umag = 140.0
    M = 16

    # aero params.
    WakeLength = 30.0*c
    Mstar = 1    # aero options.
    N = XBINPUT.NumNodesTot - 1
    VMOPTS = DerivedTypesAero.VMopts(M = M,
                                     N = N,
                                     ImageMethod = True,
                                     Mstar = Mstar,
                                     Steady = True,
                                     KJMeth = True,
                                     NewAIC = True)
    # Aero inputs.
    iMin = M - M/4
    iMax = M
    jMin = N - N/4
    jMax = N
    typeMotion = 'sin'
    betaBar = 0.0*np.pi/180.0
    omega = 30.0
    ctrlSurf = ControlSurf(iMin,
                           iMax,
                           jMin,
                           jMax,
                           typeMotion,
                           betaBar,
                           omega)
    
    VMINPUT = DerivedTypesAero.VMinput(c = c,
                                       b = XBINPUT.BeamLength,
                                       U_mag = Umag,
                                       alpha = 5.0*np.pi/180.0,
                                       theta = 0.0,
                                       WakeLength = WakeLength,
                                       ctrlSurf = ctrlSurf)
    
    # Aerolastic simulation results. # 1.02
    AELAOPTS = AeroelasticOps(ElasticAxis = ElasticAxis,
                              InertialAxis = InertialAxis,
                              AirDensity = 1.02,
                              Tight = False,
                              ImpStart = False)
    # Solve nonlinear static simulation.
    Solve_Py(XBINPUT,XBOPTS,VMOPTS,VMINPUT,AELAOPTS)
    
    
    