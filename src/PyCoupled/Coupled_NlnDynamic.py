'''@package PyCoupled.Coupled_NlnDynamic
@brief      NonlinearDynamic Beam + UVLM.
@author     Rob Simpson
@contact    r.simpson11@imperial.ac.uk 
@version    0.0
@date       11/02/2013
@pre        None
@warning    None
'''

import sys
import Main.SharPySettings as Settings
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
import PyCoupled.Coupled_NlnStatic as Static
import PyBeam.Utils.XbeamLib as xbl
from PyCoupled.Coupled_NlnStatic import AddGravityLoads
from DerivedTypesAero import ControlSurf
from collections import OrderedDict
import re
from math import pow
from PyBeam.Utils.XbeamLib import Skew
from PyAero.UVLM.Utils.DerivedTypesAero import Gust
import PyMPC.MPC as MPC
from PyCoupled.Utils.CoupledIO import OpenOutFile, WriteToOutFile

def Solve_Py(XBINPUT,XBOPTS,VMOPTS,VMINPUT,AELAOPTS,**kwords):
    """@brief Nonlinear dynamic solver using Python to solve aeroelastic
    equation.
    @details Assembly of structural matrices is carried out with 
    Fortran subroutines. Aerodynamics solved using PyAero\.UVLM.
    @warning test outstanding: test for maintaining static deflections in
    same conditions.
    TODO: Maintain static deflections in same conditions.
    @param XBINPUT Beam inputs (for initialization in Python).
    @param XBOPTS Beam solver options (for Fortran).
    @param VMOPTS UVLM solver options (for C/C++).
    @param VMINPUT UVLM solver inputs (for initialization in Python).
    @param VMUNST Unsteady input information for aero solver.
    @param AELAOPTS Options relevant to coupled aeroelastic simulations.
    
    kwords:
    @param writeDict OrderedDict of of outputs to write.
    @param mpcCont Instance of PyMPC.MPC class.
    """
        
    # Check correct solution code.
    assert XBOPTS.Solution.value == 312, ('NonlinearDynamic requested' +
                                          ' with wrong solution code')
    # Initialise static beam data.
    XBINPUT, XBOPTS, NumNodes_tot, XBELEM, PosIni, PsiIni, XBNODE, NumDof \
                = BeamInit.Static(XBINPUT,XBOPTS)
                
    # Calculate initial displacements.
    if AELAOPTS.ImpStart == False:
        XBOPTS.Solution.value = 112 # Modify options.
        VMOPTS.Steady = ct.c_bool(True)
        Rollup = VMOPTS.Rollup.value
        VMOPTS.Rollup.value = False
        # Solve Static Aeroelastic.
        PosDefor, PsiDefor, Zeta, ZetaStar, Gamma, GammaStar, Force = \
                    Static.Solve_Py(XBINPUT, XBOPTS, VMOPTS, VMINPUT, AELAOPTS)
        XBOPTS.Solution.value = 312 # Reset options.
        VMOPTS.Steady = ct.c_bool(False)
        VMOPTS.Rollup.value = Rollup
    elif AELAOPTS.ImpStart == True:
        PosDefor = PosIni.copy(order='F')
        PsiDefor = PsiIni.copy(order='F')
        Force = np.zeros((XBINPUT.NumNodesTot,6),ct.c_double,'F')
        
    # Write deformed configuration to file. TODO: tidy this away inside function.
    ofile = Settings.OutputDir + Settings.OutputFileRoot + '_SOL312_def.dat'
    if XBOPTS.PrintInfo==True:
        sys.stdout.write('Writing file %s ... ' %(ofile))
    fp = open(ofile,'w')
    fp.write('TITLE="Non-linear static solution: deformed geometry"\n')
    fp.write('VARIABLES="iElem" "iNode" "Px" "Py" "Pz" "Rx" "Ry" "Rz"\n')
    fp.close()
    if XBOPTS.PrintInfo==True:
        sys.stdout.write('done\n')
    WriteMode = 'a'
    # Write
    BeamIO.OutputElems(XBINPUT.NumElems, NumNodes_tot.value, XBELEM,
                       PosDefor, PsiDefor, ofile, WriteMode)
    
    # Initialise structural variables for dynamic analysis.
    Time, NumSteps, ForceTime, ForcedVel, ForcedVelDot,\
    PosDotDef, PsiDotDef,\
    OutGrids, PosPsiTime, VelocTime, DynOut\
        = BeamInit.Dynamic(XBINPUT,XBOPTS)
    # Delete unused variables.
    del ForceTime, OutGrids, VelocTime
        
    # Write _force file
#    ofile = Settings.OutputDir + Settings.OutputFileRoot + '_SOL312_force.dat'
#    fp = open(ofile,'w')
#    BeamIO.Write_force_File(fp, Time, ForceTime, ForcedVel, ForcedVelDot)
#    fp.close() 
    # Write _vel file   
    #TODO: write _vel file
    # Write .mrb file.
    #TODO: write .mrb file
    
    if XBOPTS.PrintInfo.value==True:
        sys.stdout.write('Solve nonlinear dynamic case in Python ... \n')
    
    # Initialise structural system tensors.
    MglobalFull = np.zeros((NumDof.value,NumDof.value), ct.c_double, 'F')
    CglobalFull = np.zeros((NumDof.value,NumDof.value), ct.c_double, 'F') 
    KglobalFull = np.zeros((NumDof.value,NumDof.value), ct.c_double, 'F') 
    FglobalFull = np.zeros((NumDof.value,NumDof.value), ct.c_double, 'F')
    Asys = np.zeros((NumDof.value,NumDof.value), ct.c_double, 'F')
    
    ms = ct.c_int()
    cs = ct.c_int()
    ks = ct.c_int()
    fs = ct.c_int()
    
    Mvel = np.zeros((NumDof.value,6), ct.c_double, 'F')
    Cvel = np.zeros((NumDof.value,6), ct.c_double, 'F')
    
#     X0    = np.zeros(NumDof.value, ct.c_double, 'F')
    X     = np.zeros(NumDof.value, ct.c_double, 'F')
    DX    = np.zeros(NumDof.value, ct.c_double, 'F')
    dXdt  = np.zeros(NumDof.value, ct.c_double, 'F')
    dXddt = np.zeros(NumDof.value, ct.c_double, 'F')
    Force_Dof = np.zeros(NumDof.value, ct.c_double, 'F')
    
    Qglobal = np.zeros(NumDof.value, ct.c_double, 'F')
    
    # Initialise rotation operators.
    Unit = np.zeros((3,3), ct.c_double, 'F')
    for i in range(3):
        Unit[i,i] = 1.0
    
    Unit4 = np.zeros((4,4), ct.c_double, 'F')
    for i in range(4):
        Unit4[i,i] = 1.0
        
    Cao = Unit.copy('F')
    Temp = Unit4.copy('F')
    
    Quat = np.zeros(4, ct.c_double, 'F')
    Quat[0] = 1.0
    
    # Extract initial displacements and velocities.
    BeamLib.Cbeam_Solv_Disp2State(NumNodes_tot, NumDof, XBINPUT, XBNODE,
                                  PosDefor, PsiDefor, PosDotDef, PsiDotDef,
                                  X, dXdt)
    
    # Approximate initial accelerations.
    PosDotDotDef = np.zeros((NumNodes_tot.value,3),ct.c_double,'F')
    PsiDotDotDef = np.zeros((XBINPUT.NumElems,Settings.MaxElNod,3),
                             ct.c_double, 'F')
    
    # Assemble matrices for dynamic analysis.
    BeamLib.Cbeam3_Asbly_Dynamic(XBINPUT, NumNodes_tot, XBELEM, XBNODE,
                                 PosIni, PsiIni, PosDefor, PsiDefor,
                                 PosDotDef, PsiDotDef, PosDotDotDef, PsiDotDotDef,
                                 Force, ForcedVel[0,:], ForcedVelDot[0,:],
                                 NumDof, Settings.DimMat,
                                 ms, MglobalFull, Mvel,
                                 cs, CglobalFull, Cvel,
                                 ks, KglobalFull, fs, FglobalFull,
                                 Qglobal, XBOPTS, Cao)
    
    # Get force vector for unconstrained nodes (Force_Dof).
    BeamLib.f_fem_m2v(ct.byref(NumNodes_tot),
                      ct.byref(ct.c_int(6)),
                      Force.ctypes.data_as(ct.POINTER(ct.c_double)),
                      ct.byref(NumDof),
                      Force_Dof.ctypes.data_as(ct.POINTER(ct.c_double)),
                      XBNODE.Vdof.ctypes.data_as(ct.POINTER(ct.c_int)) )
    
    # Get RHS at initial condition.
    Qglobal = Qglobal - np.dot(FglobalFull, Force_Dof)
    
    # Initial Accel.
    dXddt[:] = np.dot(np.linalg.inv(MglobalFull), -Qglobal)
    
    # Record position of all grid points in global FoR at initial time step.
    DynOut[0:NumNodes_tot.value,:] = PosDefor
    
    # Record state of the selected node in initial deformed configuration.
    PosPsiTime[0,:3] = PosDefor[-1,:]
    PosPsiTime[0,3:] = PsiDefor[-1,XBELEM.NumNodes[-1]-1,:]
    
    # Get gamma and beta for Newmark scheme.
    gamma = 0.5 + XBOPTS.NewmarkDamp.value
    beta = 0.25*pow((gamma + 0.5),2)
    
    # Initialize Aero       
    Section = InitSection(VMOPTS,VMINPUT,AELAOPTS.ElasticAxis)
    
    # Declare memory for Aero variables.
    ZetaDot = np.zeros((Section.shape[0],PosDefor.shape[0],3),ct.c_double,'C')
    K = VMOPTS.M.value*VMOPTS.N.value
    AIC = np.zeros((K,K),ct.c_double,'C')
    BIC = np.zeros((K,K),ct.c_double,'C')
    AeroForces = np.zeros((VMOPTS.M.value+1,VMOPTS.N.value+1,3),ct.c_double,'C')
    
    # Initialise A-frame location and orientation to be zero
    OriginA_G = np.zeros(3,ct.c_double,'C')
    PsiA_G = np.zeros(3,ct.c_double,'C')
    
    # Init external velocities.  
    Ufree = InitSteadyExternalVels(VMOPTS,VMINPUT)
    
    # Init uninit vars if an impulsive start is specified.
    if AELAOPTS.ImpStart == True:
        Zeta = np.zeros((Section.shape[0],PosDefor.shape[0],3),ct.c_double,'C')             
        Gamma = np.zeros((VMOPTS.M.value,VMOPTS.N.value),ct.c_double,'C')
        # Generate surface, wake and gamma matrices.
        CoincidentGrid(PosDefor, PsiDefor, Section, ForcedVel[0,:3], 
                       ForcedVel[0,3:], PosDotDef, PsiDotDef, XBINPUT,
                       Zeta, ZetaDot, OriginA_G, PsiA_G,
                       VMINPUT.ctrlSurf)
        # init wake grid and gamma matrix.
        ZetaStar, GammaStar = InitSteadyWake(VMOPTS,VMINPUT,Zeta,ForcedVel[0,:3])
        
    # Init GammaDot
    GammaDot = np.zeros_like(Gamma, ct.c_double, 'C')
          
    # Define tecplot stuff
    if Settings.PlotTec==True:
        FileName = Settings.OutputDir + Settings.OutputFileRoot + 'AeroGrid.dat'
        Variables = ['X', 'Y', 'Z','Gamma']        
        FileObject = PostProcess.WriteAeroTecHeader(FileName, 
                                                    'Default',
                                                    Variables)
        # Plot results of static analysis
        PostProcess.WriteUVLMtoTec(FileObject,
                                   Zeta,
                                   ZetaStar,
                                   Gamma,
                                   GammaStar,
                                   TimeStep = 0,
                                   NumTimeSteps = XBOPTS.NumLoadSteps.value,
                                   Time = 0.0,
                                   Text = False)
    
    # Open output file for writing
    if 'writeDict' in kwords and Settings.WriteOut == True:
        fp = OpenOutFile(kwords['writeDict'], XBOPTS, Settings)
        
    # Write initial outputs to file.
    if 'writeDict' in kwords and Settings.WriteOut == True:
        WriteToOutFile(kwords['writeDict'],
                       fp,
                       Time[0],
                       PosDefor,
                       PsiDefor,
                       PosIni,
                       PsiIni,
                       XBELEM,
                       ctrlSurf,
                       kwords['mpcCont'])
    # END if write

    # Time loop.
    for iStep in range(NumSteps.value):
        
        if XBOPTS.PrintInfo.value==True:
            sys.stdout.write('Time: %-10.4e\n' %(Time[iStep+1]))
            sys.stdout.write('   SubIter DeltaF     DeltaX     ResLog10\n')
        
        dt = Time[iStep+1] - Time[iStep]
        
        # Set dt for aero force calcs.
        VMOPTS.DelTime = ct.c_double(dt)
        
        # Save Gamma at iStep.
        GammaSav = Gamma.copy(order = 'C')
        # Force at current time-step
        if iStep > 0 and AELAOPTS.Tight == False:
            
            # zero aero forces.
            AeroForces[:,:,:] = 0.0
            
            # Update CRV.
            PsiA_G = xbl.quat2psi(Quat) # CRV at iStep
            
            # Update origin.
            OriginA_G[:] = OriginA_G[:] + ForcedVel[iStep-1,:3]*dt
            
            # Update control surface deflection.
            if VMINPUT.ctrlSurf != None:
                if 'mpcCont' in kwords and kwords['mpcCont'] != None:
                    uOpt = kwords['mpcCont'].getUopt(
                            getState(Gamma,GammaStar,GammaDot,X,dXdt) )
                    VMINPUT.ctrlSurf.update(Time[iStep],uOpt[0,0])
                else:
                    VMINPUT.ctrlSurf.update(Time[iStep])
            
            # Generate surface grid.
            CoincidentGrid(PosDefor, PsiDefor, Section, ForcedVel[iStep,:3], 
                           ForcedVel[iStep,3:], PosDotDef, PsiDotDef, XBINPUT,
                           Zeta, ZetaDot, OriginA_G, PsiA_G,
                           VMINPUT.ctrlSurf)
            
            # Update wake geom       
            #'roll' data.
            ZetaStar = np.roll(ZetaStar,1,axis = 0)
            GammaStar = np.roll(GammaStar,1,axis = 0)
            #overwrite grid points with TE.
            ZetaStar[0,:] = Zeta[VMOPTS.M.value,:]
            # overwrite Gamma with TE value from previous timestep.
            GammaStar[0,:] = Gamma[VMOPTS.M.value-1,:]
            
            # Apply gust velocity.
            if VMINPUT.gust != None:
                Utot = Ufree + VMINPUT.gust.Vels(Zeta)
            else:
                Utot = Ufree
            
            # Solve for AeroForces.
            UVLMLib.Cpp_Solver_VLM(Zeta, ZetaDot, Utot, ZetaStar, VMOPTS, 
                           AeroForces, Gamma, GammaStar, AIC, BIC)
            
            # Get GammaDot.
            GammaDot[:] = Gamma[:] - GammaSav[:]
            
            # Apply density scaling.
            AeroForces[:,:,:] = AELAOPTS.AirDensity*AeroForces[:,:,:]
            
            if Settings.PlotTec==True:
                PostProcess.WriteUVLMtoTec(FileObject,
                                           Zeta - OriginA_G[:],
                                           ZetaStar - OriginA_G[:],
                                           Gamma,
                                           GammaStar,
                                           TimeStep = iStep,
                                           NumTimeSteps = XBOPTS.NumLoadSteps.value,
                                           Time = Time[iStep],
                                           Text = False)

            # map AeroForces to beam.
            CoincidentGridForce(XBINPUT, PsiDefor, Section, AeroForces,
                                Force)
            
            # Add gravity loads.
            AddGravityLoads(Force, XBINPUT, XBELEM, AELAOPTS,
                            PsiDefor, VMINPUT.c)
        #END if iStep > 0
        
        # Quaternion update for orientation.
        Temp = np.linalg.inv(Unit4 + 0.25*xbl.QuadSkew(ForcedVel[iStep+1,3:])*dt)
        Quat = np.dot(Temp, np.dot(Unit4 - 0.25*xbl.QuadSkew(ForcedVel[iStep,3:])*dt, Quat))
        Quat = Quat/np.linalg.norm(Quat)
        Cao  = xbl.Rot(Quat) # transformation matrix at iStep+1
        
        if AELAOPTS.Tight == True:
            # CRV at iStep+1
            PsiA_G = xbl.quat2psi(Quat)
            # Origin at iStep+1
            OriginA_G[:] = OriginA_G[:] + ForcedVel[iStep,:3]*dt
            
            GammaSav = Gamma.copy(order = 'C')
            
        # Predictor step.
        X        = X + dt*dXdt + (0.5-beta)*dXddt*pow(dt,2.0)
        dXdt     = dXdt + (1.0-gamma)*dXddt*dt
        dXddt[:] = 0.0
        
        # Reset convergence parameters.
        Iter = 0
        ResLog10 = 0.0
        
        # Newton-Raphson loop.        
        while ( (ResLog10 > np.log10(XBOPTS.MinDelta.value)) 
                & (Iter < XBOPTS.MaxIterations.value) ):
            
            # set tensors to zero.
            Qglobal[:] = 0.0 
            Mvel[:,:] = 0.0
            Cvel[:,:] = 0.0
            MglobalFull[:,:] = 0.0
            CglobalFull[:,:] = 0.0
            KglobalFull[:,:] = 0.0
            FglobalFull[:,:] = 0.0
            
            # Update counter.
            Iter += 1
            
            if XBOPTS.PrintInfo.value==True:
                sys.stdout.write('   %-7d ' %(Iter))
            
            # nodal diplacements and velocities from state vector.
            BeamLib.Cbeam3_Solv_State2Disp(XBINPUT,
                                           NumNodes_tot,
                                           XBELEM,
                                           XBNODE,
                                           PosIni,
                                           PsiIni,
                                           NumDof,
                                           X,
                                           dXdt,
                                           PosDefor,
                                           PsiDefor,
                                           PosDotDef,
                                           PsiDotDef)
            
            # if tightly coupled is on then get new aeroforces.
            if AELAOPTS.Tight == True:
                # zero aero forces.
                AeroForces[:,:,:] = 0.0
                
                # Set gamma at t-1 to saved solution.
                Gamma[:,:] = GammaSav[:,:]
                # get new grid.
                # The rigid-body DoFs (OriginA_G,PsiA_G,ForcedVel) at time step
                # i+1 are used to converge the aeroelastic equations.
                CoincidentGrid(PosDefor, PsiDefor, Section, ForcedVel[iStep+1,:3], 
                               ForcedVel[iStep+1,3:], PosDotDef, PsiDotDef, XBINPUT,
                               Zeta, ZetaDot, OriginA_G, PsiA_G,
                               VMINPUT.ctrlSurf)
                
                # close wake.
                ZetaStar[0,:] = Zeta[VMOPTS.M.value,:]
                
                # save pereference and turn off rollup.
                Rollup = VMOPTS.Rollup.value
                VMOPTS.Rollup.value = False
                
                # Solve for AeroForces.
                UVLMLib.Cpp_Solver_VLM(Zeta, ZetaDot, Ufree, ZetaStar, VMOPTS, 
                                       AeroForces, Gamma, GammaStar, AIC, BIC)
                
                # turn rollup back to original preference
                VMOPTS.Rollup.value = Rollup
                
                # apply density scaling.
                AeroForces[:,:,:] = AELAOPTS.AirDensity*AeroForces[:,:,:]
                
                # beam forces.
                CoincidentGridForce(XBINPUT, PsiDefor, Section, AeroForces,
                                    Force)
                
                # Add gravity loads.
                AddGravityLoads(Force,XBINPUT,XBELEM,AELAOPTS,
                                PsiDefor,VMINPUT.c)
            
            #END if Tight
            
            ForcedVelLoc = ForcedVel[iStep+1,:].copy('F')
            ForcedVelDotLoc = ForcedVelDot[iStep+1,:].copy('F')
            
            # Update matrices.
            BeamLib.Cbeam3_Asbly_Dynamic(XBINPUT, NumNodes_tot, XBELEM, XBNODE,
                         PosIni, PsiIni, PosDefor, PsiDefor,
                         PosDotDef, PsiDotDef, PosDotDotDef, PsiDotDotDef,
                         Force, ForcedVelLoc, ForcedVelDotLoc,
                         NumDof, Settings.DimMat,
                         ms, MglobalFull, Mvel,
                         cs, CglobalFull, Cvel,
                         ks, KglobalFull, fs, FglobalFull,
                         Qglobal, XBOPTS, Cao)
            
            
            # Get force vector for unconstrained nodes (Force_Dof).
            BeamLib.f_fem_m2v(ct.byref(NumNodes_tot),
                              ct.byref(ct.c_int(6)),
                              Force.ctypes.data_as(ct.POINTER(ct.c_double)),
                              ct.byref(NumDof),
                              Force_Dof.ctypes.data_as(ct.POINTER(ct.c_double)),
                              XBNODE.Vdof.ctypes.data_as(ct.POINTER(ct.c_int)) )
            
            # Solve for update vector.
            # Residual.
            Qglobal = Qglobal +  np.dot(MglobalFull, dXddt) \
                              + np.dot(Mvel,ForcedVelDotLoc) \
                              - np.dot(FglobalFull, Force_Dof)
                              
            if XBOPTS.PrintInfo.value==True:                 
                sys.stdout.write('%-10.4e ' %(max(abs(Qglobal))))

            
            # Calculate system matrix for update calculation.
            Asys = KglobalFull \
                    + CglobalFull*gamma/(beta*dt) \
                    + MglobalFull/(beta*pow(dt,2.0))
            
            # Solve for update.
            DX[:] = np.dot(np.linalg.inv(Asys), -Qglobal)
            
            # Corrector step.
            X = X + DX
            dXdt = dXdt + DX*gamma/(beta*dt)
            dXddt = dXddt + DX/(beta*pow(dt,2.0))
            
            # Residual at first iteration.
            if(Iter == 1):
                Res0_Qglobal = max(abs(Qglobal)) + 1.e-16
                Res0_DeltaX  = max(abs(DX)) + 1.e-16
            
            # Update residual and compute log10.
            Res_Qglobal = max(abs(Qglobal)) + 1.e-16
            Res_DeltaX  = max(abs(DX)) + 1.e-16
            ResLog10 = max([np.log10(Res_Qglobal/Res0_Qglobal),
                            np.log10(Res_DeltaX/Res0_DeltaX)])
            
            if XBOPTS.PrintInfo.value==True:
                sys.stdout.write('%-10.4e %8.4f\n' %(max(abs(DX)),ResLog10))
            
            if ResLog10 > 2.0:
                print("Residual growing! Exit Newton-Raphson...")
                break
                
        # END Netwon-Raphson.
        
        if ResLog10 > 2.0:
                print("Residual growing! Exit time-loop...")
                debug = 'here'
                del debug
                break
        
        # Update to converged nodal displacements and velocities.
        BeamLib.Cbeam3_Solv_State2Disp(XBINPUT, NumNodes_tot, XBELEM, XBNODE,
                           PosIni, PsiIni, NumDof, X, dXdt,
                           PosDefor, PsiDefor, PosDotDef, PsiDotDef)
        
        PosPsiTime[iStep+1,:3] = PosDefor[-1,:]
        PosPsiTime[iStep+1,3:] = PsiDefor[-1,XBELEM.NumNodes[-1]-1,:]
        
        # Position of all grid points in global FoR.
        i1 = (iStep+1)*NumNodes_tot.value
        i2 = (iStep+2)*NumNodes_tot.value
        DynOut[i1:i2,:] = PosDefor
        
        # Write selected outputs
        if 'writeDict' in kwords and Settings.WriteOut == True:
            WriteToOutFile(writeDict,
                           fp,
                           Time[iStep+1],
                           PosDefor,
                           PsiDefor,
                           PosIni,
                           PsiIni,
                           XBELEM,
                           ctrlSurf,
                           kwords['mpcCont'])
        # END if write.
        
        ZetaStar[:,:] = ZetaStar[:,:] + VMINPUT.U_infty*dt
        if VMINPUT.gust != None:
            ZetaStar[:,:,:] = ZetaStar[:,:,:] + VMINPUT.gust.Vels(ZetaStar)*dt
    
    # END Time loop
    
    # Write _dyn file.
    ofile = Settings.OutputDir + Settings.OutputFileRoot + '_SOL312_dyn.dat'
    fp = open(ofile,'w')
    BeamIO.Write_dyn_File(fp, Time, PosPsiTime)
    fp.close()
    
#    "Write _shape file"
#    ofile = Settings.OutputDir + Settings.OutputFileRoot + '_SOL312_shape.dat'
#    fp = open(ofile,'w')
#    BeamIO.Write_shape_File(fp, len(Time), NumNodes_tot.value, Time, DynOut)
#    fp.close()
    
    # Close output file if it exists.
    if 'writeDict' in kwords and Settings.WriteOut == True:
        fp.close()
    
    # Close Tecplot ascii FileObject.
    if Settings.PlotTec==True:
        PostProcess.CloseAeroTecFile(FileObject)
    
    if XBOPTS.PrintInfo.value==True:
        sys.stdout.write(' ... done\n')
        
    # For interactive analysis at end of simulation set breakpoint.
    pass

def panellingFromFreq(freq,c=1.0,Umag=1.0):
    """@brief Calculate adequate spatial/temporal resolution on UVLM grid
    based on a frequency of interest.
    @param freq Frequency of interest.
    @param c chord of wing.
    @param Umag mean reltive free-stream velocity magnitude.
    @returns M Number of chordwise panels for wing.
    @returns DelTime Suggested timestep [seconds.]
    """
    k = freq*c/(2*Umag) #get expected reduced freq
    M = int(50.0*k/np.pi) #get adequate panelling
    DelTime = c/(Umag*M) #get resulting DelTime
    return M, DelTime

def getState(Gamma,GammaStar,GammaDot,X,dXdt):
    """@brief Get full (coupled aeroelastic) state vector.
    @param Bound circulation strengths.
    @param Wake circulation strengths.
    @param Bound rate of circulation change.
    @param Strutural DoF.
    @returns x state vector.
    """
    
    # Get size of each group
    nG = Gamma.size
    nGs = GammaStar.size
    nGd = GammaDot.size
    nX = X.size
    nXd = dXdt.size
    xSize = nG + nGs + nGd + nX + nXd
    # Print if you want to cross reference state sizes with MATLAB.
    check = False
    if check == True:
        print('Gamma     = ',nG)
        print('GammaStar = ',nGs)
        print('GammaDot  = ',nGd)
        print('eta+etaDot= ',nX + nXd)
        print('total     = ',xSize)
        
    x = np.zeros((xSize),ct.c_double,'C')
    x[:nG] = Gamma.flatten()
    x[nG:nG+nGs] = GammaStar.flatten()
    x[nG+nGs:nG+nGs+nGd] = GammaDot.flatten()
    x[nG+nGs+nGd:nG+nGs+nGd+nX] = X
    x[nG+nGs+nGd+nX:nG+nGs+nGd+nX+nXd] = dXdt
    return x

if __name__ == '__main__':
    # Beam options.
    XBOPTS = DerivedTypes.Xbopts(FollowerForce = ct.c_bool(False),
                                 MaxIterations = ct.c_int(50),
                                 PrintInfo = ct.c_bool(True),
                                 NumLoadSteps = ct.c_int(25),
                                 Solution = ct.c_int(312),
                                 MinDelta = ct.c_double(1e-5),
                                 NewmarkDamp = ct.c_double(5e-3))
    # beam inputs.
    XBINPUT = DerivedTypes.Xbinput(3,10)
    XBINPUT.BeamLength = 5*6.096
#     XBINPUT.BeamLength = 16.0
    XBINPUT.BeamStiffness[0,0] = 1.0e+09
    XBINPUT.BeamStiffness[1,1] = 1.0e+09
    XBINPUT.BeamStiffness[2,2] = 1.0e+09
    XBINPUT.BeamStiffness[3,3] = 0.99e+06
    XBINPUT.BeamStiffness[4,4] = 9.77e+06
    XBINPUT.BeamStiffness[5,5] = 1.0e+09
#     XBINPUT.BeamStiffness[0,0] = 1.0e+09
#     XBINPUT.BeamStiffness[1,1] = 1.0e+09
#     XBINPUT.BeamStiffness[2,2] = 1.0e+09
#     XBINPUT.BeamStiffness[3,3] = 1.0e+04
#     XBINPUT.BeamStiffness[4,4] = 2.0e+04
#     XBINPUT.BeamStiffness[5,5] = 5.0e+06
    XBINPUT.BeamStiffness[:,:] = 1.0*XBINPUT.BeamStiffness[:,:]
    XBINPUT.BeamMass[0,0] = 35.71
    XBINPUT.BeamMass[1,1] = 35.71
    XBINPUT.BeamMass[2,2] = 35.71
    XBINPUT.BeamMass[3,3] = 8.64
    XBINPUT.BeamMass[4,4] = 0.001
    XBINPUT.BeamMass[5,5] = 0.001
#     XBINPUT.BeamMass[0,0] = 0.75
#     XBINPUT.BeamMass[1,1] = 0.75
#     XBINPUT.BeamMass[2,2] = 0.75
#     XBINPUT.BeamMass[3,3] = 0.1
#     XBINPUT.BeamMass[4,4] = 0.001
#     XBINPUT.BeamMass[5,5] = 0.001
    # Off diagonal terms (in Theodorsen sectional coordinates)
    ElasticAxis = -0.34
#     ElasticAxis = -0.5
#     InertialAxis = -7.0/50.0
    InertialAxis =  -0.54
#     InertialAxis = -0.5
    x_alpha = InertialAxis - ElasticAxis
    # pitch-plunge coupling term (b-frame coordinates)
    c = 1.8288
    cgLoc = 0.5*c*np.array([0.0, x_alpha, 0.0])
    cgSkew = Skew(cgLoc)
    XBINPUT.BeamMass[:3,3:] = XBINPUT.BeamMass[0,0] * cgSkew
    XBINPUT.BeamMass[3:,:3] = XBINPUT.BeamMass[:3,3:].T
    XBINPUT.BeamMass[4,4] = 0.001 + XBINPUT.BeamMass[0,0]*pow(cgLoc[2],2.0)
    XBINPUT.BeamMass[5,5] = 0.001 + XBINPUT.BeamMass[0,0]*pow(cgLoc[1],2.0)
    
    #XBINPUT.g = 9.81
    
    # Get suggested panelling.
    Umag = 28.0
    M = 8
#     M, delTime = panellingFromFreq(70,c,Umag)
    delTime = c/(Umag*M)
    delTime = c/(Umag*M)
    # Unsteady parameters.
    XBINPUT.dt = delTime
    XBINPUT.t0 = 0.0
    XBINPUT.tfin = 5.0
    
    # Set motion of wing.
    alpha = 0.0*np.pi/180.0
    NumSteps = np.ceil( (XBINPUT.tfin + XBINPUT.dt - XBINPUT.t0) / XBINPUT.dt)
    XBINPUT.ForcedVel = np.zeros((NumSteps,6),ct.c_double,'F')
    for i in range(XBINPUT.ForcedVel.shape[0]):
        XBINPUT.ForcedVel[i,:] = [0.0, Umag*np.cos(alpha), -Umag*np.sin(alpha), 0.0, 0.0, 0.0]
    XBINPUT.ForcedVelDot = np.zeros((NumSteps,6),ct.c_double,'F')
     
     
    # aero params.
    WakeLength = 15.0*1.8 #c
#     WakeLength = 30.0
    Mstar = int(WakeLength/(delTime*Umag))
    # aero options.
    N = XBINPUT.NumNodesTot - 1
    VMOPTS = DerivedTypesAero.VMopts(M = M,
                                     N = N,
                                     ImageMethod = True,
                                     Mstar = Mstar,
                                     Steady = False,
                                     KJMeth = False,
                                     NewAIC = True,
                                     DelTime = delTime,
                                     NumCores = 4)
    # Aero inputs.
    iMin = M - M/4
    iMax = M
    jMin = 16#N - N/4
    jMax = N
    typeMotion = 'asInput'
    betaBar = 0.0*np.pi/180.0
    omega = np.pi
    ctrlSurf = ControlSurf(iMin,
                           iMax,
                           jMin,
                           jMax,
                           typeMotion,
                           betaBar,
                           omega)
    
    # Gust inputs
    gust = Gust(uMag = 17.07,
                l = 20.0,
                r = 0.0)
    
    VMINPUT = DerivedTypesAero.VMinput(c = c,
                                       b = XBINPUT.BeamLength,
                                       U_mag = 0.0,
                                       alpha = 0.0*np.pi/180.0,
                                       theta = 0.0,
                                       WakeLength = WakeLength,
                                       ctrlSurf = ctrlSurf,
                                       gust = gust)
    
    # Aerolastic simulation.
    AirDensity = 1.02
    AELAOPTS = AeroelasticOps(ElasticAxis = ElasticAxis,
                              InertialAxis = InertialAxis,
                              AirDensity = AirDensity,
                              Tight = False,
                              ImpStart = True)
    
    mpcCont = MPC.MPC('modifiedGolandControlOutput','/home/rjs10/git/SHARPy/src/PyMPC/systems/', LQR = True)
#     mpcCont = None
    
    # Live output options.
    writeDict = OrderedDict()
    writeDict['R_z_tip'] = 0
    writeDict['kappa_x_root'] = 0
    writeDict['kappa_y_root'] = 0
    writeDict['kappa_z_root'] = 0
    writeDict['u_opt_1'] = 0
    writeDict['du_opt_1'] = 0
    writeDict['contTime'] = 0
    
    Settings.OutputDir = Settings.SharPyProjectDir + "output/MPC/ModifiedGoland/testMPC/"
    Settings.OutputFileRoot = "Q28_M8N20_CS80_L20_W1707_Nmod8_LQR"
    
    # Solve nonlinear dynamic simulation.
    Solve_Py(XBINPUT, XBOPTS, VMOPTS, VMINPUT, AELAOPTS,
             writeDict =  writeDict, mpcCont = mpcCont)