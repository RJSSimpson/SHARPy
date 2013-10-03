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
import PyCoupled.Coupled_NlnStatic as Static
import PyBeam.Utils.XbeamLib as XbeamLib
from PyCoupled.Coupled_NlnStatic import AddGravityLoads
from DerivedTypesAero import ControlSurf
from scipy.linalg import expm, logm
from collections import OrderedDict
import re
from PyBeam.Utils.Misc import iNode2iElem

class VMCoupledUnstInput:
    """@brief Contains data for unsteady run of UVLM.
    @param WakeLength Length of wake in chordlengths.
    @param DelS non-dim timestep s = omega*c/U.
    @param NumChordLengths Number of chord lengths to travel 
    in prescribed simulation.
    @param VelA_G Velocity of reference frame.
    @param OmegaA_G Initial angular vel of reference frame.
    @param OriginA_G Origin of reference frame in G-frame.
    @param PsiA_G Orientation of reference frame in G-frame."""
    
    def __init__(self, VMOPTS, VMINPUT,
                  DelS, NumChordLengths,
                  VelA_G, OmegaA_G,
                  OriginA_G = np.zeros((3),ct.c_double),
                  PsiA_G = np.zeros((3),ct.c_double)):
        
        self.NumChordLengths = NumChordLengths
        self.VelA_G = VelA_G
        self.OmegaA_G = OmegaA_G
        self.OriginA_G = OriginA_G
        self.PsiA_G = PsiA_G

def Solve_Py(XBINPUT,XBOPTS,VMOPTS,VMINPUT,VMUNST,AELAOPTS,**kwords):
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
    @param writeDict OrderedDict of 'name':tuple of outputs to write."""
        
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
    Qglobal += -np.dot(FglobalFull, Force_Dof)
    
    # Initial Accel.
    dXddt[:] = np.dot(np.linalg.inv(MglobalFull), -Qglobal)
    
    # Record position of all grid points in global FoR at initial time step.
    DynOut[0:NumNodes_tot.value,:] = PosDefor
    
    # Record state of the selected node in initial deformed configuration.
    PosPsiTime[0,:3] = PosDefor[-1,:]
    PosPsiTime[0,3:] = PsiDefor[-1,XBELEM.NumNodes[-1]-1,:]
    
    # Get gamma and beta for Newmark scheme.
    gamma = 0.5 + XBOPTS.NewmarkDamp.value
    beta = 0.25*(gamma + 0.5)**2
    
    # Initialize Aero       
    Section = InitSection(VMOPTS,VMINPUT,AELAOPTS.ElasticAxis)
    
    # Declare memory for Aero variables.
    ZetaDot = np.zeros((Section.shape[0],PosDefor.shape[0],3),ct.c_double,'C')
    K = VMOPTS.M.value*VMOPTS.N.value
    AIC = np.zeros((K,K),ct.c_double,'C')
    BIC = np.zeros((K,K),ct.c_double,'C')
    AeroForces = np.zeros((VMOPTS.M.value+1,VMOPTS.N.value+1,3),ct.c_double,'C')
    
    # Init external velocities.  
    Uext = InitSteadyExternalVels(VMOPTS,VMINPUT)
    if AELAOPTS.ImpStart == True:
        Zeta = np.zeros((Section.shape[0],PosDefor.shape[0],3),ct.c_double,'C')             
        Gamma = np.zeros((VMOPTS.M.value,VMOPTS.N.value),ct.c_double,'C')
        # Generate surface, wake and gamma matrices.
        CoincidentGrid(PosDefor, PsiDefor, Section, ForcedVel[0,:3], 
                       ForcedVel[0,3:], PosDotDef, PsiDotDef, XBINPUT,
                       Zeta, ZetaDot, VMUNST.OriginA_G, VMUNST.PsiA_G,
                       VMINPUT.ctrlSurf)
        # init wake grid and gamma matrix.
        ZetaStar, GammaStar = InitSteadyWake(VMOPTS,VMINPUT,Zeta)
          
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
                                   Text = True)
    
    # Open output file for writing
    if 'writeDict' in kwords and Settings.WriteOut == True:
        writeDict = kwords['writeDict']
        ofile = Settings.OutputDir + \
                Settings.OutputFileRoot + \
                '_SOL312_out.dat'
        fp = open(ofile,'w')
        fp.write("{:<14}".format("Time"))
        for output in writeDict.keys():
            fp.write("{:<14}".format(output))
        fp.write("\n")
        fp.flush()
        
    # Write initial outputs to file.
    if 'writeDict' in kwords and Settings.WriteOut == True:
        locForces = None # Stops recalculation of forces
        fp.write("{:<14,e}".format(Time[0]))
        for myStr in writeDict.keys():
            if re.search(r'^R_.',myStr):
                if re.search(r'^R_._.', myStr):
                    index = int(myStr[4])
                elif re.search(r'root', myStr):
                    index = 0
                elif re.search(r'tip', myStr):
                    index = -1
                else:
                    raise IOError("Node index not recognised.")
                
                if myStr[2] == 'x':
                    component = 0
                elif myStr[2] == 'y':
                    component = 1
                elif myStr[2] == 'z':
                    component = 2
                else:
                    raise IOError("Displacement component not recognised.")
                
                fp.write("{:<14,e}".format(PosDefor[index,component]))
                
            elif re.search(r'^M_.',myStr):
                if re.search(r'^M_._.', myStr):
                    index = int(myStr[4])
                elif re.search(r'root', myStr):
                    index = 0
                elif re.search(r'tip', myStr):
                    index = -1
                else:
                    raise IOError("Node index not recognised.")
                
                if myStr[2] == 'x':
                    component = 0
                elif myStr[2] == 'y':
                    component = 1
                elif myStr[2] == 'z':
                    component = 2
                else:
                    raise IOError("Moment component not recognised.")
                
                if locForces == None:
                    locForces = localElasticForces(PosDefor,
                                                   PsiDefor,
                                                   XBELEM,
                                                   [index])
                
                fp.write("{:<14,e}".format(locForces[0,3+component]))
            else:
                raise IOError("writeDict key not recognised.")
        # END for myStr
        fp.write("\n")
        fp.flush()
    # END if write

    # Time loop.
    for iStep in range(NumSteps.value):
        
        if XBOPTS.PrintInfo.value==True:
            sys.stdout.write('Time: %-10.4e\n' %(Time[iStep+1]))
            sys.stdout.write('   SubIter DeltaF     DeltaX     ResLog10\n')
        
        dt = Time[iStep+1] - Time[iStep]
        
        # Set dt for aero force calcs.
        VMOPTS.DelTime = ct.c_double(dt)
        
        # Update transformation matrix for given angular velocity.
        Temp = np.linalg.inv(Unit4 + 0.25*XbeamLib.QuadSkew(ForcedVel[iStep+1,3:])*dt)
        Quat = np.dot(Temp, np.dot(Unit4 - 0.25*XbeamLib.QuadSkew(ForcedVel[iStep,3:])*dt, Quat))
        Quat = Quat/np.linalg.norm(Quat)
        Cao  = XbeamLib.Rot(Quat)
        
        # Force at current time-step
        if iStep > 0:
            
            # zero aero forces.
            AeroForces[:,:,:] = 0.0
            
            # Update geometry
            VMUNST.OriginA_G[:] += ForcedVel[iStep,:3]*dt
            VMINPUT.ctrlSurf.update(Time[iStep])
            CoincidentGrid(PosDefor, PsiDefor, Section, ForcedVel[iStep+1,:3], 
                           ForcedVel[iStep+1,3:], PosDotDef, PsiDotDef, XBINPUT,
                           Zeta, ZetaDot, VMUNST.OriginA_G, VMUNST.PsiA_G,
                           VMINPUT.ctrlSurf)
            
            # Update wake geom       
            #'roll' data.
            ZetaStar = np.roll(ZetaStar,1,axis = 0)
            GammaStar = np.roll(GammaStar,1,axis = 0)
            #overwrite grid points with TE.
            ZetaStar[0,:] = Zeta[VMOPTS.M.value,:]
            # overwrite Gamma with TE value from previous timestep.
            GammaStar[0,:] = Gamma[VMOPTS.M.value-1,:]
            
            # Solve for AeroForces
            UVLMLib.Cpp_Solver_VLM(Zeta, ZetaDot, Uext, ZetaStar, VMOPTS, 
                           AeroForces, Gamma, GammaStar, AIC, BIC)
            
            # Apply density scaling"
            AeroForces[:,:,:] = AELAOPTS.AirDensity*AeroForces[:,:,:]
            
            if Settings.PlotTec==True:
                PostProcess.WriteUVLMtoTec(FileObject,
                                           Zeta,
                                           ZetaStar,
                                           Gamma,
                                           GammaStar,
                                           TimeStep = iStep,
                                           NumTimeSteps = XBOPTS.NumLoadSteps.value,
                                           Time = Time[iStep],
                                           Text = True)

            # map AeroForces to beam.
            CoincidentGridForce(XBINPUT, PsiDefor, Section, AeroForces,
                                Force)
            
            # Add gravity loads.
            AddGravityLoads(Force, XBINPUT, XBELEM, AELAOPTS,
                            PsiDefor, VMINPUT.c)
        #END if iStep > 0
            
        # Predictor step.
        X += dt*dXdt + (0.5-beta)*dXddt*dt**2
        dXdt += (1.0-gamma)*dXddt*dt
        dXddt[:] = 0.0
        
        # Reset convergence parameters.
        Iter = 0
        ResLog10 = 0.0
            
        # Save Gamma for tightly coupled.
        if AELAOPTS.Tight == True:
            GammaSav = Gamma.copy(order = 'C')
        
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
                CoincidentGrid(PosDefor, PsiDefor, Section, ForcedVel[iStep+1,:3], 
                               ForcedVel[iStep+1,3:], PosDotDef, PsiDotDef, XBINPUT,
                               Zeta, ZetaDot, VMUNST.OriginA_G, VMUNST.PsiA_G,
                               VMINPUT.ctrlSurf)
                
                # close wake.
                ZetaStar[0,:] = Zeta[VMOPTS.M.value,:]
                
                # save pereference and turn off rollup.
                Rollup = VMOPTS.Rollup.value
                VMOPTS.Rollup.value = False
                
                # Solve for AeroForces.
                UVLMLib.Cpp_Solver_VLM(Zeta, ZetaDot, Uext, ZetaStar, VMOPTS, 
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
            
            
            # Update matrices.
            BeamLib.Cbeam3_Asbly_Dynamic(XBINPUT, NumNodes_tot, XBELEM, XBNODE,
                         PosIni, PsiIni, PosDefor, PsiDefor,
                         PosDotDef, PsiDotDef, PosDotDotDef, PsiDotDotDef,
                         Force, ForcedVel[iStep+1,:], ForcedVelDot[iStep+1,:],
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
            Qglobal += np.dot(MglobalFull, dXddt) \
                        + np.dot(Mvel,ForcedVelDot[iStep+1,:]) \
                        - np.dot(FglobalFull, Force_Dof)
                              
            if XBOPTS.PrintInfo.value==True:                 
                sys.stdout.write('%-10.4e ' %(max(abs(Qglobal))))
            
            
            # Calculate system matrix for update calculation.
            Asys = KglobalFull \
                    + CglobalFull*gamma/(beta*dt) \
                    + MglobalFull/(beta*dt**2)
                      
            
            # Solve for update.
            DX[:] = np.dot(np.linalg.inv(Asys), -Qglobal)
            
            # Corrector step"
            X += DX
            dXdt += DX*gamma/(beta*dt)
            dXddt += DX/(beta*dt**2)
            
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
            locForces = None # Stops recalculation of forces
            fp.write("{:<14,e}".format(Time[iStep+1]))
            for myStr in writeDict.keys():
                if re.search(r'^R_.',myStr):
                    if re.search(r'^R_._.', myStr):
                        index = int(myStr[4])
                    elif re.search(r'root', myStr):
                        index = 0
                    elif re.search(r'tip', myStr):
                        index = -1
                    else:
                        raise IOError("Node index not recognised.")
                    
                    if myStr[2] == 'x':
                        component = 0
                    elif myStr[2] == 'y':
                        component = 1
                    elif myStr[2] == 'z':
                        component = 2
                    else:
                        raise IOError("Displacement component not recognised.")
                    
                    fp.write("{:<14,e}".format(PosDefor[index,component]))
                    
                elif re.search(r'^M_.',myStr):
                    if re.search(r'^M_._.', myStr):
                        index = int(myStr[4])
                    elif re.search(r'root', myStr):
                        index = 0
                    elif re.search(r'tip', myStr):
                        index = -1
                    else:
                        raise IOError("Node index not recognised.")
                    
                    if myStr[2] == 'x':
                        component = 0
                    elif myStr[2] == 'y':
                        component = 1
                    elif myStr[2] == 'z':
                        component = 2
                    else:
                        raise IOError("Moment component not recognised.")
                    
                    if locForces == None:
                        locForces = localElasticForces(PosDefor,
                                                       PsiDefor,
                                                       XBELEM,
                                                       [index])
                    
                    fp.write("{:<14,e}".format(locForces[0,3+component]))
                else:
                    raise IOError("writeDict key not recognised.")
            # END for myStr
            fp.write("\n")
            fp.flush()
                    
            
        
        # 'Rollup' due to external velocities. TODO: Must add gusts here!
        ZetaStar[:,:] += VMINPUT.U_infty*dt
    
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

def localStrains(PosDefor, PsiDef, XBELEM, NodeList,
                  SO3 = False):
    """@brief Approximate strains at the midpoint of nodes.
    
    @param PosDefor Deformed nodal coordinates.
    @param PsiDef Deformed nodal rotation vectors.
    @param XBELEM Xbeam element derived type containing element data.
    @param NodeList List of node numbers for strain approximation at adjacent
            midpoint in the direction of increasing node index,
            starts from zero:
            0---x---1---x---2---x---3-- ... --(NumNodes-1)
            Strains are calculated at the x to the right of selected nodes.
    @param SO3 Flag for use of the SO(3) manifold in CRV interpolation.
    
    @details Assumes that all nodes are either two- or three- noded. 
    @warning Untested.
    """
    
    strains = np.zeros((len(NodeList),6))
    NumNodes = PosDefor.shape[0]-1
    NumNodesElem = XBELEM.NumNodes[0]
    iOut = 0 # Output index
    for iNode in NodeList:
        if iNode == NumNodes:
            raise ValueError("Midpoint strain requested beyond final node.")
        
        iElem, iiElem = iNode2iElem(iNode, NumNodes, NumNodesElem)
        
        if NumNodesElem == 2:
            
            Rdash = (PosDefor[iNode+1] - PosDefor[iNode]) / XBELEM.Length[iElem]
            # Consistent calculation of midpoint strains using SO(3) manifold.
            if SO3 == True:
#                 # Transformation at node 1
#                 CaB = XbeamLib.Psi2TransMat(PsiDef[iElem,iiElem,:])
#                 CB1a = CaB.T
#                 
#                 # Transformation at node 2
#                 CaB2 = XbeamLib.Psi2TransMat(PsiDef[iElem,iiElem+1,:])
#                 
#                 # Transformation between nodes 2 and 1
#                 CB1B2 = np.dot(CB1a,CaB2)
#                 
#                 # Crisfield and Jelenic 1999, pg.1134
#                 phiMat = logm(CB1B2) # skew symmetric matrix of rotation vector from node 1 -> node 2
#                 CaBmid = np.dot(CaB,expm(0.5*phiMat))
#                 CBaMid = CaBmid.T
#                 momentStrain = ( 1.0/XBELEM.Length[iElem] ) \
#                 			   *XbeamLib.VectofSkew(phiMat)
                raise NotImplementedError("SO(3) manifold strains.")
            else:
                # Interpolate the CRV as if it were a standard vector.
                psiDash = (PsiDef[iElem,iiElem+1,:] - PsiDef[iElem,iiElem,:]) \
                          /XBELEM.Length[iElem]
                psiMid = 0.5*(PsiDef[iElem,iiElem,:] + PsiDef[iElem,iiElem+1,:])
                CaBmid = XbeamLib.Psi2TransMat(psiMid)
                CBaMid = CaBmid.T
                momentStrain = np.dot(XbeamLib.Tangential(psiMid), psiDash) \
                			   - XBELEM.PreCurv[iElem*3:(iElem+1)*3]
            # END if SO3
            forceStrain = np.dot(CBaMid,Rdash) - np.array([1.0, 0.0 ,0.0])
            
        elif NumNodesElem == 3:
            Rdash = (PosDefor[iNode+1] - PosDefor[iNode-1]) / XBELEM.Length[iElem]
            
            # Work out what sub-element node (iiElem) we are in.
            iElem, iiElem = iNode2iElem(iNode, NumNodes, NumNodesElem)
            
            
            
            
            if SO3 == True:
                raise NotImplementedError("Only available for 2-noded just now.")
            else:
                #TODO 3-noded calc
                raise NotImplementedError("Only available for 2-noded just now.")
        else:
            raise ValueError("Only 2- or 3-noded elements are supported.")
        
        strains[iOut,:3] = np.real(forceStrain) #matrix log or exponent gives small imaginary part
        strains[iOut,3:] = np.real(momentStrain)
        iOut += 1
    # END for iElem
    return strains
    
def localElasticForces(PosDefor, PsiDef, XBELEM, ElemList):
    """@brief Approximate shear force and moments from local strain and
    stiffness matrix.
    
    @param PosDefor Deformed nodal coordinates.
    @param PsiDef Deformed nodal rotation vectors.
    @param XBELEM Xbeam element derived type containing element data.
    @param ElemList List of element numbers for strain approximation,
            starts from zero.
    @warning Untested.
    """
    
    elemF = np.zeros((len(ElemList),6))
    strains = localStrains(PosDefor, PsiDef, XBELEM, ElemList)
    iOut = 0 # Output index
    for iElem in ElemList:
        elemStrain = strains[iOut,:]
        elemK = XBELEM.Stiff[iElem*6:(iElem+1)*6,:]
        elemF[iOut,:] = np.dot(elemK,elemStrain)
        iOut += 1
    # END of iElem
    return elemF

if __name__ == '__main__':
    # Beam options.
    XBOPTS = DerivedTypes.Xbopts(FollowerForce = ct.c_bool(False),
                                 MaxIterations = ct.c_int(50),
                                 PrintInfo = ct.c_bool(True),
                                 NumLoadSteps = ct.c_int(25),
                                 Solution = ct.c_int(312),
                                 MinDelta = ct.c_double(1e-5))
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
    M, delTime = panellingFromFreq(70,c,Umag)
    print("suggested M = %d\n" % M)
    print("suggested delTime = %f\n" % delTime)
    # Unsteady parameters.
    XBINPUT.dt = delTime
    XBINPUT.t0 = 0.0
    XBINPUT.tfin = 0.5
    # aero params.
    WakeLength = 30.0*c
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
                                     DelTime = delTime)
    # Aero inputs.
    iMin = M - M/4
    iMax = M
    jMin = N - N/4
    jMax = N
    typeMotion = 'sin'
    betaBar = 1.0*np.pi/180.0
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
                                       alpha = 0.0*np.pi/180.0,
                                       theta = 0.0,
                                       ZetaDotTest = 0.0,
                                       WakeLength = WakeLength,
                                       ctrlSurf = ctrlSurf)
    # Unsteady aero inputs.
    VelA_G   = np.array(([0.0,0.0,0.0]))
    OmegaA_G = np.array(([0.0,0.0,0.0]))
    VMUNST   = VMCoupledUnstInput(VMOPTS, VMINPUT, 0.0, 0.0,
                                  VelA_G, OmegaA_G)
    # Aerolastic simulation results.
    AELAOPTS = AeroelasticOps(ElasticAxis = ElasticAxis,
                              InertialAxis = InertialAxis,
                              AirDensity = 1.02,
                              Tight = False,
                              ImpStart = False)
    # Live output options.
    writeDict = OrderedDict()
    writeDict['R_z (tip)'] = 0
    writeDict['M_x (root)'] = 0
    writeDict['M_y (root)'] = 0
    writeDict['M_z (root)'] = 0
    
    # Solve nonlinear dynamic simulation.
    Solve_Py(XBINPUT, XBOPTS, VMOPTS, VMINPUT, VMUNST, AELAOPTS,
             writeDict =  writeDict)