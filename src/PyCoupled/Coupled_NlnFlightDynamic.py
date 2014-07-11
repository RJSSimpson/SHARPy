'''@package PyCoupled.Coupled_NlnFlightDynamic
@brief      NonlinearDynamic Beam + Rigid-body dynamics + UVLM.
@author     Rob Simpson & Henrik Hesse
@contact    r.simpson11@imperial.ac.uk 
@version    0.0
@date       07/11/2013
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

def Solve_Py(XBINPUT,XBOPTS,VMOPTS,VMINPUT,AELAOPTS,**kwords):
    """@brief Nonlinear dynamic solver using Python to solve aeroelastic
    equation.
    @details Assembly of structural matrices is carried out with 
    Fortran subroutines. Aerodynamics solved using PyAero\.UVLM.
    @param XBINPUT Beam inputs (for initialization in Python).
    @param XBOPTS Beam solver options (for Fortran).
    @param VMOPTS UVLM solver options (for C/C++).
    @param VMINPUT UVLM solver inputs (for initialization in Python).
    @param VMUNST Unsteady input information for aero solver.
    @param AELAOPTS Options relevant to coupled aeroelastic simulations.
    @param writeDict OrderedDict of 'name':tuple of outputs to write.
    """
        
    # Check correct solution code.
    assert XBOPTS.Solution.value == 912, ('NonlinearFlightDynamic requested' +
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
        XBOPTS.Solution.value = 912 # Reset options.
        VMOPTS.Steady = ct.c_bool(False)
        VMOPTS.Rollup.value = Rollup
    elif AELAOPTS.ImpStart == True:
        PosDefor = PosIni.copy(order='F')
        PsiDefor = PsiIni.copy(order='F')
        Force = np.zeros((XBINPUT.NumNodesTot,6),ct.c_double,'F')
        
    # Write deformed configuration to file. TODO: tidy this away inside function.
    ofile = Settings.OutputDir + Settings.OutputFileRoot + '_SOL912_def.dat'
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
    Time, NumSteps, ForceTime, Vrel, VrelDot,\
    PosDotDef, PsiDotDef,\
    OutGrids, PosPsiTime, VelocTime, DynOut\
        = BeamInit.Dynamic(XBINPUT,XBOPTS)
    # Delete unused variables.
    del OutGrids, VelocTime
    
    if XBOPTS.PrintInfo.value==True:
        sys.stdout.write('Solve nonlinear dynamic case in Python ... \n')
    
    
    #Initialise structural system tensors
    MssFull = np.zeros((NumDof.value,NumDof.value), ct.c_double, 'F')
    CssFull = np.zeros((NumDof.value,NumDof.value), ct.c_double, 'F') 
    KssFull = np.zeros((NumDof.value,NumDof.value), ct.c_double, 'F') 
    FstrucFull = np.zeros((NumDof.value,NumDof.value), ct.c_double, 'F')
    
    ms = ct.c_int()
    cs = ct.c_int()
    ks = ct.c_int()
    fs = ct.c_int()
    
    Msr = np.zeros((NumDof.value,6), ct.c_double, 'F')
    Csr = np.zeros((NumDof.value,6), ct.c_double, 'F')
    
    X     = np.zeros(NumDof.value, ct.c_double, 'F')
    dXdt  = np.zeros(NumDof.value, ct.c_double, 'F')
    Force_Dof = np.zeros(NumDof.value, ct.c_double, 'F')
    
    Qstruc = np.zeros(NumDof.value, ct.c_double, 'F')
    
    #Initialise rigid-body system tensors
    MrsFull = np.zeros((6,NumDof.value), ct.c_double, 'F')
    CrsFull = np.zeros((6,NumDof.value), ct.c_double, 'F') 
    KrsFull = np.zeros((6,NumDof.value), ct.c_double, 'F') 
    FrigidFull = np.zeros((6,NumDof.value+6), ct.c_double, 'F')
    
    mr = ct.c_int()
    cr = ct.c_int()
    kr = ct.c_int()
    fr = ct.c_int()
    
    Mrr = np.zeros((6,6), ct.c_double, 'F')
    Crr = np.zeros((6,6), ct.c_double, 'F')
        
    Qrigid = np.zeros(6, ct.c_double, 'F')
    
    #Initialise full system tensors
    Q     = np.zeros(NumDof.value+6+4, ct.c_double, 'F')
    DQ    = np.zeros(NumDof.value+6+4, ct.c_double, 'F')
    dQdt  = np.zeros(NumDof.value+6+4, ct.c_double, 'F')
    dQddt = np.zeros(NumDof.value+6+4, ct.c_double, 'F')
    Force_All = np.zeros(NumDof.value+6, ct.c_double, 'F')

    Msys = np.zeros((NumDof.value+6+4,NumDof.value+6+4), ct.c_double, 'F')
    Csys = np.zeros((NumDof.value+6+4,NumDof.value+6+4), ct.c_double, 'F') 
    Ksys = np.zeros((NumDof.value+6+4,NumDof.value+6+4), ct.c_double, 'F') 
    Asys = np.zeros((NumDof.value+6+4,NumDof.value+6+4), ct.c_double, 'F')
    
    Qsys = np.zeros(NumDof.value+6+4, ct.c_double, 'F')
    
    #Initialise rotation operators. TODO: include initial AOA here
    currVrel=Vrel[0,:].copy('F')
    AOA  = np.arctan(currVrel[2]/-currVrel[1])
    Quat = xbl.Euler2Quat(AOA,0,0)
    Cao  = xbl.Rot(Quat)
    ACoa = np.zeros((6,6), ct.c_double, 'F')
    ACoa[:3,:3] = np.transpose(Cao)
    ACoa[3:,3:] = np.transpose(Cao)
    Cqr = np.zeros((4,6), ct.c_double, 'F')
    Cqq = np.zeros((4,4), ct.c_double, 'F')
        
    Unit4 = np.zeros((4,4), ct.c_double, 'F')
    for i in range(4):
        Unit4[i,i] = 1.0
    
    # Extract initial displacements and velocities.
    BeamLib.Cbeam_Solv_Disp2State(NumNodes_tot, NumDof, XBINPUT, XBNODE,
                                  PosDefor, PsiDefor, PosDotDef, PsiDotDef,
                                  X, dXdt)
    
    # Approximate initial accelerations.
    PosDotDotDef = np.zeros((NumNodes_tot.value,3),ct.c_double,'F')
    PsiDotDotDef = np.zeros((XBINPUT.NumElems,Settings.MaxElNod,3),
                             ct.c_double, 'F')
    
    #Populate state vector
    Q[:NumDof.value]=X.copy('F')
    dQdt[:NumDof.value]=dXdt.copy('F')
    dQdt[NumDof.value:NumDof.value+6] = Vrel[0,:].copy('F')
    dQdt[NumDof.value+6:]= Quat.copy('F')
    
    #Force at the first time-step
    Force += (XBINPUT.ForceDyn*ForceTime[0]).copy('F')
    

    #Assemble matrices and loads for structural dynamic analysis
    currVrel=Vrel[0,:].copy('F')
    tmpQuat=Quat.copy('F')
    BeamLib.Cbeam3_Asbly_Dynamic(XBINPUT, NumNodes_tot, XBELEM, XBNODE,\
                         PosIni, PsiIni, PosDefor, PsiDefor,\
                         PosDotDef, PsiDotDef, PosDotDotDef, PsiDotDotDef,\
                         Force, currVrel, 0*currVrel,\
                         NumDof, Settings.DimMat,\
                         ms, MssFull, Msr,\
                         cs, CssFull, Csr,\
                         ks, KssFull, fs, FstrucFull,\
                         Qstruc, XBOPTS, Cao)
       
    BeamLib.f_fem_m2v(ct.byref(NumNodes_tot),\
                              ct.byref(ct.c_int(6)),\
                              Force.ctypes.data_as(ct.POINTER(ct.c_double)),\
                              ct.byref(NumDof),\
                              Force_Dof.ctypes.data_as(ct.POINTER(ct.c_double)),\
                              XBNODE.Vdof.ctypes.data_as(ct.POINTER(ct.c_int)) )

    Qstruc -= np.dot(FstrucFull, Force_Dof)
    
    
    #Assemble matrices for rigid-body dynamic analysis
    BeamLib.Xbeam_Asbly_Dynamic(XBINPUT, NumNodes_tot, XBELEM, XBNODE,\
                         PosIni, PsiIni, PosDefor, PsiDefor,\
                         PosDotDef, PsiDotDef, PosDotDotDef, PsiDotDotDef,\
                         currVrel, 0*currVrel, tmpQuat,\
                         NumDof, Settings.DimMat,\
                         mr, MrsFull, Mrr,\
                         cr, CrsFull, Crr, Cqr, Cqq,\
                         kr, KrsFull, fr, FrigidFull,\
                         Qrigid, XBOPTS, Cao)
    
    BeamLib.f_fem_m2v_nofilter(ct.byref(NumNodes_tot),\
                               ct.byref(ct.c_int(6)),\
                               Force.ctypes.data_as(ct.POINTER(ct.c_double)),\
                               ct.byref(ct.c_int(NumDof.value+6)),\
                               Force_All.ctypes.data_as(ct.POINTER(ct.c_double)) )

    Qrigid -= np.dot(FrigidFull, Force_All)
        
          
#     #Separate assembly of follower and dead loads   
#     tmpForceTime=ForceTime[0].copy('F') 
#     tmpQforces,Dummy,tmpQrigid = xbl.LoadAssembly(XBINPUT, XBELEM, XBNODE, XBOPTS, NumDof, \
#                                     PosIni, PsiIni, PosDefor, PsiDefor, \
#                                     (XBINPUT.ForceStatic_foll + XBINPUT.ForceDyn_foll*tmpForceTime), \
#                                     (XBINPUT.ForceStatic_dead + XBINPUT.ForceDyn_dead*tmpForceTime), \
#                                     Cao,1)
#                            
#     Qstruc -= tmpQforces      
#     Qrigid -= tmpQrigid
    
    
    #Assemble system matrices
    Msys[:NumDof.value,:NumDof.value] = MssFull.copy('F')
    Msys[:NumDof.value,NumDof.value:NumDof.value+6] = Msr.copy('F')
    Msys[NumDof.value:NumDof.value+6,:NumDof.value] = MrsFull.copy('F')
    Msys[NumDof.value:NumDof.value+6,NumDof.value:NumDof.value+6] = Mrr.copy('F')
    Msys[NumDof.value+6:,NumDof.value+6:] = Unit4.copy('F')
       
    Qsys[:NumDof.value] = Qstruc
    Qsys[NumDof.value:NumDof.value+6] = Qrigid
    Qsys[NumDof.value+6:] = np.dot(Cqq,dQdt[NumDof.value+6:])
       

    # Initial Accel.
    dQddt[:] = np.dot(np.linalg.inv(Msys), -Qsys)
    
    
    #Record position of all grid points in global FoR at initial time step
    DynOut[0:NumNodes_tot.value,:] = PosDefor
    
    #Position/rotation of the selected node in initial deformed configuration
    PosPsiTime[0,:3] = PosDefor[-1,:]
    PosPsiTime[0,3:] = PsiDefor[-1,XBELEM.NumNodes[-1]-1,:]
    
    
    #Get gamma and beta for Newmark scheme
    gamma = 0.5 + XBOPTS.NewmarkDamp.value
    beta = 0.25*(gamma + 0.5)**2
    
    
    # Initialise Aero       
    Section = InitSection(VMOPTS,VMINPUT,AELAOPTS.ElasticAxis)
    
    # Declare memory for Aero variables.
    ZetaDot = np.zeros((Section.shape[0],PosDefor.shape[0],3),ct.c_double,'C')
    K = VMOPTS.M.value*VMOPTS.N.value
    AIC = np.zeros((K,K),ct.c_double,'C')
    BIC = np.zeros((K,K),ct.c_double,'C')
    AeroForces = np.zeros((VMOPTS.M.value+1,VMOPTS.N.value+1,3),ct.c_double,'C')
    
    # Initialise A-frame location and orientation to be zero.
    OriginA_G = np.zeros(3,ct.c_double,'C')
    PsiA_G = xbl.quat2psi(Quat) # CRV at iStep
    
    # Init external velocities.  
    Ufree = InitSteadyExternalVels(VMOPTS,VMINPUT)
    if AELAOPTS.ImpStart == True:
        Zeta = np.zeros((Section.shape[0],PosDefor.shape[0],3),ct.c_double,'C')             
        Gamma = np.zeros((VMOPTS.M.value,VMOPTS.N.value),ct.c_double,'C')
        # Generate surface, wake and gamma matrices.
        CoincidentGrid(PosDefor, PsiDefor, Section, currVrel[:3], 
                       currVrel[3:], PosDotDef, PsiDotDef, XBINPUT,
                       Zeta, ZetaDot, OriginA_G, PsiA_G,
                       VMINPUT.ctrlSurf)
        # init wake grid and gamma matrix.
        ZetaStar, GammaStar = InitSteadyWake(VMOPTS,VMINPUT,Zeta,currVrel[:3])
          
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
                '_SOL912_out.dat'
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
                    locForces = BeamIO.localElasticForces(PosDefor,
                                                          PsiDefor,
                                                          PosIni,
                                                          PsiIni,
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
        
        #calculate dt
        dt = Time[iStep+1] - Time[iStep]
        
        # Set dt for aero force calcs.
        VMOPTS.DelTime = ct.c_double(dt)
        
        #Predictor step
        Q       += dt*dQdt + (0.5-beta)*dQddt*np.power(dt,2.0)
        dQdt    += (1.0-gamma)*dQddt*dt
        dQddt[:] = 0.0
        
        # Quaternion update for orientation.
        Quat = dQdt[NumDof.value+6:].copy('F')
        Quat = Quat/np.linalg.norm(Quat)
        Cao  = xbl.Rot(Quat)
        
        #nodal diplacements and velocities from state vector
        X=Q[:NumDof.value].copy('F') 
        dXdt=dQdt[:NumDof.value].copy('F'); 
        BeamLib.Cbeam3_Solv_State2Disp(XBINPUT,NumNodes_tot,XBELEM,XBNODE,
                                       PosIni,PsiIni,NumDof,X,dXdt,
                                       PosDefor,PsiDefor,PosDotDef,PsiDotDef)
            
        # Force at current time-step. TODO: Check communication flow. 
        if iStep > 0 and AELAOPTS.Tight == False:
            
            # zero aero forces.
            AeroForces[:,:,:] = 0.0
            
            # Update CRV.
            PsiA_G = xbl.quat2psi(Quat) # CRV at iStep
            
            # Update origin.
            currVrel=Vrel[iStep-1,:].copy('F')
            OriginA_G[:] = OriginA_G[:] + currVrel[:3]*dt
            
            # Update control surface deflection.
            if VMINPUT.ctrlSurf != None:
                VMINPUT.ctrlSurf.update(Time[iStep])
            
            # Generate surface grid.
            currVrel=Vrel[iStep,:].copy('F')
            CoincidentGrid(PosDefor, PsiDefor, Section, currVrel[:3], 
                           currVrel[3:], PosDotDef, PsiDotDef, XBINPUT,
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
            
            # Solve for AeroForces
            UVLMLib.Cpp_Solver_VLM(Zeta, ZetaDot, Utot, ZetaStar, VMOPTS, 
                           AeroForces, Gamma, GammaStar, AIC, BIC)
            
            # Apply density scaling
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
            
            # Add thrust and other point loads
            Force += (XBINPUT.ForceStatic + 
                      XBINPUT.ForceDyn*ForceTime[iStep+1]).copy('F')
        #END if iStep > 0
            
            
        #Reset convergence parameters
        Iter = 0
        ResLog10 = 1.0
        
        
        #Newton-Raphson loop      
        while ( (ResLog10 > XBOPTS.MinDelta.value) \
                & (Iter < XBOPTS.MaxIterations.value) ):
                                    
            #set tensors to zero 
            MssFull[:,:] = 0.0; CssFull[:,:] = 0.0
            KssFull[:,:] = 0.0; FstrucFull[:,:] = 0.0
            Msr[:,:] = 0.0; Csr[:,:] = 0.0
            Qstruc[:] = 0.0
            
            MrsFull[:,:] = 0.0; CrsFull[:,:] = 0.0
            KrsFull[:,:] = 0.0; FrigidFull[:,:] = 0.0
            Mrr[:,:] = 0.0; Crr[:,:] = 0.0
            Qrigid[:] = 0.0
    
            Msys[:,:] = 0.0; Csys[:,:] = 0.0
            Ksys[:,:] = 0.0; Asys[:,:] = 0.0;
            Qsys[:] = 0.0
            
            # Update counter.
            Iter += 1
            
            if XBOPTS.PrintInfo.value==True:
                sys.stdout.write('   %-7d ' %(Iter))
                        
            #nodal diplacements and velocities from state vector
            X=Q[:NumDof.value].copy('F') 
            dXdt=dQdt[:NumDof.value].copy('F'); 
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


            #rigid-body velocities and orientation from state vector
            Vrel[iStep+1,:]    = dQdt[NumDof.value:NumDof.value+6].copy('F')
            VrelDot[iStep+1,:] = dQddt[NumDof.value:NumDof.value+6].copy('F')
            Quat = dQdt[NumDof.value+6:].copy('F')
            Quat = Quat/np.linalg.norm(Quat)
            Cao  = xbl.Rot(Quat)


            #Update matrices and loads for structural dynamic analysis
            tmpVrel=Vrel[iStep+1,:].copy('F')
            tmpQuat=Quat.copy('F')
            BeamLib.Cbeam3_Asbly_Dynamic(XBINPUT, NumNodes_tot, XBELEM, XBNODE,\
                                 PosIni, PsiIni, PosDefor, PsiDefor,\
                                 PosDotDef, PsiDotDef, PosDotDotDef, PsiDotDotDef,\
                                 Force, tmpVrel, 0*tmpVrel,\
                                 NumDof, Settings.DimMat,\
                                 ms, MssFull, Msr,\
                                 cs, CssFull, Csr,\
                                 ks, KssFull, fs, FstrucFull,\
                                 Qstruc, XBOPTS, Cao)
            
            BeamLib.f_fem_m2v(ct.byref(NumNodes_tot),\
                              ct.byref(ct.c_int(6)),\
                              Force.ctypes.data_as(ct.POINTER(ct.c_double)),\
                              ct.byref(NumDof),\
                              Force_Dof.ctypes.data_as(ct.POINTER(ct.c_double)),\
                              XBNODE.Vdof.ctypes.data_as(ct.POINTER(ct.c_int)) )
                    
            
            #Update matrices for rigid-body dynamic analysis
            BeamLib.Xbeam_Asbly_Dynamic(XBINPUT, NumNodes_tot, XBELEM, XBNODE,\
                                 PosIni, PsiIni, PosDefor, PsiDefor,\
                                 PosDotDef, PsiDotDef, PosDotDotDef, PsiDotDotDef,\
                                 tmpVrel, 0*tmpVrel, tmpQuat,\
                                 NumDof, Settings.DimMat,\
                                 mr, MrsFull, Mrr,\
                                 cr, CrsFull, Crr, Cqr, Cqq,\
                                 kr, KrsFull, fs, FrigidFull,\
                                 Qrigid, XBOPTS, Cao)
    
            BeamLib.f_fem_m2v_nofilter(ct.byref(NumNodes_tot),\
                                       ct.byref(ct.c_int(6)),\
                                       Force.ctypes.data_as(ct.POINTER(ct.c_double)),\
                                       ct.byref(ct.c_int(NumDof.value+6)),\
                                       Force_All.ctypes.data_as(ct.POINTER(ct.c_double)) )
        
        
            #Residual at first iteration
            if(Iter == 1):
                Res0_Qglobal = max(max(abs(Qsys)),1)
                Res0_DeltaX  = max(max(abs(DQ)),1)
              
            
            #Assemble discrete system matrices with linearised quaternion equations          
            Msys[:NumDof.value,:NumDof.value] = MssFull.copy('F')
            Msys[:NumDof.value,NumDof.value:NumDof.value+6] = Msr.copy('F')
            Msys[NumDof.value:NumDof.value+6,:NumDof.value] = MrsFull.copy('F')
            Msys[NumDof.value:NumDof.value+6,NumDof.value:NumDof.value+6] = Mrr.copy('F')
            Msys[NumDof.value+6:,NumDof.value+6:] = Unit4.copy('F')
            
            Csys[:NumDof.value,:NumDof.value] = CssFull.copy('F')
            Csys[:NumDof.value,NumDof.value:NumDof.value+6] = Csr.copy('F')
            Csys[NumDof.value:NumDof.value+6,:NumDof.value] = CrsFull.copy('F')
            Csys[NumDof.value:NumDof.value+6,NumDof.value:NumDof.value+6] = Crr.copy('F')
            
            Csys[NumDof.value+6:,NumDof.value:NumDof.value+6] = Cqr.copy('F')
            Csys[NumDof.value+6:,NumDof.value+6:] = Cqq.copy('F')
            
            Ksys[:NumDof.value,:NumDof.value] = KssFull.copy('F')
            Ksys[NumDof.value:NumDof.value+6,:NumDof.value] = KrsFull.copy('F')
            
          
#             #Separate assembly of follower and dead loads   
#             tmpForceTime=ForceTime[iStep+1].copy('F') 
#             tmpQforces,Dummy,tmpQrigid = xbl.LoadAssembly(XBINPUT, XBELEM, XBNODE, XBOPTS, NumDof, \
#                                             PosIni, PsiIni, PosDefor, PsiDefor, \
#                                             (XBINPUT.ForceStatic_foll + XBINPUT.ForceDyn_foll*tmpForceTime), \
#                                             (XBINPUT.ForceStatic_dead + XBINPUT.ForceDyn_dead*tmpForceTime), \
#                                             Cao,1)
#                                    
#             Qstruc -= tmpQforces      
#             Qrigid -= tmpQrigid
    
            
            #Compute residual to solve update vector
            Qstruc += -np.dot(FstrucFull, Force_Dof)
            Qrigid += -np.dot(FrigidFull, Force_All)
            
            Qsys[:NumDof.value] = Qstruc
            Qsys[NumDof.value:NumDof.value+6] = Qrigid
            Qsys[NumDof.value+6:] = np.dot(Cqq,dQdt[NumDof.value+6:])
            
            Qsys += np.dot(Msys,dQddt)

                
            #Calculate system matrix for update calculation
            Asys = Ksys + \
                      Csys*gamma/(beta*dt) + \
                      Msys/(beta*dt**2)
                      
            
            #Compute correction
            DQ[:] = np.dot(np.linalg.inv(Asys), -Qsys)

            Q += DQ
            dQdt += DQ*gamma/(beta*dt)
            dQddt += DQ/(beta*dt**2)
            
            
            #Update convergence criteria
            if XBOPTS.PrintInfo.value==True:                 
                sys.stdout.write('%-10.4e ' %(max(abs(Qsys))))
            
            Res_Qglobal = max(abs(Qsys))
            Res_DeltaX  = max(abs(DQ))
            
            ResLog10 = max(Res_Qglobal/Res0_Qglobal,Res_DeltaX/Res0_DeltaX)
            
            if XBOPTS.PrintInfo.value==True:
                sys.stdout.write('%-10.4e %8.4f\n' %(max(abs(DQ)),ResLog10))

        # END Netwon-Raphson.
                
                
        #update to converged nodal displacements and velocities
        X=Q[:NumDof.value].copy('F') 
        dXdt=dQdt[:NumDof.value].copy('F'); 
        BeamLib.Cbeam3_Solv_State2Disp(XBINPUT, NumNodes_tot, XBELEM, XBNODE,
                           PosIni, PsiIni, NumDof, X, dXdt,\
                           PosDefor, PsiDefor, PosDotDef, PsiDotDef)
        
        PosPsiTime[iStep+1,:3] = PosDefor[-1,:]
        PosPsiTime[iStep+1,3:] = PsiDefor[-1,XBELEM.NumNodes[-1]-1,:]
        
        #Position of all grid points in global FoR
        i1 = (iStep+1)*NumNodes_tot.value
        i2 = (iStep+2)*NumNodes_tot.value
        DynOut[i1:i2,:] = PosDefor
        
        #Export rigid-body velocities/accelerations
        if XBOPTS.OutInaframe.value==True:
            Vrel[iStep,:] = dQdt[NumDof.value:NumDof.value+6].copy('F')
            VrelDot[iStep,:] = dQddt[NumDof.value:NumDof.value+6].copy('F')
        else:
            Quat = dQdt[NumDof.value+6:].copy('F')
            Quat = Quat/np.linalg.norm(Quat)
            Cao  = xbl.Rot(Quat)
            ACoa[:3,:3] = np.transpose(Cao)
            ACoa[3:,3:] = np.transpose(Cao)
            
            Vrel[iStep,:] = np.dot(ACoa,dQdt[NumDof.value:NumDof.value+6].copy('F'))
            VrelDot[iStep,:] = np.dot(ACoa,dQddt[NumDof.value:NumDof.value+6].copy('F'))
            
        # Write selected outputs
        # tidy this away using function.
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
                        locForces = BeamIO.localElasticForces(PosDefor,
                                                              PsiDefor,
                                                              PosIni,
                                                              PsiIni,
                                                              XBELEM,
                                                              [index])
                    
                    fp.write("{:<14,e}".format(locForces[0,3+component]))
                else:
                    raise IOError("writeDict key not recognised.")
            # END for myStr
            fp.write("\n")
            fp.flush()
                    
        # 'Rollup' due to external velocities. TODO: Must add gusts here!
        ZetaStar[:,:] = ZetaStar[:,:] + VMINPUT.U_infty*dt
    
    # END Time loop
    
    # Write _dyn file.
    ofile = Settings.OutputDir + Settings.OutputFileRoot + '_SOL912_dyn.dat'
    fp = open(ofile,'w')
    BeamIO.Write_dyn_File(fp, Time, PosPsiTime)
    fp.close()
    
    #Write _shape file
    ofile = Settings.OutputDir + Settings.OutputFileRoot + '_SOL912_shape.dat'
    fp = open(ofile,'w')
    BeamIO.Write_shape_File(fp, len(Time), NumNodes_tot.value, Time, DynOut)
    fp.close()
    
    #Write rigid-body velocities
    ofile = Settings.OutputDir + Settings.OutputFileRoot + '_SOL912_rigid.dat'
    fp = open(ofile,'w')
    BeamIO.Write_rigid_File(fp, Time, Vrel, VrelDot)
    fp.close()
    
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


if __name__ == '__main__':
    # Beam options.
    XBOPTS = DerivedTypes.Xbopts(FollowerForce = ct.c_bool(False),
                                 MaxIterations = ct.c_int(50),
                                 PrintInfo = ct.c_bool(True),
                                 OutInaframe = ct.c_bool(True),
                                 NumLoadSteps = ct.c_int(1),
                                 Solution = ct.c_int(912),
                                 MinDelta = ct.c_double(1e-5),
                                 NewmarkDamp = ct.c_double(5e-3))
    # beam inputs.
    XBINPUT = DerivedTypes.Xbinput(2,4)
    XBINPUT.BeamLength = 6.096
    XBINPUT.BeamStiffness[0,0] = 1.0e+09
    XBINPUT.BeamStiffness[1,1] = 1.0e+09
    XBINPUT.BeamStiffness[2,2] = 1.0e+09
    XBINPUT.BeamStiffness[3,3] = 0.99e+06
    XBINPUT.BeamStiffness[4,4] = 9.77e+06
    XBINPUT.BeamStiffness[5,5] = 1.0e+09
    XBINPUT.BeamStiffness[:,:] = 1.0e0*XBINPUT.BeamStiffness[:,:]
    XBINPUT.BeamMass[0,0] = 35.71
    XBINPUT.BeamMass[1,1] = 35.71
    XBINPUT.BeamMass[2,2] = 35.71
    XBINPUT.BeamMass[3,3] = 8.64*10
    XBINPUT.BeamMass[4,4] = 1.0e+9
    XBINPUT.BeamMass[5,5] = 1.0e+9
    # Off diagonal terms (in Theodorsen sectional coordinates)
    ElasticAxis = -0.34
    InertialAxis = -7.0/50.0
    x_alpha = -(InertialAxis - ElasticAxis)
    # pitch-plunge coupling term (b-frame coordinates)
    c = 1.8288
    cgLoc = 0.5*c*np.array([0.0, x_alpha, 0.0])
    cgSkew = Skew(cgLoc)
#     mOff = x_alpha*(c/2)*XBINPUT.BeamMass[0,0]
#     XBINPUT.BeamMass[2,3] = -mOff
#     XBINPUT.BeamMass[0,5] = mOff
#     XBINPUT.BeamMass[3:,:3] = XBINPUT.BeamMass[:3,3:].T
    XBINPUT.BeamMass[:3,3:] = XBINPUT.BeamMass[0,0] * cgSkew
    XBINPUT.BeamMass[3:,:3] = XBINPUT.BeamMass[:3,3:].T
    XBINPUT.BeamMass[4,4] += XBINPUT.BeamMass[0,0]*pow(cgLoc[2],2.0)
    XBINPUT.BeamMass[5,5] += XBINPUT.BeamMass[0,0]*pow(cgLoc[1],2.0)
    
    # Get suggested panelling.
    Umag = 50.0
    AOA  = 5.0*np.pi/180.0
    M = 4
    delTime = c/(Umag*M)
    # Unsteady parameters.
    XBINPUT.dt = delTime
    XBINPUT.t0 = 0.0
    XBINPUT.tfin = 2
    
    # Set motion of wing.
    NumSteps = np.ceil( (XBINPUT.tfin + XBINPUT.dt - XBINPUT.t0) / XBINPUT.dt)
    XBINPUT.ForcedVel = np.zeros((NumSteps,6),ct.c_double,'F')
    for i in range(XBINPUT.ForcedVel.shape[0]):
        XBINPUT.ForcedVel[i,:] = [0.0, Umag*np.cos(AOA), -Umag*np.sin(AOA), 0.0, 0.0, 0.0]
#         XBINPUT.ForcedVel[i,:] = [0.0, Umag, 0.0, 0.0, 0.0, 0.0]
    XBINPUT.ForcedVelDot = np.zeros((NumSteps,6),ct.c_double,'F')
     
    # aero params.
    WakeLength = 30.0*c
    Mstar = int(WakeLength/(delTime*Umag))
    # aero options.
    N = XBINPUT.NumNodesTot - 1
    VMOPTS = DerivedTypesAero.VMopts(M = M,
                                     N = N,
                                     ImageMethod = False,
                                     Mstar = Mstar,
                                     Steady = False,
                                     KJMeth = False,
                                     NewAIC = True,
                                     DelTime = delTime,
                                     NumCores = 4)
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
    
    # Gust inputs
    gust = Gust(uMag = 0.0*Umag,
                h = 10.0,
                r = 0.0)
    
    VMINPUT = DerivedTypesAero.VMinput(c = c,
                                       b = XBINPUT.BeamLength,
                                       U_mag = 0.0,
                                       alpha = 0.0*np.pi/180.0,
                                       theta = 0.0,
                                       WakeLength = WakeLength,
                                       ctrlSurf = ctrlSurf,
                                       gust = None)
    
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
    
    Settings.OutputDir = Settings.SharPyProjectDir + "output/temp/"
    Settings.OutputFileRoot = "FlyingWing"
    
    # Solve nonlinear dynamic simulation.
    Solve_Py(XBINPUT, XBOPTS, VMOPTS, VMINPUT, AELAOPTS,
             writeDict =  writeDict)