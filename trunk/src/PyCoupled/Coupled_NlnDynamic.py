'''@package PyCoupled.Coupled_NlnDynamic
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
import PyCoupled.Coupled_NlnStatic as Static
import PyBeam.Utils.XbeamLib as XbeamLib
from PyCoupled.Coupled_NlnStatic import AddGravityLoads

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
    
    def __init__(self, VMOPTS, VMINPUT,\
                 DelS, NumChordLengths,\
                 VelA_G, OmegaA_G,\
                 OriginA_G = np.zeros((3),ct.c_double),\
                 PsiA_G = np.zeros((3),ct.c_double)):
        
        self.NumChordLengths = NumChordLengths
        self.VelA_G = VelA_G
        self.OmegaA_G = OmegaA_G
        self.OriginA_G = OriginA_G
        self.PsiA_G = PsiA_G

def Solve_Py(XBINPUT,XBOPTS,VMOPTS,VMINPUT,VMUNST,AELAOPTS):
    """Nonlinear dynamic solver using Python to solve aeroelastic
    equation. Assembly of structural matrices is carried out with 
    Fortran subroutines. Aerodynamics solved using PyAero\.UVLM."""
        
    """BEAM INIT ------------------------------------------------------------"""
    
    "Check correct solution code"
    assert XBOPTS.Solution.value == 312, ('NonlinearDynamic requested' +\
                                              ' with wrong solution code')
    
    "Initialise beam"
    XBINPUT, XBOPTS, NumNodes_tot, XBELEM, PosIni, PsiIni, XBNODE, NumDof \
                = BeamInit.Static(XBINPUT,XBOPTS)
    
    
    "Solve NlnStatic (PosDefor,PsiDefor,Gamma)"
    if AELAOPTS.ImpStart == False:
        "modify options"
        XBOPTS.Solution.value = 112
        VMOPTS.Steady = ct.c_bool(True)
        Rollup = VMOPTS.Rollup.value
        VMOPTS.Rollup.value = False
        
        "solve"
        PosDefor, PsiDefor, Zeta, ZetaStar, Gamma, GammaStar, Force = \
                    Static.Solve_Py(XBINPUT, XBOPTS, VMOPTS, VMINPUT, AELAOPTS)
        
        "reset options"
        XBOPTS.Solution.value = 312
        VMOPTS.Steady = ct.c_bool(False)
        VMOPTS.Rollup.value = Rollup
        
    elif AELAOPTS.ImpStart == True:
        
        "PosDefor,PsiDefor"
        PosDefor = PosIni.copy(order='F')
        PsiDefor = PsiIni.copy(order='F')
        
        "Force = 0.0"
        Force = np.zeros((XBINPUT.NumNodesTot,6),ct.c_double,'F')
    
    #END if-elif ImpStart
    
    
    "Write deformed configuration to file"
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
    
    BeamIO.OutputElems(XBINPUT.NumElems, NumNodes_tot.value, XBELEM, \
                       PosDefor, PsiDefor, ofile, WriteMode)
    
    
    "Initialise structural variables for dynamic analysis"
    Time, NumSteps, ForceTime, ForcedVel, ForcedVelDot,\
    PosDotDef, PsiDotDef,\
    OutGrids, PosPsiTime, VelocTime, DynOut\
        = BeamInit.Dynamic(XBINPUT,XBOPTS)
    
    
    "Write _force file"
#    ofile = Settings.OutputDir + Settings.OutputFileRoot + '_SOL312_force.dat'
#    fp = open(ofile,'w')
#    BeamIO.Write_force_File(fp, Time, ForceTime, ForcedVel, ForcedVelDot)
#    fp.close() 
    
    
    "Write _vel file"   
    #TODO: write _vel file
    
    
    "Write .mrb file"
    #TODO: write .mrb file
    
    
    if XBOPTS.PrintInfo.value==True:
        sys.stdout.write('Solve nonlinear dynamic case in Python ... \n')
    
    
    "Initialise structural system tensors"
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
    
    X0    = np.zeros(NumDof.value, ct.c_double, 'F')
    X     = np.zeros(NumDof.value, ct.c_double, 'F')
    DX    = np.zeros(NumDof.value, ct.c_double, 'F')
    dXdt  = np.zeros(NumDof.value, ct.c_double, 'F')
    dXddt = np.zeros(NumDof.value, ct.c_double, 'F')
    Force_Dof = np.zeros(NumDof.value, ct.c_double, 'F')
    
    Qglobal = np.zeros(NumDof.value, ct.c_double, 'F')
    
    
    "Initialise rotation operators"
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
    
    
    "Extract initial displacements and velocities"
    BeamLib.Cbeam_Solv_Disp2State(NumNodes_tot, NumDof, XBINPUT, XBNODE,\
                          PosDefor, PsiDefor, PosDotDef, PsiDotDef,
                          X, dXdt)
    
    
    "Approximate initial accelerations"
    
    "Initialise accelerations as zero arrays"
    PosDotDotDef = np.zeros((NumNodes_tot.value,3),ct.c_double,'F')
    PsiDotDotDef = np.zeros((XBINPUT.NumElems,Settings.MaxElNod,3),\
                           ct.c_double, 'F')
    

    "Assemble matrices for dynamic analysis"
    BeamLib.Cbeam3_Asbly_Dynamic(XBINPUT, NumNodes_tot, XBELEM, XBNODE,\
                         PosIni, PsiIni, PosDefor, PsiDefor,\
                         PosDotDef, PsiDotDef, PosDotDotDef, PsiDotDotDef,\
                         Force, ForcedVel[0,:], ForcedVelDot[0,:],\
                         NumDof, Settings.DimMat,\
                         ms, MglobalFull, Mvel,\
                         cs, CglobalFull, Cvel,\
                         ks, KglobalFull, fs, FglobalFull,\
                         Qglobal, XBOPTS, Cao)
    
       
    "Get force vector for unconstrained nodes (Force_Dof)"
    BeamLib.f_fem_m2v(ct.byref(NumNodes_tot),\
                              ct.byref(ct.c_int(6)),\
                              Force.ctypes.data_as(ct.POINTER(ct.c_double)),\
                              ct.byref(NumDof),\
                              Force_Dof.ctypes.data_as(ct.POINTER(ct.c_double)),\
                              XBNODE.Vdof.ctypes.data_as(ct.POINTER(ct.c_int)) )
    
    
    "Get RHS at initial condition"
    Qglobal += -np.dot(FglobalFull, Force_Dof)
    
    
    "Initial Accel"
    dXddt[:] = np.dot(np.linalg.inv(MglobalFull), -Qglobal)
    
    
    "Record position of all grid points in global FoR at initial time step"
    DynOut[0:NumNodes_tot.value,:] = PosDefor
    
    "Position/rotation of the selected node in initial deformed configuration"
    PosPsiTime[0,:3] = PosDefor[-1,:]
    PosPsiTime[0,3:] = PsiDefor[-1,XBELEM.NumNodes[-1]-1,:]
    
    
    "Get gamma and beta for Newmark scheme"
    gamma = 0.5 + XBOPTS.NewmarkDamp.value
    beta = 0.25*(gamma + 0.5)**2
    
    
    """ Initialise Aero -----------------------------------------------------"""
    
    "set AoA to zero Foo"
    if AELAOPTS.ImpStart == False:
        VMINPUT.U_infty = np.array([0.0, -np.linalg.norm(VMINPUT.U_infty), 0.0])
        
#    "set density to zero for free vibration"
#    AELAOPTS.AirDensity = 0.0
    
    "define section used to generate aero mesh"       
    Section = InitSection(VMOPTS,VMINPUT,AELAOPTS.ElasticAxis)
 
    "declare memory for Aero grid velocities"
    ZetaDot = np.zeros((Section.shape[0],PosDefor.shape[0],3),ct.c_double,'C')
    
    "Additional Aero solver variables"
    AeroForces = np.zeros((VMOPTS.M.value+1,VMOPTS.N.value+1,3),ct.c_double,'C')
    AIC = np.zeros((VMOPTS.M.value*VMOPTS.N.value,VMOPTS.M.value*VMOPTS.N.value), \
                   ct.c_double,'C')
    BIC = np.zeros((VMOPTS.M.value*VMOPTS.N.value,VMOPTS.M.value*VMOPTS.N.value), \
                   ct.c_double,'C')
    
    "init external velocities"  
    Uext = InitSteadyExternalVels(VMOPTS,VMINPUT)
    
    
    if AELAOPTS.ImpStart == True:
        "init required variables"
        Zeta = np.zeros((Section.shape[0],PosDefor.shape[0],3),ct.c_double,'C')             
        Gamma = np.zeros((VMOPTS.M.value,VMOPTS.N.value),ct.c_double,'C')
        
        "generate surface, wake and gamma matrix"
        CoincidentGrid(PosDefor, PsiDefor, Section, ForcedVel[0,:3], \
                           ForcedVel[0,3:], PosDotDef, PsiDotDef, XBINPUT,\
                           Zeta, ZetaDot, VMUNST.OriginA_G,VMUNST.PsiA_G)
        
        "init wake"
        ZetaStar, GammaStar = InitSteadyWake(VMOPTS,VMINPUT,Zeta)    
    #END if ImpStart
    
    if Settings.PlotTec==True:
        "define tecplot stuff"
        FileName = Settings.OutputDir + Settings.OutputFileRoot + 'AeroGrid.dat'
        Variables = ['X', 'Y', 'Z','Gamma']        
        FileObject = PostProcess.WriteAeroTecHeader(FileName, \
                                                    'Default',\
                                                    Variables)
        
        "plot results of static analysis"
        "surface"
        PostProcess.WriteAeroTecZone(FileObject, 'Surface', Zeta,\
                                  0,\
                                  XBOPTS.NumLoadSteps.value,\
                                  Time[0+1],\
                                  Variables, True, Gamma)
        "wake"
        PostProcess.WriteAeroTecZone(FileObject, 'Wake', ZetaStar,\
                                  0,\
                                  XBOPTS.NumLoadSteps.value,\
                                  Time[0+1],\
                                  Variables, False, GammaStar)
    #END if plot tec
    
    
    """ begin analysis ----------------------------------------------------- """
    
    "Time loop"
    for iStep in range(NumSteps.value):
        
        if XBOPTS.PrintInfo.value==True:
            sys.stdout.write('Time: %-10.4e\n' %(Time[iStep+1]))
            sys.stdout.write('   SubIter DeltaF     DeltaX     ResLog10\n')
        
        "calculate dt"
        dt = Time[iStep+1] - Time[iStep]
        
        "set for aero force calcs"
        VMOPTS.DelTime = ct.c_double(dt)
        
        
        "Update transformation matrix for given angular velocity"
        Temp = np.linalg.inv(Unit4 + 0.25*XbeamLib.QuadSkew(ForcedVel[iStep+1,3:])*dt)
        Quat = np.dot(Temp, np.dot(Unit4 - 0.25*XbeamLib.QuadSkew(ForcedVel[iStep,3:])*dt, Quat))
        Quat = Quat/np.linalg.norm(Quat)
        Cao  = XbeamLib.Rot(Quat)
        
        "Force at current time-step"
        if iStep > 0:
            
            "zero aero forces"
            AeroForces[:,:,:] = 0.0
            
            
            "update geometry"
            VMUNST.OriginA_G[:] += ForcedVel[iStep,:3]*dt
            
            CoincidentGrid(PosDefor, PsiDefor, Section, ForcedVel[iStep+1,:3], \
                           ForcedVel[iStep+1,3:], PosDotDef, PsiDotDef, XBINPUT,\
                           Zeta, ZetaDot, VMUNST.OriginA_G,VMUNST.PsiA_G)
            
            
            "update wake geom"
                   
            "'roll' data"
            ZetaStar = np.roll(ZetaStar,1,axis = 0)
            GammaStar = np.roll(GammaStar,1,axis = 0)
            
            "overwrite grid points with TE"
            ZetaStar[0,:] = Zeta[VMOPTS.M.value,:]
            
            "overwrite Gamma with TE value from previous timestep"
            GammaStar[0,:] = Gamma[VMOPTS.M.value-1,:]
            
            
            "solve for AeroForces"
            UVLMLib.Cpp_Solver_VLM(Zeta, ZetaDot, Uext, ZetaStar, VMOPTS, \
                           AeroForces, Gamma, GammaStar, AIC, BIC)
            
            
            "apply density scaling"
            AeroForces[:,:,:] = AELAOPTS.AirDensity*AeroForces[:,:,:]
            
            
            "write to tecplot file"
            if Settings.PlotTec==True:
                
                "surface"
                PostProcess.WriteAeroTecZone(FileObject, 'Surface', Zeta,\
                                  iStep,\
                                  XBOPTS.NumLoadSteps.value,\
                                  Time[iStep+1],\
                                  Variables, True, Gamma)
                "wake"
                PostProcess.WriteAeroTecZone(FileObject, 'Wake', ZetaStar,\
                                  iStep,\
                                  XBOPTS.NumLoadSteps.value,\
                                  Time[iStep+1],\
                                  Variables, False, GammaStar)
            # END if PlotTec 

            
            "map AeroForces to beam"
            CoincidentGridForce(XBINPUT, PsiDefor, Section, AeroForces,\
                        Force)
            
            
            "add gravity loads"
            AddGravityLoads(Force,XBINPUT,XBELEM,AELAOPTS.gForce)
        
        #END if iStep > 0
            
            
        
        
        "Predictor step"
        X += dt*dXdt + (0.5-beta)*dXddt*dt**2
        dXdt += (1.0-gamma)*dXddt*dt
        dXddt[:] = 0.0
        
    
        "Reset convergence parameters"
        Iter = 0
        ResLog10 = 0.0
            
        
        "save Gamma for tightly coupled"
        if AELAOPTS.Tight == True:
            GammaSav = Gamma.copy(order = 'C')
        
        "Newton-Raphson loop"        
        while ( (ResLog10 > np.log10(XBOPTS.MinDelta.value)) \
                & (Iter < XBOPTS.MaxIterations.value) ):
            
            "set tensors to zero"
            Qglobal[:] = 0.0 
            Mvel[:,:] = 0.0
            Cvel[:,:] = 0.0
            MglobalFull[:,:] = 0.0
            CglobalFull[:,:] = 0.0
            KglobalFull[:,:] = 0.0
            FglobalFull[:,:] = 0.0
            
            "Update counter"
            Iter += 1
            
            if XBOPTS.PrintInfo.value==True:
                sys.stdout.write('   %-7d ' %(Iter))
                
            
            "nodal diplacements and velocities from state vector"
            BeamLib.Cbeam3_Solv_State2Disp(XBINPUT, NumNodes_tot, XBELEM, XBNODE,
                           PosIni, PsiIni, NumDof, X, dXdt,\
                           PosDefor, PsiDefor, PosDotDef, PsiDotDef)
            
            "if tightly coupled is on then get new aeroforces"
            if AELAOPTS.Tight == True:
                "zero aero forces"
                AeroForces[:,:,:] = 0.0
                
                "set gamma at t-1 to saved solution"
                Gamma[:,:] = GammaSav[:,:]
                
                "get new grid"
                CoincidentGrid(PosDefor, PsiDefor, Section, ForcedVel[iStep+1,:3], \
                           ForcedVel[iStep+1,3:], PosDotDef, PsiDotDef, XBINPUT,\
                           Zeta, ZetaDot, VMUNST.OriginA_G,VMUNST.PsiA_G)
                
                "close wake"
                ZetaStar[0,:] = Zeta[VMOPTS.M.value,:]
                
                "save data and turn off rollup"
                Rollup = VMOPTS.Rollup.value
                VMOPTS.Rollup.value = False
                
                "solve for AeroForces"
                UVLMLib.Cpp_Solver_VLM(Zeta, ZetaDot, Uext, ZetaStar, VMOPTS, \
                               AeroForces, Gamma, GammaStar, AIC, BIC)
                
                "turn rollup back to original"
                VMOPTS.Rollup.value = Rollup
                
                "apply density scaling"
                AeroForces[:,:,:] = AELAOPTS.AirDensity*AeroForces[:,:,:]
                
                "beam forces"
                CoincidentGridForce(XBINPUT, PsiDefor, Section, AeroForces,\
                            Force)
                
                "add gravity loads"
                AddGravityLoads(Force,XBINPUT,XBELEM,AELAOPTS.gForce)
            
            #END if Tight
            
            
            "update matrices"
            BeamLib.Cbeam3_Asbly_Dynamic(XBINPUT, NumNodes_tot, XBELEM, XBNODE,\
                         PosIni, PsiIni, PosDefor, PsiDefor,\
                         PosDotDef, PsiDotDef, PosDotDotDef, PsiDotDotDef,\
                         Force, ForcedVel[iStep+1,:], ForcedVelDot[iStep+1,:],\
                         NumDof, Settings.DimMat,\
                         ms, MglobalFull, Mvel,\
                         cs, CglobalFull, Cvel,\
                         ks, KglobalFull, fs, FglobalFull,\
                         Qglobal, XBOPTS, Cao)
            
            
            "Get force vector for unconstrained nodes (Force_Dof)"
            BeamLib.f_fem_m2v(ct.byref(NumNodes_tot),\
                              ct.byref(ct.c_int(6)),\
                              Force.ctypes.data_as(ct.POINTER(ct.c_double)),\
                              ct.byref(NumDof),\
                              Force_Dof.ctypes.data_as(ct.POINTER(ct.c_double)),\
                              XBNODE.Vdof.ctypes.data_as(ct.POINTER(ct.c_int)) )
            
            
            "Solve for update vector"
            "Residual"
            Qglobal += np.dot(MglobalFull, dXddt) \
                        + np.dot(Mvel,ForcedVelDot[iStep+1,:]) \
                     - np.dot(FglobalFull, Force_Dof)
                              
            
             
            if XBOPTS.PrintInfo.value==True:                 
                sys.stdout.write('%-10.4e ' %(max(abs(Qglobal))))
            
            
            "Calculate system matrix for update calculation"
            Asys = KglobalFull + \
                      CglobalFull*gamma/(beta*dt) + \
                      MglobalFull/(beta*dt**2)
                      
            
            "Solve for update"
            DX[:] = np.dot(np.linalg.inv(Asys), -Qglobal)
            
            
            "Corrector step"
            X += DX
            dXdt += DX*gamma/(beta*dt)
            dXddt += DX/(beta*dt**2)
            
            
            "Residual at first iteration"
            if(Iter == 1):
                Res0_Qglobal = max(abs(Qglobal)) + 1.e-16
                Res0_DeltaX  = max(abs(DX)) + 1.e-16
            
            "Update residual and compute log10"
            Res_Qglobal = max(abs(Qglobal)) + 1.e-16
            Res_DeltaX  = max(abs(DX)) + 1.e-16
            
            ResLog10 = max([np.log10(Res_Qglobal/Res0_Qglobal),\
                            np.log10(Res_DeltaX/Res0_DeltaX)])
            
            if XBOPTS.PrintInfo.value==True:
                sys.stdout.write('%-10.4e %8.4f\n' %(max(abs(DX)),ResLog10))
            
            if ResLog10 > 1.0:
                print("Residual growing! Exit Newton-Raphson...")
                break
                
        "END Netwon-Raphson"
        
        if ResLog10 > 1.0:
                print("Residual growing! Exit time-loop...")
                break
        
        "update to converged nodal displacements and velocities"
        BeamLib.Cbeam3_Solv_State2Disp(XBINPUT, NumNodes_tot, XBELEM, XBNODE,
                           PosIni, PsiIni, NumDof, X, dXdt,\
                           PosDefor, PsiDefor, PosDotDef, PsiDotDef)
        
        
        PosPsiTime[iStep+1,:3] = PosDefor[-1,:]
        PosPsiTime[iStep+1,3:] = PsiDefor[-1,XBELEM.NumNodes[-1]-1,:]
        
        "Position of all grid points in global FoR"
        i1 = (iStep+1)*NumNodes_tot.value
        i2 = (iStep+2)*NumNodes_tot.value
        DynOut[i1:i2,:] = PosDefor
        
        "Rollup due to external velocities"
        ZetaStar[:,:] += VMINPUT.U_infty*dt
    
    # END Time loop
    
    "Write _dyn file"
    ofile = Settings.OutputDir + Settings.OutputFileRoot + '_SOL312_dyn.dat'
    fp = open(ofile,'w')
    BeamIO.Write_dyn_File(fp, Time, PosPsiTime)
    fp.close()
    
#    "Write _shape file"
#    ofile = Settings.OutputDir + Settings.OutputFileRoot + '_SOL312_shape.dat'
#    fp = open(ofile,'w')
#    BeamIO.Write_shape_File(fp, len(Time), NumNodes_tot.value, Time, DynOut)
#    fp.close()
    
    
    "Close Tecplot ascii FileObject"
    if Settings.PlotTec==True:
        PostProcess.CloseAeroTecFile(FileObject)
    
    if XBOPTS.PrintInfo.value==True:
        sys.stdout.write(' ... done\n')

if __name__ == '__main__':
    "solve nlnstatic problem"
    
    "beam options"
    XBOPTS = DerivedTypes.Xbopts(FollowerForce = ct.c_bool(False),\
                                 MaxIterations = ct.c_int(50),\
                                 PrintInfo = ct.c_bool(True),\
                                 NumLoadSteps = ct.c_int(25),\
                                 Solution = ct.c_int(312),\
                                 MinDelta = ct.c_double(1e-5))
    
    "beam inputs"
    XBINPUT = DerivedTypes.Xbinput(3,6)
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
    
    "off diagonal terms"
    "in Theodorsen sectional coordinates"
    ElasticAxis = -0.34
    InertialAxis = -7.0/50.0
    x_alpha = InertialAxis - ElasticAxis
    
    "pitch-plunge coupling term"
    "b-frame coordinates"
    c = 1.8288
    mOff = x_alpha*(c/2)*XBINPUT.BeamMass[0,0]
    
    XBINPUT.BeamMass[2,3] = -mOff
    XBINPUT.BeamMass[0,5] = mOff
    XBINPUT.BeamMass[3:,:3] = XBINPUT.BeamMass[:3,3:].T
    
    
    "unsteady parameters"
    OmegaF = 70.0 #approx flutter freq
    U_mag = 170.0
    k = OmegaF*c/(2*U_mag) #get expected reduced freq
    M = int(50.0*k/np.pi) #get adequate panelling
    DelTime = c/(U_mag*M) #get resulting DelTime
    
    "beam"
    XBINPUT.dt = DelTime
    XBINPUT.t0 = 0.0
    XBINPUT.tfin = 0.1
    
    "aero"
    WakeLength = 30.0*c
    Mstar = int(WakeLength/(DelTime*U_mag))
    
    
    "aero options"
    VMOPTS = DerivedTypesAero.VMopts(M = M,\
                                     N = XBINPUT.NumNodesTot - 1 ,\
                                     ImageMethod = True,\
                                     Mstar = Mstar,\
                                     Steady = False,\
                                     KJMeth = False,\
                                     NewAIC = True,\
                                     DelTime = DelTime)
    
    "aero inputs"
    VMINPUT = DerivedTypesAero.VMinput(c = c, b = XBINPUT.BeamLength,\
                                       U_mag = U_mag,\
                                       alpha = 0.05*np.pi/180.0,\
                                       theta = 0.0,\
                                       ZetaDotTest = 0.0,\
                                       WakeLength = WakeLength)
    
    "unsteady aero inputs"
    VelA_G = np.array(([0.0,0.0,0.0]))
    OmegaA_G = np.array(([0.0,0.0,0.0]))
    VMUNST = VMCoupledUnstInput(VMOPTS, VMINPUT, 123.4, WakeLength/c,\
                              VelA_G, OmegaA_G)
    
    
    
    "aeroelastic opts"
    "Density and gravity"
    AELAOPTS = AeroelasticOps(ElasticAxis,0.0,1.02)
    
    
    Solve_Py(XBINPUT,XBOPTS,VMOPTS,VMINPUT,VMUNST,AELAOPTS)