'''@package PyFSI.Beam2UVLM.test_Coincident
@brief      print out deformed aero grid.
@author     Rob Simpson
@contact    r.simpson11@imperial.ac.uk 
@version    0.0
@date       21/01/2013
@pre        None
@warning    None
'''
import Beam2UVLM
from SectionInit import FlatPlateReg
import DerivedTypes
import NonlinearStatic
import numpy as np
import SharPySettings as Settings
import subprocess
import sys
import BeamIO
import BeamLib
import BeamInit
import ctypes as ct
import XbeamLib

def WriteAeroTecHeader(FileName='AeroGrid.dat', Title='Default',\
                       Variables=['X', 'Y', 'Z']):
    """@brief Opens Filename and writes header for tecplot ascii format.
    @return open file object for writing to."""
    
    FileObject = open(FileName,'w')
    FileObject.write('TITLE="%s"\n' % (Title))
    FileObject.write('VARIABLES=')
    for Var in range(len(Variables)):
        FileObject.write('"%s" ' % (Variables[Var]))
        if Var == len(Variables):
            FileObject.write('\n')
    #END for Var
    FileObject.write('\n')
    
    return FileObject


def WriteAeroTecZone(FileObject, AeroGrid, TimeStep=-1, NumTimeSteps=0,\
                     Time=0.0, Variables=['X', 'Y', 'Z'], Text=True):
    """@brief Writes aero grid data to ascii file in tecplot ij format.
    @param FileObject Must be an open file object.
    @param Timestep -1 for static solution and >-1 thereafter."""
    
    FileObject.write('ZONE I=%d J=%d DATAPACKING=BLOCK T="Timestep %d of %d"\nSOLUTIONTIME=%f STRANDID=1\n' % \
                     (AeroGrid.shape[0],AeroGrid.shape[1],TimeStep+1,\
                      NumTimeSteps,Time))
    
    for Var in range(len(Variables)):
        for j in range(AeroGrid.shape[1]):
            for i in range(AeroGrid.shape[0]):
                FileObject.write('%f\t' % (AeroGrid[i,j,Var]))
            #END for j
        #END for i
        FileObject.write('\n')
    #END for Var
    
    if Text==True and TimeStep > -1:
        FileObject.write('TEXT T="time = %fs" CS=GRID3D AN=Headleft S=LOCAL ZN=%d X=0.0 Y=0.0 Z=-1.0\n' % (Time,TimeStep+1))
    
    return FileObject
    

def CloseAeroTecFile(FileObject):
    """TODO: could whole file writing process be encapsulated in a class?"""
    FileObject.close()
        
    
def StaticDef():
    """Set up Xbopts for nonlinear static analysis defined in 
    NonlinearStatic/testcases.pdf case 1.1"""
    XBOPTS = DerivedTypes.Xbopts()
    XBOPTS.Solution.value = 112 
    XBOPTS.NumLoadSteps.value = 10
    XBOPTS.MinDelta.value = 1e-05
    XBOPTS.FollowerForce.value = False
         
    """Set up Xbinput for nonlinear static analysis defined in 
    NonlinearStatic/testcases.pdf case 1.1 with extra out-of plane loads."""
    XBINPUT = DerivedTypes.Xbinput(2,20)
    XBINPUT.BeamLength = 5.0
    XBINPUT.BeamStiffness[0,0] = 4.8e+08
    XBINPUT.BeamStiffness[1,1] = 3.231e+08
    XBINPUT.BeamStiffness[2,2] = 3.231e+08
    XBINPUT.BeamStiffness[3,3] = 1.0e+06
    XBINPUT.BeamStiffness[4,4] = 9.346e+06
    XBINPUT.BeamStiffness[5,5] = 9.346e+06
    XBINPUT.BeamMass[0,0] = 100
    XBINPUT.BeamMass[1,1] = 100
    XBINPUT.BeamMass[2,2] = 100
    XBINPUT.BeamMass[3,3] = 10
    XBINPUT.BeamMass[4,4] = 0.0 #Neglect the cross-section bending inertia
    XBINPUT.BeamMass[5,5] = 0.0 #Neglect the cross-section bending inertia
    XBINPUT.ForceStatic[-1,2] = 6e+05
    XBINPUT.ForceStatic[-1,1] = 6e+05
    XBINPUT.ForceStatic[-1,0] = 6e+03
    
    PosDefor, PsiDefor = NonlinearStatic.Solve_Py(XBINPUT, XBOPTS)
    
    "Initialise PyFSI variables"
    M = 10
    LeadingEdge = np.array([0.0, 0.5, 0.0])
    TrailingEdge = np.array([0.0, -0.5, 0.0])
    Section = FlatPlateReg(M,LeadingEdge,TrailingEdge)
    
    VelA_A = np.zeros(3)
    OmegaA_A = np.zeros(3)
    PosDotDef = np.zeros_like(PosDefor)
    PsiDotDef = np.zeros_like(PsiDefor)
    
    AeroGrid, AeroVels = Beam2UVLM.CoincidentGrid(PosDefor, PsiDefor, Section,\
                   VelA_A, OmegaA_A, PosDotDef, PsiDotDef,
                   XBINPUT)
    
    FileObject = WriteAeroTecHeader()
    FileObject = WriteAeroTecZone(FileObject, AeroGrid)
    CloseAeroTecFile(FileObject)
    
    
def DynamicAnim():
    """Set up Xbopts for nonlinear static analysis defined in 
    NonlinearStatic/testcases.pdf case 1.2."""
    XBOPTS = DerivedTypes.Xbopts()
    XBOPTS.Solution.value = 312 
    XBOPTS.NumLoadSteps.value = 10
    XBOPTS.MinDelta.value = 1e-05
    XBOPTS.FollowerForce.value = False
    XBOPTS.NumGauss.value = 1
           
    """Set up Xbinput for nonlinear static analysis defined in 
    NonlinearStatic/testcases.pdf case 1.2."""
    XBINPUT = DerivedTypes.Xbinput(2,20)
    XBINPUT.BeamLength = 5.0
    XBINPUT.BeamStiffness[0,0] = 4.8e+08
    XBINPUT.BeamStiffness[1,1] = 3.231e+08
    XBINPUT.BeamStiffness[2,2] = 3.231e+08
    XBINPUT.BeamStiffness[3,3] = 1.0e+06
    XBINPUT.BeamStiffness[4,4] = 9.346e+06
    XBINPUT.BeamStiffness[5,5] = 9.346e+06
    XBINPUT.BeamMass[0,0] = 100
    XBINPUT.BeamMass[1,1] = 100
    XBINPUT.BeamMass[2,2] = 100
    XBINPUT.BeamMass[3,3] = 10
    XBINPUT.BeamMass[4,4] = 0.0001
    XBINPUT.BeamMass[5,5] = 0.0001
    
    "Dynamic parameters"
    XBINPUT.t0 = 0.0
    XBINPUT.tfin = 1.0
    XBINPUT.dt = 0.001
    XBINPUT.Omega = 40.0
    XBINPUT.ForceDyn[-1,2] = 180e+03
    XBINPUT.ForceDyn[-1,1] = 240e+03
    XBINPUT.ForceDyn[-1,3] = 4800.0
    XBINPUT.ForcingType = 'RampSin'
    XBINPUT.RampTime = 1.0
    
    """Nonlinear dynamic structural solver in Python. Assembly of matrices 
    is carried out with Fortran subroutines."""
    
    "Check correct solution code"
    assert XBOPTS.Solution.value == 312, ('NonlinearDynamic requested' +\
                                              ' with wrong solution code')
    
    "Initialise beam"
    XBINPUT, XBOPTS, NumNodes_tot, XBELEM, PosIni, PsiIni, XBNODE, NumDof \
                = BeamInit.Static(XBINPUT,XBOPTS)
    
    
    "Set initial conditions as undef config"
    PosDefor = PosIni.copy(order='F')
    PsiDefor = PsiIni.copy(order='F')
    
    
    "Solve static"
    BeamLib.Cbeam3_Solv_NonlinearStatic(XBINPUT, XBOPTS, NumNodes_tot, XBELEM,\
                            PosIni, PsiIni, XBNODE, NumDof,\
                            PosDefor, PsiDefor)
    
    
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
    
    "Initialise variables for dynamic analysis"
    Time, NumSteps, ForceTime, ForcedVel, ForcedVelDot,\
    PosDotDef, PsiDotDef,\
    OutGrids, PosPsiTime, VelocTime, DynOut\
        = BeamInit.Dynamic(XBINPUT,XBOPTS)
        
        
    "Calculate Aerodynamic grid"
    "Init Variables"
    M = 10
    LeadingEdge = np.array([0.0, 0.5, 0.0])
    TrailingEdge = np.array([0.0, -0.5, 0.0])
    Section = FlatPlateReg(M,LeadingEdge,TrailingEdge)
    
    AeroGrid, AeroVels = Beam2UVLM.CoincidentGrid(PosDefor, PsiDefor, Section,\
                   ForcedVel[0,0:3], ForcedVel[0,3:], PosDotDef, PsiDotDef,
                   XBINPUT)
    
    "Create output and plot initial deformed shape"
    FileObject = WriteAeroTecHeader('AeroGridAnim.tec','Aerodynamic Grid')
    FileObject = WriteAeroTecZone(FileObject, AeroGrid)
    
    
    "Write _force file"
    ofile = Settings.OutputDir + Settings.OutputFileRoot + '_SOL312_force.dat'
    fp = open(ofile,'w')
    BeamIO.Write_force_File(fp, Time, ForceTime, ForcedVel, ForcedVelDot)
    fp.close() 
    
    
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
    
    "Force at the first time-step"
    Force = (XBINPUT.ForceStatic + XBINPUT.ForceDyn*ForceTime[0]).copy('F')
    

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
    
    
    "Time loop"
    for iStep in range(NumSteps.value):
        
        if XBOPTS.PrintInfo.value==True:
            sys.stdout.write('  Time: %-10.4e\n' %(Time[iStep+1]))
            sys.stdout.write('   SubIter DeltaF     DeltaX     ResLog10\n')
        
        
        "calculate dt"
        dt = Time[iStep+1] - Time[iStep]
        
        
        "Update transformation matrix for given angular velocity"
        Temp = np.linalg.inv(Unit4 + 0.25*XbeamLib.QuadSkew(ForcedVel[iStep+1,3:])*dt)
        Quat = np.dot(Temp, np.dot(Unit4 - 0.25*XbeamLib.QuadSkew(ForcedVel[iStep,3:])*dt, Quat))
        Quat = Quat/np.linalg.norm(Quat)
        Cao  = XbeamLib.Rot(Quat)
        
        
        "Predictor step"
        X += dt*dXdt + (0.5-beta)*dXddt*dt**2
        dXdt += (1.0-gamma)*dXddt*dt
        dXddt[:] = 0.0
        
        
        "Force at current time-step"
        Force = (XBINPUT.ForceStatic + \
                 XBINPUT.ForceDyn*ForceTime[iStep+1]).copy('F')
        
    
        "Reset convergence parameters"
        Iter = 0
        ResLog10 = 0.0
        
        
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
                
        "END Netwon-Raphson"
        
        
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
        
        "Create aero grid"
        AeroGrid, AeroVels = Beam2UVLM.CoincidentGrid(PosDefor, PsiDefor, Section,\
                   ForcedVel[iStep,0:3], ForcedVel[iStep,3:], PosDotDef, PsiDotDef,
                   XBINPUT)
    
        "Plot deformed shape"
        FileObject = WriteAeroTecZone(FileObject, AeroGrid,\
                                       iStep, NumSteps.value, Time[iStep+1])
    
    "END Time loop"
    
    "Write _dyn file"
    ofile = Settings.OutputDir + Settings.OutputFileRoot + '_SOL312_dyn.dat'
    fp = open(ofile,'w')
    BeamIO.Write_dyn_File(fp, Time, PosPsiTime)
    fp.close()
    
    "Write _shape file"
    ofile = Settings.OutputDir + Settings.OutputFileRoot + '_SOL312_shape.dat'
    fp = open(ofile,'w')
    BeamIO.Write_shape_File(fp, len(Time), NumNodes_tot.value, Time, DynOut)
    fp.close()
    
    "Close Tecplot ascii FileObject"
    CloseAeroTecFile(FileObject)
    
    if XBOPTS.PrintInfo.value==True:
        sys.stdout.write(' ... done\n')

if __name__ == '__main__':
    "Print Aero grid after static solution"
    StaticDef()
    "Create Dynamic Animation"
    DynamicAnim()
    
    