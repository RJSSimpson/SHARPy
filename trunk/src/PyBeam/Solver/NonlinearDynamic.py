'''@package PyBeam.Solver.NonlinearDynamic
@brief      Nonlinear dynamic solvers.
@author     Rob Simpson
@contact    r.simpson11@imperial.ac.uk
@version    0.0
@date       10/12/2012
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

def Solve_F90(XBINPUT,XBOPTS):
    """@brief Nonlinear dynamic structural solver using f90 solve routine."""
    
    "Check correct solution code"
    assert XBOPTS.Solution.value == 312, ('NonlinearDynamic (F90) requested' +\
                                              ' with wrong solution code')
    
    
    "Initialise variables for static analysis"
    XBINPUT, XBOPTS, NumNodes_tot, XBELEM, PosIni, PsiIni, XBNODE, NumDof \
                = BeamInit.Static(XBINPUT, XBOPTS)
    
    
    "Change solution code to NonlinearStatic"
    XBOPTS.Solution.value = 112
    
    
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
    
    
    "Change solution code to NonlinearDynamic"
    XBOPTS.Solution.value = 312
    
    
    "Initialise variables for dynamic analysis"
    Time, NumSteps, ForceTime, ForcedVel, ForcedVelDot,\
    PosDotDef, PsiDotDef,\
    OutGrids, PosPsiTime, VelocTime, DynOut\
        = BeamInit.Dynamic(XBINPUT,XBOPTS)
    
    """TODO: move time vector init to input"""
    
    "Write _force file"
    ofile = Settings.OutputDir + Settings.OutputFileRoot + '_SOL312_force.dat'
    fp = open(ofile,'w')
    BeamIO.Write_force_File(fp, Time, ForceTime, ForcedVel, ForcedVelDot)
    fp.close() 
    
    
    "Write _vel file"   
    #TODO: write _vel file
    
    
    "Write .mrb file"
    #TODO: write .mrb file
    
    "Solve dynamic using f90 solver"
    BeamLib.Cbeam3_Solv_NonlinearDynamic(XBINPUT, XBOPTS, NumNodes_tot, XBELEM,\
            PosIni, PsiIni, XBNODE, NumDof,\
            PosDefor, PsiDefor,\
            NumSteps, Time, ForceTime, ForcedVel, ForcedVelDot,\
            PosDotDef, PsiDotDef,\
            PosPsiTime, VelocTime, DynOut, OutGrids)
    
    
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
    

def Solve_Py(XBINPUT,XBOPTS):
    """Nonlinear dynamic structural solver using python to solve residual
    equation. Assembly of matrices is carried out with fortran subroutines."""
    
    "Check correct solution code"
    assert XBOPTS.Solution.value == 312, ('NonlinearDynamic requested' +\
                                              ' with wrong solution code')
    
    "Initialise beam"
    XBINPUT, XBOPTS, NumNodes_tot, XBELEM, PosIni, PsiIni, XBNODE, NumDof \
                = BeamInit.Static(XBINPUT,XBOPTS)
    
    
    "Set initial conditions as undef config"
    PosDefor = PosIni.copy(order='F')
    PsiDefor = PsiIni.copy(order='F')
    
    "TODO: Solve static"
    
    
    if XBOPTS.PrintInfo.value==True:
        sys.stdout.write('Solve nonlinear dynamic case in Python ... \n')
    
    
    "Initialise structural system tensors"
    MglobalFull = np.zeros((NumDof,NumDof), ct.c_double, 'F'); ms = ct.c_int()
    CglobalFull = np.zeros((NumDof,NumDof), ct.c_double, 'F'); cs = ct.c_int()
    KglobalFull = np.zeros((NumDof,NumDof), ct.c_double, 'F'); ks = ct.c_int()
    FglobalFull = np.zeros((NumDof,NumDof), ct.c_double, 'F'); fs = ct.c_int()
    
    X0    = np.zeros(NumDof, ct.c_double, 'F')
    X     = np.zeros(NumDof, ct.c_double, 'F')
    DX    = np.zeros(NumDof, ct.c_double, 'F')
    dXdt  = np.zeros(NumDof, ct.c_double, 'F')
    dXddt = np.zeros(NumDof, ct.c_double, 'F')
    
    Qglobal = np.zeros(NumDof, ct.c_double, 'F')
    
    
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
    
    
    "TODO: Extract initial displacements and velocities"
    """f_cbeam3_solv_disp2state( \
        byref(NumNodes_tot), \
        byref(NumDof), \
        byref(c_int(NumElems)), \
        byref(Nod_Master), \
        byref(Nod_Vdof), \
        byref(Nod_Fdof), \
        PosDefor.ctypes.data_as(POINTER(c_double)), \
        PsiDefor.ctypes.data_as(POINTER(c_double)), \
        PosDotDefor.ctypes.data_as(POINTER(c_double)), \
        PsiDotDefor.ctypes.data_as(POINTER(c_double)), \
        X.ctypes.data_as(POINTER(c_double)), \
        dXdt.ctypes.data_as(POINTER(c_double)) )"""
    
    
    if XBOPTS.PrintInfo.value==True:
        sys.stdout.write(' ... done\n')
    
    
    "Write deformed configuration to file"
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
    
    BeamIO.OutputElems(XBINPUT.NumElems, NumNodes_tot.value, XBELEM, \
                       PosDefor, PsiDefor, ofile, WriteMode)
    
    
    "Print deformed configuration"
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
        
    
    "Return solution as optional output argument"
    return PosDefor, PsiDefor
    
    
    

if __name__ == '__main__':
    """Main"""
    
    """Set up Xbopts for nonlinear static analysis defined in 
    NonlinearStatic/testcases.pdf case 1.1."""
    XBOPTS = DerivedTypes.Xbopts()
    XBOPTS.Solution.value = 312 
    XBOPTS.NumLoadSteps.value = 10
    XBOPTS.MinDelta.value = 1e-05
    XBOPTS.FollowerForce.value = False
         
    """Set up Xbinput for nonlinear static analysis defined in 
    NonlinearStatic/testcases.pdf case 1.1."""
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
    XBINPUT.BeamMass[4,4] = 0.0001 #Must be non-zero
    XBINPUT.BeamMass[5,5] = 0.0001 #Must be non-zero
    XBINPUT.ForceStatic[-1,2] = 6e+05
    XBINPUT.tfin = 1.0
    XBINPUT.dt = 0.1
    
    Solve_F90(XBINPUT,XBOPTS)
    