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
    ofile = Settings.OutputFileRoot + '_SOL312_def.dat'
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
    ofile = Settings.OutputFileRoot + '_SOL312_force.dat'
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
    ofile = Settings.OutputFileRoot + '_SOL312_dyn.dat'
    fp = open(ofile,'w')
    BeamIO.Write_dyn_File(fp, Time, PosPsiTime)
    fp.close()
    
    
    "Write _shape file"
    ofile = Settings.OutputFileRoot + '_SOL312_shape.dat'
    fp = open(ofile,'w')
    BeamIO.Write_shape_File(fp, len(Time), NumNodes_tot.value, Time, DynOut)
    fp.close()
    
    
    

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
    