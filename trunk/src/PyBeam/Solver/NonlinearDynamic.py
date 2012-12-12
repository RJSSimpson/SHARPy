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
import ctypes as ct #http://docs.python.org/3.2/library/ctypes.html
import SharPySettings as Settings
import DerivedTypes
import BeamIO
import BeamLib
import BeamInit
import NonlinearStatic

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
    
    
    

    """subroutine wrap_cbeam3_solv_nlndyn(iOut,NumDof,NumSteps,Time,&
&                    NumElems, NumNodes, MemNo, Conn,        &!for do_xbelem_var
&                    Master_Array,                             &!for do_xbelem_var
&                    Length, PreCurv,                        &!for do_xbelem_var
&                    Psi, Vector, Mass_Array,                &!for do_xbelem_var
&                    Stiff_Array,                            &!for do_xbelem_var
&                    InvStiff_Array, RBMass_Array,            &!for do_xbelem_var
&                    NumNodes_tot, Master, Vdof, Fdof,        &!for pack_xbnode
&                    F0_Vec,Fa_Vec,Ftime,                            &
&                   Vrel_Vec, VrelDot_Vec, Coords_Vec, Psi0_Vec, PosDefor_Vec,    &
&                    PsiDefor_Vec, PosDotDefor_Vec, PsiDotDefor_Vec,        &
&                   PosPsiTime_Vec,VelocTime_Vec,DynOut_Vec,            &
&                   OutGrids,                                &
&                   FollowerForce, FollowerForceRig,        &!for pack_xbopts
&                    PrintInfo, OutInBframe, OutInaframe,    &!for pack_xbopts
&                    ElemProj, MaxIterations, NumLoadSteps,    &!for pack_xbopts
&                    NumGauss, Solution, DeltaCurved,         &!for pack_xbopts
&                    MinDelta, NewmarkDamp)                     !for pack_xbopts"""

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
    XBINPUT = DerivedTypes.Xbinput(2,20,10)
    XBINPUT.NumElems = 20
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
    XBINPUT.ForceStatic[2] = 6e+05
    
    Solve_F90(XBINPUT,XBOPTS)
    