'''@package PyBeam.Main.SharPySettings
@brief      Nonlinear static solvers and related tests.
@author     Rob Simpson
@contact    r.simpson11@imperial.ac.uk
@version    0.0
@date       25/10/2012
@pre        None
@warning    None
'''

import sys
import ctypes as ct #http://docs.python.org/3.2/library/ctypes.html
import SharPySettings as Settings
import DerivedTypes
import BeamIO
import Input
import BeamLib
import BeamInit

def Solve_F90(XBINPUT,XBOPTS):
    """ Nonlinear static structural solver using f90 solve routine."""
    
    "Check correct solution code"
    assert XBOPTS.Solution.value == 112, ('NonlinearStatic_F90 requested' +\
                                              ' with wrong solution code')
    
    "Initialise beam"
    XBINPUT, XBOPTS, NumNodes_tot, XBELEM, PosIni, PsiIni,\
            ForceStatic, XBNODE, NumDof \
                = BeamInit.Initialise(XBINPUT,XBOPTS)
    
    
    "Set initial conditions as undef config"
    PosDefor = PosIni.copy(order='F')
    PsiDefor = PsiIni.copy(order='F')
    
    
    if XBOPTS.PrintInfo==True:
        sys.stdout.write('Solve nonlinear static case (using .f90 routines) ... \n')
    
    BeamLib.f_cbeam3_solv_nlnstatic(ct.byref(NumDof),\
                            ct.byref(ct.c_int(XBINPUT.NumElems)),\
                            XBELEM.NumNodes.ctypes.data_as(ct.POINTER(ct.c_int)),\
                            XBELEM.MemNo.ctypes.data_as(ct.POINTER(ct.c_int)),\
                            XBELEM.Conn.ctypes.data_as(ct.POINTER(ct.c_int)),\
                            XBELEM.Master.ctypes.data_as(ct.POINTER(ct.c_int)),\
                            XBELEM.Length.ctypes.data_as(ct.POINTER(ct.c_double)),\
                            XBELEM.PreCurv.ctypes.data_as(ct.POINTER(ct.c_double)),\
                            XBELEM.Psi.ctypes.data_as(ct.POINTER(ct.c_double)),\
                            XBELEM.Vector.ctypes.data_as(ct.POINTER(ct.c_double)),\
                            XBELEM.Mass.ctypes.data_as(ct.POINTER(ct.c_double)),\
                            XBELEM.Stiff.ctypes.data_as(ct.POINTER(ct.c_double)),\
                            XBELEM.InvStiff.ctypes.data_as(ct.POINTER(ct.c_double)),\
                            XBELEM.RBMass.ctypes.data_as(ct.POINTER(ct.c_double)), \
                            ct.byref(NumNodes_tot),\
                            XBNODE.Master.ctypes.data_as(ct.POINTER(ct.c_int)),\
                            XBNODE.Vdof.ctypes.data_as(ct.POINTER(ct.c_int)),\
                            XBNODE.Fdof.ctypes.data_as(ct.POINTER(ct.c_int)),\
                            ForceStatic.ctypes.data_as(ct.POINTER(ct.c_double)),\
                            PosIni.ctypes.data_as(ct.POINTER(ct.c_double)),\
                            PsiIni.ctypes.data_as(ct.POINTER(ct.c_double)),\
                            PosDefor.ctypes.data_as(ct.POINTER(ct.c_double)),\
                            PsiDefor.ctypes.data_as(ct.POINTER(ct.c_double)),\
                            ct.byref(XBOPTS.FollowerForce),\
                            ct.byref(XBOPTS.FollowerForceRig),\
                            ct.byref(XBOPTS.PrintInfo),\
                            ct.byref(XBOPTS.OutInBframe),\
                            ct.byref(XBOPTS.OutInaframe),\
                            ct.byref(XBOPTS.ElemProj),\
                            ct.byref(XBOPTS.MaxIterations),\
                            ct.byref(XBOPTS.NumLoadSteps),\
                            ct.byref(XBOPTS.NumGauss),\
                            ct.byref(XBOPTS.Solution),\
                            ct.byref(XBOPTS.DeltaCurved),\
                            ct.byref(XBOPTS.MinDelta),\
                            ct.byref(XBOPTS.NewmarkDamp) )
    
    if XBOPTS.PrintInfo==True:
        sys.stdout.write(' ... done\n')
    
    
    "Write deformed configuration to file"
    ofile = Settings.OutputFileRoot + '_SOL112_def.dat'
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
    
    "Print deformed configuration"
    if XBOPTS.PrintInfo==True:
        sys.stdout.write('--------------------------------------\n')
        sys.stdout.write('NONLINEAR STATIC SOLUTION\n')
        sys.stdout.write('%10s %10s %10s\n' %('X','Y','Z'))
        for inodi in range(NumNodes_tot.value):
            sys.stdout.write(' ')
            for inodj in range(3):
                sys.stdout.write('%12.5e' %(PosDefor[inodi,inodj]))
                sys.stdout.write('\n')
        sys.stdout.write('--------------------------------------\n')
    

if __name__ == '__main__':
    """Set up Xbopts for nonlinear static analysis defined in input_rob.f90
    TPY0 test case"""
    XBOPTS = DerivedTypes.Xbopts()
    XBOPTS.Solution.value = 112 
    XBOPTS.NumLoadSteps.value = 10
    XBOPTS.MinDelta.value = 1e-04       
    """Set up Xbinput for nonlinear static analysis defined in input_rob.f90
    TPY0 test case"""
    XBINPUT = DerivedTypes.Xbinput()
    XBINPUT.NumElems = 8
    XBINPUT.BeamLength = 16.0
    XBINPUT.BeamStiffness[0,0] = 1.0e+09
    XBINPUT.BeamStiffness[1,1] = 1.0e+09
    XBINPUT.BeamStiffness[2,2] = 1.0e+09
    XBINPUT.BeamStiffness[3,3] = 1.0e+04
    XBINPUT.BeamStiffness[4,4] = 2.0e+04
    XBINPUT.BeamStiffness[5,5] = 4.0e+06
    XBINPUT.BeamMass[0,0] = 0.75
    XBINPUT.BeamMass[1,1] = 0.75
    XBINPUT.BeamMass[2,2] = 0.75
    XBINPUT.BeamMass[3,3] = 0.1
    XBINPUT.BeamMass[4,4] = 0.001
    XBINPUT.BeamMass[5,5] = 0.001
    XBINPUT.ForceStatic[2] = 800

    Solve_F90(XBINPUT,XBOPTS)