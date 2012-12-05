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
import numpy as np
import SharPySettings as Settings
import DerivedTypes
import BeamIO
import Input

BeamPath = Settings.BeamLibDir + Settings.BeamLibName
BeamLib = ct.cdll.LoadLibrary(BeamPath)

f_input_setup               = BeamLib.__test_MOD_wrap_input_setup
f_input_elem                = BeamLib.__test_MOD_wrap_input_elem
f_input_node                = BeamLib.__test_MOD_wrap_input_node
f_xbeam_undef_geom          = BeamLib.__test_MOD_wrap_xbeam_undef_geom
f_xbeam_undef_dofs          = BeamLib.__test_MOD_wrap_xbeam_undef_dofs
f_cbeam3_solv_nlnstatic      = BeamLib.__test_MOD_wrap_cbeam3_solv_nlnstatic

"""ctypes does not check whether the correct number OR type of input arguments
are passed to each of these functions - great care must be taken to ensure the
number of arguments and the argument types are correct.
TODO: library loading within a separate module which also contains
developer provided argtypes which must be a sequence of C data types.
http://docs.python.org/2/library/ctypes.html"""

f_input_setup.restype               = None
f_input_elem.restype                = None
f_input_node.restype                = None
f_xbeam_undef_geom.restype          = None
f_xbeam_undef_dofs.restype          = None
f_cbeam3_solv_nlnstatic.restype     = None


def Solve_F90(XBINPUT,XBOPTS):
    """ Nonlinear static structural solver using f90 solve routine."""
    
    "Check correct solution code"
    assert XBOPTS.Solution.value == 112, ('NonlinearStatic_F90 requested' +\
                                              ' with wrong solution code')
    
    "Declare variables not dependent on NumNodes_tot"
    XBELEM      = DerivedTypes.Xbelem(XBINPUT.NumElems,Settings.MaxElNod)
    NumDof      = ct.c_int()
    PsiIni      = np.zeros((XBINPUT.NumElems,Settings.MaxElNod,3),\
                         dtype=ct.c_double, order='F')
    
    sys.stdout.write('Setup testcase ... ')
    
    "Check inputs"
    XBINPUT, XBOPTS = Input.Setup(XBINPUT, XBOPTS)
    
    "Set-up element properties"
    NumNodes_tot, XBELEM = Input.Elem(XBINPUT, XBOPTS, XBELEM)
    
    "Set-up nodal properties"
    PosIni, PhiNodes, ForceStatic, BoundConds = \
        Input.Node(XBINPUT, XBOPTS, NumNodes_tot, XBELEM)
        
    "Compute initial (undeformed) geometry"
    f_xbeam_undef_geom( \
        ct.byref(ct.c_int(XBINPUT.NumElems)), \
        ct.byref(NumNodes_tot), \
        XBELEM.NumNodes.ctypes.data_as(ct.POINTER(ct.c_int)), \
        XBELEM.MemNo.ctypes.data_as(ct.POINTER(ct.c_int)), \
        XBELEM.Conn.ctypes.data_as(ct.POINTER(ct.c_int)), \
        XBELEM.Master.ctypes.data_as(ct.POINTER(ct.c_int)), \
        XBELEM.Length.ctypes.data_as(ct.POINTER(ct.c_double)), \
        XBELEM.PreCurv.ctypes.data_as(ct.POINTER(ct.c_double)), \
        XBELEM.Psi.ctypes.data_as(ct.POINTER(ct.c_double)), \
        XBELEM.Vector.ctypes.data_as(ct.POINTER(ct.c_double)), \
        XBELEM.Mass.ctypes.data_as(ct.POINTER(ct.c_double)), \
        XBELEM.Stiff.ctypes.data_as(ct.POINTER(ct.c_double)), \
        XBELEM.InvStiff.ctypes.data_as(ct.POINTER(ct.c_double)), \
        XBELEM.RBMass.ctypes.data_as(ct.POINTER(ct.c_double)),\
        PosIni.ctypes.data_as(ct.POINTER(ct.c_double)), \
        PhiNodes.ctypes.data_as(ct.POINTER(ct.c_double)), \
        PsiIni.ctypes.data_as(ct.POINTER(ct.c_double)),\
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

    "Write to undeformed geometry to file"
    WriteMode = 'a'
    BeamIO.WriteUndefGeometry(XBINPUT.NumElems,NumNodes_tot.value,XBELEM,\
                              PosIni,PsiIni,\
                              Settings.OutputFileRoot + '_SOL112',WriteMode)
    
    
    
    
    sys.stdout.write('Identify nodal degrees of freedom ... ')
    
    XBNODE = DerivedTypes.Xbnode(NumNodes_tot.value)
    
    f_xbeam_undef_dofs( \
        ct.byref(ct.c_int(XBINPUT.NumElems)), \
        ct.byref(NumNodes_tot), \
        XBELEM.NumNodes.ctypes.data_as(ct.POINTER(ct.c_int)), \
        XBELEM.MemNo.ctypes.data_as(ct.POINTER(ct.c_int)), \
        XBELEM.Conn.ctypes.data_as(ct.POINTER(ct.c_int)), \
        XBELEM.Master.ctypes.data_as(ct.POINTER(ct.c_int)), \
        XBELEM.Length.ctypes.data_as(ct.POINTER(ct.c_double)), \
        XBELEM.PreCurv.ctypes.data_as(ct.POINTER(ct.c_double)), \
        XBELEM.Psi.ctypes.data_as(ct.POINTER(ct.c_double)), \
        XBELEM.Vector.ctypes.data_as(ct.POINTER(ct.c_double)), \
        XBELEM.Mass.ctypes.data_as(ct.POINTER(ct.c_double)), \
        XBELEM.Stiff.ctypes.data_as(ct.POINTER(ct.c_double)), \
        XBELEM.InvStiff.ctypes.data_as(ct.POINTER(ct.c_double)), \
        XBELEM.RBMass.ctypes.data_as(ct.POINTER(ct.c_double)), \
        BoundConds.ctypes.data_as(ct.POINTER(ct.c_int)), \
        XBNODE.Master.ctypes.data_as(ct.POINTER(ct.c_int)), \
        XBNODE.Vdof.ctypes.data_as(ct.POINTER(ct.c_int)), \
        XBNODE.Fdof.ctypes.data_as(ct.POINTER(ct.c_int)), \
        ct.byref(NumDof) )
    
    sys.stdout.write('done\n')
    
    
    PosDefor = PosIni.copy(order='F') #Set initial conditions to undef config
    PsiDefor = PsiIni.copy(order='F')
    
    
    sys.stdout.write('Solve nonlinear static case (using .f90 routines) ... \n')
    
    f_cbeam3_solv_nlnstatic(ct.byref(NumDof),\
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
    
    sys.stdout.write(' ... done\n')
    
    
    ofile = Settings.OutputFileRoot + '_SOL112_def.dat'
    sys.stdout.write('Writing file %s ... ' %(ofile))
    fp = open(ofile,'w')
    fp.write('TITLE="Non-linear static solution: deformed geometry"\n')
    fp.write('VARIABLES="iElem" "iNode" "Px" "Py" "Pz" "Rx" "Ry" "Rz"\n')
    fp.close()
    sys.stdout.write('done\n')
    WriteMode = 'a'
    
    BeamIO.OutputElems(XBINPUT.NumElems, NumNodes_tot.value, XBELEM, \
                       PosDefor, PsiDefor, ofile, WriteMode)
    
    
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