'''@package PyBeam.Utils.BeamLib
@brief      Loads the f90 subroutines.
@author     Rob Simpson
@contact    r.simpson11@imperial.ac.uk 
@version    0.0
@date       10/12/2012
@pre        None
@warning    None
'''

import ctypes as ct #http://docs.python.org/3.2/library/ctypes.html
import SharPySettings as Settings

BeamPath = Settings.BeamLibDir + Settings.BeamLibName
BeamLib = ct.cdll.LoadLibrary(BeamPath)

f_input_setup               = BeamLib.__test_MOD_wrap_input_setup
f_input_elem                = BeamLib.__test_MOD_wrap_input_elem
f_input_node                = BeamLib.__test_MOD_wrap_input_node
f_xbeam_undef_geom          = BeamLib.__test_MOD_wrap_xbeam_undef_geom
f_xbeam_undef_dofs          = BeamLib.__test_MOD_wrap_xbeam_undef_dofs
f_cbeam3_solv_nlnstatic     = BeamLib.__test_MOD_wrap_cbeam3_solv_nlnstatic
f_cbeam3_solv_nlndyn        = BeamLib.__test_MOD_wrap_cbeam3_solv_nlndyn

"""ctypes does not check whether the correct number OR type of input arguments
are passed to each of these functions - great care must be taken to ensure the
number of arguments and the argument types are correct.
TODO: developer provided argtypes which must be a sequence of C data types.
http://docs.python.org/2/library/ctypes.html"""

f_input_setup.restype               = None
f_input_elem.restype                = None
f_input_node.restype                = None
f_xbeam_undef_geom.restype          = None
f_xbeam_undef_dofs.restype          = None
f_cbeam3_solv_nlnstatic.restype     = None
f_cbeam3_solv_nlndyn.restype        = None

def Cbeam3_Solv_NonlinearStatic(XBINPUT, XBOPTS, NumNodes_tot, XBELEM, PosIni, \
            PsiIni, XBNODE, NumDof, PosDefor, PsiDefor):
    """@brief Python wrapper for f_cbeam3_solv_nlnstatic
    
    @details Numpy arrays are mutable so the changes (solution) made here are
     reflected in the data of the calling script after execution."""
    
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
                XBINPUT.ForceStatic.ctypes.data_as(ct.POINTER(ct.c_double)),\
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
    

def Cbeam3_Solv_NonlinearDynamic(XBINPUT, XBOPTS, NumNodes_tot, XBELEM, PosIni,\
            PsiIni, XBNODE, NumDof, PosDefor, PsiDefor, NumSteps, Time,\
            ForceTime, ForcedVel, ForcedVelDot, PosDotDef, PsiDotDef,\
            PosPsiTime, VelocTime, DynOut, OutGrids):
    """@brief Python wrapper for f_cbeam3_solv_nlndyn"""
    
    f_cbeam3_solv_nlndyn(ct.byref(ct.c_int(XBINPUT.iOut)),\
                ct.byref(NumDof),\
                ct.byref(NumSteps),\
                Time.ctypes.data_as(ct.POINTER(ct.c_double)),\
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
                XBINPUT.ForceStatic.ctypes.data_as(ct.POINTER(ct.c_double)),\
                XBINPUT.ForceDyn.ctypes.data_as(ct.POINTER(ct.c_double)),\
                ForceTime.ctypes.data_as(ct.POINTER(ct.c_double)),\
                ForcedVel.ctypes.data_as(ct.POINTER(ct.c_double)),\
                ForcedVelDot.ctypes.data_as(ct.POINTER(ct.c_double)),\
                PosIni.ctypes.data_as(ct.POINTER(ct.c_double)),\
                PsiIni.ctypes.data_as(ct.POINTER(ct.c_double)),\
                PosDefor.ctypes.data_as(ct.POINTER(ct.c_double)),\
                PsiDefor.ctypes.data_as(ct.POINTER(ct.c_double)),\
                PosDotDef.ctypes.data_as(ct.POINTER(ct.c_double)),\
                PsiDotDef.ctypes.data_as(ct.POINTER(ct.c_double)),\
                PosPsiTime.ctypes.data_as(ct.POINTER(ct.c_double)),\
                VelocTime.ctypes.data_as(ct.POINTER(ct.c_double)),\
                DynOut.ctypes.data_as(ct.POINTER(ct.c_double)),\
                OutGrids.ctypes.data_as(ct.POINTER(ct.c_bool)),\
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

    
if __name__ == '__main__':
    pass