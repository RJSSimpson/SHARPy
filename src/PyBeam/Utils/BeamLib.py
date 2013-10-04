'''@package PyBeam.Utils.BeamLib
@brief      Loads the f90 subroutines with Python wrappers.
@author     Rob Simpson
@contact    r.simpson11@imperial.ac.uk 
@version    0.0
@date       10/12/2012
@pre        Compile the beam code into wrapped dynamic library.
@warning    None
'''

import ctypes as ct #http://docs.python.org/3.2/library/ctypes.html
import SharPySettings as Settings
import lib_fem
import lib_cbeam3

BeamPath = Settings.BeamLibDir + Settings.BeamLibName
BeamLib = ct.cdll.LoadLibrary(BeamPath)

f_input_setup               = BeamLib.__test_MOD_wrap_input_setup
f_input_elem                = BeamLib.__test_MOD_wrap_input_elem
f_input_node                = BeamLib.__test_MOD_wrap_input_node
f_xbeam_undef_geom          = BeamLib.__test_MOD_wrap_xbeam_undef_geom
f_xbeam_undef_dofs          = BeamLib.__test_MOD_wrap_xbeam_undef_dofs
f_cbeam3_solv_nlnstatic     = BeamLib.__test_MOD_wrap_cbeam3_solv_nlnstatic
f_cbeam3_solv_nlndyn        = BeamLib.__test_MOD_wrap_cbeam3_solv_nlndyn
f_cbeam3_asbly_static       = BeamLib.__test_MOD_wrap_cbeam3_asbly_static
f_cbeam3_solv_disp2state    = BeamLib.__test_MOD_wrap_cbeam3_solv_disp2state
f_fem_m2v                   = BeamLib.__test_MOD_wrap_fem_m2v
f_cbeam3_solv_update_static = BeamLib.__test_MOD_wrap_cbeam3_solv_update_static
f_cbeam3_solv_disp2state    = BeamLib.__test_MOD_wrap_cbeam3_solv_disp2state
f_cbeam3_solv_nlndyn_accel  = BeamLib.__test_MOD_wrap_cbeam3_solv_nlndyn_accel
f_cbeam3_asbly_dynamic = BeamLib.__test_MOD_wrap_cbeam3_asbly_dynamic
f_cbeam3_solv_state2disp = BeamLib.__test_MOD_wrap_cbeam3_solv_state2disp

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
f_cbeam3_asbly_static.restype       = None
f_cbeam3_solv_disp2state.restype    = None
f_fem_m2v.restype                   = None
f_cbeam3_solv_nlndyn_accel.restype  = None
f_cbeam3_asbly_dynamic.restype      = None
f_cbeam3_solv_state2disp.restype    = None

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
    """@brief Python wrapper for f_cbeam3_solv_nlndyn
    
    @details Numpy arrays are mutable so the changes (solution) made here are
     reflected in the data of the calling script after execution."""
    
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


def Cbeam3_Asbly_Static(XBINPUT, NumNodes_tot, XBELEM, XBNODE,\
                        PosIni, PsiIni, PosDefor, PsiDefor,\
                        ForceStatic, NumDof,\
                        ks, KglobalFull, fs, FglobalFull, Qglobal,
                        XBOPTS):
    """@brief Python wrapper for f_cbeam3_asbly_static.
    
    @details Numpy arrays are mutable so the changes made here are 
    reflected in the data of the calling script after execution."""
    
    f_cbeam3_asbly_static( \
                ct.byref(ct.c_int(XBINPUT.NumElems)),\
                ct.byref(NumNodes_tot),\
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
                XBNODE.Master.ctypes.data_as(ct.POINTER(ct.c_int)),\
                XBNODE.Vdof.ctypes.data_as(ct.POINTER(ct.c_int)),\
                XBNODE.Fdof.ctypes.data_as(ct.POINTER(ct.c_int)),\
                PosIni.ctypes.data_as(ct.POINTER(ct.c_double)),\
                PsiIni.ctypes.data_as(ct.POINTER(ct.c_double)),\
                PosDefor.ctypes.data_as(ct.POINTER(ct.c_double)),\
                PsiDefor.ctypes.data_as(ct.POINTER(ct.c_double)),\
                ForceStatic.ctypes.data_as(ct.POINTER(ct.c_double)),\
                ct.byref(ct.c_int(Settings.DimMat)), \
                ct.byref(NumDof), \
                ct.byref(ks), \
                KglobalFull.ctypes.data_as(ct.POINTER(ct.c_double)), \
                ct.byref(fs), \
                FglobalFull.ctypes.data_as(ct.POINTER(ct.c_double)), \
                Qglobal.ctypes.data_as(ct.POINTER(ct.c_double)),\
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
    
    
def Cbeam_Solv_Disp2State(NumNodes_tot, NumDof, XBINPUT, XBNODE,\
                          PosDefor, PsiDefor, PosDotDef, PsiDotDef,
                          x, dxdt):
    """@brief Python wrapper for f_cbeam3_solv_disp2state.
    
    @details Numpy arrays are mutable so the changes made here are 
    reflected in the data of the calling script after execution."""
    
    f_cbeam3_solv_disp2state( \
                ct.byref(NumNodes_tot), \
                ct.byref(NumDof), \
                ct.byref(ct.c_int(XBINPUT.NumElems)), \
                XBNODE.Master.ctypes.data_as(ct.POINTER(ct.c_int)),\
                XBNODE.Vdof.ctypes.data_as(ct.POINTER(ct.c_int)),\
                XBNODE.Fdof.ctypes.data_as(ct.POINTER(ct.c_int)),\
                PosDefor.ctypes.data_as(ct.POINTER(ct.c_double)),\
                PsiDefor.ctypes.data_as(ct.POINTER(ct.c_double)),\
                PosDotDef.ctypes.data_as(ct.POINTER(ct.c_double)), \
                PsiDotDef.ctypes.data_as(ct.POINTER(ct.c_double)), \
                x.ctypes.data_as(ct.POINTER(ct.c_double)), \
                dxdt.ctypes.data_as(ct.POINTER(ct.c_double)) )
    

def Cbeam3_Solv_Update_Static(XBINPUT, NumNodes_tot, XBELEM, XBNODE, NumDof,\
                              DeltaS, PosIni, PsiIni, PosDefor, PsiDefor):
    """@brief Wrapper for f_cbeam3_solv_update_static.
    
    @details Numpy arrays are mutable so the changes made here are 
    reflected in the data of the calling script after execution."""
    
    f_cbeam3_solv_update_static( \
                ct.byref(ct.c_int(XBINPUT.NumElems)), \
                ct.byref(NumNodes_tot), \
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
                XBNODE.Master.ctypes.data_as(ct.POINTER(ct.c_int)),\
                XBNODE.Vdof.ctypes.data_as(ct.POINTER(ct.c_int)),\
                XBNODE.Fdof.ctypes.data_as(ct.POINTER(ct.c_int)),\
                ct.byref(NumDof), \
                DeltaS.ctypes.data_as(ct.POINTER(ct.c_double)), \
                PosIni.ctypes.data_as(ct.POINTER(ct.c_double)), \
                PsiIni.ctypes.data_as(ct.POINTER(ct.c_double)), \
                PosDefor.ctypes.data_as(ct.POINTER(ct.c_double)), \
                PsiDefor.ctypes.data_as(ct.POINTER(ct.c_double)) )



def Cbeam3_Solv_NonlinearDynamicAccel(XBINPUT, XBOPTS, NumNodes_tot, XBELEM, PosIni,\
            PsiIni, XBNODE, NumDof, PosDefor, PsiDefor, NumSteps, Time,\
            ForceTime, ForcedVel, ForcedVelDot, PosDotDef, PsiDotDef,\
            PosDotDotDef, PsiDotDotDef,\
            PosPsiTime, VelocTime, DynOut, OutGrids):
    """@brief Python wrapper for f_cbeam3_solv_nlndyn
    
    @details Numpy arrays are mutable so the changes (solution) made here are
     reflected in the data of the calling script after execution."""
    
    f_cbeam3_solv_nlndyn_accel(ct.byref(ct.c_int(XBINPUT.iOut)),\
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
                PosDotDotDef.ctypes.data_as(ct.POINTER(ct.c_double)),\
                PsiDotDotDef.ctypes.data_as(ct.POINTER(ct.c_double)),\
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
    

def Cbeam3_Asbly_Dynamic(XBINPUT, NumNodes_tot, XBELEM, XBNODE,\
                         PosIni, PsiIni, PosDefor, PsiDefor,\
                         PosDotDef, PsiDotDef, PosDotDotDef, PsiDotDotDef,\
                         Force, Vrel, VrelDot,\
                         NumDof, DimMat,\
                         ms, MglobalFull, Mvel,\
                         cs, CglobalFull, Cvel,\
                         ks, KglobalFull, fs, FglobalFull,\
                         Qglobal, XBOPTS, Cao):
    """@brief Wrapper for f_cbeam3_asbly_dynamic."""

    f_cbeam3_asbly_dynamic( \
            ct.byref(ct.c_int(XBINPUT.NumElems)), \
            ct.byref(NumNodes_tot), \
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
            XBNODE.Master.ctypes.data_as(ct.POINTER(ct.c_int)),\
            XBNODE.Vdof.ctypes.data_as(ct.POINTER(ct.c_int)),\
            XBNODE.Fdof.ctypes.data_as(ct.POINTER(ct.c_int)), \
            PosIni.ctypes.data_as(ct.POINTER(ct.c_double)),\
            PsiIni.ctypes.data_as(ct.POINTER(ct.c_double)),\
            PosDefor.ctypes.data_as(ct.POINTER(ct.c_double)),\
            PsiDefor.ctypes.data_as(ct.POINTER(ct.c_double)),\
            PosDotDef.ctypes.data_as(ct.POINTER(ct.c_double)),\
            PsiDotDef.ctypes.data_as(ct.POINTER(ct.c_double)),\
            PosDotDotDef.ctypes.data_as(ct.POINTER(ct.c_double)), \
            PsiDotDotDef.ctypes.data_as(ct.POINTER(ct.c_double)), \
            Force.ctypes.data_as(ct.POINTER(ct.c_double)), \
            Vrel.ctypes.data_as(ct.POINTER(ct.c_double)), \
            VrelDot.ctypes.data_as(ct.POINTER(ct.c_double)), \
            ct.byref(NumDof), \
            ct.byref(ct.c_int(DimMat)), \
            ct.byref(ms), \
            MglobalFull.ctypes.data_as(ct.POINTER(ct.c_double)), \
            Mvel.ctypes.data_as(ct.POINTER(ct.c_double)), \
            ct.byref(cs), \
            CglobalFull.ctypes.data_as(ct.POINTER(ct.c_double)), \
            Cvel.ctypes.data_as(ct.POINTER(ct.c_double)), \
            ct.byref(ks), \
            KglobalFull.ctypes.data_as(ct.POINTER(ct.c_double)), \
            ct.byref(fs), \
            FglobalFull.ctypes.data_as(ct.POINTER(ct.c_double)), \
            Qglobal.ctypes.data_as(ct.POINTER(ct.c_double)), \
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
            ct.byref(XBOPTS.NewmarkDamp),\
            Cao.ctypes.data_as(ct.POINTER(ct.c_double)) )
    
    
def Cbeam3_Solv_State2Disp(XBINPUT, NumNodes_tot, XBELEM, XBNODE,
                           PosIni, PsiIni, NumDof, X, dXdt,\
                           PosDefor, PsiDefor, PosDotDefor, PsiDotDefor):
    """@brief Wrapper for f_cbeam3_solv_state2disp."""
    
    f_cbeam3_solv_state2disp( \
                ct.byref(ct.c_int(XBINPUT.NumElems)), \
                ct.byref(NumNodes_tot), \
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
                XBNODE.Master.ctypes.data_as(ct.POINTER(ct.c_int)),\
                XBNODE.Vdof.ctypes.data_as(ct.POINTER(ct.c_int)),\
                XBNODE.Fdof.ctypes.data_as(ct.POINTER(ct.c_int)),\
                PosIni.ctypes.data_as(ct.POINTER(ct.c_double)),\
                PsiIni.ctypes.data_as(ct.POINTER(ct.c_double)),\
                ct.byref(NumDof), \
                X.ctypes.data_as(ct.POINTER(ct.c_double)), \
                dXdt.ctypes.data_as(ct.POINTER(ct.c_double)), \
                PosDefor.ctypes.data_as(ct.POINTER(ct.c_double)), \
                PsiDefor.ctypes.data_as(ct.POINTER(ct.c_double)), \
                PosDotDefor.ctypes.data_as(ct.POINTER(ct.c_double)), \
                PsiDotDefor.ctypes.data_as(ct.POINTER(ct.c_double)) )

    
if __name__ == '__main__':
    print(lib_fem.__doc__)
    print(lib_cbeam3.__doc__)