'''@package PyBeam.Utils.BeamInit
@brief      For initialising variables required for beam simulations.
@author     Rob Simpson
@contact    r.simpson11@imperial.ac.uk
@version    0.0
@date       10/12/2012
@pre        None
@warning    None'''

import sys
import DerivedTypes
import SharPySettings as Settings
import numpy as np
import ctypes as ct
import BeamIO
import BeamLib
import Input

def Static(XBINPUT,XBOPTS):
    """@brief Initialise everything needed for Static beam simulation."""
    
    "Declare variables not dependent on NumNodes_tot"
    XBELEM      = DerivedTypes.Xbelem(XBINPUT.NumElems,Settings.MaxElNod)
    NumDof      = ct.c_int()
    PsiIni      = np.zeros((XBINPUT.NumElems,Settings.MaxElNod,3),\
                         dtype=ct.c_double, order='F')
    
    "Check inputs"
    XBINPUT, XBOPTS = Input.Setup(XBINPUT, XBOPTS)
    
    "Set-up element properties"
    NumNodes_tot, XBELEM = Input.Elem(XBINPUT, XBOPTS, XBELEM)
    
    "Set-up nodal properties"
    PosIni, PhiNodes, BoundConds = \
        Input.Node(XBINPUT, XBOPTS, NumNodes_tot, XBELEM)
        
    "Compute initial (undeformed) geometry"
    BeamLib.f_xbeam_undef_geom( \
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
    WriteMode = 'w'
    if XBOPTS.PrintInfo==True:
        sys.stdout.write('Writing file (mode: %s) %s' %(WriteMode,\
                            Settings.OutputFileRoot + '_SOL112'\
                            + '_und.dat\n'))
        
    BeamIO.WriteUndefGeometry(XBINPUT.NumElems,NumNodes_tot.value,XBELEM,\
                              PosIni,PsiIni,\
                              Settings.OutputFileRoot + '_SOL112',WriteMode)
    
    
    if XBOPTS.PrintInfo==True:
        sys.stdout.write('Identify nodal degrees of freedom ... ')
    
    XBNODE = DerivedTypes.Xbnode(NumNodes_tot.value)
    
    BeamLib.f_xbeam_undef_dofs( \
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
    
    if XBOPTS.PrintInfo==True:
        sys.stdout.write('done\n')
        
    return XBINPUT, XBOPTS, NumNodes_tot, XBELEM, PosIni, PsiIni,\
            XBNODE, NumDof

if __name__ == '__main__':
    XBINPUT = DerivedTypes.Xbinput()
    XBOPTS = DerivedTypes.Xbopts()
    
    XBINPUT.BeamStiffness[0,0] = 1.0e+09 #otherwise inverse fails
    XBINPUT.BeamStiffness[1,1] = 1.0e+09
    XBINPUT.BeamStiffness[2,2] = 1.0e+09
    XBINPUT.BeamStiffness[3,3] = 1.0e+04
    XBINPUT.BeamStiffness[4,4] = 2.0e+04
    XBINPUT.BeamStiffness[5,5] = 4.0e+06
    
    XBINPUT, XBOPTS, NumNodes_tot, XBELEM, PosIni, PsiIni,\
            ForceStatic, XBNODE, NumDof \
                = Static(XBINPUT,XBOPTS)               
