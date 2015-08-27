'''@package PyAero.UVLM.Utils.Vorticity
@brief      Low level functions for evaluating the biot-savart law.
@author     Rob Simpson
@contact    r.simpson11@imperial.ac.uk 
@version    0.0
@date       29/01/2013
@pre        None
@warning    None
'''

import ctypes as ct
import SharPySettings as Settings
import numpy as np

UVLMPath = Settings.UVLMLibDir + Settings.UVLMLibName
UVLMLib = ct.cdll.LoadLibrary(UVLMPath)

cpp_vorticity_biotsegment = UVLMLib.cpp_wrap_vorticity_biotsegment
cpp_vorticity_biotsegment_map = UVLMLib.cpp_wrap_vorticity_biotsegment_map
c_vorticity_biotsegment = UVLMLib.c_wrap_vorticity_biotsegment
cpp_test_biotsegment = UVLMLib.cpp_wrap_test_biotsegment
c_test_biotsegment = UVLMLib.c_wrap_test_biotsegment
cpp_solver_vlm = UVLMLib.cpp_wrap_solver_vlm
cpp_AIC = UVLMLib.cpp_wrap_AIC
cpp_dAgamma0_dZeta = UVLMLib.cpp_wrap_dAgamma0_dZeta
cpp_dWzetaPri0_dZeta = UVLMLib.cpp_wrap_dWzetaPri0_dZeta
cpp_genW = UVLMLib.cpp_wrap_genW
cpp_genH = UVLMLib.cpp_wrap_genH
cpp_genXi = UVLMLib.cpp_wrap_genXi
cpp_AIC3 = UVLMLib.cpp_wrap_AIC3
cpp_AIC3s = UVLMLib.cpp_wrap_AIC3s
cpp_AIC3s_noTE = UVLMLib.cpp_wrap_AIC3s_noTE
cpp_AIC3noTE = UVLMLib.cpp_wrap_AIC3noTE
cpp_dA3gamma0_dZeta = UVLMLib.cpp_wrap_dA3gamma0_dZeta
cpp_Y1 = UVLMLib.cpp_wrap_Y1
cpp_Y2 = UVLMLib.cpp_wrap_Y2
cpp_Y3 = UVLMLib.cpp_wrap_Y3
cpp_Y4 = UVLMLib.cpp_wrap_Y4
cpp_Y5 = UVLMLib.cpp_wrap_Y5
cpp_KJforces = UVLMLib.cpp_wrap_KJMethodForces
cpp_KJforces_vC = UVLMLib.cpp_wrap_KJMethodForces_vC
cpp_KJforces_vC_mod = UVLMLib.cpp_wrap_KJMethodForces_vC_mod
cpp_q_k = UVLMLib.cpp_wrap_q_k
cpp_dAs3gam0_dZeta_num = UVLMLib.cpp_wrap_dAs3gam0_dZeta_numerical

"""ctypes does not check whether the correct number OR type of input arguments
are passed to each of these functions - great care must be taken to ensure the
number of arguments and the argument types are correct.
TODO: developer provided argtypes which must be a sequence of C data types.
http://docs.python.org/2/library/ctypes.html"""

cpp_vorticity_biotsegment.restype = None
cpp_vorticity_biotsegment_map.restype = None
c_vorticity_biotsegment.restype = None
cpp_test_biotsegment.restype = None
c_test_biotsegment.restype = None
cpp_solver_vlm.restype = None
cpp_AIC.restype = None
cpp_dAgamma0_dZeta.restype = None
cpp_dWzetaPri0_dZeta.restype = None
cpp_genW.restype = None
cpp_genH.restype = None
cpp_genXi.restype = None
cpp_AIC3.restype = None
cpp_AIC3s.restype = None
cpp_AIC3s_noTE.restype = None
cpp_AIC3noTE.restype = None
cpp_dA3gamma0_dZeta.restype = None
cpp_Y1.restype = None
cpp_Y2.restype = None
cpp_Y3.restype = None
cpp_Y4.restype = None
cpp_Y5.restype = None
cpp_KJforces.restype = None
cpp_KJforces_vC.restype = None
cpp_KJforces_vC_mod.restype = None
cpp_q_k.restype = ct.c_uint
cpp_dAs3gam0_dZeta_num.restype = None

# arg types
cpp_q_k.argtypes = [ct.c_uint, ct.c_uint, ct.c_uint]

def AreEqual(arg1,arg2):
    """@brief Returns true if argumants are equal."""
    if arg1 == arg2:
        return True
    else:
        return False

def Cpp_Vorticity_BiotSegment(xP,x1,x2,Gamma,Uind):
    """@brief Wrapper for c_vorticity_biotsegment."""
    
    cpp_vorticity_biotsegment(xP.ctypes.data_as(ct.POINTER(ct.c_double)), \
                            x1.ctypes.data_as(ct.POINTER(ct.c_double)), \
                            x2.ctypes.data_as(ct.POINTER(ct.c_double)), \
                            ct.byref(ct.c_double(Gamma)), \
                            Uind.ctypes.data_as(ct.POINTER(ct.c_double)) )
    
    
def Cpp_Vorticity_BiotSegment_Map(xP,x1,x2,Gamma,Uind):
    """@brief Wrapper for c_vorticity_biotsegment."""
    
    cpp_vorticity_biotsegment_map(xP.ctypes.data_as(ct.POINTER(ct.c_double)), \
                            x1.ctypes.data_as(ct.POINTER(ct.c_double)), \
                            x2.ctypes.data_as(ct.POINTER(ct.c_double)), \
                            ct.byref(ct.c_double(Gamma)), \
                            Uind.ctypes.data_as(ct.POINTER(ct.c_double)) )
    
    
def C_Vorticity_BiotSegment(xP,x1,x2,Gamma,Uind):
    """@brief Wrapper for c_vorticity_biotsegment."""
    
    c_vorticity_biotsegment(xP.ctypes.data_as(ct.POINTER(ct.c_double)), \
                            x1.ctypes.data_as(ct.POINTER(ct.c_double)), \
                            x2.ctypes.data_as(ct.POINTER(ct.c_double)), \
                            ct.byref(ct.c_double(Gamma)), \
                            Uind.ctypes.data_as(ct.POINTER(ct.c_double)) )
    
    
def Cpp_test_biotsegment(NumTests):
    """@brief Wrapper for cpp_test_biotsegment."""
    
    cpp_test_biotsegment(ct.byref(ct.c_int(NumTests)))
    
    
def C_test_biotsegment(NumTests):
    """@brief Wrapper for cpp_test_biotsegment."""
    
    c_test_biotsegment(ct.byref(ct.c_int(NumTests)))


def Cpp_Solver_VLM(Zeta, ZetaDot, Uext, ZetaStar, VMOPTS, Forces, \
                   Gamma, GammaStar, AIC = None, BIC = None):
    """@details wrapper for cpp_solver_vlm."""
    
    "If memory for AIC and BIC is not allocated, allocate here."
    if AIC == None and BIC == None:
        AIC = np.zeros((VMOPTS.M.value*VMOPTS.N.value,\
                        VMOPTS.M.value*VMOPTS.N.value), \
                        ct.c_double,'C')
        BIC = np.zeros((VMOPTS.M.value*VMOPTS.N.value,\
                        VMOPTS.M.value*VMOPTS.N.value), \
                        ct.c_double,'C')
    
    cpp_solver_vlm(Zeta.ctypes.data_as(ct.POINTER(ct.c_double)),\
                   ZetaDot.ctypes.data_as(ct.POINTER(ct.c_double)),\
                   Uext.ctypes.data_as(ct.POINTER(ct.c_double)),\
                   ZetaStar.ctypes.data_as(ct.POINTER(ct.c_double)),\
                   ct.byref(VMOPTS.M),\
                   ct.byref(VMOPTS.N),\
                   ct.byref(VMOPTS.ImageMethod), \
                   ct.byref(VMOPTS.Mstar),\
                   ct.byref(VMOPTS.Steady),\
                   ct.byref(VMOPTS.KJMeth),\
                   ct.byref(VMOPTS.NewAIC),\
                   ct.byref(VMOPTS.DelTime),\
                   ct.byref(VMOPTS.Rollup),\
                   ct.byref(VMOPTS.NumCores),\
                   Forces.ctypes.data_as(ct.POINTER(ct.c_double)), \
                   Gamma.ctypes.data_as(ct.POINTER(ct.c_double)), \
                   GammaStar.ctypes.data_as(ct.POINTER(ct.c_double)), \
                   AIC.ctypes.data_as(ct.POINTER(ct.c_double)), \
                   BIC.ctypes.data_as(ct.POINTER(ct.c_double)))

def Colloc(zeta_G_panel):
    """@brief Collocation from grid position.
        
    @param zeta_G_panel A (2,2,3) array of (i,j) panel corner points.
    @return colloc panel collocation point.
    """
    
    if zeta_G_panel.shape[0] > 2 or zeta_G_panel.shape[1] > 2:
        raise TypeError("Extent of corner point array greater than 2x2(x3).")
    elif zeta_G_panel.shape[0] < 2 or zeta_G_panel.shape[1] < 2:
        raise TypeError("Extent of corner point array less than 2x2(x3).")
    
    return ( 0.25 * (zeta_G_panel[0,0,:] +
                     zeta_G_panel[0,1,:] +
                     zeta_G_panel[1,1,:] +
                     zeta_G_panel[1,0,:]  ) )
    
def Cpp_AIC(zetaSrc,mSrc,nSrc,zetaTgt,mTgt,nTgt,AIC,imageMeth=False):
    """@details Wrapper for cpp_wrap_AIC.
    @param zetaSrc Source grid points.
    @param mSrc Chordwise panels.
    @param nSrc Spanwise.
    @param zetaTgt Target grid points.
    @param mTgt Chordwise panels.
    @param nTgt Spanwise.
    @param imageMeth Flag, default False. Image across xy-plane. 
    @return AIC Influence coefficient matrix."""
    
    # make sure y>=0 for all zetas
    if imageMeth:
        assert np.all( zetaSrc[1::3] >= 0.0 ), "source grid defintion in negative y"
        assert np.all( zetaTgt[1::3] >= 0.0 ), "target grid defintion in negative y"
    # end if
    
    cpp_AIC(zetaSrc.ctypes.data_as(ct.POINTER(ct.c_double)),
            ct.byref(ct.c_uint(mSrc)),
            ct.byref(ct.c_uint(nSrc)),
            zetaTgt.ctypes.data_as(ct.POINTER(ct.c_double)),
            ct.byref(ct.c_uint(mTgt)),
            ct.byref(ct.c_uint(nTgt)),
            ct.byref(ct.c_bool(imageMeth)),
            AIC.ctypes.data_as(ct.POINTER(ct.c_double)))
    

def Cpp_dAgamma0_dZeta(zetaSrc,mSrc,nSrc,gamma0,zetaTgt,mTgt,nTgt,dAgam0):
    """@details Wrapper for cpp_wrap_dAgamma0_dZeta.
    @param zetaSrc Source grid points.
    @param mSrc Chordwise panels.
    @param nSrc Spanwise.
    @param gamma0 Reference circulation distribution.
    @param zetaTgt Target grid points.
    @param mTgt Chordwise panels.
    @param nTgt Spanwise.
    @return dAgam0 variation of influence coefficient matrix times gamma."""
    cpp_dAgamma0_dZeta(zetaSrc.ctypes.data_as(ct.POINTER(ct.c_double)),
                       ct.byref(ct.c_uint(mSrc)),
                       ct.byref(ct.c_uint(nSrc)),
                       gamma0.ctypes.data_as(ct.POINTER(ct.c_double)),
                       zetaTgt.ctypes.data_as(ct.POINTER(ct.c_double)),
                       ct.byref(ct.c_uint(mTgt)),
                       ct.byref(ct.c_uint(nTgt)),
                       dAgam0.ctypes.data_as(ct.POINTER(ct.c_double)))
    
def Cpp_dWzetaPri0_dZeta(zeta,m,n,zetaPri,dWzetaPri):
    """@details Wrapper for cpp_wrap_dWzetaPri0_dZeta.
    @param zeta Grid points.
    @param m Chordwise panels.
    @param n Spanwise.
    @param zetaPri0 Reference grid velocities.
    @param dWzetaPri Variation of downwash matrix times zetaPri."""
    cpp_dWzetaPri0_dZeta(zeta.ctypes.data_as(ct.POINTER(ct.c_double)),
                         ct.byref(ct.c_uint(m)),
                         ct.byref(ct.c_uint(n)),
                         zetaPri.ctypes.data_as(ct.POINTER(ct.c_double)),
                         dWzetaPri.ctypes.data_as(ct.POINTER(ct.c_double)))
    
def Cpp_genW(zeta,m,n,W):
    """@details Wrapper for cpp_wrap_genW.
    @param zeta Source grid points.
    @param m Chordwise panels.
    @param n Spanwise.
    @return W Interpolation and projection matrix."""
    cpp_genW(zeta.ctypes.data_as(ct.POINTER(ct.c_double)),
             ct.byref(ct.c_uint(m)),
             ct.byref(ct.c_uint(n)),
             W.ctypes.data_as(ct.POINTER(ct.c_double)));
             
def Cpp_genH(m,n,H):
    """@details Wrapper for cpp_wrap_genH.
    @param m Chordwise panels.
    @param n Spanwise.
    @return H sgement to lattice vertex distribution matrix."""
    cpp_genH(ct.byref(ct.c_uint(m)),
             ct.byref(ct.c_uint(n)),
             H.ctypes.data_as(ct.POINTER(ct.c_double)));
             
def Cpp_genXi(m,n,eta1,eta2,Xi):
    """@details Wrapper for cpp_wrap_genXi.
    @param m Chordwise panels.
    @param n Spanwise.
    @param eta1 Computational coordinate in axis 1.
    @param eta2 Computational coordinate in axis 2.
    @return Xi Lattice vertex to collocation interpolation matrix (3Kx3K_zeta)."""
    cpp_genXi(ct.byref(ct.c_uint(m)),
             ct.byref(ct.c_uint(n)),
             ct.byref(ct.c_double(eta1)),
             ct.byref(ct.c_double(eta2)),
             Xi.ctypes.data_as(ct.POINTER(ct.c_double)));

def Cpp_AIC3(zetaSrc,mSrc,nSrc,zetaTgt,mTgt,nTgt,AIC3):
    """@details Wrapper for cpp_wrap_AICnoTE.
    @param zetaSrc Source grid points.
    @param mSrc Chordwise panels.
    @param nSrc Spanwise.
    @param zetaTgt Target grid points.
    @param mTgt Chordwise panels.
    @param nTgt Spanwise.
    @return AIC3 3-component influence coefficient matrix."""
    cpp_AIC3(zetaSrc.ctypes.data_as(ct.POINTER(ct.c_double)),
            ct.byref(ct.c_uint(mSrc)),
            ct.byref(ct.c_uint(nSrc)),
            zetaTgt.ctypes.data_as(ct.POINTER(ct.c_double)),
            ct.byref(ct.c_uint(mTgt)),
            ct.byref(ct.c_uint(nTgt)),
            AIC3.ctypes.data_as(ct.POINTER(ct.c_double)))

def Cpp_AIC3s(zetaSrc,mSrc,nSrc,zetaTgt,mTgt,nTgt,AIC3):
    """@details Wrapper for cpp_wrap_AICnoTE.
    @param zetaSrc Source grid points.
    @param mSrc Chordwise panels.
    @param nSrc Spanwise.
    @param zetaTgt Target grid points.
    @param mTgt Chordwise panels.
    @param nTgt Spanwise.
    @return AIC3 3-component influence coefficient matrix at segment midpoints."""
    cpp_AIC3s(zetaSrc.ctypes.data_as(ct.POINTER(ct.c_double)),
            ct.byref(ct.c_uint(mSrc)),
            ct.byref(ct.c_uint(nSrc)),
            zetaTgt.ctypes.data_as(ct.POINTER(ct.c_double)),
            ct.byref(ct.c_uint(mTgt)),
            ct.byref(ct.c_uint(nTgt)),
            AIC3.ctypes.data_as(ct.POINTER(ct.c_double)))
 
def Cpp_AIC3s_noTE(zetaSrc,mSrc,nSrc,zetaTgt,mTgt,nTgt,wakeSrc,AIC3):
    """@details Wrapper for cpp_wrap_AICnoTE.
    @param zetaSrc Source grid points.
    @param mSrc Chordwise panels.
    @param nSrc Spanwise.
    @param zetaTgt Target grid points.
    @param mTgt Chordwise panels.
    @param nTgt Spanwise.
    @return AIC3 3-component influence coefficient matrix at segment midpoints."""
    cpp_AIC3s_noTE(zetaSrc.ctypes.data_as(ct.POINTER(ct.c_double)),
            ct.byref(ct.c_uint(mSrc)),
            ct.byref(ct.c_uint(nSrc)),
            zetaTgt.ctypes.data_as(ct.POINTER(ct.c_double)),
            ct.byref(ct.c_uint(mTgt)),
            ct.byref(ct.c_uint(nTgt)),
            ct.byref(ct.c_bool(wakeSrc)),
            AIC3.ctypes.data_as(ct.POINTER(ct.c_double)))
       
def Cpp_AIC3noTE(zetaSrc,mSrc,nSrc,zetaTgt,mTgt,nTgt,isWakeSrc,AIC3):
    """@details Wrapper for cpp_wrap_AICnoTE.
    @param zetaSrc Source grid points.
    @param mSrc Chordwise panels.
    @param nSrc Spanwise.
    @param zetaTgt Target grid points.
    @param mTgt Chordwise panels.
    @param nTgt Spanwise.
    @param isWakeSrc.
    @return AIC3 3-component influence coefficient matrix."""
    cpp_AIC3noTE(zetaSrc.ctypes.data_as(ct.POINTER(ct.c_double)),
            ct.byref(ct.c_uint(mSrc)),
            ct.byref(ct.c_uint(nSrc)),
            zetaTgt.ctypes.data_as(ct.POINTER(ct.c_double)),
            ct.byref(ct.c_uint(mTgt)),
            ct.byref(ct.c_uint(nTgt)),
            ct.byref(ct.c_bool(isWakeSrc)),
            AIC3.ctypes.data_as(ct.POINTER(ct.c_double)))
    
def Cpp_dA3gamma0_dZeta(zetaSrc,mSrc,nSrc,gamma0,zetaTgt,mTgt,nTgt,dA3gam0):
    """@details Wrapper for cpp_wrap_dAgamma0_dZeta.
    @param zetaSrc Source grid points.
    @param mSrc Chordwise panels.
    @param nSrc Spanwise.
    @param gamma0 Reference circulation distribution.
    @param zetaTgt Target grid points.
    @param mTgt Chordwise panels.
    @param nTgt Spanwise.
    @return dAgam0 variation of influence coefficient matrix times gamma."""
    cpp_dA3gamma0_dZeta(zetaSrc.ctypes.data_as(ct.POINTER(ct.c_double)),
                       ct.byref(ct.c_uint(mSrc)),
                       ct.byref(ct.c_uint(nSrc)),
                       gamma0.ctypes.data_as(ct.POINTER(ct.c_double)),
                       zetaTgt.ctypes.data_as(ct.POINTER(ct.c_double)),
                       ct.byref(ct.c_uint(mTgt)),
                       ct.byref(ct.c_uint(nTgt)),
                       dA3gam0.ctypes.data_as(ct.POINTER(ct.c_double)))
    
def Cpp_Y1(vC,zeta,m,n,Y1):
    """@details Wrapper for cpp_wrap_Y1.
    @param vC Collocation point fluid-grid relative velocities.
    @param zeta Lattice vertices.
    @param m Chordwise panels.
    @param n Spanwise.
    @return Y1 Matrix (12K x K) transforming dGamma."""
    cpp_Y1(vC.ctypes.data_as(ct.POINTER(ct.c_double)),
           zeta.ctypes.data_as(ct.POINTER(ct.c_double)),
           ct.byref(ct.c_uint(m)),
           ct.byref(ct.c_uint(n)),
           Y1.ctypes.data_as(ct.POINTER(ct.c_double)))
    
def Cpp_Y2(gamma,vC,m,n,Y2):
    """@details Wrapper for cpp_wrap_Y2.
    @param gamma Panel circulation strengths.
    @param vC Collocation point fluid-grid relative velocities.
    @param m Chordwise panels.
    @param n Spanwise.
    @return Y2 Matrix (12K x 3K_zeta) transforming dZeta."""
    cpp_Y2(gamma.ctypes.data_as(ct.POINTER(ct.c_double)),
           vC.ctypes.data_as(ct.POINTER(ct.c_double)),
           ct.byref(ct.c_uint(m)),
           ct.byref(ct.c_uint(n)),
           Y2.ctypes.data_as(ct.POINTER(ct.c_double)))
    
def Cpp_Y3(gamma,zeta,m,n,Y3):
    """@details Wrapper for cpp_wrap_Y3.
    @param gamma Panel circulation strengths.
    @param zeta Lattice vertices.
    @param m Chordwise panels.
    @param n Spanwise.
    @return Y3 Matrix (12K x 3K) transforming coolocation velocities."""
    cpp_Y3(gamma.ctypes.data_as(ct.POINTER(ct.c_double)),
           zeta.ctypes.data_as(ct.POINTER(ct.c_double)),
           ct.byref(ct.c_uint(m)),
           ct.byref(ct.c_uint(n)),
           Y3.ctypes.data_as(ct.POINTER(ct.c_double)))
    
def Cpp_Y4(zeta,m,n,Y4):
    """@details Wrapper for cpp_wrap_Y4.
    @param zeta Lattice vertices.
    @param m Chordwise panels.
    @param n Spanwise.
    @return Y4 Matrix (3K x K) transforming dGammaPrime."""
    cpp_Y4(zeta.ctypes.data_as(ct.POINTER(ct.c_double)),
           ct.byref(ct.c_uint(m)),
           ct.byref(ct.c_uint(n)),
           Y4.ctypes.data_as(ct.POINTER(ct.c_double)))
    
def Cpp_Y5(gammaPri,zeta,m,n,Y5):
    """@details Wrapper for cpp_wrap_Y5.
    @param gammaPri Rate of change of circulation strengths.
    @param zeta Lattice vertices.
    @param m Chordwise panels.
    @param n Spanwise.
    @return Y5 Matrix (3K x 3K_zeta) transforming dZeta."""
    cpp_Y5(gammaPri.ctypes.data_as(ct.POINTER(ct.c_double)),
           zeta.ctypes.data_as(ct.POINTER(ct.c_double)),
           ct.byref(ct.c_uint(m)),
           ct.byref(ct.c_uint(n)),
           Y5.ctypes.data_as(ct.POINTER(ct.c_double)))
    
def Cpp_KJForces(zeta,gamma,zetaW,gammaW,zetaDot,nu,VMOPTS,gamma_tm1,f):
    """@details Wrapper force Joukowski method force calculations.
    @param zeta Lattice vertices.
    @param gamma Bound circulations.
    @param zetaW Wake lattice vertices.
    @param gammaW Wake circulations.
    @param zetaDot Bound lattice vertex velocities.
    @param nu External velocities.
    @param VMOPTS Derived type discretization and solver info.
    @return f forces at the lattice vertices.
    @note Dimensional variables used here. 
    """
    
    cpp_KJforces(zeta.ctypes.data_as(ct.POINTER(ct.c_double)),
                 gamma.ctypes.data_as(ct.POINTER(ct.c_double)),
                 zetaW.ctypes.data_as(ct.POINTER(ct.c_double)),
                 gammaW.ctypes.data_as(ct.POINTER(ct.c_double)),
                 zetaDot.ctypes.data_as(ct.POINTER(ct.c_double)),
                 nu.ctypes.data_as(ct.POINTER(ct.c_double)),
                 ct.byref(VMOPTS.M),
                 ct.byref(VMOPTS.N),
                 ct.byref(VMOPTS.ImageMethod),
                 ct.byref(VMOPTS.Mstar),
                 ct.byref(VMOPTS.Steady),
                 ct.byref(VMOPTS.KJMeth),
                 ct.byref(VMOPTS.NewAIC),
                 ct.byref(VMOPTS.DelTime),
                 ct.byref(VMOPTS.Rollup),
                 ct.byref(VMOPTS.NumCores),
                 gamma_tm1.ctypes.data_as(ct.POINTER(ct.c_double)),
                 f.ctypes.data_as(ct.POINTER(ct.c_double)) )
    
def Cpp_KJForces_vC(zeta,gamma,zetaW,gammaW,zetaDot,nu,VMOPTS,gamma_tm1,f):
    """@details Wrapper force Joukowski method force calculations.
    @param zeta Lattice vertices.
    @param gamma Bound circulations.
    @param zetaW Wake lattice vertices.
    @param gammaW Wake circulations.
    @param zetaDot Bound lattice vertex velocities.
    @param nu External velocities.
    @param VMOPTS Derived type discretization and solver info.
    @return f forces at the lattice vertices.
    @note Dimensional variables used here. 
    """
    
    cpp_KJforces_vC(zeta.ctypes.data_as(ct.POINTER(ct.c_double)),
                 gamma.ctypes.data_as(ct.POINTER(ct.c_double)),
                 zetaW.ctypes.data_as(ct.POINTER(ct.c_double)),
                 gammaW.ctypes.data_as(ct.POINTER(ct.c_double)),
                 zetaDot.ctypes.data_as(ct.POINTER(ct.c_double)),
                 nu.ctypes.data_as(ct.POINTER(ct.c_double)),
                 ct.byref(VMOPTS.M),
                 ct.byref(VMOPTS.N),
                 ct.byref(VMOPTS.ImageMethod),
                 ct.byref(VMOPTS.Mstar),
                 ct.byref(VMOPTS.Steady),
                 ct.byref(VMOPTS.KJMeth),
                 ct.byref(VMOPTS.NewAIC),
                 ct.byref(VMOPTS.DelTime),
                 ct.byref(VMOPTS.Rollup),
                 ct.byref(VMOPTS.NumCores),
                 gamma_tm1.ctypes.data_as(ct.POINTER(ct.c_double)),
                 f.ctypes.data_as(ct.POINTER(ct.c_double)) )
    
    
def Cpp_KJForces_vC_mod(zeta,gamma,zetaW,gammaW,zetaDot,nu,VMOPTS,gamma_tm1,f):
    """@details Wrapper force Joukowski method force calculations.
    @param zeta Lattice vertices.
    @param gamma Bound circulations.
    @param zetaW Wake lattice vertices.
    @param gammaW Wake circulations.
    @param zetaDot Bound lattice vertex velocities.
    @param nu External velocities.
    @param VMOPTS Derived type discretization and solver info.
    @return f forces at the lattice vertices.
    @note Dimensional variables used here. 
    """
    
    cpp_KJforces_vC_mod(zeta.ctypes.data_as(ct.POINTER(ct.c_double)),
                 gamma.ctypes.data_as(ct.POINTER(ct.c_double)),
                 zetaW.ctypes.data_as(ct.POINTER(ct.c_double)),
                 gammaW.ctypes.data_as(ct.POINTER(ct.c_double)),
                 zetaDot.ctypes.data_as(ct.POINTER(ct.c_double)),
                 nu.ctypes.data_as(ct.POINTER(ct.c_double)),
                 ct.byref(VMOPTS.M),
                 ct.byref(VMOPTS.N),
                 ct.byref(VMOPTS.ImageMethod),
                 ct.byref(VMOPTS.Mstar),
                 ct.byref(VMOPTS.Steady),
                 ct.byref(VMOPTS.KJMeth),
                 ct.byref(VMOPTS.NewAIC),
                 ct.byref(VMOPTS.DelTime),
                 ct.byref(VMOPTS.Rollup),
                 ct.byref(VMOPTS.NumCores),
                 gamma_tm1.ctypes.data_as(ct.POINTER(ct.c_double)),
                 f.ctypes.data_as(ct.POINTER(ct.c_double)) )
    
def Cpp_q_k(k,N,no):
    """@details Get vertex index from panel number and corner number.
    @param k Panel index (starts from zero).
    @param N Spanwise panels.
    @param no Corner number, 1-4.
    """
    return cpp_q_k(k,N,no)

def Cpp_dAs3gamma0_dZeta_num(zetaSrc,mSrc,nSrc,gamma0,zetaTgt,mTgt,nTgt,dAgam0):
    """@details Wrapper for cpp_wrap_dAs3gam0_dZeta.
    @param zetaSrc Source grid points.
    @param mSrc Chordwise panels.
    @param nSrc Spanwise.
    @param gamma0 Reference circulation distribution.
    @param zetaTgt Target grid points.
    @param mTgt Chordwise panels.
    @param nTgt Spanwise.
    @return dAgam0 variation of segment midpoint coefficient matrix times gamma."""
    cpp_dAs3gam0_dZeta_num(zetaSrc.ctypes.data_as(ct.POINTER(ct.c_double)),
                       ct.byref(ct.c_uint(mSrc)),
                       ct.byref(ct.c_uint(nSrc)),
                       gamma0.ctypes.data_as(ct.POINTER(ct.c_double)),
                       zetaTgt.ctypes.data_as(ct.POINTER(ct.c_double)),
                       ct.byref(ct.c_uint(mTgt)),
                       ct.byref(ct.c_uint(nTgt)),
                       dAgam0.ctypes.data_as(ct.POINTER(ct.c_double)))
    

if __name__ == '__main__':
    xP = np.array([0.0,0.0,-1.0], ct.c_double, order='C')
    x1 = np.array([-0.5,0.0,0.0], ct.c_double, order='C')
    x2 = np.array([0.5,0.0,0.0], ct.c_double, order='C')
    Gamma = 1.0
    Uind = np.array([0.0,0.0,0.0], ct.c_double, order='C')
    Cpp_Vorticity_BiotSegment_Map(xP, x1, x2, Gamma, Uind)
    print(Uind)