'''@package PyFSI.Beam2Aero_Lin
@brief      Generate matrices that transform linearized structural DoFs into
             linearized grid velocities.
@author     Rob Simpson
@contact    r.simpson11@imperial.ac.uk 
@version    0.0
@date       03/07/2014
@pre        None
@warning    None
'''

import numpy as np
import PyBeam.Utils.XbeamLib as xbl

def zetaDotSubMat(r,psi,rDot,psiDot,v_a,omega_a,psi_G2a,xi,xiDot,iLin,jLin,
                    nBeam,nSection):
    """@brief create submatrix(3,N_states) for the linearized velocity at any
    point on the surface, zeta, defined by beam DoFs R and Psi, and the
    cross-sectional coordinate xi.
    @param r Position vector on beam, defined in the a-frame.
    @param psi CRV of orientation of beam cross-section.
    @param rDot Velocity of beam at r, defined in the a-frame.
    @param psiDot Time rate of change of psi.
    @param v_a Velocity of a-frame, defined in the a-frame.
    @param omega_a Angular vel of a-frame, defined in a-frame.
    @param psi_G2a CRV of orientation of a-frame relative to earth.
    @param xi Cross-sectional coordinate, defined in B-frame.
    @param xiDot Time rate of change of xi, defined in B-frame.
    @param iLin index of xi on the sectional DoF used for the linear analysis.
    @param jLin index of R on the beam DoF used for the linear analysis.
    @param nBeam Number of points on the beam axis.
    @param nSection Number of points in the section definition.
    
    @returns subMat Linear transformation from FBD + sectional DoFs to surface
              velocity (G-frame) at point r_0 + C_Ga*(R + C(Psi)*xi).
    @details This routine calculates the matrix corresponding to any point
              on the beam axis and corresponding cross-section.
              Interpolation of the primary beam and section DoFs may pre-
              or post- multiply this matrix in an aeroelastic analysis.
    """
    # Total states on RHS of matrix
    nStates = 2*6*nBeam + 9 + 2*3*nSection
    # Initialise subMat
    subMat = np.zeros((3,nStates))
    
    # term 1: \delta(C^{Ga}v_a)
    cGa = xbl.Psi2TransMat(psi_G2a)
    skewVa = xbl.Skew(v_a)
    tangGa = xbl.Tangential(psi_G2a);
    # term for \delta \psi_{Ga}
    subMat[:,2*6*nBeam+6:2*6*nBeam+9:1] += -np.dot(cGa,np.dot(skewVa,tangGa))
    # term for \delta v_a
    subMat[:,2*6*nBeam:2*6*nBeam+3:1] += cGa
    
    # term 2: \delta(C^{Ga} \tilde{\Omega_a_a} R_a)
    skewOmAa = xbl.Skew(omega_a)
    skewRa = xbl.Skew(r)
    # term for \delta \psi_{Ga}
    subMat[:,2*6*nBeam+6:2*6*nBeam+9:1] += -np.dot(cGa,
                                            np.dot(xbl.Skew(np.dot(skewOmAa,r)),
                                            tangGa))
    # term for \delta \Omega_a_a
    subMat[:,2*6*nBeam+3:2*6*nBeam+6:1] += -np.dot(cGa,skewRa)
    # term for \delta R_a
    subMat[:,6*jLin:6*jLin+3:1] += np.dot(cGa,skewOmAa)
    
    # term 3: \delta(C^{Ga}\tilde{\Omega_a_a}C^{aB}\xi_B)
    cAb = xbl.Psi2TransMat(psi)
    skewXi = xbl.Skew(xi)
    tangAb = xbl.Tangential(psi)
    # term for \delta \psi_{Ga}
    subMat[:,2*6*nBeam+6:2*6*nBeam+9:1] += -np.dot(cGa,
                                            np.dot(xbl.Skew(np.dot(skewOmAa,np.dot(cAb,xi))),
                                            tangGa))
    # term for \delta \Omega_a_a
    subMat[:,2*6*nBeam+3:2*6*nBeam+6:1] += -np.dot(cGa,xbl.Skew(np.dot(cAb,xi)))
    # term for \delta \psi_{aB}
    subMat[:,6*jLin+3:6*jLin+6:1] += -np.dot(cGa,
                                      np.dot(skewOmAa,
                                      np.dot(cAb,
                                      np.dot(skewXi,
                                      tangAb))))
    # term for \delta xi_B
    subMat[:,2*6*nBeam+9+3*iLin:2*6*nBeam+9+3*iLin+3:1] += np.dot(cGa,
                                                           np.dot(skewOmAa,
                                                           cAb))
    
    # term 4: \delta(C^{Ga}\dot{R}_a)
    skewRdot = xbl.Skew(rDot)
    # term for \delta \psi_{Ga}
    subMat[:,2*6*nBeam+6:2*6*nBeam+9:1] += -np.dot(cGa,
                                            np.dot(skewRdot,
                                            tangGa))
    # term for \delta \dot{R}_a
    subMat[:,6*nBeam+6*jLin:6*nBeam+6*jLin+3:1] += cGa
    
    # term 5: \delta(C^{Ga}C^{aB}\dot{\xi}_B)
    skewXiDot = xbl.Skew(xiDot)
    # term for \delta \psi_{Ga}
    subMat[:,2*6*nBeam+6:2*6*nBeam+9:1] += -np.dot(cGa,
                                            np.dot(xbl.Skew(np.dot(cAb,xiDot)),
                                            tangGa))
    # term for \delta \psi_{aB}
    subMat[:,6*jLin+3:6*jLin+6:1] += -np.dot(cGa,
                                      np.dot(cAb,
                                      np.dot(skewXiDot,
                                      tangAb)))
    # term for \delta \dot{\xi}_B
    subMat[:,2*6*nBeam+9+3*nSection+3*iLin:2*6*nBeam+9+3*nSection+3*iLin+3] += np.dot(cGa,cAb)
    
    # term 6: \delta(C^{Ga}C^{aB}\tilde{T(\psi_{aB})\dot{\psi}_{aB}\xi_B)
    # Note: eqn. rearranged to \delta(-C^{Ga}C^{aB}\tilde{\xi}T(\psi_{aB})\dot{\psi}_{aB})
    # term for  \delta \psi_{Ga}
    subMat[:,2*6*nBeam+6:2*6*nBeam+9:1] += np.dot(cGa,
                                           np.dot(xbl.Skew(
                                                  np.dot(cAb,
                                                  np.dot(skewXi,
                                                  np.dot(tangAb,
                                                  psiDot)))),
                                           tangGa))
    # term for \delta \psi_{aB}
    subMat[:,6*jLin+3:6*jLin+6:1] += np.dot(cGa,
                                     np.dot(cAb,
                                     np.dot(xbl.Skew(
                                            np.dot(skewXi,
                                            np.dot(tangAb,
                                            psiDot))),
                                     tangAb)))
    # term for \delta \xi_B
    subMat[:,2*6*nBeam+9+3*iLin:2*6*nBeam+9+3*iLin+3:1] += np.dot(cGa,
                                                           np.dot(cAb,
                                                           xbl.Skew(
                                                               np.dot(tangAb,
                                                               psiDot))))
    # term for \delta T(\psi_{aB})\dot{\psi_{aB}}
    subMat[:,6*jLin+3:6*jLin+6:1] += -np.dot(cGa,
                                      np.dot(cAb,
                                      np.dot(skewXi,
                                      xbl.a1(psi,psiDot))))
    # term for \delta \dot{\psi}_{aB}
    subMat[:,6*nBeam+6*jLin+3:6*nBeam+6*jLin+6:1] += -np.dot(cGa,
                                                      np.dot(cAb,
                                                      np.dot(skewXi,
                                                      tangAb)))
    return subMat

def printBlocks(mat,iBlock,jBlock):
    """@brief Prints matrices in specified block format.
    @param mat Matrix to print, numpy.array type.
    @param iBlock Size of block in i index.
    @param jBlock Size of block in j index.
    """
    
    nI,nJ = np.shape(mat)
    i=0
    j=0
    while i < nI:
        while j < nJ:
            print('block [%d:%d,%d:%d]:'%(i,i+iBlock,j,j+jBlock))
            print(mat[i:i+iBlock,j:j+jBlock],'\n')
            j+=jBlock
        # end while j
        i+=iBlock
    # end while i
      

if __name__ == '__main__':
    # test zetaDotSubMat
    r=np.array([1, 0, 0])
    psi=np.array([np.pi/4.0, np.pi/4.0, 0.0])
    rDot=np.array([1, 0, 0])
    psiDot=np.array([np.pi/4.0, np.pi/4.0, 0.0])
    v_a=np.array([0, -10.0, 0])
    omega_a=np.array([0, 0, np.pi/4.0])
    psi_G2a=np.array([np.pi/4.0, 0, 0])
    xi=np.array([0, -1, 0])
    xiDot=np.array([0, 0, -0.1])
    iLin=0
    jLin=0
    nBeam=1
    nSection=1
    zetaDotSubMat(r,psi,rDot,psiDot,v_a,omega_a,psi_G2a,xi,xiDot,iLin,jLin,
                    nBeam,nSection)
    