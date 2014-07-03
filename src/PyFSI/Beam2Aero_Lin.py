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
              velocity at point r_0 + C_Ga*(R + C(Psi)*xi).
    @details This routine calculates the matrix corresponding to any point
              on the beam axis and corresponding cross-section.
              Interpolation of the primary beam and section DoFs may pre-
              or post- multiply this matrix in an aeroelastic analysis.
    """
    
    nStates = 2*6*nBeam + 9 + 2*3*nSection
    # Initialise subMat
    subMat = np.zeros(3,nStates)
    
    return subMat

if __name__ == '__main__':
    pass