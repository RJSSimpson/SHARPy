'''PyAero.UVLM.Utils.Linear
@brief      Generating linear UVLM models.
@author     Rob Simpson
@contact    r.simpson11@imperial.ac.uk 
@version    0.0
@date       09/07/2015
@pre        None
@warning    None
'''

import numpy as np
from UVLMLib import Cpp_AIC, Cpp_dAgamma0_dZeta, Cpp_genW, Cpp_dWzetaPri0_dZeta
from UVLMLib import Cpp_genH, Cpp_genXi, Cpp_AIC3, Cpp_dA3gamma0_dZeta, Cpp_Y1
from UVLMLib import Cpp_Y2, Cpp_Y3, Cpp_Y4, Cpp_Y5

def genSSuvlm(gam,gamW,gamPri,zeta,zetaW,zetaPri,m,n,mW,delS):
    """@details generate state-space matrices for linear UVLM.
    @param gam Reference circulation distribution on body.
    @param gamW Reference circulation distribution in the wake.
    @param gamPri Reference Rate of change of circulation on body.
    @param zeta Lattice vertices.
    @param zetaW Wake lattice vertices.
    @param zetaPri Vertex velocities.
    @param m Chordwise panels.
    @param n Spanwise.
    @param mW Chordwise panels in wake.
    @param delS Nondimensional timestep.
    @return E LHS state transfer matrix.
    @return F RHS state transfer matrix.
    @return G RHS Input matrix.
    @return C Output matrix.
    @return D Feedthrough matrix.
    @note All of the above should be nondimensional."""
    
    # init matrices
    E = np.zeros((2*m*n+mW*n,2*m*n+mW*n))
    F = np.zeros((2*m*n+mW*n,2*m*n+mW*n))
    G = np.zeros((2*m*n+mW*n,9*(m+1)*(n+1)))
    C = np.zeros((3*(m+1)*(n+1),2*m*n+mW*n))
    D = np.zeros((3*(m+1)*(n+1),9*(m+1)*(n+1)))
    
    # populate E
    AIC = np.zeros((m*n,m*n))
    Cpp_AIC(zeta, m, n, zeta, m, n, AIC)
    AICw= np.zeros((m*n,mW*n))
    Cpp_AIC(zetaW, mW, n, zeta, m, n, AICw)
    eye1 = np.eye(mW*n)
    eye2 = np.eye(m*n)
    E[0:m*n,0:m*n] = AIC
    E[0:m*n,m*n:m*n+mW*n] = AICw
    E[m*n:m*n+mW*n,m*n:m*n+mW*n] = eye1
    E[m*n+mW*n:,0:m*n] = -eye2
    E[m*n+mW*n:,m*n+mW*n:] = delS*eye1
    
    #populate F
    #TODO
    
    return E,F,G,C,D

if __name__ == '__main__':
    m=1
    n=1
    mW=1
    delS=1
    gam=np.ones((m*n))
    gamW=np.ones((mW*n))
    gamPri=np.ones((m*n))
    chords = np.linspace(0.0, 1.0, m+1, True)
    chordsW = np.linspace(1.0, 2.0, m+1, True)
    spans = np.linspace(0.0, 1.0, n+1, True)
    zeta=np.zeros(3*len(chords)*len(spans))
    zetaW=np.zeros(3*len(chordsW)*len(spans))
    kk=0
    for c in chords:
        for s in spans:
            zeta[3*kk]=c
            zeta[3*kk+1]=s
            kk=kk+1
        # end for s
    # end for c
    kk=0
    for c in chordsW:
        for s in spans:
            zetaW[3*kk]=c
            zetaW[3*kk+1]=s
            kk=kk+1
        # end for s
    # end for c
    zetaPri = np.ones((3*len(chords)*len(spans)))
    E,F,G,C,D = genSSuvlm(gam, gamW, gamPri, zeta, zetaW, zetaPri, m, n, mW, delS)
    print(E)