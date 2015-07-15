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

np.set_printoptions(precision = 3)

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
    E[0:m*n,0:m*n] = AIC
    E[0:m*n,m*n:m*n+mW*n] = AICw
    E[m*n:m*n+mW*n,m*n:m*n+mW*n] = np.eye(mW*n)
    E[m*n+mW*n:,0:m*n] = -np.eye(m*n)
    E[m*n+mW*n:,m*n+mW*n:] = delS*np.eye(m*n)
    
    # populate F
    Cgam = np.zeros((mW*n,m*n))
    Cgam[0:n,m*n-n:m*n] = np.eye(n)
    CgamW = np.zeros((mW*n,mW*n))
    CgamW[n:,0:mW*n-n] = np.eye(n*(mW-1))
    F[m*n:m*n+mW*n,0:m*n] = Cgam
    F[m*n:m*n+mW*n,m*n:m*n+mW*n] = CgamW
    F[m*n+mW*n:,0:m*n] = np.eye(m*n)
    
    # populate G
    W = np.zeros((m*n,3*(m+1)*(n+1)))
    Cpp_genW(zeta, m, n, W)
    dAgam_dZeta = np.zeros((m*n,3*(m+1)*(n+1)))
    Cpp_dAgamma0_dZeta(zeta, m, n, gam, zeta, m, n, dAgam_dZeta)
    dAwGamW_dZeta = np.zeros((m*n,3*(m+1)*(n+1)))
    Cpp_dAgamma0_dZeta(zetaW, mW, n, gamW, zeta, m, n, dAwGamW_dZeta)
    dWzetaPri_dZeta = np.zeros((m*n,3*(m+1)*(n+1)))
    Cpp_dWzetaPri0_dZeta(zeta, m, n, zetaPri, dWzetaPri_dZeta)
    G[0:m*n,0:3*(m+1)*(n+1)] = W
    G[0:m*n,3*(m+1)*(n+1):6*(m+1)*(n+1)] = - dAgam_dZeta - dAwGamW_dZeta + dWzetaPri_dZeta
    G[0:m*n,6*(m+1)*(n+1):] = -W
    
    # populate output matrices C and D
    Xi = np.zeros((3*m*n,3*(m+1)*(n+1)))
    H = np.zeros((3*(m+1)*(n+1),12*m*n))
    Y1 = np.zeros((12*m*n,m*n))
    Y2 = np.zeros((12*m*n,3*(m+1)*(n+1)))
    Y3 = np.zeros((12*m*n,3*m*n))
    Y4 = np.zeros((3*m*n,m*n))
    Y5 = np.zeros((3*m*n,3*(m+1)*(n+1)))
    AIC3 = np.zeros((3*m*n,m*n))
    AIC3w = np.zeros((3*m*n,mW*n))
    dA3gam_dZeta = np.zeros((3*m*n,3*(m+1)*(n+1)))
    dAw3gamW_dZeta = np.zeros((3*m*n,3*(m+1)*(n+1)))
    
    Cpp_genXi(m, n, 0.5, 0.5, Xi)
    Cpp_genH(m, n, H)
    Cpp_AIC3(zeta, m, n, zeta, m, n, AIC3)
    Cpp_AIC3(zetaW, mW, n, zeta, m, n, AIC3w)
    Cpp_dA3gamma0_dZeta(zeta, m, n, gam, zeta, m, n, dA3gam_dZeta)
    Cpp_dA3gamma0_dZeta(zetaW, mW, n, gamW, zeta, m, n, dAw3gamW_dZeta)
    
    # generate Y matrices
    vC0 = np.zeros((3*m*n)) # collocation fluid-grid relative velocities
    vC0[:] = np.dot(AIC3,gam) + np.dot(AIC3w,gamW) -np.dot(Xi,zetaPri)
    Cpp_Y1(vC0, zeta, m, n, Y1)
    Cpp_Y2(gam, vC0, m, n, Y2)
    Cpp_Y3(gam, zeta, m, n, Y3)
    Cpp_Y4(zeta, m, n, Y4)
    Cpp_Y5(gamPri, zeta, m, n, Y5)
    
    # Matrix C
    C[:,0:m*n] = np.dot(H, (Y1 - np.dot(Y3,AIC3)))
    C[:,m*n:m*n+mW*n] = np.dot(H, np.dot(Y3,AIC3w))
    C[:,m*n+mW*n:] = np.dot(np.transpose(Xi),Y4)
    
    # Matrix D
    D[:,0:3*(m+1)*(n+1)] = np.dot(H, np.dot(Y3, Xi))
    D[:,3*(m+1)*(n+1):6*(m+1)*(n+1)] = ( np.dot(H,Y2) 
        - np.dot(H, np.dot(Y3, dA3gam_dZeta + dAw3gamW_dZeta))
        + np.dot(np.transpose(Xi),Y5)                        )
    D[:,6*(m+1)*(n+1):] = -np.dot(H, np.dot(Y3, Xi))
    
    return E,F,G,C,D

if __name__ == '__main__':
    m=10
    n=1
    mW=100
    delS=1
    gam=np.ones((m*n))
    gamW=np.ones((mW*n))
    gamPri=np.ones((m*n))
    chords = np.linspace(0.0, 1.0, m+1, True)
    chordsW = np.linspace(1.0, 3.0, mW+1, True)
    spans = np.linspace(0.0, 1.0, n+1, True)
    zeta=np.zeros(3*len(chords)*len(spans))
    zetaW=np.zeros(3*len(chordsW)*len(spans))
    zetaPri = np.ones((3*len(chords)*len(spans)))
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
    E,F,G,C,D = genSSuvlm(gam, gamW, gamPri, zeta, zetaW, zetaPri, m, n, mW, delS)
    print("\n E matrix:\n",E)
    print("\n F matrix:\n",F)
    print("\n G matrix:\n",G)
    print("\n C matrix:\n",C)
    print("\n D[:,0:3K_\zeta] matrix:\n",D[:,0:3*(m+1)*(n+1)])
    print("\n D[:,3K_\zeta:6K_\zeta] matrix:\n",D[:,3*(m+1)*(n+1):6*(m+1)*(n+1)])
    print("\n D[:,6*K_\zeta:] matrix:\n",D[:,6*(m+1)*(n+1):])