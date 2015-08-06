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
from UVLMLib import Cpp_AIC, Cpp_dAgamma0_dZeta, Cpp_genW, Cpp_dWzetaPri0_dZeta,\
    Cpp_dAs3gamma0_dZeta_num
from UVLMLib import Cpp_genH, Cpp_genXi, Cpp_AIC3, Cpp_dA3gamma0_dZeta, Cpp_Y1
from UVLMLib import Cpp_Y2, Cpp_Y3, Cpp_Y4, Cpp_Y5, Cpp_AIC3noTE, Cpp_AIC3s
from scipy.io.matlab.mio import savemat
import SharPySettings as Settings
import getpass

np.set_printoptions(precision = 4)

def genSSuvlm(gam,gamW,gamPri,zeta,zetaW,zetaPri,nu,m,n,mW,delS):
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
    F[m*n+mW*n:,0:m*n] = -np.eye(m*n)
    
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
    Y3 = np.zeros((12*m*n,12*m*n))
    Y4 = np.zeros((3*m*n,m*n))
    Y5 = np.zeros((3*m*n,3*(m+1)*(n+1)))
    AICs3 = np.zeros((12*m*n,m*n))
    AICs3w = np.zeros((12*m*n,mW*n))
    dAs3gam_dZeta = np.zeros((12*m*n,3*(m+1)*(n+1)))
    dAs3wGamW_dZeta = np.zeros((12*m*n,3*(m+1)*(n+1)))
    
    # gen interpolation and AIC matrices
    Cpp_genXi(m, n, 0.5, 0.5, Xi)
    Cpp_genH(m, n, H)
    Cpp_AIC3s(zeta, m, n, zeta, m, n, AICs3)
    Cpp_AIC3s(zetaW, mW, n, zeta, m, n, AICs3w)
    Cpp_dAs3gamma0_dZeta_num(zeta, m, n, gam, zeta, m, n, dAs3gam_dZeta)
    Cpp_dAs3gamma0_dZeta_num(zetaW, mW, n, gamW, zeta, m, n, dAs3wGamW_dZeta)
     
    # generate Y matrices
    vM0 = np.zeros((12*m*n)) # collocation fluid-grid relative velocities
    vM0[:] = np.dot(AICs3,gam) + np.dot(AICs3w,gamW) + np.dot(H.transpose(),nu) - np.dot(H.transpose(),zetaPri)
    Cpp_Y1(vM0, zeta, m, n, Y1)
    Cpp_Y2(gam, vM0, m, n, Y2)
    Cpp_Y3(gam, zeta, m, n, Y3)
    Cpp_Y4(zeta, m, n, Y4)
    Cpp_Y5(gamPri, zeta, m, n, Y5)
     
    # Matrix C
    C[:,0:m*n] = np.dot(H, Y1) - np.dot(H , np.dot(Y3,AICs3))
    C[:,m*n:m*n+mW*n] = -np.dot(H, np.dot(Y3,AICs3w))
    C[:,m*n+mW*n:] = np.dot(np.transpose(Xi),Y4)
     
    # Matrix D
    D[:,0:3*(m+1)*(n+1)] = np.dot(H, np.dot(Y3, H.transpose()))
    D[:,3*(m+1)*(n+1):6*(m+1)*(n+1)] = ( np.dot(H,Y2) 
        - np.dot(H, np.dot(Y3, dAs3gam_dZeta + dAs3wGamW_dZeta))
        + np.dot(np.transpose(Xi),Y5)                        )
    D[:,6*(m+1)*(n+1):] = -np.dot(H, np.dot(Y3, H.transpose()))
    
    return E,F,G,C,D

def genLinearAerofoil(m,mW,e,f,writeToMat = False):
    """@brief Generate linear model of aerofoil.
    @param m Chordwise panels.
    @param mW Chordwise panels in wake.
    @param delS Non-dim time step.
    @param e location of pitch axis aft of LE [0,1].
    @param f location of flap hinge aft of LE [0,1].
    @param writeToMat write to file in Settings.OutputDir.
    """
    
    # infer delS from body discretization
    delS = 2.0/m
    
    # initialise states and inputs
    gam=np.zeros((m))
    gamW=np.zeros((mW))
    gamPri=np.zeros((m))
    chords = np.linspace(0.0, 1.0, m+1, True)
    chordsW = np.linspace(1.0, 1.0+mW*delS/2.0, mW+1, True)
    spans = np.linspace(-1000.0, 1000.0, 2, True)
    zeta=np.zeros(3*len(chords)*len(spans))
    zetaW=np.zeros(3*len(chordsW)*len(spans))
    zetaPri = np.zeros((3*len(chords)*len(spans)))
    zetaPri[0::3] = -1.0;
    nu = np.zeros_like(zetaPri)
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
    
    # generate model
    E,F,G,C,D = genSSuvlm(gam, gamW, gamPri, zeta, zetaW, zetaPri, nu, m, 1, mW, delS)
    
    # convert inputs from general kinematics to aerofil DoFs
    T = np.zeros((9*(m+1)*2,5))
    for i in range(m+1):
        # alpha, alphaPrime
        T[3*(m+1)*2+6*i+2,0] = -(zeta[6*i]+0.25/m-e)
        T[3*(m+1)*2+6*i+5,0] = -(zeta[6*i+3]+0.25/m-e)
        T[6*i+2,1] = -(zeta[6*i]+0.25/m-e)*4.0
        T[6*i+5,1] = -(zeta[6*i+3]+0.25/m-e)*4.0
        # plunge
        T[6*i+2,2] = -1
        T[6*i+5,2] = -1
        # beta, betaPrime
        if (zeta[6*i] + 0.25/m) > f:
            T[3*(m+1)*2+6*i+2,3] = -(zeta[6*i]+0.25/m-f)
            T[3*(m+1)*2+6*i+5,3] = -(zeta[6*i+3]+0.25/m-f)
            T[6*i+2,4] = -(zeta[6*i]+0.25/m-f)*4.0
            T[6*i+5,4] = -(zeta[6*i+3]+0.25/m-f)*4.0
    
    G_s = np.dot(G,T)
    D_s = np.dot(D,T)
    
    # get coefficients as output
    T_coeff = np.zeros((2,3*(m+1)*2))
    T_coeff[0,0::3] = 1.0 #drag
    T_coeff[1,2::3] = 1.0 #lift
    C_coeff = np.dot(T_coeff,C)
    D_coeff = np.dot(T_coeff,D)
    D_s_coeff = np.dot(T_coeff,D_s)

    if writeToMat == True:
        fileName = Settings.OutputDir + 'aerofoil_m' + str(m) + 'mW' + str(mW) + 'delS' + str(delS)
        savemat(fileName,
                {'E':E, 'F':F, 'G':G, 'C':C, 'D':D, 'm':m, 'mW':mW, 'delS':delS,
                 'G_s':G_s, 'D_s':D_s,
                 'C_coeff':C_coeff, 'D_coeff':D_coeff, 'D_s_coeff':D_s_coeff},
                True)
    # end if
    
    return E,F,G,C,D

if __name__ == '__main__':
    Settings.OutputDir = '/home/' + getpass.getuser() + '/Documents/MATLAB/newUVLM/aerofoil/'
    m=10
    e=0.25
    f=0.75
    for mW in (10*m,):
        genLinearAerofoil(m,mW,e,f,writeToMat = True)