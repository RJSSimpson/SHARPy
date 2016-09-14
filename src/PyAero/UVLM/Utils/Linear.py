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
from XbeamLib import Skew
from PyMPC.ssdiscrete import StateSpace
from scipy.signal import dlsim

np.set_printoptions(precision = 4)

def genSSuvlm(gam,gamW,gamPri,zeta,zetaW,zetaPri,nu,m,n,mW,delS,imageMeth=False):
    """@details generate state-space matrices for linear UVLM.
    @param gam Reference circulation distribution on body.
    @param gamW Reference circulation distribution in the wake.
    @param gamPri Reference Rate of change of circulation on body.
    @param zeta Lattice vertices.
    @param zetaW Wake lattice vertices.
    @param zetaPri Vertex velocities.
    @param nu Atmospheric velocities.
    @param m Chordwise panels.
    @param n Spanwise.
    @param mW Chordwise panels in wake.
    @param delS Nondimensional timestep.
    @param imageMeth Use image method across xz-plane.
    @return E LHS state transfer matrix.
    @return F RHS state transfer matrix.
    @return G RHS Input matrix.
    @return C Output matrix.
    @return D Feedthrough matrix.
    @note All of the above should be nondimensional."""
    
    midPoint = False# use midpoint approximation for derivative of gamma 
    
    # init matrices
    E = np.zeros((2*m*n+mW*n,2*m*n+mW*n))
    F = np.zeros((2*m*n+mW*n,2*m*n+mW*n))
    G = np.zeros((2*m*n+mW*n,9*(m+1)*(n+1)))
    C = np.zeros((3*(m+1)*(n+1),2*m*n+mW*n))
    D = np.zeros((3*(m+1)*(n+1),9*(m+1)*(n+1)))
    
    # populate E
    AIC = np.zeros((m*n,m*n))
    Cpp_AIC(zeta, m, n, zeta, m, n, AIC, imageMeth)
    AICw= np.zeros((m*n,mW*n))
    Cpp_AIC(zetaW, mW, n, zeta, m, n, AICw, imageMeth)
    E[0:m*n,0:m*n] = AIC
    E[0:m*n,m*n:m*n+mW*n] = AICw
    E[m*n:m*n+mW*n,m*n:m*n+mW*n] = np.eye(mW*n)
    if midPoint == True:
        E[m*n+mW*n:,0:m*n] = np.eye(m*n)
        E[m*n+mW*n:,m*n+mW*n:] = -0.5*delS*np.eye(m*n)
    else:
        E[m*n+mW*n:,0:m*n] = -np.eye(m*n)
        E[m*n+mW*n:,m*n+mW*n:] = delS*np.eye(m*n)
    
    # populate F
    Cgam = np.zeros((mW*n,m*n))
    Cgam[0:n,m*n-n:m*n] = np.eye(n)
    CgamW = np.zeros((mW*n,mW*n))
    CgamW[n:,0:mW*n-n] = np.eye(n*(mW-1))
    #CgamW[-1:,-1]=0.99975
    F[m*n:m*n+mW*n,0:m*n] = Cgam
    F[m*n:m*n+mW*n,m*n:m*n+mW*n] = CgamW
    if midPoint == True:
        F[m*n+mW*n:,0:m*n] = np.eye(m*n)
        F[m*n+mW*n:,m*n+mW*n:] = 0.5*delS*np.eye(m*n)
    else:
        F[m*n+mW*n:,0:m*n] = -np.eye(m*n)
    
    # populate G
    W = np.zeros((m*n,3*(m+1)*(n+1)))
    Cpp_genW(zeta, m, n, W)
    dAgam_dZeta = np.zeros((m*n,3*(m+1)*(n+1)))
    Cpp_dAgamma0_dZeta(zeta, m, n, gam, zeta, m, n, dAgam_dZeta, imageMeth)
    dAwGamW_dZeta = np.zeros((m*n,3*(m+1)*(n+1)))
    Cpp_dAgamma0_dZeta(zetaW, mW, n, gamW, zeta, m, n, dAwGamW_dZeta, imageMeth)
    dWzetaPri_dZeta = np.zeros((m*n,3*(m+1)*(n+1)))
    dWnuPri_dZeta = np.zeros((m*n,3*(m+1)*(n+1)))
    Cpp_dWzetaPri0_dZeta(zeta, m, n, zetaPri, dWzetaPri_dZeta)
    Cpp_dWzetaPri0_dZeta(zeta, m, n, nu, dWnuPri_dZeta)
    G[0:m*n,0:3*(m+1)*(n+1)] = 2*W
    G[0:m*n,3*(m+1)*(n+1):6*(m+1)*(n+1)] = - dAgam_dZeta - dAwGamW_dZeta + 2*dWzetaPri_dZeta - 2*dWnuPri_dZeta
    G[0:m*n,6*(m+1)*(n+1):] = -2*W
    
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
    Cpp_AIC3s(zeta, m, n, zeta, m, n, AICs3, imageMeth)
    Cpp_AIC3s(zetaW, mW, n, zeta, m, n, AICs3w, imageMeth)
    Cpp_dAs3gamma0_dZeta_num(zeta, m, n, gam, zeta, m, n, dAs3gam_dZeta, imageMeth)
    Cpp_dAs3gamma0_dZeta_num(zetaW, mW, n, gamW, zeta, m, n, dAs3wGamW_dZeta, imageMeth)
     
    # generate Y matrices
    vM0 = np.zeros((12*m*n)) # collocation fluid-grid relative velocities
    vM0[:] = np.dot(AICs3,gam) + np.dot(AICs3w,gamW) + 2.0*np.dot(H.transpose(),nu) - 2.0*np.dot(H.transpose(),zetaPri)
    Cpp_Y1(vM0, zeta, m, n, Y1)
    Cpp_Y2(gam, vM0, m, n, Y2)
    Cpp_Y3(gam, zeta, m, n, Y3)
    Cpp_Y4(zeta, m, n, Y4)
    Cpp_Y5(gamPri, zeta, m, n, Y5)
     
    # Matrix C
    C[:,0:m*n] = np.dot(H, Y1) - np.dot(H , np.dot(Y3,AICs3))
    C[:,m*n:m*n+mW*n] = -np.dot(H, np.dot(Y3,AICs3w))
    #C[:,m*n+mW*n:] = 2.0*np.dot(np.transpose(Xi),Y4)
    
    # Matrix D
    D[:,0:3*(m+1)*(n+1)] = 2.0*np.dot(H, np.dot(Y3, H.transpose()))
    D[:,3*(m+1)*(n+1):6*(m+1)*(n+1)] = ( np.dot(H,Y2) 
        - np.dot(H, np.dot(Y3, dAs3gam_dZeta + dAs3wGamW_dZeta))
        + np.dot(np.transpose(Xi),Y5)                          )
    D[:,6*(m+1)*(n+1):] = -2.0*np.dot(H, np.dot(Y3, H.transpose()))
    
    return E,F,G,C,D

def genLinearAerofoil(m,mW,writeToMat = False,e=0.25,f=0.75):
    """@brief Generate linear model of aerofoil.
    @param m Chordwise panels.
    @param mW Chordwise panels in wake.
    @param delS Non-dim time step.
    @param e location of pitch axis aft of LE [0,1], default 0.25.
    @param f location of flap hinge aft of LE [0,1], default 0.75.
    @param writeToMat write to file in Settings.OutputDir.
    """
    
    n=1 # number of spanwise panels
    
    # infer delS from body discretization
    delS = 2.0/m
    
    # hack
#     factor = 4 
#     mW=factor*mW
#     delS = delS/float(factor)
    
    # initialise states and inputs
    gam=np.zeros((m*n))
    gamW=np.zeros((mW*n))
    gamPri=np.zeros((m*n))
    chords = np.linspace(0.0, 1.0, m+1, True)
    chordsW = np.linspace(1.0, 1.0+mW*delS/2.0, mW+1, True)
    spans = np.linspace(0.0, 2000.0, n+1, True)
    zeta=np.zeros(3*len(chords)*len(spans))
    zetaW=np.zeros(3*len(chordsW)*len(spans))
    zetaPri = np.zeros((3*len(chords)*len(spans)))
    #zetaPri[0::3] = -1.0
    nu = np.zeros_like(zetaPri)
    nu[0::3] = 0.5
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
    E,F,G,C,D = genSSuvlm(gam, gamW, gamPri, zeta, zetaW, zetaPri, nu, m, n, mW, delS, True)
    
    # convert inputs from general kinematics to aerofil DoFs
    T = np.zeros((9*(m+1)*(n+1),5))
    for i in range(m+1):
        for j in range(n+1):
            q=i*(n+1)+j
            # alpha, alphaPrime
            T[3*(m+1)*(n+1)+3*q+2,0] = -(zeta[3*q]+0.25/m-e)
            T[3*q+2,1] = -(zeta[3*q]+0.25/m-e)
            # plunge
            T[3*q+2,2] = -1
            # beta, betaPrime
            if zeta[3*q]+0.25/m > f:
                T[3*(m+1)*(n+1)+3*q+2,3] = -(zeta[3*q]+0.25/m-f)
                T[3*q+2,4] = -(zeta[3*q]+0.25/m-f)
                
                
    G_s = np.dot(G,T)
    D_s = np.dot(D,T)
    
    # get coefficients as output
    T_coeff = np.zeros((3,3*(m+1)*(n+1)))
    T_coeff[0,0::3] = 1.0 #drag
    T_coeff[1,2::3] = 1.0 #lift
    for i in range(m+1):
        for j in range(n+1):
            q=i*(n+1)+j
            # moment = r*L, +ve nose-up, at quarter chord
            T_coeff[2,3*q+2] = -(zeta[3*q]+0.25/m-0.25)
    
    C_coeff = np.dot(T_coeff,C)
    D_coeff = np.dot(T_coeff,D)
    D_s_coeff = np.dot(T_coeff,D_s)

    if writeToMat == True:
        fileName = Settings.OutputDir + 'TESTaerofoil_m' + str(m) + 'mW' + str(mW) + 'delS' + str(delS)
        if e != 0.25:
            fileName += 'e'+str(e)
        if f != 0.75:
            fileName += 'f'+str(f)
        savemat(fileName,
                {'E':E, 'F':F, 'G':G, 'C':C, 'D':D, 'm':m, 'mW':mW, 'delS':delS,
                 'G_s':G_s, 'D_s':D_s,
                 'C_coeff':C_coeff, 'D_coeff':D_coeff, 'D_s_coeff':D_s_coeff,
                 'T_coeff':T_coeff},
                True)
    # end if
    
    return E,F,G,C,D,delS

def genLinearRectWing(AR,m,mW,n,e=0.25,f=0.75,writeToMat = False, imageMeth = False):
    """@brief Generate linear model of rectangular wing.
    @param AR Aspect ration
    @param m Chordwise panels.
    @param mW Chordwise panels in wake.
    @param n Spanwise panels.
    @param delS Non-dim time step.
    @param e location of pitch axis aft of LE [0,1], default 0.25.
    @param f location of flap hinge aft of LE [0,1], default 0.75.
    @param writeToMat write to file in Settings.OutputDir.
    @param imageMeth Use image method across xz-plane.
    """
    
    # infer delS from body discretization
    delS = 2.0/m
    
    # initialise states and inputs
    gam=np.zeros((m*n))
    gamW=np.zeros((mW*n))
    gamPri=np.zeros((m*n))
    chords = np.linspace(0.0, 1.0, m+1, True)
    chordsW = np.linspace(1.0, 1.0+mW*delS/2.0, mW+1, True)
    if imageMeth:
        spans = np.linspace(0.0, AR/2.0, n+1, True)
    else:
        spans = np.linspace(-AR/2.0, AR/2.0, n+1, True)
    # end if
    zeta=np.zeros(3*len(chords)*len(spans))
    zetaW=np.zeros(3*len(chordsW)*len(spans))
    zetaPri = np.zeros((3*len(chords)*len(spans)))
    zetaPri[0::3] = -0.5
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
    E,F,G,C,D = genSSuvlm(gam, gamW, gamPri, zeta, zetaW, zetaPri, nu, m, n, mW, delS, imageMeth)
    
    # convert inputs from general kinematics to aerofil DoFs
    T = np.zeros((9*(m+1)*(n+1),5))
    for i in range(m+1):
        for j in range(n+1):
            q=i*(n+1)+j
            # alpha, alphaPrime
            T[3*(m+1)*(n+1)+3*q+2,0] = -(zeta[3*q]+0.25/m-e)
            T[3*q+2,1] = -(zeta[3*q]+0.25/m-e)
            # plunge
            T[3*q+2,2] = -1
            # beta, betaPrime
            if zeta[3*q]+0.25/m > f:
                T[3*(m+1)*(n+1)+3*q+2,3] = -(zeta[3*q]+0.25/m-f)
                T[3*q+2,4] = -(zeta[3*q]+0.25/m-f)
                
                
    G_s = np.dot(G,T)
    D_s = np.dot(D,T)
    
    # get coefficients as output
    T_coeff = np.zeros((3,3*(m+1)*(n+1)))
    T_coeff[0,0::3] = 1.0 #drag
    T_coeff[1,2::3] = 1.0 #lift
    for i in range(m+1):
        for j in range(n+1):
            q=i*(n+1)+j
            # moment = r*L, +ve nose-up, at quarter chord
            T_coeff[2,3*q+2] = -(zeta[3*q]+0.25/m-0.25)
    
    C_coeff = np.dot(T_coeff,C)
    D_coeff = np.dot(T_coeff,D)
    D_s_coeff = np.dot(T_coeff,D_s)
    
    # spanwise lift distribution as output
    T_span=np.zeros((n+1,3*(m+1)*(n+1)))
    for jj in range(n+1):
        T_span[jj,3*jj+2::3*(n+1)] = 1.0

    if writeToMat == True:
        fileName = Settings.OutputDir + 'rectWingAR' + str(AR) + '_m' + str(m) + 'mW' + str(mW) + 'n' + str(n) + 'delS' + str(delS)
        if e != 0.25:
            fileName += 'e'+str(e)
        if f != 0.75:
            fileName += 'f'+str(f)
        if imageMeth != False:
            fileName += 'half'
        savemat(fileName,
                {'E':E, 'F':F, 'G':G, 'C':C, 'D':D, 'm':m, 'mW':mW, 'delS':delS,
                 'G_s':G_s, 'D_s':D_s,
                 'C_coeff':C_coeff, 'D_coeff':D_coeff, 'D_s_coeff':D_s_coeff,
                 'T_coeff':T_coeff, 'T_span':T_span,
                 'AR':AR, 'm':m, 'mW':mW, 'n':n, 'zeta':zeta},
                True)
    # end if
    
    return E,F,G,C,D

def nln2linStates(Zeta, ZetaStar, Gamma, GammaStar, Uext, M, N, mW, chord, vRef = 1.0):
    """@details generate state-space matrices for linear UVLM.
    @param Zeta UVLM surface coordinates in a-frame (dimensional).
    @param ZetaStar UVLM wake coordinates in a-frame (dimensional).
    @param Gamma UVLM surface circulation (dimensional).
    @param GammaStar UVLM wake circulation (dimensional).
    @param mW number of wake panels.
    @return gam
    @return gamW
    @return gamPri
    @return zeta
    @return zetaW
    @return zetaPri
    @return nu
    @note All of the return values are nondimensional and in the coordinate system for the linearized aerodynamics."""
    
    # circulation dist
    gam = Gamma.flatten()/(chord*vRef)
    gamW = np.zeros(mW*N)
    for i in range(mW):
        gamW[i*N:(i+1)*N]=GammaStar/(chord*vRef)
        
    # geometry/external vels
    zeta=np.zeros((M+1,N+1,3))
    nu=np.zeros((M+1,N+1,3))
    beam2aero = np.zeros((3,3))
    beam2aero[0,1]=-1.0
    beam2aero[1,0]=1.0
    beam2aero[2,2]=1.0
    for i in range(M+1):
        for j in range(N+1):
            zeta[i,j,:]=np.dot(beam2aero,Zeta[i,j,:])
            nu[i,j,:]=np.dot(beam2aero,Uext[i,j,:])
    zeta = zeta.flatten()/chord
    nu=nu.flatten()/2.0/vRef
    # wake geom
    zetaW = np.zeros((mW+1,N+1,3))
    # get wake direction from root TE and far wake point
    delW=(ZetaStar[1,0,:]-ZetaStar[0,0,:])/np.linalg.norm(ZetaStar[1,0,:]-ZetaStar[0,0,:])*chord/M
    for i in range(mW+1):
        for j in range(N+1):
            zetaW[i,j,:]=Zeta[M,j,:]+i*delW
            zetaW[i,j,:]=np.dot(beam2aero,zetaW[i,j,:])
    zetaW = zetaW.flatten()/chord
    
    # init uninit for SS uvlm description (zero reference conditions)
    gamPri = (np.zeros_like(Gamma)).flatten() #check these they might have to be flat
    
    # set relative velocity
    zetaPri = (np.zeros_like(Zeta)).flatten() #remember to nondim if nonzero
    
    return gam, gamW, gamPri, zeta, zetaW, zetaPri, nu, beam2aero

def runLinearAero(E,F,G,C,D,delS,nT,u,x0 = None):
    """@details run time-domain simulation of linear aerodynamics.
    @param E Discrete-time state-space matrix.
    @param F Discrete-time state-space matrix.
    @param G Discrete-time state-space matrix.
    @param C Discrete-time state-space matrix.
    @param D Discrete-time state-space matrix.
    @param delS Non-dimensional time step of model.
    @param nT Number of time-steps to run simulation.
    @param u Inputs, an nT x m array.
    @return x State history, an nT x n array.
    @return y Output history, an nT x l array.
    """
    
    invE = np.linalg.inv(E)    
    # run simulation
    tOut, yOut, xOut = dlsim((np.dot(invE,F),np.dot(invE,G),C,D,delS),
                             u,
                             None,
                             x0)
    
    return tOut, yOut, xOut
    
                        

if __name__ == '__main__':
    if True:
        Settings.OutputDir = '/home/' + getpass.getuser() + '/Desktop/'
        # AR=1e10
        for m in (30,):
            for mW in (30*m,):
                #genLinearRectWing(AR,m,mW,40,e=0.33,writeToMat = True,imageMeth = True)
                genLinearAerofoil(m,mW,writeToMat = True)
                
    if False:
        # aerofoil model
        m=20
        mW=20*m
        E,F,G,C,D,delS = genLinearAerofoil(m,mW)
        # simulation inputs
        nT=10 # number of time-steps
        u = np.zeros((nT,D.shape[1]))
        tOut,yOut,xOut = runLinearAero(E,F,G,C,D,delS,nT,u)
        