'''@package PyBeam.Utils.XbeamLib
@brief      Functions found in lib_xbeam.f90, lib_rotvect.f90
@author     Rob Simpson
@contact    r.simpson11@imperial.ac.uk 
@version    0.0
@date       15/01/2013
@pre        None
@warning    None
'''

import numpy as np
import ctypes as ct
import SharPySettings as Settings
from scipy import linalg

def QuadSkew(Omega):
    """@brief Calculate incremental update operator for Quaternion based on Omega.
    Quaternion ODE to compute orientation of body-fixed frame a.
    Duplicated from lib_xbeam.f90
    See Shearer and Cesnik, 'Nonlinear Flight Dynamics of Very Flexible Aircraft',
    Journal of Aircraft 2007; 44 (5): 1528-1545
    """
    
    QS = np.zeros((4,4), ct.c_double, 'F')
    
    QS[0,1], QS[3,2], QS[1,0], QS[2,3] = \
        Omega[0], Omega[0], -Omega[0], -Omega[0]
        
    QS[0,2], QS[1,3], QS[2,0], QS[3,1] = \
        Omega[1], Omega[1], -Omega[1], -Omega[1]
        
    QS[0,3], QS[2,1], QS[1,2], QS[3,0] = \
        Omega[2], Omega[2], -Omega[2], -Omega[2]
        
    return QS


def Rot(q1):
    """@brief Calculate rotation matrix based on quaternions.
    See Aircraft Control and Simulation, pag. 31, by Stevens, Lewis.
    """
    
    q = q1.copy(order='F')
    q = q/np.linalg.norm(q)
    
    RotMat = np.zeros((3,3), ct.c_double, 'F')
    
    RotMat[0,0] = q[0]**2 + q[1]**2 - q[2]**2 - q[3]**2
    RotMat[1,1] = q[0]**2 - q[1]**2 + q[2]**2 - q[3]**2
    RotMat[2,2] = q[0]**2 - q[1]**2 - q[2]**2 + q[3]**2
    
    RotMat[0,1] = 2.*(q[1]*q[2] + q[0]*q[3])
    RotMat[1,0] = 2.*(q[1]*q[2] - q[0]*q[3])
    
    RotMat[0,2] = 2.*(q[1]*q[3] - q[0]*q[2])
    RotMat[2,0] = 2.*(q[1]*q[3] + q[0]*q[2])
    
    RotMat[1,2] = 2.*(q[2]*q[3] + q[0]*q[1])
    RotMat[2,1] = 2.*(q[2]*q[3] - q[0]*q[1])
    
    return RotMat


def Psi2TransMat(Psi):
    """@brief Calculates the transformation matrix associated with CRV Psi."""
    TransMat = np.zeros((3,3))
    Norm = TriadNorm(Psi)
    SkewPsi=Skew(Psi)
    "Check norm of Psi for linearised OR fully-nonlinear matrix"
    if Norm < Settings.RotEpsilon:
        TransMat[:,:] = np.identity(3) + SkewPsi + 0.5*np.dot(SkewPsi,SkewPsi) 
    else:
        TransMat[:,:] = np.identity(3) \
                      + np.sin(Norm)/Norm*SkewPsi \
                      + (1 - np.cos(Norm))/Norm**2 * np.dot(SkewPsi,SkewPsi)
    
    return TransMat


def TriadNorm(Triad):
    """@brief Returns the norm of a 3x1 tensor (Triad)."""
    return np.sqrt(Triad[0]**2+Triad[1]**2+Triad[2]**2)


def Skew(Vector):
    """@brief Returns the skew-symmetric Matrix associated to Vector"""
    
    SkewMat = np.zeros((3,3))
    
    SkewMat[0,1]=-Vector[2]
    SkewMat[0,2]= Vector[1]
    SkewMat[1,0]= Vector[2]
    SkewMat[1,2]=-Vector[0]
    SkewMat[2,0]=-Vector[1]
    SkewMat[2,1]= Vector[0]
    
    return SkewMat

def VectofSkew(Skew):
    """@brief Returns the vector defining skew-symmetric matrix Skew."""
    
    Vect = np.array([Skew[2,1],Skew[0,2],Skew[1,0]])
    
    return Vect


def Tangential(Psi):
    """@brief Calculates the tangential operator related to Psi, a cartesian
     rotation vector."""
    Tang = np.zeros((3,3))
    Norm = TriadNorm(Psi)
    SkewPsi=Skew(Psi)
    "Check norm of Psi for linearised OR fully-nonlinear matrix"
    if Norm < Settings.RotEpsilon:
        Tang[:,:] = np.identity(3) - 0.5* SkewPsi + 1.0/6.0*np.dot(SkewPsi,SkewPsi) 
    else:
        Tang[:,:] = np.identity(3) \
                    - (1-np.cos(Norm))/Norm**2 * SkewPsi \
                    + (1 - np.sin(Norm)/Norm) * np.dot(SkewPsi,SkewPsi)/Norm**2
    
    return Tang
     


if __name__ == '__main__':
    
    "create matrices that correspond to Psi of nodes 1 and 2"
    Psi1 = np.array([0.1, 0.1, 0.1])
    Psi2 = np.array([0.11, 0.12, 0.13])

    CaB1 = Psi2TransMat(Psi1)
    CaB2 = Psi2TransMat(Psi2)
    
    "Calc local rotation in frame B1"
    CB1B2 = np.dot(CaB1.T,CaB2)
    
    "Get skew-symmetric matrix"
    phitilda12 = linalg.logm(CB1B2)
    
    "get vector"
    phi12 = VectofSkew(phitilda12)
    
    "get C matrix at midpoint"
    CaBh = np.dot(CaB1,linalg.expm(0.5 * phitilda12))
    print('Objectively interpolated matrix = \n',CaBh,'\n')
    
    "now try by direct interp of CRV"
    Psih_err = 0.5*(Psi2 + Psi1)
    CaBh_err = Psi2TransMat(Psih_err)
    print('Non-objective matrix = \n',CaBh_err,'\n')
    
    "diff"
    print('diff = \n',CaBh_err - CaBh,'\n')
    
    "get objectively interped Psi"
    Psih = VectofSkew(linalg.logm(CaBh))
    print('Objectively interpolated Psi = \n',Psih,'\n')
    
    "diff"
    print('diff = \n',Psih_err - Psih,'\n')
    
    
    
    