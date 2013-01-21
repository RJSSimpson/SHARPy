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
    Psi = np.zeros(3,ct.c_double,'F')
    Psi[0] = np.pi
    Psi2TransMat(Psi)