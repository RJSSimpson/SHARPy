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
from PyBeam.Utils.Misc import isodd

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
    # Check norm of Psi for linearised OR fully-nonlinear matrix.
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
    # Check norm of Psi for linearised OR fully-nonlinear matrix.
    if Norm < Settings.RotEpsilon:
        Tang[:,:] = np.identity(3) - 0.5* SkewPsi + 1.0/6.0*np.dot(SkewPsi,SkewPsi) 
    else:
        Tang[:,:] = np.identity(3) \
                    - (1-np.cos(Norm))/Norm**2 * SkewPsi \
                    + (1 - np.sin(Norm)/Norm) * np.dot(SkewPsi,SkewPsi)/Norm**2
    
    return Tang

def AddGravityLoads(BeamForces,XBINPUT,XBELEM,AELAOPTS,PsiDefor,
                      chord = 1.0):
    """@brief Apply nodal gravity loads.
    @param BeamForces Nodal forces to update.
    @param XBINPUT Xbeam inputs.
    @param XBELEM Xbeam element information.
    @param AELAOPTS Aeroelastic options.
    @param PsiA_G Cartesian rotation vector describing orientation of a-frame
           with respect to earth.
    @param chord Aerofoil chord, assumed to be 1.0 if ommited. 
    @details Offset of mass centroid from elastic axis is currently calculated
             using the AELAOPTS.ElasticAxis and .InertialAxis parameters which 
             are analogous to that used by Theodorsen. It therefore assumes the 
             aerofoil section is defined on the y-axis and hence the moment arm
             points in the B-frame y-direction.
    @warning Assumes even distribution of nodes along beam.
    @warning Not valid for twisted aerofoil sections, i.e those that are defined
             with some theta angle.
    """
    
    MassPerLength = XBINPUT.BeamMass[0,0]
    
    if XBINPUT.NumNodesElem == 2:
        NodeSpacing = XBELEM.Length[0]
    elif XBINPUT.NumNodesElem == 3:
        NodeSpacing = 0.5*XBELEM.Length[0]
        
    ForcePerNode = -NodeSpacing * MassPerLength * XBINPUT.g 
    
    # Obtain transformation from Earth to a-frame.
    CGa = Psi2TransMat(XBINPUT.PsiA_G)
    CaG = CGa.T
    
    # Force in a-frame.
    Force_a = np.dot(CaG,np.array([0.0,0.0,ForcePerNode]))
    
    # Indices for boundary nodes.
    Root = 0
    Tip = BeamForces.shape[0]-1
    
    # Apply forces.
    BeamForces[Root+1:Tip,:3] += Force_a
    BeamForces[Root,:3] += 0.5*Force_a
    BeamForces[Tip,:3] += 0.5*Force_a
    
    # Get number of nodes per beam element.
    NumNodesElem = XBINPUT.NumNodesElem
    
    # Loop through nodes to get moment arm at each.
    for iNode in range(XBINPUT.NumNodesTot):
        
        # Work out what element we are in (works for 2 and 3-noded).
        if iNode == 0:
            iElem = 0
        elif iNode < XBINPUT.NumNodesTot-1:
            iElem = int(iNode/(NumNodesElem-1))
        elif iNode == XBINPUT.NumNodesTot-1:
            iElem = int((iNode-1)/(NumNodesElem-1))
            
        # Work out what sub-element node (iiElem) we are in.
        if NumNodesElem == 2:
            if iNode < XBINPUT.NumNodesTot-1:
                iiElem = 0
            elif iNode == XBINPUT.NumNodesTot-1:
                iiElem = 1
                
        elif NumNodesElem == 3:
            iiElem = 0
            if iNode == XBINPUT.NumNodesTot-1:
                iiElem = 2 
            elif isodd(iNode):
                iiElem = 1
        
        # Calculate transformation matrix for each node.
        CaB = Psi2TransMat(PsiDefor[iElem,iiElem,:])
        
        # Define moment arm in B-frame
        # Moment arm to elastic axis defined using EA and IA of Theodorsen.
        armY = -(AELAOPTS.InertialAxis - AELAOPTS.ElasticAxis)*chord/2.0
        armY_a = np.dot(CaB,np.array([0.0, armY, 0.0]))
        
        # Calculate moment
        if (iNode == Root or iNode == Tip):
            BeamForces[iNode,3:] += np.cross(armY_a, 0.5*Force_a)
        else:
            BeamForces[iNode,3:] += np.cross(armY_a, Force_a)


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
    
    
    
    