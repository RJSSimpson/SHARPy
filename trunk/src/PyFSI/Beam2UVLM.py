'''@package PyFSI.Beam2UVLM
@brief      Displacement extrapolation from beam FEM to UVLM grid.
@author     Rob Simpson
@contact    r.simpson11@imperial.ac.uk 
@version    0.0
@date       18/01/2013
@pre        None
@warning    None
'''

import numpy as np
import ctypes as ct
import SharPySettings as Settings
from PyBeam.Utils.XbeamLib import Psi2TransMat
from PyBeam.Utils.XbeamLib import Skew
from PyBeam.Utils.XbeamLib import Tangential

def isodd(num):
    """@brief returns True if Num is odd"""
    return bool(num & 1)

def CoincidentGrid(PosDefor, PsiDefor, Section,\
                   VelA_A, OmegaA_A, PosDotDef, PsiDotDef,
                   XBINPUT):
    """@brief Creates aero grid and velocities 
    centred on beam nodes.
    @param PosDefor Array of beam nodal displacements.
    @param PsiDefor Array of beam nodal rotations.
    @param Section Array containing sectional camberline coordinates.
    @param VelA_A Velocity of A-frame projected in A-frame.
    @param OmegaA_A Angular vel of A-frame projected in A-frame.
    @param PosDotDef Array of beam nodal velocities.
    @param PsiDotDef Array of beam nodal angular velocities.
    
    @details All displacements and velocities are projected in A-frame. 
    All Omegas are 'inertial' angular velocities, so their magnitudes are not 
    affected by what frame they are projected in (objectivity).
    TODO: check 2- and 3-noded elements.
    """
    
    #Initialise data for UVLM
    AeroGrid = np.zeros((PosDefor.shape[0],Section.shape[0],3),ct.c_double,'C')
    AeroVels = np.zeros((PosDefor.shape[0],Section.shape[0],3),ct.c_double,'C')
    
    NumNodesElem = XBINPUT.NumNodesElem
    NumNE_array = np.zeros(5,ct.c_int)
    NumNE_array[:] = 3
            
    for iNode in range(PosDefor.shape[0]):
        "Work out what element we are in (works for 2 and 3-noded)"
        if iNode == 0:
            iElem = 0
        elif iNode < PosDefor.shape[0]-1:
            iElem = int(iNode/(NumNodesElem-1))
        elif iNode == PosDefor.shape[0]-1:
            iElem = int((iNode-1)/(NumNodesElem-1))
            
        "Work out what sub-element node we are in"
        if NumNodesElem == 2:
            if iNode < PosDefor.shape[0]-1:
                iiElem = 0
            elif iNode == PosDefor.shape[0]-1:
                iiElem = 1
                
        elif NumNodesElem == 3:
            iiElem = 0
            if iNode == PosDefor.shape[0]-1:
                iiElem = 2 
            elif isodd(iNode):
                iiElem = 1
        
        "Calculate transformation matrix for each node"
        CaB = Psi2TransMat(PsiDefor[iElem,iiElem,:])
        
        "loop through section coordinates"
        for jSection in range(Section.shape[0]):
            "Calculate instantaneous aero grid points"
            AeroGrid[iNode,jSection,:] = PosDefor[iNode,:] \
                                         + np.dot(CaB,Section[jSection,:])
                                         
            "Calc inertial angular velocity of B-frame projected in B-frame"
            Omega_B_B = np.dot( Tangential(PsiDefor[iElem,iiElem,:]) \
                        , PsiDotDef[iElem,iiElem,:] ) \
                        + np.dot(CaB.T, OmegaA_A)
                        
            "Calc inertial velocity at grid points (projected in A-frame.)"                             
            AeroVels[iNode,jSection,:] = VelA_A \
                        + np.dot(Skew(OmegaA_A),PosDefor[iNode,:]) \
                        + PosDotDef[iNode,:] \
                        + np.dot(CaB, np.dot(\
                        Skew(Omega_B_B), Section[jSection,:] ))
            
        #END for jSection
    #END for iNode
    return AeroGrid, AeroVels    
    
    
     
 

if __name__ == '__main__':
    import DerivedTypes
    XBINPUT = DerivedTypes.Xbinput(3,10)
    AeroM = 5
    
    PosDefor = np.zeros((XBINPUT.NumNodesTot,3),ct.c_double,'F')
    for i in range(PosDefor.shape[0]):
        PosDefor[i,0] = i
         
    PosDotDef = np.zeros((XBINPUT.NumNodesTot,3),ct.c_double,'F')
    
    Section = np.zeros((AeroM,3),ct.c_double,'C')
    for j in range(Section.shape[0]):
        Section[j,1] = j*(1/AeroM)
    
    
    PsiDefor = np.zeros((XBINPUT.NumElems,Settings.MaxElNod,3),\
                         dtype=ct.c_double, order='F')
    
    PsiDotDef = np.zeros((XBINPUT.NumElems,Settings.MaxElNod,3),\
                         dtype=ct.c_double, order='F')
    
    VelA_A = np.zeros(3)
    OmegaA_A = np.zeros(3)
    
    AeroGrid, AeroVels = CoincidentGrid(PosDefor, PsiDefor, Section, VelA_A, \
                                        OmegaA_A, PosDotDef, PsiDotDef, XBINPUT)
    
    print(AeroGrid, '\n')
    print(AeroVels, '\n')