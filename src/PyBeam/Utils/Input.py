'''@package PyBeam.Utils.Input
@brief      Input functions to populate DerivedTypes for PyBeam simulations.
@author     Rob Simpson
@contact    r.simpson11@imperial.ac.uk
@version    0.0
@date       28/11/2012
@pre        None
@warning    None
'''

import sys
import DerivedTypes
import SharPySettings as Settings
import numpy as np
import ctypes as ct

def Setup(XBINPUT, XBOPTS):        
    """@brief Set up beam input and check options for inconsistencies.
    @details Initial functionality adapted from FxInput_PyFx.f90."""
    
    # Check number of nodes per element and Gauss points.
    if XBINPUT.NumNodesElem != 2 and XBINPUT.NumNodesElem != 3:
        sys.stderr.write('Invalid number of nodes per element\n')
    elif XBINPUT.NumNodesElem == 2 and XBOPTS.NumGauss.value != 1:
        sys.stderr.write('Number of nodes per element (2) inconsistent with' +\
                         ' number of gauss points (%d)'\
                         %(XBOPTS.NumGauss.value) +\
                         ' Gauss points set to (1)\n')
        XBOPTS.NumGauss.value = 1
    elif XBINPUT.NumNodesElem == 3 and XBOPTS.NumGauss.value != 2:
        sys.stderr.write('Number of nodes per element (3) inconsistent with' +\
                         ' number of gauss points(%d)'\
                         %(XBOPTS.NumGauss.value)  +\
                         ' Gauss points set to (2)\n')
        XBOPTS.NumGauss.value = 2    
        
    # General set-up of test case
#     XBINPUT = DerivedTypes.Xbinput(2,XBINPUT.NumElems)
    XBINPUT.BeamLength = 10.0

    SectWidth =0.1
    SectHeight=0.05
    E  = 70.e9
    G  = E/(2.0*(1.0+0.3))
    rho= 2700.0
    
    XBINPUT.BeamMass[0,0] = rho*SectWidth*SectHeight
    XBINPUT.BeamMass[1,1] = XBINPUT.BeamMass[0,0]
    XBINPUT.BeamMass[2,2] = XBINPUT.BeamMass[0,0]
    XBINPUT.BeamMass[3,3] = rho*SectWidth*SectHeight*(SectHeight**2.0+SectWidth**2.0)/12.0
    XBINPUT.BeamMass[4,4] = rho*(SectWidth*SectHeight**3.0)/12.0
    XBINPUT.BeamMass[5,5] = rho*(SectHeight*SectWidth**3.0)/12.0

    XBINPUT.BeamStiffness[0,0] = SectWidth*SectHeight*E
    XBINPUT.BeamStiffness[1,1] = (5.0/6.0)*SectWidth*SectHeight*G
    XBINPUT.BeamStiffness[2,2] = (5.0/6.0)*SectWidth*SectHeight*G
    XBINPUT.BeamStiffness[3,3] = G*0.22887745119762*SectWidth*(SectHeight**3.0)
    XBINPUT.BeamStiffness[4,4] = E*(SectWidth*(SectHeight**3.0)/12.0)
    XBINPUT.BeamStiffness[5,5] = E*(SectHeight*(SectWidth**3.0)/12.0)

    XBOPTS.FollowerForce = ct.c_bool(False)
    XBOPTS.FollowerForceRig = ct.c_bool(False)
    XBOPTS.OutInaframe = ct.c_bool(False)
    XBOPTS.OutInBframe = ct.c_bool(False)
    XBOPTS.NumLoadSteps = ct.c_int(1)
    
    # Ensure that follower forces are consistently defined
    if XBOPTS.FollowerForce == ct.c_bool(True):
        XBOPTS.FollowerForceRig = ct.c_bool(True)
    
    sys.stderr.flush()    
    return XBINPUT, XBOPTS


def Elem(XBINPUT, XBOPTS, XBELEM):
    """@brief Set up Xbelem derived type.
    @details Initial functionality adapted from 
    FxInput_PyFx.f90."""
    
    # Check number of nodes per element and set default connectivity.
    assert (XBINPUT.NumNodesElem == 2 or XBINPUT.NumNodesElem == 3),\
            'Invalid number of nodes per element'
    if XBINPUT.NumNodesElem == 2:
        # connectivities of right arm
        for ElemNo in range(0,int(XBINPUT.NumElems/2)):
            i = ElemNo*Settings.MaxElNod
            XBELEM.Conn[i]=ElemNo+1
            XBELEM.Conn[i+1]=ElemNo+2
            XBELEM.NumNodes[ElemNo] = 2
        # connectivities of left arm    
        for ElemNo in range(int(XBINPUT.NumElems/2),XBINPUT.NumElems):
            i = ElemNo*Settings.MaxElNod
            XBELEM.Conn[i]=ElemNo+2
            XBELEM.Conn[i+1]=ElemNo+3
            XBELEM.NumNodes[ElemNo] = 2
        XBELEM.Conn[i+1]=1
            
        NumNodes_tot = ct.c_int(XBINPUT.NumElems + 1) 
            
    elif XBINPUT.NumNodesElem == 3:
        # connectivities of right arm
        for ElemNo in range(0,int(XBINPUT.NumElems/2)):
            i = ElemNo*Settings.MaxElNod
            XBELEM.Conn[i]=2*ElemNo+1
            XBELEM.Conn[i+1]=2*ElemNo+3
            XBELEM.Conn[i+2]=2*ElemNo+2
            XBELEM.NumNodes[ElemNo] = 3
        # connectivities of left arm 
        for ElemNo in range(int(XBINPUT.NumElems/2),XBINPUT.NumElems):
            i = ElemNo*Settings.MaxElNod
            XBELEM.Conn[i]=2*ElemNo+2
            XBELEM.Conn[i+1]=2*ElemNo+4
            XBELEM.Conn[i+2]=2*ElemNo+3
            XBELEM.NumNodes[ElemNo] = 3
        XBELEM.Conn[i+1]=1
        
        NumNodes_tot = ct.c_int(2*XBINPUT.NumElems + 1)
    
    # Calculate Beam Inverse Stiffness Matrix.        
    InverseStiffness = np.linalg.inv(XBINPUT.BeamStiffness)
    
    # Populate XBELEM
    for ElemNo in range(0,XBINPUT.NumElems):
        i_start = ElemNo*6
        i_end = ElemNo*6+6
        XBELEM.Mass[i_start:i_end,:] = XBINPUT.BeamMass
        XBELEM.Stiff[i_start:i_end,:] = XBINPUT.BeamStiffness
        XBELEM.InvStiff[i_start:i_end,:] = InverseStiffness
    
    # Element orientation and member number.
    for ElemNo in range(0,XBINPUT.NumElems):
        i = ElemNo*3
        XBELEM.Vector[i+1] = 1.0
        XBELEM.MemNo[ElemNo] = 0
    
    return NumNodes_tot, XBELEM


def Node(XBINPUT, XBOPTS, NumNodes_tot, XBELEM):
    """@brief Define nodal properties."""
    
    # Custom double L frame starting at . : |_._|
    PosIni = np.zeros((NumNodes_tot.value,3), dtype=ct.c_double, order='F')
    for NodeNo in range(0,NumNodes_tot.value):
        if NodeNo <= (NumNodes_tot.value-1)/3:
            PosIni[NodeNo,0] = XBINPUT.BeamLength*NodeNo/((NumNodes_tot.value-1)/3)
        elif NodeNo <= (NumNodes_tot.value-1)/2:
            PosIni[NodeNo,0] = XBINPUT.BeamLength
            PosIni[NodeNo,2] = XBINPUT.BeamLength/2*\
                            (NodeNo-(NumNodes_tot.value-1)/3)/((NumNodes_tot.value-1)/6)
        elif NodeNo <= 4*(NumNodes_tot.value-1)/6:
            PosIni[NodeNo,0] =-XBINPUT.BeamLength
            PosIni[NodeNo,2] = XBINPUT.BeamLength/2*\
                            (1-(NodeNo-1-((NumNodes_tot.value-1)/2))/((NumNodes_tot.value-1)/6))
        else:
            PosIni[NodeNo,0] = XBINPUT.BeamLength*\
                            (-1+(NodeNo-1-(4*(NumNodes_tot.value-1)/6))/((NumNodes_tot.value-1)/3))

    # Define pre-twist.
    PhiNodes = np.zeros(NumNodes_tot.value, dtype=ct.c_double, order='F')
    
    # Define external forces at nodes
#     XBINPUT.ForceStatic[(NumNodes_tot.value-1)/2,0:3]   = [0,0, 1000]
#     XBINPUT.ForceStatic[(NumNodes_tot.value-1)/2+1,0:3] = [0,0,-1000]
    
#     XBINPUT.ForceStatic_dead[0,0:3]   = [0,125,0]
#       
#     XBINPUT.ForceStatic_dead[(NumNodes_tot.value-1)/2,0:3]   = [0,0, 1000]
#     XBINPUT.ForceStatic_dead[(NumNodes_tot.value-1)/2+1,0:3] = [0,0,-1000]
#         
#     XBINPUT.ForceStatic_foll[(NumNodes_tot.value-1)/2,3:]   = [0, 2000,0]
#     XBINPUT.ForceStatic_foll[(NumNodes_tot.value-1)/2+1,3:] = [0,-2000,0]
#          
#     XBINPUT.ForceStatic_dead[2*(NumNodes_tot.value-1)/6,0:3]   = [0, 100,0]
#     XBINPUT.ForceStatic_dead[4*(NumNodes_tot.value-1)/6+1,0:3] = [0,-100,0]
    
    # Declare and populate boundary conditions.
    BoundConds = np.zeros(NumNodes_tot.value,dtype=ct.c_int,order='F')
    BoundConds[0] = 1
    BoundConds[(NumNodes_tot.value-1)/2]   = -1
    BoundConds[(NumNodes_tot.value-1)/2+1] = -1
        
    return PosIni, PhiNodes, BoundConds
    

def DynSetup(XBINPUT, XBOPTS):        
    """@brief Setup of dynamic simulation properties.
    @details Initial functionality adapted from FxInput_PyFx.f90."""
    
    XBINPUT.tfin = 10.0
    XBINPUT.dt = 0.01
    XBOPTS.NewmarkDamp = ct.c_double(1.0e-2)
    XBOPTS.MinDelta = ct.c_double(1.0e-3)
        
    # Define external forces at nodes
    NumNodes_tot = np.shape(XBINPUT.ForceStatic)[0]
    
#     XBINPUT.ForceDyn=1*XBINPUT.ForceStatic
#     XBINPUT.ForceStatic=0*XBINPUT.ForceStatic
    
    XBINPUT.ForceDyn_dead[0,0:3]   = [0,250,0]
      
    XBINPUT.ForceDyn_dead[(NumNodes_tot-1)/2,0:3]   = [0,0, 1000]
    XBINPUT.ForceDyn_dead[(NumNodes_tot-1)/2+1,0:3] = [0,0,-1000]
        
    XBINPUT.ForceDyn_foll[(NumNodes_tot-1)/2,3:]   = [0, 2000,0]
    XBINPUT.ForceDyn_foll[(NumNodes_tot-1)/2+1,3:] = [0,-2000,0]
         
    XBINPUT.ForceDyn_dead[2*(NumNodes_tot-1)/6,0:3]   = [0, 100,0]
    XBINPUT.ForceDyn_dead[4*(NumNodes_tot-1)/6+1,0:3] = [0,-100,0]
    
    XBINPUT.ForcingType = 'Ramp'
    XBINPUT.RampTime = 10.0      #XBINPUT.tfin
    
    
    return XBINPUT, XBOPTS


if __name__ == '__main__':
    XBINPUT = DerivedTypes.Xbinput()
    XBINPUT.NumNodesElem = 3
    XBOPTS = DerivedTypes.Xbopts()
    XBINPUT, XBOPTS = Setup(XBINPUT, XBOPTS)
    XBINPUT, XBOPTS = Setup(XBINPUT, XBOPTS)
    
    XBINPUT.NumNodesElem = 3
    XBINPUT.NumElems = 8
    XBINPUT.BeamStiffness[0,0] = 1.0e+09
    XBINPUT.BeamStiffness[1,1] = 1.0e+09
    XBINPUT.BeamStiffness[2,2] = 1.0e+09
    XBINPUT.BeamStiffness[3,3] = 1.0e+04
    XBINPUT.BeamStiffness[4,4] = 2.0e+04
    XBINPUT.BeamStiffness[5,5] = 4.0e+06
    XBELEM = DerivedTypes.Xbelem(XBINPUT.NumElems, Settings.MaxElNod)
    NumNodes_tot, XBELEM = Elem(XBINPUT, XBOPTS, XBELEM)
    
    XBINPUT.ForceStatic[2] = 12345.6
    PosIni, PhiNodes, ForceStatic, BoundConds = Node(XBINPUT, XBOPTS,\
                                                      NumNodes_tot, XBELEM)
    
    