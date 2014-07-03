'''@package SHARPy.input.DancingFrame.PyBeemInput
@brief      Input functions to populate DerivedTypes for PyBeam simulations.
@author     Rob Simpson
@author     Henrik Hesse
@contact    r.simpson11@imperial.ac.uk
@version    0.0
@date       03/12/13
@pre        None
@warning    None
'''

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
    
    XBINPUT.ForceDyn_dead[0,0:3]   = [0,250,0]
      
    XBINPUT.ForceDyn_dead[(NumNodes_tot-1)/2,0:3]   = [0,0, 1000]
    XBINPUT.ForceDyn_dead[(NumNodes_tot-1)/2+1,0:3] = [0,0,-1000]
        
    XBINPUT.ForceDyn_foll[(NumNodes_tot-1)/2,3:]   = [0, 2000,0]
    XBINPUT.ForceDyn_foll[(NumNodes_tot-1)/2+1,3:] = [0,-2000,0]
         
    XBINPUT.ForceDyn_dead[2*(NumNodes_tot-1)/6,0:3]   = [0, 100,0]
    XBINPUT.ForceDyn_dead[4*(NumNodes_tot-1)/6+1,0:3] = [0,-100,0]
    
    XBINPUT.ForcingType = 'Ramp'
    XBINPUT.RampTime = 10.0
    
    return XBINPUT, XBOPTS