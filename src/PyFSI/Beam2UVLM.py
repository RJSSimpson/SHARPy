'''@package PyFSI.Beam2UVLM
@brief      Displacement and velocity extrapolation from beam FEM to UVLM grid.
@author     Rob Simpson
@contact    r.simpson11@imperial.ac.uk 
@version    1.0
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
from PyAero.UVLM.Utils import DerivedTypesAero
from PyBeam.Utils.Misc import iNode2iElem

def CoincidentGrid(PosDefor, PsiDefor, Section,
                   VelA_A, OmegaA_A, PosDotDef, PsiDotDef,
                   XBINPUT, AeroGrid, AeroVels,
                   OriginA_G = None,
                   PsiA_G = None,
                   ctrlSurf = None):
    """@brief Creates aero grid and velocities 
    based on beam nodal information.
    @param PosDefor Array of beam nodal displacements.
    @param PsiDefor Array of beam nodal rotations.
    @param Section Array containing sectional camberline coordinates.
    @param VelA_A Velocity of a-frame projected in a-frame.
    @param OmegaA_A Angular vel of a-frame projected in a-frame.
    @param PosDotDef Array of beam nodal velocities.
    @param PsiDotDef Array of beam nodal angular velocities.
    @param XBINPUT Beam input options.
    @param AeroGrid Aerodynamic grid to be updated.
    @param AeroVels Aerodynamic grid velocities to be updated
    @param OriginA_G Origin of a-frame.
    @param PsiA_G Attitude of a-frame w.r.t earth.
    @param ctrlSurf Control surface definition. 
    
    @details AeroGrid and AeroVels are projected in G-frame.
    """
    
    NumNodes = PosDefor.shape[0]
    NumNodesElem = XBINPUT.NumNodesElem
            
    for iNode in range(PosDefor.shape[0]):
        # Work out what element we are in.
        iElem, iiElem = iNode2iElem(iNode, NumNodes, NumNodesElem)
        
        # Calculate transformation matrix for each node.
        CaB = Psi2TransMat(PsiDefor[iElem,iiElem,:])
        
        # Loop through section coordinates.
        for jSection in range(Section.shape[0]):
            # Calculate instantaneous aero grid points
            AeroGrid[jSection,iNode,:] = PosDefor[iNode,:] \
                                         + np.dot(CaB,Section[jSection,:])
                                         
            # Calculate inertial angular velocity of B-frame projected in 
            # the B-frame.
            Omega_B_B = (np.dot(Tangential(PsiDefor[iElem,iiElem,:]),
                                PsiDotDef[iElem,iiElem,:])
                         + np.dot(CaB.T, OmegaA_A))
                        
            # Calculate inertial velocity at grid points (projected in A-frame).                             
            AeroVels[jSection,iNode,:] = (VelA_A
                        + np.dot(Skew(OmegaA_A),PosDefor[iNode,:])
                        + PosDotDef[iNode,:]
                        + np.dot(CaB, np.dot(
                        Skew(Omega_B_B), Section[jSection,:])))
            
            if (ctrlSurf != None and 
                jSection > ctrlSurf.iMin and jSection <= ctrlSurf.iMax and
                iNode >= ctrlSurf.jMin and iNode <= ctrlSurf.jMax):
                # Apply changes to Zeta and ZetaDot in a-frame due to control
                # surface movement.
                
                # Check if control surface indices are within range of grid
                assert (ctrlSurf.iMin >= 0), 'CtrlSurf iMin less then zero'
                assert (ctrlSurf.iMax < Section.shape[0]), ('CtrlSurf iMax '
                        + 'greater than section iMax')
                assert (ctrlSurf.jMin >= 0), 'CtrlSurf jMin less then zero'
                assert (ctrlSurf.jMax < PosDefor.shape[0]), ('CtrlSurf jMax '
                        + 'greater than section jMax')
                
                # Determine hinge point
                hingePoint = Section[ctrlSurf.iMin,:]
                
                # Use a CRV to define rotation in sectional frame
                hingeRot = np.array([ctrlSurf.beta, 0.0, 0.0])
                hingeRotDot = np.array([ctrlSurf.betaDot, 0.0, 0.0])
                
                # Find lever arm from hinge point
                leverArm = Section[jSection,:] - hingePoint
                
                # Deformed position of node in B-frame
                CBBnew = Psi2TransMat(hingeRot)
                newSectionPoint = hingePoint + np.dot(CBBnew,leverArm)
                
                # Overwrite corresponding element in Aerogrid
                AeroGrid[jSection,iNode,:] = (PosDefor[iNode,:]
                                              + np.dot(CaB,newSectionPoint))
                
                # Overwrite corresponding element in AeroVels
                # There is an extra velocity contribution due to the rotation 
                # of the hinge, with reference frame Bnew. Hence, 
                # CBBnew * ( Omega_Bnew_Bnew X Point__Bnew), which is the
                # velocity at a point on the flap due to the rotation
                # of the hinge projected in the B-frame. is added.
                # Note: Skew(hingeRotDot) is an acceptable definition of the
                # angular rate as it corresponds to planar motion within the
                # section. In general Skew(T(\psi)\dot{\psi}) is required.
                # Warning: should omega_Bnew_Bnew  (Skew(hingeRotDot))
                # be an inertial velocity?
                AeroVels[jSection,iNode,:] = (VelA_A
                        + np.dot(Skew(OmegaA_A),PosDefor[iNode,:])
                        + PosDotDef[iNode,:]
                        + np.dot(CaB, 
                                 np.dot(Skew(Omega_B_B), newSectionPoint)
                                 + np.dot(CBBnew,
                                          np.dot(Skew(hingeRotDot),
                                                 leverArm))))
                
        #END for jSection
    #END for iNode
    
    if ( (OriginA_G != None) and (PsiA_G != None) ):
        # Get transformation from a-frame to earth frame.
        CaG = Psi2TransMat(PsiA_G)
        CGa = CaG.T
        
        # Add the origin to grids in earth frame and transform velocities.
        for iNode in range(PosDefor.shape[0]):
            for jSection in range(Section.shape[0]):
                AeroGrid[jSection,iNode,:] = OriginA_G + \
                                             np.dot(
                                             CGa, AeroGrid[jSection,iNode,:])
                                             
                AeroVels[jSection,iNode,:] = np.dot(
                                             CGa, AeroVels[jSection,iNode,:])
            #END for jSection
        #END for iNode
    #END if    

def CoincidentGridForce(XBINPUT, PsiDefor, Section, AeroForces,
                           BeamForces):
    """@brief Calculates aero forces and moments on the beam nodes from those on
    the aerodynamic grid.
    @param XBINPUT Beam input options.
    @param PsiDefor Array of beam nodal rotations.
    @param Section Array containing sectional camberline coordinates.
    @param AeroForces Aerodynamic grid forces to be mapped to beam nodes.
    @param BeamForces BeamForces to overwrite.
    
    @details All Beam forces calculated in in a-frame, while aero forces 
    are defined in earth frame."""
    
    # Zero beam forces.
    BeamForces[:,:] = 0.0
    
    # Get transformation matrix between a-frame and earth.
    CGa = Psi2TransMat(XBINPUT.PsiA_G)
    CaG = CGa.T
    
    # Num Nodes
    NumNodes = XBINPUT.NumNodesTot
    
    # Get number of nodes per beam element.
    NumNodesElem = XBINPUT.NumNodesElem
    
    # Loop along beam length.
    for iNode in range(NumNodes):
        
        iElem, iiElem = iNode2iElem(iNode, NumNodes, NumNodesElem)
        
        # Calculate transformation matrix for each node.
        CaB = Psi2TransMat(PsiDefor[iElem,iiElem,:])
        
        # Loop through each sectional coordinate.
        for jSection in range(Section.shape[0]):
            
            # Get Section coord and AeroForce in a-frame.
            Section_A = np.dot(CaB,Section[jSection,:])
            AeroForce_A = np.dot(CaG,AeroForces[jSection][iNode][:])
            
            # Calc moment.
            BeamForces[iNode,3:] += np.cross(Section_A,
                                             AeroForce_A)
            
            # Calc force.
            BeamForces[iNode,:3] += AeroForce_A
        
        # END for jSection
    #END for iNode

def InitSection(VMOPTS,VMINPUT,ElasticAxis):
    """@brief Initialise section based on aero & aeroelastic elastic options.
    @param VMOPTS UVLM solver options.
    @param VMINPUT UVLM solver inputs.
    @param ElasticAxis Position of elastic axis, defined as Theodorsen's a 
    parameter - i.e the number of semi-chords aft of midchord.
    @return Section cross section coordinates."""
    
    Section = np.zeros((VMOPTS.M.value+1,3),ct.c_double,'C')
    
    # Calculate rotation due to twist.
    Psi = np.zeros((3))
    Psi[0] = VMINPUT.theta
    
    # Get transformation matrix.
    R = Psi2TransMat(Psi)
    
    # Get delta chord.
    DeltaC = VMINPUT.c/VMOPTS.M.value
    
    # Get UVLM discretisation offset.
    Offset = 0.25*DeltaC
    
    # based on elastic axis at (z= 0, y = 0) get leading edge position.
    LeadingEdge = 0.5*VMINPUT.c + ElasticAxis*0.5*VMINPUT.c - Offset
    
    # calculate section coordinates.
    for j in range(VMOPTS.M.value+1):
        Section[j,:] = np.dot(R,[0.0,LeadingEdge-j*DeltaC,0.0])
        
    return Section
 

if __name__ == '__main__':
    import DerivedTypes
    XBINPUT = DerivedTypes.Xbinput(3,1)
    AeroM = 1
    
    #Initialise for Beam
    PosDefor = np.zeros((XBINPUT.NumNodesTot,3),ct.c_double,'F')
    for i in range(PosDefor.shape[0]):
        PosDefor[i,0] = i
         
    PosDotDef = np.zeros((XBINPUT.NumNodesTot,3),ct.c_double,'F')
    
    PsiDefor = np.zeros((XBINPUT.NumElems,Settings.MaxElNod,3),\
                         dtype=ct.c_double, order='F')
    
    PsiDotDef = np.zeros((XBINPUT.NumElems,Settings.MaxElNod,3),\
                         dtype=ct.c_double, order='F')
    
    VelA_A = np.zeros(3)
    OmegaA_A = np.zeros(3)
    
    BeamForces = np.zeros((XBINPUT.NumNodesTot,6),ct.c_double,'F')
    
    #Initialise data for UVLM
    VMOPTS = DerivedTypesAero.VMopts(1,1,True)
    VMINPUT = DerivedTypesAero.VMinput(1.0, 1.0, 25.0, 2.0*np.pi/180.0, \
                                             0.0*np.pi/180.0)
    
    Section = InitSection(VMOPTS,VMINPUT,ElasticAxis=-0.5)
    
    AeroGrid = np.zeros((Section.shape[0],PosDefor.shape[0],3),ct.c_double,'C')
    AeroVels = np.zeros((Section.shape[0],PosDefor.shape[0],3),ct.c_double,'C')
    
    CoincidentGrid(PosDefor, PsiDefor, Section, VelA_A, \
                   OmegaA_A, PosDotDef, PsiDotDef, XBINPUT,\
                   AeroGrid, AeroVels)
    
    print(AeroGrid, '\n')
    
    #Create forces
    AeroForces = np.zeros((AeroM+1,XBINPUT.NumNodesTot,3),ct.c_double,'C')
    AeroForces[:,:,2] = 1.0
    
    #map to beam forces
    CoincidentGridForce(XBINPUT, PsiDefor, Section, AeroForces,\
                        BeamForces)
    
    print(BeamForces, '\n')
    
    
    #forces must be reset
    BeamForces[:,:] = 0.0
    
    #rotate by 90 degrees and check
    PsiDefor[:,:,0] = np.pi/2.0
    
    #create updated grid
    CoincidentGrid(PosDefor, PsiDefor, Section, VelA_A, \
                   OmegaA_A, PosDotDef, PsiDotDef, XBINPUT,\
                   AeroGrid, AeroVels)
    
    print(AeroGrid, '\n')
    
    #map to beam forces
    CoincidentGridForce(XBINPUT, PsiDefor, Section, AeroForces,\
                        BeamForces)
    
    print(BeamForces, '\n')
    
    #check applied forces
    #declare temp traids
    Psi = np.zeros((3))
    
    #loop through all beam nodes
    BeamForceCheck = np.zeros((3))
    for iNode in range(XBINPUT.NumNodesTot):
        BeamForceCheck[:] += BeamForces[iNode,:3]
    
    #account for rotation of aerodynamic FoR (freestream velocity)
    Psi[0] = VMINPUT.alpha
    
    #get transformation matrix
    CalphaG = Psi2TransMat(Psi)
    
    print(np.dot(CalphaG,BeamForceCheck))
    
