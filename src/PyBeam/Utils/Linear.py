'''PyBeam.Utils.Linear
@brief      Generating linear beam models.
@author     Rob Simpson
@contact    r.simpson11@imperial.ac.uk 
@version    0.0
@date       09/07/2015
@pre        None
@warning    None
'''

import BeamLib
import numpy as np
import ctypes as ct
import SharPySettings as Settings
import PyBeam.Utils.DerivedTypes as DerivedTypes
import BeamInit
from PyFSI.Beam2Aero_Lin import printBlocks

np.set_printoptions(precision = 4)

def genSSbeam(XBINPUT,
               NumNodes_tot,
               NumDof,
               XBELEM,
               XBNODE,
               PosIni,
               PsiIni,
               PosDefor,
               PsiDefor,
               StaticForces,
               XBOPTS):
    """@details Generate state-space matrices for linear beam dynamics.
    @param XBINPUT Beam inputs.
    @param NumNodes_tot Total number of nodes in the model.
    @param NumDof Total number of degrees of freedom.
    @param XBELEM Beam FE node information.
    @param XBNODE Beam FE node information.
    @param PosIni Undeformed beam displacements.
    @param PsiIni Undeformed beam FE rotations.
    @param PosDefor Deformed beam displacements.
    @param PsiDefor Deformed beam FE rotations.
    @param StaticForces Total static forces.
    @param XBOPTS Beam options.
    @return A State transfer matrix.
    @return B Input matrix (assumed to be nodal forces and moments).
    @return C Output matrix (nodal velocities and displacements/rotations).
    """
    
    # Init zero-valued vectors
    PosDotDef  = np.zeros((XBINPUT.NumNodesTot,3),ct.c_double,'F')
    PsiDotDef  = np.zeros((XBINPUT.NumElems,Settings.MaxElNod,3),
                           ct.c_double, 'F')
    PosDotDotDef  = np.zeros((XBINPUT.NumNodesTot,3),ct.c_double,'F')
    PsiDotDotDef  = np.zeros((XBINPUT.NumElems,Settings.MaxElNod,3),
                           ct.c_double, 'F')
    ForcedVel = np.zeros((1,6),ct.c_double,'F') # check this definition as 2d array is OK within interface
    ForcedVelDot = np.zeros((1,6),ct.c_double,'F')
    
    # Init matrices and sparse counter for fortran
    MglobalFull = np.zeros((NumDof.value,NumDof.value), ct.c_double, 'F')
    CglobalFull = np.zeros((NumDof.value,NumDof.value), ct.c_double, 'F') 
    KglobalFull = np.zeros((NumDof.value,NumDof.value), ct.c_double, 'F') 
    FglobalFull = np.zeros((NumDof.value,NumDof.value), ct.c_double, 'F')
    ms = ct.c_int()
    cs = ct.c_int()
    ks = ct.c_int()
    fs = ct.c_int()
    Mvel = np.zeros((NumDof.value,6), ct.c_double, 'F')
    Cvel = np.zeros((NumDof.value,6), ct.c_double, 'F')
    Qglobal = np.zeros(NumDof.value, ct.c_double, 'F')
    
    # Init under assumption the the inertial frame and a-frame are coincident
    Cao = np.zeros((3,3),ct.c_double, 'F')
    for i in range(3):
        Cao[i,i] = 1.0
    
    BeamLib.Cbeam3_Asbly_Dynamic(XBINPUT,
                                 NumNodes_tot,
                                 XBELEM,
                                 XBNODE,
                                 PosIni,
                                 PsiIni,
                                 PosDefor,
                                 PsiDefor,
                                 PosDotDef, #zero
                                 PsiDotDef, #zero
                                 PosDotDotDef, #zero
                                 PsiDotDotDef, #zero
                                 StaticForces, # check if I affect matrices
                                 ForcedVel[0,:], # check if I make a difference
                                 ForcedVelDot[0,:], # zero
                                 NumDof,
                                 Settings.DimMat,
                                 ms,
                                 MglobalFull,
                                 Mvel,
                                 cs,
                                 CglobalFull, # check if non-zero?
                                 Cvel,
                                 ks,
                                 KglobalFull,
                                 fs,
                                 FglobalFull, # make sure forces are applied in B-frame
                                 Qglobal, # check how this is defined but I think it should be zero
                                 XBOPTS,
                                 Cao)
    
    printBlocks(KglobalFull,6,6)
    printBlocks(MglobalFull,6,6)

if __name__ == '__main__':
    # Define Patil's wing FE model
    
    # beam options.
    XBOPTS = DerivedTypes.Xbopts(FollowerForce = ct.c_bool(False),
                                 MaxIterations = ct.c_int(50),
                                 PrintInfo = ct.c_bool(True),
                                 NumLoadSteps = ct.c_int(25),
                                 Solution = ct.c_int(112),
                                 MinDelta = ct.c_double(1e-4))
     
    # beam inputs.
    XBINPUT = DerivedTypes.Xbinput(3,2)
    XBINPUT.BeamLength = 16.0
    XBINPUT.BeamStiffness[0,0] = 1.0e+09
    XBINPUT.BeamStiffness[1,1] = 1.0e+09
    XBINPUT.BeamStiffness[2,2] = 1.0e+09
    XBINPUT.BeamStiffness[3,3] = 1.0e+04
    XBINPUT.BeamStiffness[4,4] = 2.0e+04
    XBINPUT.BeamStiffness[5,5] = 5.0e+06
     
    XBINPUT.BeamStiffness[:,:] = 1.0*XBINPUT.BeamStiffness[:,:]
     
    XBINPUT.BeamMass[0,0] = 0.75
    XBINPUT.BeamMass[1,1] = 0.75
    XBINPUT.BeamMass[2,2] = 0.75
    XBINPUT.BeamMass[3,3] = 0.1
    XBINPUT.BeamMass[4,4] = 0.001
    XBINPUT.BeamMass[5,5] = 0.001
    
    # Initialize beam.
    XBINPUT, XBOPTS, NumNodes_tot, XBELEM, PosIni, PsiIni, XBNODE, NumDof \
                    = BeamInit.Static(XBINPUT,XBOPTS)
    # Set initial conditions as undef config.
    PosDefor = PosIni.copy(order='F')
    PsiDefor = PsiIni.copy(order='F')
    
    genSSbeam(XBINPUT,
              NumNodes_tot,
              NumDof,
              XBELEM,
              XBNODE,
              PosIni,
              PsiIni,
              PosDefor,
              PsiDefor,
              XBINPUT.ForceStatic, # initialised as zero
              XBOPTS)