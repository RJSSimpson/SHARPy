'''@package PyBeam.Utils.DerivedTypes
@brief      A collection of derived types that mirror those in 
xbeam_shared.f90 plus extras.
@author     Rob Simpson
@contact    r.simpson11@imperial.ac.uk
@author     A.Da-Ronch (Liverpool) 
@version    2.0
@date       14/11/2012
@pre        None
@warning    Stiffness parameter Sigma not implemented.'''

import sys
import ctypes
import numpy as np #http://docs.scipy.org/doc/numpy/reference/
import ctypes as ct

def init_XBELEM(NumElems,MaxElNod):
    """@brief Initialize XBELEM dictionary
    @author A.Da Ronch
    20121001 creation
    """
    
    XBELEM = {'NumNodes': (ctypes.c_int *NumElems)(),
              'MemNo'   : (ctypes.c_int *NumElems)(),
              'Conn'    : (ctypes.c_int *(NumElems*MaxElNod))(),
              'Master'  : (ctypes.c_int *(NumElems*MaxElNod*2))(),
              'Length'  : np.zeros(NumElems, dtype=np.double, order="Fortran"), #order doesn't matter it's 1-D!
              'PreCurv' : np.zeros(3*NumElems, dtype=np.double, order="Fortran"),
              'Psi'     : np.zeros(3*NumElems, dtype=np.double, order="Fortran"),
              'Vector'  : np.zeros(3*NumElems, dtype=np.double, order="Fortran"),
              'Mass'    : np.zeros((6*NumElems,6), dtype=np.double, order="Fortran"),
              'Stiff'   : np.zeros((6*NumElems,6), dtype=np.double, order="Fortran"),
              'InvStiff': np.zeros((6*NumElems,6), dtype=np.double, order="Fortran"),
              'RBMass'  : np.zeros((MaxElNod*NumElems,6,6), dtype=np.double, \
                                                               order="Fortran")
              }
    
    return XBELEM



def init_XBNODE(NumNodes_tot):
    """
    Function   : init_XBNODE
    Description: Initialize XBNODE dictionary
    Input      :
    
    Date        Author      Note
    20121001    A.Da Ronch  Creation
    """
    
    XBNODE = {'Nod_Master': (ctypes.c_int *(2*NumNodes_tot.value))(),
              'Nod_Vdof'  : (ctypes.c_int *NumNodes_tot.value)(),
              'Nod_Fdof'  : (ctypes.c_int *NumNodes_tot.value)()
             }

    return XBNODE



def pack_XBELEM(NumNodes,MemNo,Conn,Master,Length,PreCurv,Psi,Vector,Mass,
        Stiff,InvStiff,RBMass,INFO):
    """
    Function   : pack_XBELEM
    Description: Update XBELEM dictionary
    Input      : 
    
    Date        Author      Note
    20121001    A.Da Ronch  Creation
    """
    
    if( (INFO.lower() != 'f') and
        (INFO.lower() != 'fortran') and
        (INFO.lower() != 'c') ):
        sys.stderr.write('\n\n  **** \n')
        exit(-1)

    XBELEM = {'NumNodes': NumNodes,
              'MemNo'   : MemNo,
              'Conn'    : Conn,
              'Master'  : Master,
              'Length'  : Length.copy(order=INFO),
              'PreCurv' : PreCurv.copy(order=INFO),
              'Psi'     : Psi.copy(order=INFO),
              'Vector'  : Vector.copy(order=INFO),
              'Mass'    : Mass.copy(order=INFO),
              'Stiff'   : Stiff.copy(order=INFO),
              'InvStiff': InvStiff.copy(order=INFO),
              'RBMass'  : RBMass.copy(order=INFO)
              }
    
    return XBELEM



def pack_XBNODE(Nod_Master,Nod_Vdof,Nod_Fdof):
    """
    Function   : pack_XBNODE
    Description: Update XBELEM dictionary
    Input      : 
    
    Date        Author      Note
    20121001    A.Da Ronch  Creation
    """
    
    XBNODE = {'Nod_Master': Nod_Master,
              'Nod_Vdof'  : Nod_Vdof,
              'Nod_Fdof'  : Nod_Fdof
             }
    
    return XBNODE



def unpack_XBELEM(XBELEM,INFO):
    """
    Function   : unpack_XBELEM
    Description: Update XBELEM dictionary
    Input      : 
    
    Date        Author      Note
    20121001    A.Da Ronch  Creation
    """
    
    if( (INFO.lower() != 'f') and
        (INFO.lower() != 'fortran') and
        (INFO.lower() != 'c') ):
        sys.stderr.write('\n\n  **** \n')
        exit(-1)
    
    NumNodes = XBELEM['NumNodes']
    MemNo    = XBELEM['MemNo']
    Conn     = XBELEM['Conn']
    Master   = XBELEM['Master']
    Length   = XBELEM['Length'].copy(order=INFO)
    PreCurv  = XBELEM['PreCurv'].copy(order=INFO)
    Psi      = XBELEM['Psi'].copy(order=INFO)
    Vector   = XBELEM['Vector'].copy(order=INFO)
    Mass     = XBELEM['Mass'].copy(order=INFO)
    Stiff    = XBELEM['Stiff'].copy(order=INFO)
    InvStiff = XBELEM['InvStiff'].copy(order=INFO)
    RBMass   = XBELEM['RBMass'].copy(order=INFO)
    
    #http://docs.scipy.org/doc/numpy/reference/generated/numpy.ndarray.copy.html

    return NumNodes,MemNo,Conn,Master,Length,PreCurv,Psi,Vector,Mass,Stiff, \
        InvStiff,RBMass




def unpack_XBNODE(XBNODE):
    """
    Function   : unpack_XBNODE
    Description: Update XBNODE dictionary
    Input      : 
    
    Date        Author      Note
    20121001    A.Da Ronch  Creation
    """
    
    Nod_Master = XBNODE['Nod_Master']
    Nod_Vdof   = XBNODE['Nod_Vdof']
    Nod_Fdof   = XBNODE['Nod_Fdof']
    
    return Nod_Master,Nod_Vdof,Nod_Fdof

class Xbopts:
    """@brief Simple struct-like class containing xbeam options from 
    xbeam_shared.f90.
    
    @param FollowerForce Follower force flag, bool.
    @param FollowerForceRigid Follower force in body-fixed frame flag, bool.
    @param PrintInfo Print information to stdout, bool.
    @param OutInBframe Print velocities in B-frame (if not, use a-frame)
    @param OutInaframe Print velocities in a-frame (if not, use inertial frame)
    @param ElemProj Element info computed in (1) global frame (2) fixed element
                    frame (3) moving element frame.
    @param MaxIterations Maximum number of iterations, int.
    @param NumLoadSteps Number of load increments, int.
    @param NumGauss Number of Gauss points in element integration
    @param Solution Solution process 
                = 102/112: cbeam3 linear/nonlinear static
                = 202/212: cbeam3 linear/nonlinear structural dynamic
                = 302/312: cbeam3 linear/nonlinear static + structural dynamic
                = 900/910:        linear/nonlinear rigid-body dynamic
                = 902/912: cbeam3 linear/nonlinear flexible-body dynamic
                =     922: cbeam3 nonlinear static + flexible-body dynamic
                =     952: cbeam3 linear flexible with nonlinear rigid-body dyn
    @param DeltaCurved Min. angle for two vectors to be parallel, double.
    @param MinDelta Convergence param for Newton-Rhaphson iterations, double.
    @param NewmarkDamp Numerical damping in Newmark integration scheme, double.
    """

    def __init__(self, FollowerForce = ct.c_bool(False), \
                 FollowerForceRig = ct.c_bool(True), \
                 PrintInfo = ct.c_bool(False), OutInBframe = ct.c_bool(True), \
                 OutInaframe = ct.c_bool(False), ElemProj = ct.c_int(0), \
                 MaxIterations = ct.c_int(99), NumLoadSteps = ct.c_int(5), \
                 NumGauss = ct.c_int(1), Solution = ct.c_int(111), \
                 DeltaCurved = ct.c_double(1.0e-5), \
                 MinDelta = ct.c_double(1.0e-8), \
                 NewmarkDamp = ct.c_double(1.0e-4) ):
        """@brief Default initialisation (as defined in .f90).
        
        Non-default initialisation requires everything up to the argument in
        the list is specified (except from self)."""
        self.FollowerForce = FollowerForce
        self.FollowerForceRig = FollowerForceRig
        self.PrintInfo = PrintInfo
        self.OutInBframe = OutInBframe
        self.OutInaframe = OutInaframe
        self.ElemProj = ElemProj
        self.MaxIterations = MaxIterations
        self.NumLoadSteps = NumLoadSteps
        self.NumGauss = NumGauss
        self.Solution = Solution
        self.DeltaCurved = DeltaCurved
        self.MinDelta = MinDelta
        self.NewmarkDamp = NewmarkDamp
        
        
class Xbinput:
    """@brief Contains inputs for PyBeam functions.
    
    @param NumNodesElem Number of nodes in each beam element, either 2 or 3.
    @param NumElems Number of elements in the beam model.
    @param NumSteps Number of timesteps.
    @param BeamLength Length of beam.
    @param BeamStiffness 6-by-6 sectional stiffness matrix.
    @param BeamMass 6-by-6 sectional mass matrix.
    @param BConds 2-char string with boundary condition, 'CF' = Clamped-Free.
    @param Sigma Stiffness parameter. Not Implemented.
    @param iOut Output file.
    @param t0 Initial time.
    @param tfin Final time.
    @param dt Timestep.
    @param Omega Angular velocity of oscillatory motions.
    @param NumNodesTot total number of nodes in the model.
    @param ForceStatic NumNodes-6 static force vector at nodes.
    @param ForceDyn Numnodes-6 dynamic forces at nodes.
    @param ForcingType Type of dynamic forcing."""

    def __init__(self, NumNodesElem, NumElems,\
                 BeamLength = 1.0,\
                 BeamStiffness = np.zeros((6,6),ct.c_double,'F'),\
                 BeamMass = np.zeros((6,6),ct.c_double,'F'),\
                 BConds = 'CF',\
                 Sigma = 1.0,\
                 iOut = 1,\
                 t0 = 0.0,\
                 tfin = 0.0,\
                 dt = 0.0,\
                 Omega = 0.0,\
                 ForcingType = 'Const',
                 RampTime = 0.0):
        """@brief NumNodesElem and NumElems must be specified for initialisation
        of force arrays.
        
        @param NumNodes Local variable storing number of nodes in model for
         array init.
        """
        self.NumNodesElem = NumNodesElem
        self.NumElems = NumElems
        self.BeamLength = BeamLength
        self.BeamStiffness = BeamStiffness
        self.BeamMass = BeamMass
        self.BConds = BConds
        self.Sigma = Sigma
        self.iOut = iOut
        self.t0 = t0
        self.tfin = tfin
        self.dt = dt
        self.Omega = Omega
        self.ForcingType = ForcingType
        self.RampTime = RampTime
        
        "Check number of nodes per element"
        if self.NumNodesElem != 2 and self.NumNodesElem != 3:
            sys.stderr.write('Invalid number of nodes per element\n')
        elif self.NumNodesElem == 2:
            NumNodesTot = NumElems + 1
        elif self.NumNodesElem == 3:
            NumNodesTot = 2*NumElems + 1
            
        self.NumNodesTot = NumNodesTot
        
        "Init nodal arrays"
        self.ForceStatic = np.zeros((NumNodesTot,6),ct.c_double,'F')
        self.ForceDyn = np.zeros((NumNodesTot,6),ct.c_double,'F')
        
        
def dump(obj):
    """@brief Prints all attributes of an object"""
    for attr in dir(obj):
        print("obj.%s = %s" % (attr, getattr(obj, attr)))

class Xbelem:
    """@brief Pythonic version of fortran arrays containing derived type 'Elem'
    data, one set of arrays for each element, as defined in xbeam_shared.f90."""
    def __init__(self, NumElems, MaxElNod):
        """@brief Create 'NumElems' arrays with zero entries.
        
        @param NumNodes Number of nodes in each element.
        @param MemNo Member to which the element belongs.
        @param Conn Connectivities.
        @param Master Master node for each node j in each element:
        (j,m): Node belongs to master elem m (or 0 if current is master).
        (j,n): Node n within master element (or 0 if current is master).
        @param Length Undeformed length of each element.
        @param PreCurv Undeformed curvature of the element.
        @param Psi Undeformed rotation vector of element frame.
        @param Vector Element orientation vector - goes along the local Y axis.
        @param Mass Mass matrix (constant along the element).
        @param Stiff Stiffness matrix (constant along the element).
        @param InvStiff Inverse of stiffness matrix.
        @param RBMass Non-Structural (lumped) mass at element nodes.
        @details Memory mapping in line with f90:Elem(i)%Array(j,k,l) reference 
        requirement using existing f90:do_xbelem_var protocol in 
        Fx_Wrapper_PyFx.f90 by using the 
        np.Array((j*i,k,l),dtype=ct.*,order='F') syntax."""
        
        self.NumNodes = np.zeros(NumElems, dtype=ct.c_int, order='F')
        self.MemNo = np.zeros(NumElems, dtype=ct.c_int, order='F')
        self.Conn = np.zeros((MaxElNod*NumElems), dtype=ct.c_int, order='F')
        self.Master = np.zeros((MaxElNod*NumElems*2), dtype=ct.c_int, order='F')
        self.Length = np.zeros(NumElems, dtype=ct.c_double, order='F')
        self.PreCurv = np.zeros((3*NumElems), dtype=ct.c_double, order='F')
        self.Psi = np.zeros((3*NumElems), dtype=ct.c_double, order='F')
        self.Vector = np.zeros((3*NumElems), dtype=ct.c_double, order='f')
        self.Mass = np.zeros((6*NumElems,6), dtype=ct.c_double, order='F')
        self.Stiff = np.zeros((6*NumElems,6), dtype=ct.c_double, order='F')
        self.InvStiff = np.zeros((6*NumElems,6), dtype=ct.c_double, order='F')
        self.RBMass = np.zeros((MaxElNod*NumElems,6,6), dtype=ct.c_double, \
                                                               order='F')
        
class Xbnode:
    """@brief Pythonic nodal information as defined in xbeam_shared.f90
    
    @param Master Master node info for each node:
     (master elem, node within master elem).
    @param Vdof Degree of freedom in velocity vector.
    @param Fdof Degree of freedom in force vector."""
    def __init__(self, NumNodesTot):
        self.Master = np.zeros(2*NumNodesTot,dtype=ct.c_int,order='F')       
        self.Vdof = np.zeros(NumNodesTot,dtype=ct.c_int,order='F') 
        self.Fdof = np.zeros(NumNodesTot,dtype=ct.c_int,order='F')

if(__name__ == '__main__'):
    print('Batch run of PyAsblyStruct.py...')
    
    print('Test: Xbopts class (see script)...')

    print('Default initialisation:')
    XBOPTS = Xbopts()
    print('FollowerForce = %r' %(XBOPTS.FollowerForce))

    print('Custom initialisation:')
    XBOPTS2 = Xbopts(False)
    print('FollowerForce = %r' %(XBOPTS2.FollowerForce))
    print('Solution = %r' %(XBOPTS.Solution))

    print('Default arguments:')
    print(Xbopts.__init__.__defaults__)
    print('FINISHED')
    
    print('Test: Xbinput class (see script)...')
    XBINPUT = Xbinput()
    print(XBINPUT.__dict__)
    print('FINISHED')
    
    print('Test init_XBELEM...')
    NumElems = 3
    MaxElNod = 3
    init_XBELEM(NumElems,MaxElNod)
    print('FINISHED')
    
    print('Test Xbelem class memory')
    XBELEM=Xbelem(2,3)
    print(XBELEM.Mass)