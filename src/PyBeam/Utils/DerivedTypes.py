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
import numpy as np #http://docs.scipy.org/doc/numpy/reference/
import ctypes as ct

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
    
    @warning If FollowerForce = ct.c_bool(True), beam forces must be applied
    in the local FoR, B.
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
        """@brief Default initialisation is as defined in the original fortran
        derived type."""
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
    @param ForcingType Type of dynamic forcing.
    @param g acceleration due to gravity.
    """

    def __init__(self, NumNodesElem, NumElems,
                 BeamLength = 1.0,
                 BeamStiffness = np.zeros((6,6),ct.c_double,'F'),
                 BeamMass = np.zeros((6,6),ct.c_double,'F'),
                 BConds = 'CF',
                 Sigma = 1.0,
                 iOut = 1,
                 t0 = 0.0,
                 tfin = 0.0,
                 dt = 0.0,
                 Omega = 0.0,
                 ForcingType = 'Const',
                 RampTime = 0.0,
                 g = 0.0,
                 PsiA_G = np.array([0.0,0.0,0.0])):
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
        self.g = g
        self.PsiA_G = PsiA_G
        
        # Check number of nodes per element.
        if self.NumNodesElem != 2 and self.NumNodesElem != 3:
            sys.stderr.write('Invalid number of nodes per element\n')
        elif self.NumNodesElem == 2:
            NumNodesTot = NumElems + 1
        elif self.NumNodesElem == 3:
            NumNodesTot = 2*NumElems + 1
            
        self.NumNodesTot = NumNodesTot
        
        # Initialise nodal arrays.
        self.ForceStatic = np.zeros((NumNodesTot,6),ct.c_double,'F')
        self.ForceDyn = np.zeros((NumNodesTot,6),ct.c_double,'F')
        
        
def dump(obj):
    """@brief Prints all attributes of an object"""
    for attr in dir(obj):
        print("obj.%s = %s" % (attr, getattr(obj, attr)))

class Xbelem:
    """@brief Pythonic version of fortran arrays containing derived type 'Elem'
    data, one set of arrays for each element, as defined in xbeam_shared.f90."""
    def __init__(self, NumElems, MaxElNod = 3):
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
    XBINPUT = Xbinput(3,8)
    print(XBINPUT.__dict__)
    print('FINISHED')
    
    print('Test Xbelem class memory')
    XBELEM=Xbelem(2,3)
    print(XBELEM.Mass)