!-> Copyright by Imperial College London, 2008
!
!-> Module.- XBEAM_SHARED Rafa Palacios. 15Jul2008 - Last Update 07Jan2011 Henrik Hesse
!
!-> Description.-
!
!  This module defines shared parameters and type definitions shared by all
!  XBeam routines.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module xbeam_shared
 implicit none


! Problem constants.
 integer,parameter:: MaxElNod=3                ! Max number of nodes per element.
 real(8),parameter:: Pi=3.14159265358979
 real(8),parameter,dimension(3,3):: Unit= &    ! Unit matrix.
&         reshape((/1.d0,0.d0,0.d0,0.d0,1.d0,0.d0,0.d0,0.d0,1.d0/),(/3,3/))


! Define derived types with default values:

! Element information.
 type xbelem
  integer:: NumNodes               ! Number of nodes in the element.
  integer:: MemNo                  ! Member to which the element belongs.
  integer:: Conn   (MaxElNod)      ! Connectivities.
  integer:: Master (MaxElNod,2)    ! Master node for each node j in the element.
                                   ! (j,1): # master elem  (or 0 if current is master).
                                   ! (j,2): Node within master element (or 0 if current is master).
  real(8):: Length                 ! Length in the undeformed configuration.
  real(8):: PreCurv (3)            ! Initial curvature of the element.
  real(8):: Psi     (3)            ! Rotation vector of (undeformed) element frame.
  real(8):: Vector  (3)            ! Element orientation vector. It goes along the local Y axis.
  real(8):: Mass    (6,6)          ! Mass matrix (constant along the element).
  real(8):: Stiff   (6,6)          ! Element stiffness (constant along the element).
  real(8):: InvStiff(6,6)          ! Inverse of the element stiffness (constant along the element).
  real(8):: RBMass  (MaxElNod,6,6) ! Non-Structural (lumped) mass at element nodes.
 end type xbelem

! Nodal information.
 type xbnode
  integer:: Master (2)             ! Master node for current node.
                                   ! (1): # master elem.
                                   ! (2): # node within master element.
  integer:: Vdof                   ! Corresponding degree of freedom in the velocity vector.
  integer:: Fdof                   ! Corresponding degree of freedom in the force vector.
 end type xbnode

! Simulation options (with default values).
 type xbopts
  logical:: FollowerForce   =.true.   ! =T: Follower force.
  logical:: FollowerForceRig=.true.   ! =T: Follower force in the body-fixed frame.
  logical:: PrintInfo    =.true.      ! =T: Print information on screen.
  logical:: OutInBframe  =.true.      ! =T print velocities in B-frame (if not, use a-frame)
  logical:: OutInaframe  =.false.     ! =T print velocities in a-frame (if not, Inertial frame)
  integer:: ElemProj     = 0          ! =0: Element info computed in the global frame.
                                      ! =1: Element info computed in a fixed element frame.
                                      ! =2: Element info computed in a moving element frame.
  integer:: MaxIterations=99          ! Maximum number of iterations.
  integer:: NumLoadSteps=5            ! Number of load increments.
  integer:: NumGauss=1                ! Number of Gauss points in the integration.
  integer:: Solution=111              ! Solution process:
                                      ! =102/112: cbeam3 linear/nonlinear static.
                                      ! =202/212: cbeam3 linear/nonlinear structural dynamic
                                      ! =302/312: cbeam3 linear/nonlinear static + structural dynamic
                                      ! =900/910:        linear/nonlinear rigid-body dynamic
                                      ! =902/912: cbeam3 linear/nonlinear flexible-body dynamic
                                      ! =    922: cbeam3 nonlinear static + flexible-body dynamic
                                      ! =    952: cbeam3 linear flexible with nonlinear rigid-body dynamic
  real(8):: DeltaCurved=1d-5          ! Minimum angle for two unit vectors to be parallel.
  real(8):: MinDelta=1d-8             ! Convergence parameter for Newton-Raphson iterations.
  real(8):: NewmarkDamp=1.d-4         ! Numerical damping in the Newmark integration scheme.
 end type

end module xbeam_shared
