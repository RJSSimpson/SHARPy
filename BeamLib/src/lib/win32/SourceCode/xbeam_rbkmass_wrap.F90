subroutine mexFunction(nlhs, plhs, nrhs, prhs)
 
  !I/O for gateway routine
  integer,intent(in)	:: nrhs,nlhs
  integer,intent(in)	:: prhs(*)
  integer,intent(inout)	:: plhs(*)
  !External functions
  integer	:: mxCreateDoubleMatrix, mxGetPr, mxCreateNumericArray, mxClassIDFromClassName, mxIsNumeric
  integer	:: mxGetM, mxGetN
  integer	:: m, n, size 
  !Pointers
  integer	:: NumNodesElem_pr, r0_pr, Ri_pr, RiDDot_pr, VrelDot_pr, NodalMass_pr, Kmass_pr
  !I/O variables for computational routine  
  integer	:: NumNodesElem, NumGauss, MaxElNod
  real(8)	:: r0(3,6), Ri(3,6), RiDDot(3,6), VrelDot(6), NodalMass(3,6,6), Kmass(6,18)

 
!!!!!!!!!!!!!!!!!!!!!!
!     Read input     !
!!!!!!!!!!!!!!!!!!!!!!

  !Check for proper number of arguments
  if(nrhs .ne. 7) then
    call mexErrMsgTxt('Seven inputs required.')
  elseif(nlhs .ne. 1) then
    call mexErrMsgTxt('One output required.')
  endif

  !Number of nodes in the element
  NumNodesElem_pr = mxGetPr(prhs(1))
  m = mxGetM(prhs(1))
  n = mxGetN(prhs(1))
  size = m*n
  call mxCopyPtrToInteger4(NumNodesElem_pr,NumNodesElem,size)

  !Initial coordinates/CRV of nodes in the element
  r0_pr = mxGetPr(prhs(2)) 
  m = mxGetM(prhs(2))
  n = mxGetN(prhs(2))
  size = m*n
  call mxCopyPtrToReal8(r0_pr,r0,size)
  
  !Current coordinates/CRV of nodes in the element
  Ri_pr = mxGetPr(prhs(3)) 
  m = mxGetM(prhs(3))
  n = mxGetN(prhs(3))
  size = m*n
  call mxCopyPtrToReal8(Ri_pr,Ri,size)

  !Second time derivative of Ri
  RiDDot_pr = mxGetPr(prhs(4)) 
  m = mxGetM(prhs(4))
  n = mxGetN(prhs(4))
  size = m*n
  call mxCopyPtrToReal8(RiDDot_pr,RiDDot,size)
  
  !Current time derivative of Vrel
  VrelDot_pr = mxGetPr(prhs(5)) 
  m = mxGetM(prhs(5))
  n = mxGetN(prhs(5))
  size = m*n
  call mxCopyPtrToReal8(VrelDot_pr,VrelDot,size)  

  !Element inertial properties
  NodalMass_pr = mxGetPr(prhs(6))
  m = mxGetM(prhs(6))
  n = mxGetN(prhs(6))
  size = m*n   
  call mxCopyPtrToReal8(NodalMass_pr,NodalMass,size)
  
  !Tangent gyroscopic stiffness matrix
  Kmass_pr = mxGetPr(prhs(7))
  m = mxGetM(prhs(7))
  n = mxGetN(prhs(7))
  size = m*n
  call mxCopyPtrToReal8(Kmass_pr,Kmass,size) 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     Call computational routine     !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !Create matrices for the return argument
  !Tangent stiffness matrix of element
  MaxElNod = 3
  m = 6
  n = MaxElNod*6
  size = m*n
  plhs(1) = mxCreateDoubleMatrix(m,n,0)  
  Kmass_pr = mxGetPr(plhs(1))
  
  !Call the computational subroutine.
  call xbeam_rbkmass_wrap (NumNodesElem,r0,Ri,RiDDot,VrelDot,NodalMass,Kmass)


!!!!!!!!!!!!!!!!!!!!!!!!
!     Write output     !
!!!!!!!!!!!!!!!!!!!!!!!!

  MaxElNod = 3
  m = 6
  n = MaxElNod*6
  size = m*n
  call mxCopyReal8ToPtr(Kmass,Kmass_pr,size)

   
!!!!!!!!!!!!!!!!!!
!     Finish     !
!!!!!!!!!!!!!!!!!! 
  

  return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-> Subroutine XBEAM_RBKMASS
!
!-> Description:
!
!    Compute the lumped mass contributions to the tangent damping matrix
!    for motions of reference system.
!
!-> Remarks.-
!
!   1) New values are added to those already in the matrix.
!   2) The lumped masses are defined by the inertia tensor of the nonstructural masses
!      in the section of the node, and expressed with respect to the local B frame.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine xbeam_rbkmass_wrap (NumNodesElem,r0,Ri,RiDDot,VrelDot,NodalMass,Kmass)
  use lib_xbeam
  use lib_fem
  use lib_rot
  use lib_rotvect

! I/O Variables.
  integer,intent(in)   :: NumNodesElem     ! Number of nodes in the element.
  real(8),intent(in)   :: r0       (3,6)   ! Initial position/orientation of grid points.
  real(8),intent(in)   :: Ri       (3,6)   ! Current position/orientation of grid points.
  real(8),intent(in)   :: RiDDot   (3,6)   ! Current time derivative of position/CRV of points.
  real(8),intent(in)   :: VrelDot  (6)     ! Linear/angular acceleration of reference frame.
  real(8),intent(in)   :: NodalMass(3,6,6) ! Inertia tensor of lumped masses at element nodes.
  real(8),intent(inout):: Kmass     (6,18)  ! Tangent gyroscopic damping matrix.

  call xbeam_rbkmass (NumNodesElem,r0,Ri,RiDDot,VrelDot,NodalMass,Kmass)
 
  return
 end subroutine xbeam_rbkmass_wrap
