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
  integer	:: NumNodesElem_pr, r0_pr, Ri_pr, RiDDot_pr, VrelDot_pr, ElemMass_pr, Kmass_pr, NumGauss_pr
  !I/O variables for computational routine  
  integer	:: NumNodesElem, NumGauss, MaxElNod
  real(8)	:: r0(3,6), Ri(3,6), RiDDot(3,6), VrelDot(6), ElemMass(6,6), Kmass(6,18)  
  
  
!!!!!!!!!!!!!!!!!!!!!!
!     Read input     !
!!!!!!!!!!!!!!!!!!!!!!

  !Check for proper number of arguments
  if(nrhs .ne. 8) then
    call mexErrMsgTxt('Eight inputs required.')
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
  ElemMass_pr = mxGetPr(prhs(6))
  m = mxGetM(prhs(6))
  n = mxGetN(prhs(6))
  size = m*n   
  call mxCopyPtrToReal8(ElemMass_pr,ElemMass,size)
  
  !Tangent gyroscopic stiffness matrix
  Kmass_pr = mxGetPr(prhs(7))
  m = mxGetM(prhs(7))
  n = mxGetN(prhs(7))
  size = m*n
  call mxCopyPtrToReal8(Kmass_pr,Kmass,size)

  !Number of Gauss points
  NumGauss_pr = mxGetPr(prhs(8))
  m = mxGetM(prhs(8))
  n = mxGetN(prhs(8))
  size = m*n
  call mxCopyPtrToInteger4(NumGauss_pr,NumGauss,size)  
  

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
  call xbeam_kmass_wrap (NumNodesElem,r0,Ri,RiDDot,VrelDot,ElemMass,Kmass,NumGauss)


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
!-> Subroutine xbeam_KMASS
!
!-> Description:
!
!    Compute tangent damping matrix for motions of reference system.
!
!-> Remarks.-
!
!   1) New values are added to those already in the matrix.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine xbeam_kmass_wrap (NumNodesElem,r0,Ri,RiDDot,VrelDot,ElemMass,Kmass,NumGauss)
  use lib_xbeam
  
! I/O Variables.
  integer,intent(in)    :: NumNodesElem     ! Number of nodes in the element.
  real(8),intent(in)    :: r0       (3,6)   ! Initial position/orientation of grid points.
  real(8),intent(in)    :: Ri       (3,6)   ! Current position/orientation of grid points.
  real(8),intent(in)    :: RiDDot   (3,6)   ! Second time derivative of Ri.
  real(8),intent(in)    :: VrelDot  (6)     ! Currnet time derivative of Vrel.
  real(8),intent(in)    :: ElemMass (6,6)   ! Inertia properties in the element.
  real(8),intent(inout) :: Kmass    (6,18)  ! Tangent gyroscopic stiffness matrix.
  integer,intent(in)    :: NumGauss         ! Number of Gauss points in the element.

  call xbeam_kmass (NumNodesElem,r0,Ri,RiDDot,VrelDot,ElemMass,Kmass,NumGauss)
  
  return
 end subroutine xbeam_kmass_wrap
