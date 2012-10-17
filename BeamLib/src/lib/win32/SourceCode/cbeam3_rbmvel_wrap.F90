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
  integer	:: NumNodesElem_pr, r0_pr, Ri_pr, NodalMass_pr, Mvel_pr
  !I/O variables for computational routine  
  integer	:: NumNodesElem, NumGauss, MaxElNod
  real(8)	:: r0(3,6), Ri(3,6), NodalMass(3,6,6), Mvel(18,6)

  
!!!!!!!!!!!!!!!!!!!!!!
!     Read input     !
!!!!!!!!!!!!!!!!!!!!!!

  !Check for proper number of arguments
  if(nrhs .ne. 5) then
    call mexErrMsgTxt('Five inputs required.')
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

  !Element inertial properties
  NodalMass_pr = mxGetPr(prhs(4))
  m = mxGetM(prhs(4))
  n = mxGetN(prhs(4))
  size = m*n   
  call mxCopyPtrToReal8(NodalMass_pr,NodalMass,size)
  
  !Element reference-system mass matrix
  Mvel_pr = mxGetPr(prhs(5))
  m = mxGetM(prhs(5))
  n = mxGetN(prhs(5))
  size = m*n
  call mxCopyPtrToReal8(Mvel_pr,Mvel,size)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     Call computational routine     !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !Create matrices for the return argument
  !Tangent stiffness matrix of element
  MaxElNod = 3
  m = MaxElNod*6
  n = 6
  size = m*n
  plhs(1) = mxCreateDoubleMatrix(m,n,0)  
  Mvel_pr = mxGetPr(plhs(1))
  
  !Call the computational subroutine.
  call cbeam3_rbmvel_wrap (NumNodesElem,r0,Ri,NodalMass,Mvel)


!!!!!!!!!!!!!!!!!!!!!!!!
!     Write output     !
!!!!!!!!!!!!!!!!!!!!!!!!

  MaxElNod = 3
  m = MaxElNod*6
  n = 6
  size = m*n
  call mxCopyReal8ToPtr(Mvel,Mvel_pr,size)

   
!!!!!!!!!!!!!!!!!!
!     Finish     !
!!!!!!!!!!!!!!!!!! 
  
  return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-> Subroutine CBEAM3_RBMVEL
!
!-> Description:
!
!    Compute the contribution of lumped masses to the mass matrix due to
!    the reference system accelerations.
!
!-> Remarks.-
!
!   1) New values are added to those already in the matrix.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine cbeam3_rbmvel_wrap (NumNodesElem,r0,Ri,NodalMass,Mvel)
  use lib_cbeam3
  use lib_rot
  use lib_rotvect

! I/O Variables.
  integer,intent(in)   :: NumNodesElem      ! Number of nodes in the element.
  real(8),intent(in)   :: r0       (3,6)    ! Initial position/orientation of grid points.
  real(8),intent(in)   :: Ri       (3,6)    ! Current position/orientation of grid points.
  real(8),intent(in)   :: NodalMass(3,6,6)  ! Inertia tensor of lumped masses at element nodes.
  real(8),intent(inout):: Mvel     (18,6)    ! Tangent mass matrix.

  call cbeam3_rbmvel (NumNodesElem,r0,Ri,NodalMass,Mvel)

  return
 end subroutine cbeam3_rbmvel_wrap