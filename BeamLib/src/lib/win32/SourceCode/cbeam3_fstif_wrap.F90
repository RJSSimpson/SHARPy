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
  integer	:: NE_pr, r0_pr, Ri_pr, ElemStiff_pr, Qstiff_pr, NumGauss_pr, Qstiff_out_pr
  !I/O variables for computational routine  
  integer	:: NE, NumGauss, MaxElNod
  real(8)	:: r0(3,6), Ri(3,6), ElemStiff(6,6), Qstiff(18)
  
  
!!!!!!!!!!!!!!!!!!!!!!
!     Read input     !
!!!!!!!!!!!!!!!!!!!!!!

  !Check for proper number of arguments
  if(nrhs .ne. 6) then
    call mexErrMsgTxt('Six inputs required.')
  elseif(nlhs .ne. 1) then
    call mexErrMsgTxt('One output required.')
  endif

  !Number of nodes in the element
  NE_pr = mxGetPr(prhs(1))
  m = mxGetM(prhs(1))
  n = mxGetN(prhs(1))
  size = m*n
  call mxCopyPtrToInteger4(NE_pr,NE,size)

  !Initial coordinates/CRV of nodes in the element
  r0_pr = mxGetPr(prhs(2)) 
  m = mxGetM(prhs(2))
  n = mxGetN(prhs(2))
  size = m*n
  !allocate(r0(m,n))
  call mxCopyPtrToReal8(r0_pr,r0,size)
  
  !Current coordinates/CRV of nodes in the element
  Ri_pr = mxGetPr(prhs(3)) 
  m = mxGetM(prhs(3))
  n = mxGetN(prhs(3))
  size = m*n
  !allocate(Ri(m,n))
  call mxCopyPtrToReal8(Ri_pr,Ri,size)

  !Element stiffness
  ElemStiff_pr = mxGetPr(prhs(4))
  m = mxGetM(prhs(4))
  n = mxGetN(prhs(4))
  size = m*n   
  !allocate(ElemStiff(m,n))
  call mxCopyPtrToReal8(ElemStiff_pr,ElemStiff,size)
  
  !Geometric tangent stiffness matrix
  Qstiff_pr = mxGetPr(prhs(5))
  m = mxGetM(prhs(5))
  n = mxGetN(prhs(5))
  size = m*n
  !allocate(Qstiff(m))
  call mxCopyPtrToReal8(Qstiff_pr,Qstiff,size)

  !Number of Gauss points
  NumGauss_pr = mxGetPr(prhs(6))
  m = mxGetM(prhs(6))
  n = mxGetN(prhs(6))
  size = m*n
  call mxCopyPtrToInteger4(NumGauss_pr,NumGauss,size)  
  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     Call computational routine     !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !Create matrices for the return argument
  !Tangent stiffness matrix of element
  MaxElNod = 3
  m = MaxElNod*6
  n = 1
  size = m*n
  plhs(1) = mxCreateDoubleMatrix(m,n,0)  
  Qstiff_out_pr = mxGetPr(plhs(1))
  
  !Call the computational subroutine.
  call cbeam3_fstif_wrap(NE,r0,Ri,ElemStiff,Qstiff,NumGauss)


!!!!!!!!!!!!!!!!!!!!!!!!
!     Write output     !
!!!!!!!!!!!!!!!!!!!!!!!!

  MaxElNod = 3
  m = MaxElNod*6
  n = 1
  size = m*n
  call mxCopyReal8ToPtr(Qstiff,Qstiff_out_pr,size)

   
!!!!!!!!!!!!!!!!!!
!     Finish     !
!!!!!!!!!!!!!!!!!! 
  
  !deallocate(r0,Ri,ElemStiff,Qstiff)   

  return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-> Subroutine CBEAM3_FSTIF
!
!-> Description:
!
!    Compute element stiffness forces.
!
!-> Remarks.-
!
!   1) New values are added to those already in the vector.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine cbeam3_fstif_wrap (NumNodesElem,r0,Ri,ElemStiff,Qstiff,NumGauss)
  use lib_cbeam3
  
! I/O Variables.
  integer,intent(in)   :: NumNodesElem      ! Number of nodes in the element.
  real(8),intent(in)   :: r0       (3,6)    ! Initial position/orientation of grid points.
  real(8),intent(in)   :: Ri       (3,6)    ! Current position/orientation of grid points.
  real(8),intent(in)   :: ElemStiff(6,6)    ! Stiffness properties in the element.
  real(8),intent(inout):: Qstiff   (18)     ! Discrete generalized stiffness forces.
  integer,intent(in)   :: NumGauss          ! Number of Gauss points in the element.
  
  call cbeam3_fstif (NumNodesElem,r0,Ri,ElemStiff,Qstiff,NumGauss)

  return
 end subroutine cbeam3_fstif_wrap