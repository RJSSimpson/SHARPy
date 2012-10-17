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
  integer	:: r_pr, V_pr, Psi_pr, DeltaCurved_pr
  !I/O variables for computational routine
  real(8)   :: DeltaCurved
  real(8)   :: r(2,3), V(3), Psi(3)
  
  
!!!!!!!!!!!!!!!!!!!!!!
!     Read input     !
!!!!!!!!!!!!!!!!!!!!!!

  !Check for proper number of arguments. 
  if(nrhs .ne. 3) then
    call mexErrMsgTxt('Three inputs required.')
  elseif(nlhs .ne. 1) then
    call mexErrMsgTxt('One output required.')
  endif

  !Position vector of initial and final node
  r_pr = mxGetPr(prhs(1))
  m = mxGetM(prhs(1))
  n = mxGetN(prhs(1))
  size = m*n
  call mxCopyPtrToReal8(r_pr,r,size)

  !Cross-sectional orientation vector
  V_pr = mxGetPr(prhs(2)) 
  m = mxGetM(prhs(2))
  n = mxGetN(prhs(2))
  size = m*n
  call mxCopyPtrToReal8(V_pr,V,size)
  
  !Minimum angle for two unit vectors to be parallel
  DeltaCurved_pr = mxGetPr(prhs(3)) 
  m = mxGetM(prhs(3))
  n = mxGetN(prhs(3))
  size = m*n
  call mxCopyPtrToReal8(DeltaCurved_pr,DeltaCurved,size)
   
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     Call computational routine     !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !Create matrices for the return argument
  Psi = 0.d0
  
  !Call the computational subroutine.
  call bgeom_elemframe_wrap (r,V,Psi,DeltaCurved)


!!!!!!!!!!!!!!!!!!!!!!!!
!     Write output     !
!!!!!!!!!!!!!!!!!!!!!!!!

  m = 1
  n = 3
  size = m*n
  plhs(1) = mxCreateDoubleMatrix(m,n,0)  
  Psi_pr = mxGetPr(plhs(1))
  call mxCopyReal8ToPtr(Psi,Psi_pr,size)

   
!!!!!!!!!!!!!!!!!!
!     Finish     !
!!!!!!!!!!!!!!!!!! 
  
  return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-> Subroutine BGEOM_ELEMFRAME
!
!-> Description:
!
!    Compute element frame as defined by MSC conventions.
!
!-> Remarks.-
!
!  1) Element frame is the 'a' frame in the co-rotational approach.
!
!  2) X axis is defined by the two end nodes of the elements.  The X-Z plane
!     is orthogonal to the given V vector.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine bgeom_elemframe_wrap (r,V,Psi,Delta)
  use lib_bgeom

! I/O Variables.
  real(8),intent(in) :: r(2,3)         ! Position vector of initial and final node.
  real(8),intent(in) :: V(3)           ! Cross-sectional orientation vector.
  real(8),intent(out):: Psi(3)         ! CRV of the element orientation.
  real(8),intent(in) :: Delta          ! Minimum angle for two unit vectors to be parallel.

  call bgeom_elemframe (r,V,Psi,Delta)

  return
 end subroutine bgeom_elemframe_wrap