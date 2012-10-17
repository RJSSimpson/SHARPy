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
  integer	:: PosGrids_pr, PsiGrids_pr, kb_pr
  !I/O variables for computational routine
  real(8)   :: PosGrids(2,3), PsiGrids(3,3), kb(3)
  
  
!!!!!!!!!!!!!!!!!!!!!!
!     Read input     !
!!!!!!!!!!!!!!!!!!!!!!

  !Check for proper number of arguments. 
  if(nrhs .ne. 2) then
    call mexErrMsgTxt('Two inputs required.')
  elseif(nlhs .ne. 1) then
    call mexErrMsgTxt('One output required.')
  endif

  !Position vector at nodes
  PosGrids_pr = mxGetPr(prhs(1))
  m = mxGetM(prhs(1))
  n = mxGetN(prhs(1))
  size = m*n
  call mxCopyPtrToReal8(PosGrids_pr,PosGrids,size)

  !CRV at nodes
  PsiGrids_pr = mxGetPr(prhs(2)) 
  m = mxGetM(prhs(2))
  n = mxGetN(prhs(2))
  size = m*n
  call mxCopyPtrToReal8(PsiGrids_pr,PsiGrids,size)
  
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     Call computational routine     !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !Create matrices for the return argument
  kb = 0.d0
  
  !Call the computational subroutine
  call bgeom_elemcurv2_wrap (PosGrids,PsiGrids,kb)

!!!!!!!!!!!!!!!!!!!!!!!!
!     Write output     !
!!!!!!!!!!!!!!!!!!!!!!!!

  m = 1
  n = 3
  size = m*n
  plhs(1) = mxCreateDoubleMatrix(m,n,0)  
  kb_pr = mxGetPr(plhs(1))
  call mxCopyReal8ToPtr(kb,kb_pr,size)

   
!!!!!!!!!!!!!!!!!!
!     Finish     !
!!!!!!!!!!!!!!!!!! 
  
  return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-> Subroutine BGEOM_ELEMCURV2
!
!-> Description:
!
!    Compute curvature at the center of a 2-noded beam element for given 
!    cartesian rotation vector at its nodes.
!
!-> Remarks.-
!
!  1) This is vector kb in the undeformed frame, b.
!
!  2) Element goes from node A to node B.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine bgeom_elemcurv2_wrap (PosGrids,PsiGrids,kb)
  use lib_bgeom

! I/O Variables.
  real(8),intent(in) :: PosGrids (2,3) ! Position vector at nodes.
  real(8),intent(in) :: PsiGrids (3,3) ! CRV at nodes
  real(8),intent(out):: kb     (3)     ! Curvature of the element.

  call bgeom_elemcurv2 (PosGrids,PsiGrids,kb)

  return
 end subroutine bgeom_elemcurv2_wrap