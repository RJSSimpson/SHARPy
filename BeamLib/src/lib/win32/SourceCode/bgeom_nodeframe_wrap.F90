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
  integer	:: NumNodesElem_pr, Coords_pr, Phi_pr, V_pr, Psi_pr, Delta_pr
  !I/O variables for computational routine
  integer   :: NumNodesElem
  real(8)   :: Delta
  real(8)   :: Coords(3,3), Phi(3), V(3), Psi(3,3)
  
  
!!!!!!!!!!!!!!!!!!!!!!
!     Read input     !
!!!!!!!!!!!!!!!!!!!!!!

  !Check for proper number of arguments. 
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

  !Coordinates of the grid points
  Coords_pr = mxGetPr(prhs(2)) 
  m = mxGetM(prhs(2))
  n = mxGetN(prhs(2))
  size = m*n
  call mxCopyPtrToReal8(Coords_pr,Coords,size)
  
  !Pretwist at grid points
  Phi_pr = mxGetPr(prhs(3)) 
  m = mxGetM(prhs(3))
  n = mxGetN(prhs(3))
  size = m*n
  call mxCopyPtrToReal8(Phi_pr,Phi,size)
  
  !Element cross-sectional orientation
  V_pr = mxGetPr(prhs(4)) 
  m = mxGetM(prhs(4))
  n = mxGetN(prhs(4))
  size = m*n
  call mxCopyPtrToReal8(V_pr,V,size)
  
  !Minimum angle for two unit vectors to be parallel
  Delta_pr = mxGetPr(prhs(5)) 
  m = mxGetM(prhs(5))
  n = mxGetN(prhs(5))
  size = m*n
  call mxCopyPtrToReal8(Delta_pr,Delta,size)
   
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     Call computational routine     !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !Create matrices for the return argument
  !Psi = 0.d0
  
  !Call the computational subroutine.
  call bgeom_nodeframe_wrap (NumNodesElem,Coords,Phi,V,Psi,Delta)


!!!!!!!!!!!!!!!!!!!!!!!!
!     Write output     !
!!!!!!!!!!!!!!!!!!!!!!!!

  m = 3
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
!-> Subroutine BGEOM_NODEFRAME
!
!-> Description:
!
!    For a given node, compute undeformed coordinate system ('b' frame) at 
!    each of the element nodes.
!
!-> Remarks.-
!
!  1) V is only used for collinear points in an element. This is defined
!     by the convergence parameter (Delta).
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine bgeom_nodeframe_wrap (NumNodesElem,Coords,Phi,V,Psi,Delta)
  use lib_bgeom

! I/O Variables.
  integer,intent(in) :: NumNodesElem   ! Number of nodes in the element.
  real(8),intent(in) :: Coords  (3,3)  ! Coordinates of the grid points.
  real(8),intent(in) :: Phi     (3)    ! Pretwist at grid points.
  real(8),intent(in) :: V       (3)    ! Element cross-sectional orientation.
  real(8),intent(out):: Psi     (3,3)  ! CRV at each node.
  real(8),intent(in) :: Delta          ! Minimum angle for two unit vectors to be parallel.

  call bgeom_nodeframe (NumNodesElem,Coords,Phi,V,Psi,Delta)

  return
 end subroutine bgeom_nodeframe_wrap