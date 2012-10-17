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
  integer	:: NumNodesElem_pr, TreeConn_pr, Psi0_pr, AllPsi_pr, Psi_pr, PsiOut_pr
  !I/O variables for computational routine  
  integer	:: NE, NumNodesElem
  integer	:: TreeConn(3,2)
  real(8)	:: Psi0(3,3), Psi(3,3)
  real(8), allocatable :: AllPsi(:,:,:)
  
  
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

  !Connectivity tree
  TreeConn_pr = mxGetPr(prhs(2)) 
  m = mxGetM(prhs(2))
  n = mxGetN(prhs(2))
  size = m*n
  call mxCopyPtrToInteger4(TreeConn_pr,TreeConn,size)
  
  !Initial CRV from global to nodes in element
  Psi0_pr = mxGetPr(prhs(3)) 
  m = mxGetM(prhs(3))
  n = mxGetN(prhs(3))
  size = m*n
  call mxCopyPtrToReal8(Psi0_pr,Psi0,size)

  !Initial CRV from global to nodes in all elements
  AllPsi_pr = mxGetPr(prhs(4))
  m = mxGetM(prhs(4))  
  size = m*3*3
  NE = m
  allocate(AllPsi(m,3,3))
  call mxCopyPtrToReal8(AllPsi_pr,AllPsi,size)
  
  !Current CRV in element
  Psi_pr = mxGetPr(prhs(5))
  m = mxGetM(prhs(5))  
  n = mxGetN(prhs(5)) 
  size = m*n
  call mxCopyPtrToReal8(Psi_pr,Psi,size)
  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     Call computational routine     !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !Create matrices for the return argument
  !Element transformation matrix
  m = 3
  n = 3
  size = m*n
  plhs(1) = mxCreateDoubleMatrix(m,n,0)  
  PsiOut_pr = mxGetPr(plhs(1))
  
  !Call the computational subroutine.
  call cbeam3_projm2s_wrap (NumNodesElem,TreeConn,Psi0,AllPsi,Psi,NE)

!!!!!!!!!!!!!!!!!!!!!!!!
!     Write output     !
!!!!!!!!!!!!!!!!!!!!!!!!

  m = 3
  n = 3
  size = m*n
  call mxCopyReal8ToPtr(Psi,PsiOut_pr,size)

   
!!!!!!!!!!!!!!!!!!
!     Finish     !
!!!!!!!!!!!!!!!!!! 
  
  deallocate(AllPsi)   

  return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-> Subroutine CBEAM3_PROJM2S
!
!-> Description:
!
!     Compute transformation for common element equations. This is done to
!     include discontinuities between adjacent elements
!
!     Project rotation from master to slave grid node.
!
!-> Remarks.-
!      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine cbeam3_projm2s_wrap (NumNodesElem,TreeConn,Psi0,AllPsi,Psi,NE)
  use lib_cbeam3

! I/O Variables.
  integer            :: NE
  integer,intent(in) :: NumNodesElem       ! Number of nodes in the element.
  integer,intent(in) :: TreeConn  (3,2)    ! Connectivity tree.
  real(8),intent(in) :: Psi0      (3,3)    ! Initial CRV from global to nodes in element.
  real(8),intent(in) :: AllPsi    (NE,3,3) ! Initial CRV from global to nodes in all elements.
  real(8),intent(out):: Psi       (3,3)    ! Current CRV in element.

  call cbeam3_projm2s (NumNodesElem,TreeConn,Psi0,AllPsi,Psi)

  return
 end subroutine cbeam3_projm2s_wrap