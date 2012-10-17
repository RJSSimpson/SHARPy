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
  integer	:: NumNodesElem_pr, TreeConn_pr, Psi02_pr, AllPsi_pr, SB2B1_pr, Telem_pr
  !I/O variables for computational routine  
  integer	:: NE, NumNodesElem
  integer	:: TreeConn(3,2)
  real(8)	:: Psi02(3,3),  SB2B1(18,18)
  real(8),allocatable	:: AllPsi(:,:,:)
  
  
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
  Psi02_pr = mxGetPr(prhs(3)) 
  m = mxGetM(prhs(3))
  n = mxGetN(prhs(3))
  size = m*n
  call mxCopyPtrToReal8(Psi02_pr,Psi02,size)

  !Initial CRV from global to nodes in all elements
  AllPsi_pr = mxGetPr(prhs(4))
  m = mxGetM(prhs(4))  
  size = m*3*3
  NE = m
  allocate(AllPsi(m,3,3))
  call mxCopyPtrToReal8(AllPsi_pr,AllPsi,size)
  
  !Element transformation matrix
  SB2B1_pr = mxGetPr(prhs(5))
  m = mxGetM(prhs(5))
  n = mxGetN(prhs(5))
  size = m*n
  call mxCopyPtrToReal8(SB2B1_pr,SB2B1,size)
  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     Call computational routine     !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !Create matrices for the return argument
  !Element transformation matrix
  MaxElNod = 3
  m = MaxElNod*6
  n = MaxElNod*6
  size = m*n
  plhs(1) = mxCreateDoubleMatrix(m,n,0)  
  Telem_pr = mxGetPr(plhs(1))
  
  !Call the computational subroutine.
  call cbeam3_projs2m_wrap (NumNodesElem,TreeConn,Psi02,AllPsi,SB2B1,NE)


!!!!!!!!!!!!!!!!!!!!!!!!
!     Write output     !
!!!!!!!!!!!!!!!!!!!!!!!!

  MaxElNod = 3
  m = MaxElNod*6
  n = MaxElNod*6
  size = m*n
  call mxCopyReal8ToPtr(SB2B1,Telem_pr,size)

   
!!!!!!!!!!!!!!!!!!
!     Finish     !
!!!!!!!!!!!!!!!!!! 
  
  deallocate(AllPsi)   

  return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-> Subroutine CBEAM3_PROJS2M
!
!-> Description:
!
!     Compute transformation from slave to master nodal equations. This is done to
!     include discontinuities between adjacent elements
!
!-> Remarks.-
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine cbeam3_projs2m_wrap (NumNodesElem,TreeConn,Psi02,AllPsi,Telem,NE)
  use lib_cbeam3

! I/O Variables.
  integer,intent(in) :: NE                 ! Total number of elements.
  integer,intent(in) :: NumNodesElem       ! Number of nodes in the element.
  integer,intent(in) :: TreeConn  (3,2)    ! Connectivity tree.
  real(8),intent(in) :: Psi02     (3,3)    ! Initial CRV from global to nodes in element.
  real(8),intent(in) :: AllPsi    (NE,3,3)  ! Initial CRV from global to nodes in all elements.
  real(8),intent(out):: Telem     (18,18)  ! Element transformation matrix.

  call cbeam3_projs2m (NumNodesElem,TreeConn,Psi02,AllPsi,Telem)

  return
 end subroutine cbeam3_projs2m_wrap