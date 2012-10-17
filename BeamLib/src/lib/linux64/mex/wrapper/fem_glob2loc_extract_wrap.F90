#include "fintrf.h"
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
  integer	:: ElemNodes_pr, GlobVector_pr, LocVector_pr, NumNodeElem_pr
  !I/O variables for computational routine  
  integer	:: NumNodes, NumNodeElem
  integer   :: ElemNodes(3)
  real(8)	:: LocVector(3,3)
  real(8), allocatable	:: GlobVector(:,:)
    
    
!!!!!!!!!!!!!!!!!!!!!!
!     Read input     !
!!!!!!!!!!!!!!!!!!!!!!

  !Check for proper number of arguments. 


  !Nodes in the current element
  ElemNodes_pr = mxGetPr(prhs(1))
  m = mxGetM(prhs(1))
  n = mxGetN(prhs(1))
  size = m*n
  call mxCopyPtrToInteger8(ElemNodes_pr,ElemNodes,size)

  !Vector of global values
  GlobVector_pr = mxGetPr(prhs(2)) 
  m = mxGetM(prhs(2))
  n = mxGetN(prhs(2))
  size = m*n
  allocate(GlobVector(m,n))
  NumNodes = m
  call mxCopyPtrToReal8(GlobVector_pr,GlobVector,size)  

  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     Call computational routine     !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !Create matrices for the return argument
  !Vector of local values
  MaxElNod = 3
  
  !Call the computational subroutine.
  call fem_glob2loc_extract_wrap(ElemNodes,GlobVector,LocVector,NumNodeElem,NumNodes,MaxElNod)
  


!!!!!!!!!!!!!!!!!!!!!!!!
!     Write output     !
!!!!!!!!!!!!!!!!!!!!!!!!

  plhs(1) = mxCreateDoubleMatrix(MaxElNod,3,0)  
  LocVector_pr = mxGetPr(plhs(1))

  size = MaxElNod*3
  call mxCopyReal8ToPtr(LocVector,LocVector_pr,size)
  
  plhs(2) = mxCreateNumericArray(1,1,mxClassIDFromClassName('int32'),0)
  NumNodeElem_pr = mxGetPr(plhs(2))
  
  call mxCopyInteger8ToPtr(NumNodeElem,NumNodeElem_pr,1)
  
  
!!!!!!!!!!!!!!!!!!
!     Finish     !
!!!!!!!!!!!!!!!!!! 
  
  deallocate(GlobVector)   

  return
end


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!->Subroutine FEM_GLOB2LOC_EXTRACT.
!
!->Description.-
!
!    Extract the values corresponding to given nodes from a vector of global
!    data.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine fem_glob2loc_extract_wrap (ElemNodes,GlobVector,LocVector,NumNodeElem,NumNodes,MaxElNod)
  use lib_fem

! I/O Variables
  integer,intent(in)    :: MaxElNod                 !Maximum number of nodes per element (3)
  integer,intent(in)    :: NumNodes                 !Total number of nodes
  integer,intent(in)    :: ElemNodes(MaxElNod)      !Nodes in the current element
  real(8),intent(in)    :: GlobVector(NumNodes,3)   !Vector of global values   
  real(8),intent(inout) :: LocVector(MaxElNod,3)    !Vector of local values     
  integer,intent(inout) :: NumNodeElem              !Actual number of nodes in the element

  call fem_glob2loc_extract (ElemNodes,GlobVector,LocVector,NumNodeElem)

  return
 end subroutine fem_glob2loc_extract_wrap
