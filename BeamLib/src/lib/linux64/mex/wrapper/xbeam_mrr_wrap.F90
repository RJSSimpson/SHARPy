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
  integer	:: NumNodesElem_pr, r0_pr, Ri_pr, ElemMass_pr, MRR_pr, NumGauss_pr
  !I/O variables for computational routine  
  integer	:: NumNodesElem, NumGauss, MaxElNod
  real(8)	:: r0(3,6), Ri(3,6), ElemMass(6,6), MRR(6,6)
   
  
!!!!!!!!!!!!!!!!!!!!!!
!     Read input     !
!!!!!!!!!!!!!!!!!!!!!!

  !Check for proper number of arguments


  !Number of nodes in the element
  NumNodesElem_pr = mxGetPr(prhs(1))
  m = mxGetM(prhs(1))
  n = mxGetN(prhs(1))
  size = m*n
  call mxCopyPtrToInteger8(NumNodesElem_pr,NumNodesElem,size)

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
  ElemMass_pr = mxGetPr(prhs(4))
  m = mxGetM(prhs(4))
  n = mxGetN(prhs(4))
  size = m*n   
  call mxCopyPtrToReal8(ElemMass_pr,ElemMass,size)
  
  !Tangent mass matrix
  MRR_pr = mxGetPr(prhs(5))
  m = mxGetM(prhs(5))
  n = mxGetN(prhs(5))
  size = m*n
  call mxCopyPtrToReal8(MRR_pr,MRR,size)

  !Number of Gauss points
  NumGauss_pr = mxGetPr(prhs(6))
  m = mxGetM(prhs(6))
  n = mxGetN(prhs(6))
  size = m*n
  call mxCopyPtrToInteger8(NumGauss_pr,NumGauss,size)  
  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!     Call computational routine     !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !Create matrices for the return argument
  !Tangent stiffness matrix of element
  MaxElNod = 3
  m = 6
  n = 6
  size = m*n
  plhs(1) = mxCreateDoubleMatrix(m,n,0)  
  MRR_pr = mxGetPr(plhs(1))
  
  !Call the computational subroutine.
  call xbeam_mrr_wrap (NumNodesElem,r0,Ri,ElemMass,MRR,NumGauss)


!!!!!!!!!!!!!!!!!!!!!!!!
!     Write output     !
!!!!!!!!!!!!!!!!!!!!!!!!

  MaxElNod = 3
  m = 6
  n = 6
  size = m*n
  call mxCopyReal8ToPtr(MRR,MRR_pr,size)

!!!!!!!!!!!!!!!!!!
!     Finish     !
!!!!!!!!!!!!!!!!!! 

  return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-> Subroutine xbeam_mrr
!
!-> Description:
!
!    Compute the mass matrix corresponding to the reference system accelerations.
!
!-> Remarks.-
!
!   1) New values are added to those already in the matrix.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine xbeam_mrr_wrap (NumNodesElem,r0,Ri,ElemMass,MRR,NumGauss)
  use lib_xbeam
  
! I/O Variables.
  integer,intent(in)   :: NumNodesElem      ! Number of nodes in the element.
  real(8),intent(in)   :: r0       (3,6)    ! Initial position/orientation of grid points.
  real(8),intent(in)   :: Ri       (3,6)    ! Current position/orientation of grid points.
  real(8),intent(in)   :: ElemMass (6,6)    ! Inertial properties in the element.
  real(8),intent(inout):: MRR     (6,6)    ! Tangent mass matrix.
  integer,intent(in)   :: NumGauss          ! Number of Gauss points in the element.

  call xbeam_mrr (NumNodesElem,r0,Ri,ElemMass,MRR,NumGauss)

  return
 end subroutine xbeam_mrr_wrap
