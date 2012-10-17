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
  integer	:: NumGauss_pr, NumNodesElem_pr, Coords_pr, Arc_pr
  !I/O variables for computational routine
  integer   :: NumGauss, NumNodesElem
  real(8)   :: Arc
  real(8)   :: Coords(3,3)
  
  
!!!!!!!!!!!!!!!!!!!!!!
!     Read input     !
!!!!!!!!!!!!!!!!!!!!!!

  !Check for proper number of arguments. 
  if(nrhs .ne. 3) then
    call mexErrMsgTxt('Three inputs required.')
  elseif(nlhs .ne. 1) then
    call mexErrMsgTxt('One output required.')
  endif

  !Number of Gauss integration points
  NumGauss_pr = mxGetPr(prhs(1))
  m = mxGetM(prhs(1))
  n = mxGetN(prhs(1))
  size = m*n
  call mxCopyPtrToInteger4(NumGauss_pr,NumGauss,size)

  !Number of nodes in the element
  NumNodesElem_pr = mxGetPr(prhs(2)) 
  m = mxGetM(prhs(2))
  n = mxGetN(prhs(2))
  size = m*n
  call mxCopyPtrToInteger4(NumNodesElem_pr,NumNodesElem,size)
  
  !Coordinates of the grid points
  Coords_pr = mxGetPr(prhs(3)) 
  m = mxGetM(prhs(3))
  n = mxGetN(prhs(3))
  size = m*n
  call mxCopyPtrToReal8(Coords_pr,Coords,size)
   
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     Call computational routine     !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !Total arc length
  Arc = 0.d0
  
  !Call the computational subroutine.
  call bgeom_elemlength_wrap (NumGauss,NumNodesElem,Coords,Arc)


!!!!!!!!!!!!!!!!!!!!!!!!
!     Write output     !
!!!!!!!!!!!!!!!!!!!!!!!!

  m = 1
  n = 1
  size = m*n
  plhs(1) = mxCreateDoubleMatrix(m,n,0)  
  Arc_pr = mxGetPr(plhs(1))
  call mxCopyReal8ToPtr(Arc,Arc_pr,size)

   
!!!!!!!!!!!!!!!!!!
!     Finish     !
!!!!!!!!!!!!!!!!!! 
  
  return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-> Subroutine BGEOM_ELEMLENGTH
!
!-> Description:
!
!    Compute length of an element using Gaussian integration.
!
!-> Remarks.-
!
!   1) The element can have an arbitrary number of nodes. Isoparametric 
!      interpolation is used.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine bgeom_elemlength_wrap (NumGauss,NumNodesElem,Coords,Arclength)
  use lib_bgeom

! I/O Variables.
  integer,intent(in) :: NumGauss          ! Number of Gauss integration points.
  integer,intent(in) :: NumNodesElem      ! Number of nodes in the element.
  real(8),intent(in) :: Coords    (3,3)   ! Coordinates of the grid points.
  real(8),intent(out):: ArcLength         ! Element total arclength.

  call bgeom_elemlength (NumGauss,NumNodesElem,Coords,Arclength)

  return
 end subroutine bgeom_elemlength_wrap