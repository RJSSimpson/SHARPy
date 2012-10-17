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
  integer	:: NumNodesElem_pr, Ri_pr, Fi_pr, Flags_pr, Kmat_pr, FollowerForce_pr
  !I/O variables for computational routine  
  integer	:: NE, NumNodesElem, FF, i
  integer   :: Flags_aux(3)
  real(8)	:: Ri(3,6), Fi(3,6), Kmat(18,18)
  logical   :: FollowerForce
  logical   :: Flags(3)
   
  
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
  NumNodesElem_pr = mxGetPr(prhs(1))
  m = mxGetM(prhs(1))
  n = mxGetN(prhs(1))
  size = m*n
  call mxCopyPtrToInteger4(NumNodesElem_pr,NumNodesElem,size)

  !Current position/orientation of grid points
  Ri_pr = mxGetPr(prhs(2)) 
  m = mxGetM(prhs(2))
  n = mxGetN(prhs(2))
  size = m*n
  call mxCopyPtrToReal8(Ri_pr,Ri,size)
    
  !Current position/orientation of grid points
  Fi_pr = mxGetPr(prhs(3)) 
  m = mxGetM(prhs(3))
  n = mxGetN(prhs(3))
  size = m*n
  call mxCopyPtrToReal8(Fi_pr,Fi,size)
  
  !Identify master nodes
  Flags_pr = mxGetPr(prhs(4)) 
  m = mxGetM(prhs(4))
  n = mxGetN(prhs(4))
  size = m*n
  call mxCopyPtrToInteger4(Flags_pr,Flags_aux,size)
  
  !Material tangent stiffness matrix
  Kmat_pr = mxGetPr(prhs(5))
  m = mxGetM(prhs(5))
  n = mxGetN(prhs(5))
  size = m*n
  call mxCopyPtrToReal8(Kmat_pr,Kmat,size)

  !Follower force
  FollowerForce_pr = mxGetPr(prhs(6))
  call mxCopyPtrToInteger4(FollowerForce_pr,FF,1)
  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     Call computational routine     !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !Create matrices for the return argument
  !Tangent stiffness matrix of element
  MaxElNod = 3
  m = MaxElNod*6
  n = MaxElNod*6
  size = m*n
  plhs(1) = mxCreateDoubleMatrix(m,n,0)  
  Kmat_pr = mxGetPr(plhs(1))
  
  if (FF.eq.1) then
    FollowerForce = .true.
  else if (FF.eq.0) then
    FollowerForce = .false.
  end if
  
  do i=1,3
    if (Flags_aux(i).eq.1) then
      Flags(i) = .true.
    else if (Flags_aux(i).eq.0) then
      Flags(i) = .false.
    end if
  end do
  
  !Call the computational subroutine.
  call cbeam3_dqext_wrap (NumNodesElem,Ri,Fi,Flags,Kmat,FollowerForce)

!!!!!!!!!!!!!!!!!!!!!!!!
!     Write output     !
!!!!!!!!!!!!!!!!!!!!!!!!

  MaxElNod = 3
  m = MaxElNod*6
  n = MaxElNod*6
  size = m*n
  call mxCopyReal8ToPtr(Kmat,Kmat_pr,size)

   
!!!!!!!!!!!!!!!!!!
!     Finish     !
!!!!!!!!!!!!!!!!!! 
  
  return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-> Subroutine CBEAM3_DQEXT
!
!-> Description:
!
!    Matrix of derivatives of the influence coefficients for external forces.
!
!-> Remarks.-
!
!  1) Follower forces are assumed to be given in the local deformed frame (B frame).
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine cbeam3_dqext_wrap (NumNodesElem,Ri,Fi,Flags,Kmat,FollowerForce)
  use lib_cbeam3

! I/O Variables.
  integer,intent(in)   :: NumNodesElem    ! Number of nodes in the element.
  real(8),intent(in)   :: Ri    (3,6)     ! Current position/orientation of grid points.
  real(8),intent(in)   :: Fi    (3,6)     ! Current forces/moments at grid points.
  logical,intent(in)   :: Flags (3)       ! Identify master nodes.
  real(8),intent(inout):: Kmat  (18,18)   ! Derivative of the external forces.
  logical,intent(in)   :: FollowerForce   ! =T if follower forces.
 
  call cbeam3_dqext (NumNodesElem,Ri,Fi,Flags,Kmat,FollowerForce)

  return
 end subroutine cbeam3_dqext_wrap