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
  integer	:: NumNodesElem_pr, Ri_pr, Flags_pr, FollowerForce_pr, FollowerForceRig_pr, Cao_pr, Fmat_pr
  !I/O variables for computational routine  
  integer	:: NE, NumNodesElem, FF, FFRig, i
  integer   :: Flags_aux(3)
  real(8)	:: Ri(3,6), Cao(3,3), Fmat(6,18)
  logical   :: FollowerForce, FollowerForceRig
  logical   :: Flags(3)
   
  
!!!!!!!!!!!!!!!!!!!!!!
!     Read input     !
!!!!!!!!!!!!!!!!!!!!!!

  !Check for proper number of arguments
  if(nrhs .ne. 6) then
    call mexErrMsgTxt('Four inputs required.')
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
  
  !Identify master nodes
  Flags_pr = mxGetPr(prhs(3)) 
  m = mxGetM(prhs(3))
  n = mxGetN(prhs(3))
  size = m*n
  call mxCopyPtrToInteger4(Flags_pr,Flags_aux,size)

  !Follower force
  FollowerForce_pr = mxGetPr(prhs(4))
  call mxCopyPtrToInteger4(FollowerForce_pr,FF,1)
    
  !Follower force in a frame
  FollowerForceRig_pr = mxGetPr(prhs(5))
  call mxCopyPtrToInteger4(FollowerForceRig_pr,FFRig,1)
  
  !Rotation matrix Cao
  Cao_pr = mxGetPr(prhs(6)) 
  m = mxGetM(prhs(6))
  n = mxGetN(prhs(6))
  size = m*n
  call mxCopyPtrToReal8(Cao_pr,Cao,size)
  
  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     Call computational routine     !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !Create matrices for the return argument
  !Forces/moments on the element nodes
  MaxElNod = 3
  m = 6
  n = MaxElNod*6
  size = m*n
  plhs(1) = mxCreateDoubleMatrix(m,n,0)  
  Fmat_pr = mxGetPr(plhs(1))
  
  if (FF.eq.1) then
    FollowerForce = .true.
  else if (FF.eq.0) then
    FollowerForce = .false.
  end if
  if (FFRig.eq.1) then
    FollowerForceRig = .true.
  else if (FFRig.eq.0) then
    FollowerForceRig = .false.
  end if
  
  do i=1,3
    if (Flags_aux(i).eq.1) then
      Flags(i) = .true.
    else if (Flags_aux(i).eq.0) then
      Flags(i) = .false.
    end if
  end do
  
  !Call the computational subroutine.
  call xbeam_fext_wrap (NumNodesElem,Ri,Flags,Fmat,FollowerForce,FollowerForceRig,Cao)

!!!!!!!!!!!!!!!!!!!!!!!!
!     Write output     !
!!!!!!!!!!!!!!!!!!!!!!!!

  MaxElNod = 3
  m = 6
  n = MaxElNod*6
  size = m*n
  call mxCopyReal8ToPtr(Fmat,Fmat_pr,size)

   
!!!!!!!!!!!!!!!!!!
!     Finish     !
!!!!!!!!!!!!!!!!!! 
  
  return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-> Subroutine xbeam_FEXT
!
!-> Description:
!
!    Matrix of influence coefficients for external forces.
!
!-> Remarks.-
!
!  1) Follower forces are assumed to be given in the local deformed frame (B frame).
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine xbeam_fext_wrap (NumNodesElem,Ri,Flags,Fmat,FollowerForce,FollowerForceRig,Cao)
  use lib_xbeam

! I/O Variables.
  integer,intent(in)   :: NumNodesElem      ! Number of nodes in the element.
  real(8),intent(in)   :: Ri    (3,6)       ! Current position/orientation of grid points.
  logical,intent(in)   :: Flags (3)         ! Identify master nodes.
  real(8),intent(inout):: Fmat  (6,18)      ! Forces/moments on the element nodes.
  logical,intent(in)   :: FollowerForce     ! =T if follower forces.
  logical,intent(in)   :: FollowerForceRig  ! =T if follower force in the body-fixed frame.
  real(8),intent(in)   :: Cao   (3,3)       ! Rotation operator from inertial frame to reference frame
  
  Fmat=0.d0
  call xbeam_fext (NumNodesElem,Ri,Flags,Fmat,FollowerForce,FollowerForceRig,Cao)

  return
 end subroutine xbeam_fext_wrap