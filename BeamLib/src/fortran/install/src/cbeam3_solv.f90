!-> Copyright by Imperial College London, 2009
!
!-> Module.- CBEAM3_SOLV Rafa Palacios. 27Aug2009.
!
!-> Language: FORTRAN90, Free format.
!
!-> Description.-
!
!  Solve beam equations.
!
!-> Subroutines.-
!
!    -cbeam3_solv_nlnstatic     : Nonlinear static solution.
!    -cbeam3_solv_linstatic     : Linear static solution.
!    -cbeam3_solv_modal         : Find natural modes.
!    -cbeam3_solv_nlndyn        : Solve nonlinear dynamic problem.
!    -cbeam3_solv_lindyn        : Solve linear dynamic problem.
!    -cbeam3_solv_update_static : Update solution variables in static problem.
!    -cbeam3_solv_update_lindyn : Update solution variables in linear dynamic problem.
!    -cbeam3_solv_state2disp    : Write state vector for displacements/rotations.
!    -cbeam3_solv_disp2state    : Write state vector for displacements/rotations.
!
!-> Remarks.-
!
!-> Modifications.-
! 20120317 A.Da Ronch Conditional compilation added (NOLAPACK)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module cbeam3_solv
  use xbeam_shared
  implicit none

! Private variables.
  integer,private,parameter:: DimMat=24    ! Memory index for sparse matrices.

  real(8),private,parameter,dimension(4,4):: Unit4= &       ! 4x4 Unit matrix.
&         reshape((/1.d0,0.d0,0.d0,0.d0,0.d0,1.d0,0.d0,0.d0,0.d0,0.d0,1.d0,0.d0,0.d0,0.d0,0.d0,1.d0/),(/4,4/))

 contains


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-> Subroutine CBEAM3_SOLV_NLNSTATIC
!
!-> Description:
!
!    Steady-state solution of multibeam problem under applied forces.
!
!-> Remarks.-
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine cbeam3_solv_nlnstatic (NumDof,Elem,Node,AppForces,Coords,Psi0, &
&                                  PosDefor,PsiDefor,Options)
  use lib_fem
  use lib_sparse
#ifdef NOLAPACK
  use lib_lu
#else
  use interface_lapack
#endif
  use cbeam3_asbly

! I/O Variables.
  integer,      intent(in)   :: NumDof            ! Number of independent DoFs.
  type(xbelem),intent(in)    :: Elem(:)           ! Element information.
  type(xbnode),intent(in)    :: Node(:)           ! Nodal information.
  real(8),      intent(in)   :: AppForces (:,:)   ! Applied nodal forces.
  real(8),      intent(in)   :: Coords   (:,:)    ! Initial coordinates of the grid points.
  real(8),      intent(in)   :: Psi0     (:,:,:)  ! Initial CRV of the nodes in the elements.
  real(8),      intent(inout):: PosDefor (:,:)    ! Current coordinates of the grid points
  real(8),      intent(inout):: PsiDefor (:,:,:)  ! Current CRV of the nodes in the elements.
  type(xbopts),intent(in)    :: Options           ! Solver parameters.

! Local variables.
  real(8):: Delta                          ! Flag for convergence in iterations.
  integer:: fs                             ! Current storage size of force matrix.
  integer:: iLoadStep                      ! Counter in the load steps.
  integer:: Iter                           ! Counter on iterations.
  integer:: k                              ! Auxiliary integer.
  integer:: ks                             ! Current storage size of stiffness matrix.

  integer,allocatable::      ListIN (:)    ! List of independent nodes.
  real(8),allocatable::      Qglobal(:)    ! Global vector of discrete generalize forces.
  real(8),allocatable::      DeltaX (:)    ! Unknown in the linearized system.
  type(sparse),allocatable:: Kglobal(:)    ! Global stiffness matrix in sparse storage.
  type(sparse),allocatable:: Fglobal(:)    ! Influence coefficients matrix for applied forces.


! Initialize geometric constants and system state.
  allocate (ListIN (size(Node)))
  do k=1,size(Node)
    ListIN(k)=Node(k)%Vdof
  end do

! Allocate memory for solver (Use a conservative estimate of the size of the matrix Kglobal).
  allocate (Kglobal(DimMat*NumDof)); call sparse_zero (ks,Kglobal)
  allocate (Qglobal(NumDof));    Qglobal= 0.d0
  allocate (DeltaX (NumDof));    DeltaX = 0.d0
  allocate (Fglobal(DimMat*NumDof)); call sparse_zero (fs,Fglobal)

! Apply loads in substeps.
  do iLoadStep=1,Options%NumLoadSteps
    Iter  = 0
    Delta = Options%MinDelta+1.d0

! Iteration until convergence.
    do while (Delta.gt.Options%MinDelta)
      Iter= Iter+1
      if (Iter.gt.Options%MaxIterations) STOP 'Static equations did not converge (17235)'
      if (Options%PrintInfo) write (*,'(5X,A,I4,A,I4,$)') 'Load Step',iLoadStep,' Subiteration',Iter

! Assembly matrices and functional.
      Qglobal=0.d0
      call sparse_zero (ks,Kglobal)
      call sparse_zero (fs,Fglobal)

      call cbeam3_asbly_static (Elem,Node,Coords,Psi0,PosDefor,PsiDefor,AppForces*dble(iLoadStep)/dble(Options%NumLoadSteps), &
&                               ks,Kglobal,fs,Fglobal,Qglobal,Options)

! Add forces on the unconstrained nodes.
      Qglobal= Qglobal - dble(iLoadStep)/dble(Options%NumLoadSteps) * &
&              sparse_matvmul(fs,Fglobal,NumDof,fem_m2v(AppForces,NumDof,Filter=ListIN))

! Solve equation and update the global vectors.
!      call sparse_precond(ks,Kglobal,Qglobal,Options%MinDelta)

#ifdef NOLAPACK
      call lu_sparse(ks,Kglobal,-Qglobal,DeltaX)
#else
      call lapack_sparse (ks,Kglobal,-Qglobal,DeltaX)
#endif
      call cbeam3_solv_update_static (Elem,Node,Psi0,DeltaX,PosDefor,PsiDefor)

! Convergence parameter.
      Delta=maxval(abs(DeltaX))+maxval(abs(Qglobal))
      if (Options%PrintInfo) write (*,'(2X,A,1PE10.3,A,1PE10.3)') &
&                           'DeltaF=',maxval(abs(Qglobal)), ' DeltaX=',maxval(abs(DeltaX))
    end do
  end do

  call flush(6) !Flush stdout buffer so SharPy output looks nice

  deallocate (Kglobal,Qglobal,DeltaX)

  return
 end subroutine cbeam3_solv_nlnstatic


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-> Subroutine CBEAM3_SOLV_LINSTATIC
!
!-> Description:
!
!    Linear steady-state solution of multibeam problem under applied forces.
!
!-> Remarks.-
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine cbeam3_solv_linstatic (NumDof,Elem,Node,AppForces,Coords,Psi0, &
&                                  PosDefor,PsiDefor,Options)

  use lib_fem
  use lib_sparse
#ifdef NOLAPACK
  use lib_lu
#else
  use interface_lapack
#endif
  use cbeam3_asbly

! I/O Variables.
  integer,      intent(in)   :: NumDof            ! Number of independent DoFs.
  type(xbelem),intent(in)    :: Elem(:)           ! Element information.
  type(xbnode),intent(in)    :: Node(:)           ! Nodal information.
  real(8),      intent(in)   :: AppForces (:,:)   ! Applied nodal forces.
  real(8),      intent(in)   :: Coords   (:,:)    ! Initial coordinates of the grid points.
  real(8),      intent(in)   :: Psi0     (:,:,:)  ! Initial CRV of the nodes in the elements.
  real(8),      intent(inout):: PosDefor (:,:)    ! Current coordinates of the grid points
  real(8),      intent(inout):: PsiDefor (:,:,:)  ! Current CRV of the nodes in the elements.
  type(xbopts),intent(in)    :: Options           ! Solver parameters.

! Local variables.
  integer:: fs                             ! Current storage size of force matrix.
  integer:: k                              ! Auxiliary integer.
  integer:: ks                             ! Current storage size of stiffness matrix.
  integer:: NumN                           ! Number of nodes in the model.

  integer,allocatable::      ListIN (:)    ! List of independent nodes.
  real(8),allocatable::      DeltaX (:)    ! Unknown in the linearized system.
  real(8),allocatable::      Qglobal(:)    ! Global vector of discrete generalize forces.
  type(sparse),allocatable:: Fglobal(:)    ! Influence coefficients matrix for applied forces.
  type(sparse),allocatable:: Kglobal(:)    ! Global stiffness matrix in sparse storage.

! Initialize.
  NumN=size(Node)
  allocate (ListIN (NumN));
  do k=1,NumN
    ListIN(k)=Node(k)%Vdof
  end do

! Allocate memory for solver (Use a conservative estimate of the size of the matrix Kglobal).
  allocate (Kglobal(DimMat*NumDof)); call sparse_zero (ks,Kglobal)
  allocate (Fglobal(DimMat*NumDof)); call sparse_zero (fs,Fglobal)
  allocate (Qglobal(NumDof));        Qglobal= 0.d0
  allocate (DeltaX (NumDof));        DeltaX = 0.d0

! Assembly matrices and functional.
  Qglobal=0.d0
  call sparse_zero (ks,Kglobal)
  call sparse_zero (fs,Fglobal)

  call cbeam3_asbly_static (Elem,Node,Coords,Psi0,PosDefor,PsiDefor,AppForces, &
&                           ks,Kglobal,fs,Fglobal,Qglobal,Options)

! Forces on the unconstrained nodes.
  Qglobal= sparse_matvmul(fs,Fglobal,NumDof,fem_m2v(AppForces,NumDof,Filter=ListIN))

! Solve equation and update the global vectors.
#ifdef NOLAPACK
  call lu_sparse(ks,Kglobal,Qglobal,DeltaX)
#else
  call lapack_sparse (ks,Kglobal,Qglobal,DeltaX)
#endif
  call cbeam3_solv_update_static (Elem,Node,Psi0,DeltaX,PosDefor,PsiDefor)

  deallocate (Kglobal,Qglobal,DeltaX)
  return
 end subroutine cbeam3_solv_linstatic




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-> Subroutine CBEAM3_SOLV_MODAL
!
!-> Description:
!
!   Find linear vibration modes.
!
!-> Remarks.-
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine cbeam3_solv_modal  (iOut,NumDof,Elem,Node,Vrel,Coords,Psi0, &
&                               PosDefor,PsiDefor,Options)
  use lib_fem
  use lib_sparse
#ifndef NOLAPACK
  use interface_lapack
#endif
  use cbeam3_asbly

! I/O Variables.
  integer,      intent(in) :: iOut              ! Output file.
  integer,      intent(in) :: NumDof            ! Number of independent DoFs.
  type(xbelem), intent(in) :: Elem      (:)     ! Element information.
  type(xbnode), intent(in) :: Node      (:)     ! Nodal information.
  real(8),      intent(in) :: Vrel      (:,:)   ! Time history of the velocities of the reference frame.
  real(8),      intent(in) :: Coords    (:,:)   ! Initial coordinates of the grid points.
  real(8),      intent(in) :: Psi0      (:,:,:) ! Initial CRV of the nodes in the elements.
  real(8),      intent(in) :: PosDefor  (:,:)   ! Current coordinates of the grid points
  real(8),      intent(in) :: PsiDefor  (:,:,:) ! Current CRV of the nodes in the elements.
  type(xbopts), intent(in) :: Options           ! Solver parameters.

! Local variables.
  integer:: k                            ! Counters.
  integer:: cs,ks,ms
  type(sparse),allocatable:: Cglobal(:)     ! Sparse damping matrix.
  type(sparse),allocatable:: Kglobal(:)     ! Global stiffness matrix in sparse storage.
  type(sparse),allocatable:: Mglobal(:)     ! Global mass matrix in sparse storage.


! Allocate memory for solver (Use a conservative estimate of the size of the matrices).
  allocate (Mglobal(DimMat*NumDof)); call sparse_zero (ms,Mglobal)
  allocate (Cglobal(DimMat*NumDof)); call sparse_zero (cs,Cglobal)
  allocate (Kglobal(DimMat*NumDof)); call sparse_zero (ks,Kglobal)

! Compute tangent matrices at initial time.
  call cbeam3_asbly_modal (Elem,Node,Coords,Psi0,PosDefor,PsiDefor,Vrel(1,:), &
&                          ms,Mglobal,cs,Cglobal,ks,Kglobal,Options)

! Write matrices in files.
  open (unit=71,file='Msparse',status='replace')
  do k=1,ms
    write (71,'(2I4,1PE18.10)') Mglobal(k)%i,Mglobal(k)%j,Mglobal(k)%a
  end do
  close (71)

  open (unit=72,file='Csparse',status='replace')
  do k=1,cs
    write (72,'(2I4,1PE18.10)') Cglobal(k)%i,Cglobal(k)%j,Cglobal(k)%a
  end do
  close (72)

  open (unit=73,file='Ksparse',status='replace')
  do k=1,ks
    write (73,'(2I4,1PE18.10)') Kglobal(k)%i,Kglobal(k)%j,Kglobal(k)%a
  end do
  close (73)

! End of routine.
  deallocate (Mglobal,Kglobal,Cglobal)
  return
 end subroutine cbeam3_solv_modal



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-> Subroutine CBEAM3_SOLV_NLNDYN
!
!-> Description:
!
!    Linear dynamic solution of multibeam problem under applied forces.
!
!-> Remarks.-
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine cbeam3_solv_nlndyn (iOut,NumDof,Time,Elem,Node,F0,Fa,Ftime,              &
&                               Vrel,VrelDot,Coords,Psi0,PosDefor,PsiDefor,          &
&                               PosDotDefor,PsiDotDefor,PosPsiTime,VelocTime,DynOut, &
&                               OutGrids,Options)
  use lib_fem
  use lib_rot
  use lib_rotvect
  use lib_lu
  use lib_out
  use lib_sparse
#ifdef NOLAPACK
  use lib_lu
#else
  use interface_lapack
#endif
  use cbeam3_asbly
  use lib_xbeam

! I/O Variables.
  integer,      intent(in)   :: iOut              ! Output file.
  integer,      intent(in)   :: NumDof            ! Number of independent DoFs.
  real(8),      intent(in)   :: Time      (:)     ! Time steps.
  type(xbelem), intent(in)   :: Elem      (:)     ! Element information.
  type(xbnode), intent(in)   :: Node      (:)     ! Nodal information.
  real(8),      intent(in)   :: F0        (:,:)   ! Applied static nodal forces.
  real(8),      intent(in)   :: Fa        (:,:)   ! Amplitude of the dynamic nodal forces.
  real(8),      intent(in)   :: Ftime     (:)     ! Time history of the applied forces.
  real(8),      intent(in)   :: Vrel      (:,:)   ! Time history of the velocities of the reference frame.
  real(8),      intent(in)   :: VrelDot   (:,:)   ! Time history of the accelerations of the reference frame.
  real(8),      intent(in)   :: Coords    (:,:)   ! Initial coordinates of the grid points.
  real(8),      intent(in)   :: Psi0      (:,:,:) ! Initial CRV of the nodes in the elements.
  real(8),      intent(inout):: PosDefor  (:,:)   ! Current coordinates of the grid points
  real(8),      intent(inout):: PsiDefor  (:,:,:) ! Current CRV of the nodes in the elements.
  real(8),      intent(inout):: PosDotDefor(:,:)  ! Current time derivatives of the coordinates of the grid points
  real(8),      intent(inout):: PsiDotDefor(:,:,:)! Current time derivatives of the CRV of the nodes in the elements.
  real(8),      intent(out)  :: PosPsiTime(:,:)   ! Time-history of the position/rotation at selected nodes.
  real(8),      intent(out)  :: VelocTime (:,:)   ! Time-history of the time derivatives at selected nodes.
  real(8),      intent(out)  :: DynOut    (:,:)   ! Time-history of displacement of all nodes.
  logical,      intent(inout):: OutGrids  (:)     ! Output grids.
  type(xbopts),intent(in)    :: Options           ! Solver parameters.

! Local variables.
  real(8):: beta,gamma                     ! Newmark coefficients.
  real(8):: dt                             ! Time step
  integer:: k                              ! Counters.
  integer:: iStep,Iter                     ! Counters on time steps and subiterations.
  real(8):: MinDelta                       ! Value of Delta for convergence.
  integer:: NumE(1)                        ! Number of elements in the model.
  integer:: NumN                           ! Number of nodes in the model.

  real(8),allocatable:: dXdt(:),dXddt(:)  ! Generalized coordinates and derivatives.
  real(8),allocatable:: X(:), DX(:)

  integer,allocatable::  ListIN     (:)    ! List of independent nodes.

  ! Define variables for system matrices.
  integer:: as,cs,ks,ms,fs
  type(sparse),allocatable:: Asys   (:)    ! System matrix for implicit Newmark method.
  type(sparse),allocatable:: Cglobal(:)    ! Sparse damping matrix.
  type(sparse),allocatable:: Kglobal(:)    ! Global stiffness matrix in sparse storage.
  type(sparse),allocatable:: Mglobal(:)    ! Global mass matrix in sparse storage.
  type(sparse),allocatable:: Fglobal(:)
  real(8),allocatable::      Qglobal(:)    ! Global vector of discrete generalize forces.
  real(8),allocatable::      Mvel(:,:)     ! Mass and damping from the motions of reference system.
  real(8),allocatable::      Cvel(:,:)     ! Mass and damping from the motions of reference system.

  ! Define vaiables for output information.
  character(len=80)  ::  Text          ! Text with current time information to print out.
  type(outopts)      ::  OutOptions    ! Output options.
  real(8),allocatable::  Displ(:,:)    ! Current nodal displacements/rotations.
  real(8),allocatable::  Veloc(:,:)    ! Current nodal velocities.

  ! Rotation operator for body-fixed frame using quaternions
  real(8) :: Cao(3,3)
  real(8) :: Quat(4)
  real(8) :: Temp(4,4)

! Initialize.
  NumN=size(Node)
  NumE(1)=size(Elem)
  allocate (ListIN (NumN));
  do k=1,NumN
    ListIN(k)=Node(k)%Vdof
  end do
  gamma=1.d0/2.d0+Options%NewmarkDamp
  beta =1.d0/4.d0*(gamma+0.5d0)*(gamma+0.5d0)

! Allocate memory for solver (Use a conservative estimate of the size of the matrices).
  allocate (Asys   (DimMat*NumDof)); call sparse_zero (as,Asys)
  allocate (Mglobal(DimMat*NumDof)); call sparse_zero (ms,Mglobal)
  allocate (Cglobal(DimMat*NumDof)); call sparse_zero (cs,Cglobal)
  allocate (Kglobal(DimMat*NumDof)); call sparse_zero (ks,Kglobal)
  allocate (Fglobal(DimMat*NumDof)); call sparse_zero (fs,Fglobal)
  allocate (Qglobal(NumDof));   Qglobal= 0.d0
  allocate (Mvel   (NumDof,6)); Mvel   = 0.d0
  allocate (Cvel   (NumDof,6)); Cvel   = 0.d0

  allocate (X     (NumDof)); X      = 0.d0
  allocate (DX    (NumDof)); DX     = 0.d0
  allocate (dXdt  (NumDof)); dXdt   = 0.d0
  allocate (dXddt (NumDof)); dXddt  = 0.d0

! Compute system information at initial condition.
  allocate (Veloc(NumN,6)); Veloc=0.d0
  allocate (Displ(NumN,6)); Displ=0.d0

! Allocate quaternions and rotation operator for initially undeformed system
  Quat = (/1.d0,0.d0,0.d0,0.d0/); Cao = Unit; Temp = Unit4

! Extract initial displacements and velocities.
  call cbeam3_solv_disp2state (Node,PosDefor,PsiDefor,PosDotDefor,PsiDotDefor,X,dXdt)

! Compute initial acceleration (we are neglecting qdotdot in Kmass).
  call cbeam3_asbly_dynamic (Elem,Node,Coords,Psi0,PosDefor,PsiDefor,PosDotDefor,PsiDotDefor,           &
&                            0.d0*PosDefor,0.d0*PsiDefor,F0+Ftime(1)*Fa,Vrel(1,:),VrelDot(1,:),         &
&                            ms,Mglobal,Mvel,cs,Cglobal,Cvel,ks,Kglobal,fs,Fglobal,Qglobal,Options,Cao)

  Qglobal= Qglobal - sparse_matvmul(fs,Fglobal,NumDof,fem_m2v(F0+Ftime(1)*Fa,NumDof,Filter=ListIN))

  call sparse_addsparse(0,0,ms,Mglobal,as,Asys)
#ifdef NOLAPACK
  call lu_sparse(as,Asys,-Qglobal,dXddt)
#else
  call lapack_sparse (as,Asys,-Qglobal,dXddt)
#endif

! Loop in the time steps.
  do iStep=1,size(Time)-1
    dt= Time(iStep+1)-Time(iStep)
    if (Options%PrintInfo) then
      call out_time(iStep,Time(iStep+1),Text)
      write (*,'(5X,A,$)') trim(Text)
    end if
! Update transformation matrix for given angular velocity
    call lu_invers ((Unit4+0.25d0*xbeam_QuadSkew(Vrel(iStep+1,4:6))*(Time(iStep+1)-Time(iStep))),Temp)
    Quat=matmul(Temp,matmul((Unit4-0.25d0*xbeam_QuadSkew(Vrel(iStep,4:6))*(Time(iStep+1)-Time(iStep))),Quat))
    Cao = xbeam_Rot(Quat)

! Predictor step.
    X    = X + dt*dXdt + (0.5d0-beta)*dt*dt*dXddt
    dXdt = dXdt + (1.d0-gamma)*dt*dXddt
    dXddt= 0.d0

! Iteration until convergence.
    do Iter=1,Options%MaxIterations+1
      if (Iter.gt.Options%MaxIterations) STOP 'Solution did not converge (18235)'

! Update nodal positions and velocities .
      call cbeam3_solv_state2disp (Elem,Node,Coords,Psi0,X,dXdt,PosDefor,PsiDefor,PosDotDefor,PsiDotDefor)

! Compute system functionals and matrices. (Use initial accelerations for Kgyr).
      Qglobal= 0.d0
      Mvel   = 0.d0
      Cvel   = 0.d0
      call sparse_zero (ms,Mglobal)
      call sparse_zero (cs,Cglobal)
      call sparse_zero (ks,Kglobal)
      call sparse_zero (fs,Fglobal)

      call cbeam3_asbly_dynamic (Elem,Node,Coords,Psi0,PosDefor,PsiDefor,PosDotDefor,PsiDotDefor,                  &
&                                0.d0*PosDefor,0.d0*PsiDefor,F0+Ftime(iStep+1)*Fa,Vrel(iStep+1,:),VrelDot(iStep+1,:),    &
&                                ms,Mglobal,Mvel,cs,Cglobal,Cvel,ks,Kglobal,fs,Fglobal,Qglobal,Options,Cao)

! Compute admissible error.
      MinDelta=Options%MinDelta*max(1.d0,maxval(abs(Qglobal)))

! Compute the residual.
      Qglobal= Qglobal + sparse_matvmul(ms,Mglobal,NumDof,dXddt) + matmul(Mvel,Vreldot(iStep+1,:))
      Qglobal= Qglobal - sparse_matvmul(fs,Fglobal,NumDof,fem_m2v(F0+Ftime(iStep+1)*Fa,NumDof,Filter=ListIN))

! Check convergence.
      if (maxval(abs(DX)+abs(Qglobal)).lt.MinDelta) then
        if (Options%PrintInfo) then
          write (*,'(5X,A,I4,A,1PE12.3)') 'Subiteration',Iter, '  Delta=', maxval(abs(Qglobal))
        end if
        exit
      end if

! Compute Jacobian
      call sparse_zero (as,Asys)
      call sparse_addsparse(0,0,ks,Kglobal,as,Asys,Factor=1.d0)
      call sparse_addsparse(0,0,cs,Cglobal,as,Asys,Factor=gamma/(beta*dt))
      call sparse_addsparse(0,0,ms,Mglobal,as,Asys,Factor=1.d0/(beta*dt*dt))

! Calculation of the correction.
#ifdef NOLAPACK
      call lu_sparse(as,Asys,-Qglobal,DX)
#else
      call lapack_sparse (as,Asys,-Qglobal,DX)
#endif
      X    = X     + DX
      dXdt = dXdt  + gamma/(beta*dt)*DX
      dXddt= dXddt + 1.d0/(beta*dt*dt)*DX
    end do

! Update nodal positions and velocities on the current converged time step.
    call cbeam3_solv_state2disp (Elem,Node,Coords,Psi0,X,dXdt,PosDefor,PsiDefor,PosDotDefor,PsiDotDefor)


!!! Postprocesing (for single cantilever beams) !!!
! Store data in output variables (V_B, Omega_B).
!    if (any(OutGrids)) then
!      Veloc=0.d0
!      Displ=0.d0
!
!      do k=1,NumN
!        Displ(k,1:3)= PosDefor(k,1:3)-Coords(k,1:3)
!      end do
!
!      Veloc(1,1:3)= PosDotDefor(1,:)+rot_cross(Vrel(iStep,4:6), PosDefor(1,:))+Vrel(iStep,1:3)
!      do k=2,NumN
!        Displ(k,4:6)=PsiDefor(k-1,2,:)
!        CBa=rotvect_psi2mat(PsiDefor(k-1,2,:))
!        Veloc(k,1:3)= matmul(CBa, PosDotDefor(k,:) + Vrel(iStep,1:3)      &
!&                               + rot_cross(Vrel(iStep,4:6),PosDefor(k,:)))
!        Veloc(k,4:6)= matmul(rotvect_psi2rot(PsiDefor(k-1,2,:)),PsiDotDefor(k-1,2,:)) &
!&                   + matmul(CBa,Vrel(iStep,4:6))
!      end do
!
!!  Write output information in output file.
!      OutOptions%PrintDispl=.true.
!      OutOptions%PrintVeloc=.true.
!      call out_title   (iOut,trim(Text))
!      call out_outgrid (iOut,'NODE',OutOptions,1,NumE,6,OutGrids,DISPL=Displ,VELOC=Veloc)
!
!! Write output to export to main program (obsolete!).
!      VelocTime (iStep+1,:)= Veloc (1:NumN,3)
!    end if

! Write output to export to main program (obsolete!).
    PosPsiTime(iStep+1,1:3)= PosDefor(NumN,:)
    PosPsiTime(iStep+1,4:6)= PsiDefor(NumE(1),Elem(NumE(1))%NumNodes,:)

    do k=1,NumN
        DynOut(iStep*NumN+k,:) = PosDefor(k,:)
    end do

  end do

  deallocate (ListIN,Mvel,Cvel)
  deallocate (Asys,Fglobal,Mglobal)
  deallocate (Kglobal,Cglobal,Qglobal)
  deallocate (X,DX,dXdt,dXddt)
  deallocate (Displ,Veloc)
  return
 end subroutine cbeam3_solv_nlndyn


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-> Subroutine CBEAM3_SOLV_NLNDYN_ACCEL
!
!-> Description:
!
!    Nonlinear dynamic solution of multibeam problem under applied forces
!    with accelerations added as input/output variable.
!
!-> Remarks: Initial accelerations should be calculated outside this routine.
!-> Author: Rob Simpson
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine cbeam3_solv_nlndyn_accel (iOut,NumDof,Time,Elem,Node,F0,Fa,Ftime,        &
&                               Vrel,VrelDot,Coords,Psi0,PosDefor,PsiDefor,          &
&                               PosDotDefor,PsiDotDefor,PosDDot,PsiDDot,             &
&                               PosPsiTime,VelocTime,DynOut, &
&                               OutGrids,Options)
  use lib_fem
  use lib_rot
  use lib_rotvect
  use lib_lu
  use lib_out
  use lib_sparse
#ifdef NOLAPACK
  use lib_lu
#else
  use interface_lapack
#endif
  use cbeam3_asbly
  use lib_xbeam

! I/O Variables.
  integer,      intent(in)   :: iOut              ! Output file.
  integer,      intent(in)   :: NumDof            ! Number of independent DoFs.
  real(8),      intent(in)   :: Time      (:)     ! Time steps.
  type(xbelem), intent(in)   :: Elem      (:)     ! Element information.
  type(xbnode), intent(in)   :: Node      (:)     ! Nodal information.
  real(8),      intent(in)   :: F0        (:,:)   ! Applied static nodal forces.
  real(8),      intent(in)   :: Fa        (:,:)   ! Amplitude of the dynamic nodal forces.
  real(8),      intent(in)   :: Ftime     (:)     ! Time history of the applied forces.
  real(8),      intent(in)   :: Vrel      (:,:)   ! Time history of the velocities of the reference frame.
  real(8),      intent(in)   :: VrelDot   (:,:)   ! Time history of the accelerations of the reference frame.
  real(8),      intent(in)   :: Coords    (:,:)   ! Initial coordinates of the grid points.
  real(8),      intent(in)   :: Psi0      (:,:,:) ! Initial CRV of the nodes in the elements.
  real(8),      intent(inout):: PosDefor  (:,:)   ! Current coordinates of the grid points
  real(8),      intent(inout):: PsiDefor  (:,:,:) ! Current CRV of the nodes in the elements.
  real(8),      intent(inout):: PosDotDefor(:,:)  ! Current time derivatives of the coordinates of the grid points
  real(8),      intent(inout):: PsiDotDefor(:,:,:)! Current time derivatives of the CRV of the nodes in the elements.
  real(8),      intent(inout):: PosDDot(:,:)  ! Current accelerations of the coordinates of the grid points
  real(8),      intent(inout):: PsiDDot(:,:,:)! Current accelerations of the CRV of the nodes in the elements.
  real(8),      intent(out)  :: PosPsiTime(:,:)   ! Time-history of the position/rotation at selected nodes.
  real(8),      intent(out)  :: VelocTime (:,:)   ! Time-history of the time derivatives at selected nodes.
  real(8),      intent(out)  :: DynOut    (:,:)   ! Time-history of displacement of all nodes.
  logical,      intent(inout):: OutGrids  (:)     ! Output grids.
  type(xbopts),intent(in)    :: Options           ! Solver parameters.

! Local variables.
  real(8):: beta,gamma                     ! Newmark coefficients.
  real(8):: dt                             ! Time step
  integer:: k                              ! Counters.
  integer:: iStep,Iter                     ! Counters on time steps and subiterations.
  real(8):: MinDelta                       ! Value of Delta for convergence.
  integer:: NumE(1)                        ! Number of elements in the model.
  integer:: NumN                           ! Number of nodes in the model.

  real(8),allocatable:: dXdt(:),dXddt(:)  ! Generalized coordinates and derivatives.
  real(8),allocatable:: X(:), DX(:)

  integer,allocatable::  ListIN     (:)    ! List of independent nodes.

  ! Define variables for system matrices.
  integer:: as,cs,ks,ms,fs
  type(sparse),allocatable:: Asys   (:)    ! System matrix for implicit Newmark method.
  type(sparse),allocatable:: Cglobal(:)    ! Sparse damping matrix.
  type(sparse),allocatable:: Kglobal(:)    ! Global stiffness matrix in sparse storage.
  type(sparse),allocatable:: Mglobal(:)    ! Global mass matrix in sparse storage.
  type(sparse),allocatable:: Fglobal(:)
  real(8),allocatable::      Qglobal(:)    ! Global vector of discrete generalize forces.
  real(8),allocatable::      Mvel(:,:)     ! Mass and damping from the motions of reference system.
  real(8),allocatable::      Cvel(:,:)     ! Mass and damping from the motions of reference system.

  ! Define vaiables for output information.
  character(len=80)  ::  Text          ! Text with current time information to print out.
  type(outopts)      ::  OutOptions    ! Output options.
  real(8),allocatable::  Displ(:,:)    ! Current nodal displacements/rotations.
  real(8),allocatable::  Veloc(:,:)    ! Current nodal velocities.

  ! Rotation operator for body-fixed frame using quaternions
  real(8) :: Cao(3,3)
  real(8) :: Quat(4)
  real(8) :: Temp(4,4)

! Initialize.
  NumN=size(Node)
  NumE(1)=size(Elem)
  allocate (ListIN (NumN));
  do k=1,NumN
    ListIN(k)=Node(k)%Vdof
  end do
  gamma=1.d0/2.d0+Options%NewmarkDamp
  beta =1.d0/4.d0*(gamma+0.5d0)*(gamma+0.5d0)

! Allocate memory for solver (Use a conservative estimate of the size of the matrices).
  allocate (Asys   (DimMat*NumDof)); call sparse_zero (as,Asys)
  allocate (Mglobal(DimMat*NumDof)); call sparse_zero (ms,Mglobal)
  allocate (Cglobal(DimMat*NumDof)); call sparse_zero (cs,Cglobal)
  allocate (Kglobal(DimMat*NumDof)); call sparse_zero (ks,Kglobal)
  allocate (Fglobal(DimMat*NumDof)); call sparse_zero (fs,Fglobal)
  allocate (Qglobal(NumDof));   Qglobal= 0.d0
  allocate (Mvel   (NumDof,6)); Mvel   = 0.d0
  allocate (Cvel   (NumDof,6)); Cvel   = 0.d0

  allocate (X     (NumDof)); X      = 0.d0
  allocate (DX    (NumDof)); DX     = 0.d0
  allocate (dXdt  (NumDof)); dXdt   = 0.d0
  allocate (dXddt (NumDof)); dXddt  = 0.d0

! Compute system information at initial condition.
  allocate (Veloc(NumN,6)); Veloc=0.d0
  allocate (Displ(NumN,6)); Displ=0.d0

! Allocate quaternions and rotation operator for initially undeformed system
  Quat = (/1.d0,0.d0,0.d0,0.d0/); Cao = Unit; Temp = Unit4

! Extract initial displacements, velocities AND accelerations.
  call cbeam3_solv_disp2state (Node,PosDefor,PsiDefor,PosDotDefor,PsiDotDefor,X,dXdt)
  call cbeam3_solv_accel2state (Node,PosDDot,PsiDDot,dXddt)

! Loop in the time steps.
  do iStep=1,size(Time)-1
    dt= Time(iStep+1)-Time(iStep)
    if (Options%PrintInfo) then
      call out_time(iStep,Time(iStep+1),Text)
      write (*,'(5X,A,$)') trim(Text)
    end if
! Update transformation matrix for given angular velocity
    call lu_invers ((Unit4+0.25d0*xbeam_QuadSkew(Vrel(iStep+1,4:6))*(Time(iStep+1)-Time(iStep))),Temp)
    Quat=matmul(Temp,matmul((Unit4-0.25d0*xbeam_QuadSkew(Vrel(iStep,4:6))*(Time(iStep+1)-Time(iStep))),Quat))
    Cao = xbeam_Rot(Quat)

! Predictor step.
    X    = X + dt*dXdt + (0.5d0-beta)*dt*dt*dXddt
    dXdt = dXdt + (1.d0-gamma)*dt*dXddt
    dXddt= 0.d0

! Iteration until convergence.
    do Iter=1,Options%MaxIterations+1
      if (Iter.gt.Options%MaxIterations) STOP 'Solution did not converge (18235)'

! Update nodal positions and velocities .
      call cbeam3_solv_state2disp (Elem,Node,Coords,Psi0,X,dXdt,PosDefor,PsiDefor,PosDotDefor,PsiDotDefor)

! Compute system functionals and matrices. (Use initial accelerations for Kgyr).
      Qglobal= 0.d0
      Mvel   = 0.d0
      Cvel   = 0.d0
      call sparse_zero (ms,Mglobal)
      call sparse_zero (cs,Cglobal)
      call sparse_zero (ks,Kglobal)
      call sparse_zero (fs,Fglobal)

      call cbeam3_asbly_dynamic (Elem,Node,Coords,Psi0,PosDefor,PsiDefor,PosDotDefor,PsiDotDefor,                  &
&                                0.d0*PosDefor,0.d0*PsiDefor,F0+Ftime(iStep+1)*Fa,Vrel(iStep+1,:),VrelDot(iStep+1,:),    &
&                                ms,Mglobal,Mvel,cs,Cglobal,Cvel,ks,Kglobal,fs,Fglobal,Qglobal,Options,Cao)

! Compute admissible error.
      MinDelta=Options%MinDelta*max(1.d0,maxval(abs(Qglobal)))

! Compute the residual.
      Qglobal= Qglobal + sparse_matvmul(ms,Mglobal,NumDof,dXddt) + matmul(Mvel,Vreldot(iStep+1,:))
      Qglobal= Qglobal - sparse_matvmul(fs,Fglobal,NumDof,fem_m2v(F0+Ftime(iStep+1)*Fa,NumDof,Filter=ListIN))

! Check convergence.
      if (maxval(abs(DX)+abs(Qglobal)).lt.MinDelta) then
        if (Options%PrintInfo) then
          write (*,'(5X,A,I4,A,1PE12.3)') 'Subiteration',Iter, '  Delta=', maxval(abs(Qglobal))
        end if
        exit
      end if

! Compute Jacobian
      call sparse_zero (as,Asys)
      call sparse_addsparse(0,0,ks,Kglobal,as,Asys,Factor=1.d0)
      call sparse_addsparse(0,0,cs,Cglobal,as,Asys,Factor=gamma/(beta*dt))
      call sparse_addsparse(0,0,ms,Mglobal,as,Asys,Factor=1.d0/(beta*dt*dt))

! Calculation of the correction.
#ifdef NOLAPACK
      call lu_sparse(as,Asys,-Qglobal,DX)
#else
      call lapack_sparse (as,Asys,-Qglobal,DX)
#endif
      X    = X     + DX
      dXdt = dXdt  + gamma/(beta*dt)*DX
      dXddt= dXddt + 1.d0/(beta*dt*dt)*DX
    end do

! Update nodal positions and velocities on the current converged time step.
    call cbeam3_solv_state2disp (Elem,Node,Coords,Psi0,X,dXdt,PosDefor,PsiDefor,PosDotDefor,PsiDotDefor)
    call cbeam3_solv_state2accel(Elem,Node,dXddt,PosDDot,PsiDDot)

!!! Postprocesing (for single cantilever beams) !!!
! Store data in output variables (V_B, Omega_B).
!    if (any(OutGrids)) then
!      Veloc=0.d0
!      Displ=0.d0
!
!      do k=1,NumN
!        Displ(k,1:3)= PosDefor(k,1:3)-Coords(k,1:3)
!      end do
!
!      Veloc(1,1:3)= PosDotDefor(1,:)+rot_cross(Vrel(iStep,4:6), PosDefor(1,:))+Vrel(iStep,1:3)
!      do k=2,NumN
!        Displ(k,4:6)=PsiDefor(k-1,2,:)
!        CBa=rotvect_psi2mat(PsiDefor(k-1,2,:))
!        Veloc(k,1:3)= matmul(CBa, PosDotDefor(k,:) + Vrel(iStep,1:3)      &
!&                               + rot_cross(Vrel(iStep,4:6),PosDefor(k,:)))
!        Veloc(k,4:6)= matmul(rotvect_psi2rot(PsiDefor(k-1,2,:)),PsiDotDefor(k-1,2,:)) &
!&                   + matmul(CBa,Vrel(iStep,4:6))
!      end do
!
!!  Write output information in output file.
!      OutOptions%PrintDispl=.true.
!      OutOptions%PrintVeloc=.true.
!      call out_title   (iOut,trim(Text))
!      call out_outgrid (iOut,'NODE',OutOptions,1,NumE,6,OutGrids,DISPL=Displ,VELOC=Veloc)
!
!! Write output to export to main program (obsolete!).
!      VelocTime (iStep+1,:)= Veloc (1:NumN,3)
!    end if

! Write output to export to main program (obsolete!).
    PosPsiTime(iStep+1,1:3)= PosDefor(NumN,:)
    PosPsiTime(iStep+1,4:6)= PsiDefor(NumE(1),Elem(NumE(1))%NumNodes,:)

    do k=1,NumN
        DynOut(iStep*NumN+k,:) = PosDefor(k,:)
    end do

  end do

  deallocate (ListIN,Mvel,Cvel)
  deallocate (Asys,Fglobal,Mglobal)
  deallocate (Kglobal,Cglobal,Qglobal)
  deallocate (X,DX,dXdt,dXddt)
  deallocate (Displ,Veloc)
  return
 end subroutine cbeam3_solv_nlndyn_accel

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-> Subroutine CBEAM3_SOLV_LINDYN
!
!-> Description:
!
!    Linear dynamic solution of multibeam problem for given applied forces.
!
!-> Remarks.-
!
!   1) This routine assumes that the motion of the reference system is
!      prescribed.
!
!   2) Coords/Psi0 defines the unloaded geometry, while PosDefor/PsiDefor
!      brings the initial geometry (after static equilibrium). F0 should
!      include the static loading for that equilibrium.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine cbeam3_solv_lindyn (iOut,NumDof,Time,Elem,Node,F0,Fa,Ftime,              &
&                               Vrel,VrelDot,Coords,Psi0,PosDefor,PsiDefor,          &
&                               PosDotDefor,PsiDotDefor,PosPsiTime,VelocTime,DynOut, &
&                               OutGrids,Options)
  use lib_rot
  use lib_fem
  use lib_sparse
  use lib_out
  use lib_lu
#ifndef NOLAPACK
  use interface_lapack
#endif
  use cbeam3_asbly
  use lib_xbeam

! I/O Variables.
  integer,      intent(in)   :: iOut              ! Output file.
  integer,      intent(in)   :: NumDof            ! Number of independent DoFs.
  real(8),      intent(in)   :: Time      (:)     ! Time steps.
  type(xbelem),intent(in)    :: Elem      (:)     ! Element information.
  type(xbnode),intent(in)    :: Node      (:)     ! Nodal information.
  real(8),      intent(in)   :: F0        (:,:)   ! Applied static nodal forces.
  real(8),      intent(in)   :: Fa        (:,:)   ! Amplitude of the dynamic nodal forces.
  real(8),      intent(in)   :: Ftime     (:)     ! Time history of the applied forces.
  real(8),      intent(in)   :: Vrel      (:,:)   ! Time history of the velocities of the reference frame.
  real(8),      intent(in)   :: VrelDot   (:,:)   ! Time history of the accelerations of the reference frame.
  real(8),      intent(in)   :: Coords    (:,:)   ! Undeformed coordinates of the grid points.
  real(8),      intent(in)   :: Psi0      (:,:,:) ! Undeformed CRV of the nodes in the elements.
  real(8),      intent(inout):: PosDefor  (:,:)   ! Initial/final position vector of grid points
  real(8),      intent(inout):: PsiDefor  (:,:,:) ! Initial/final CRV of the nodes in the elements.
  real(8),      intent(inout):: PosDotDefor(:,:)  ! Current time derivatives of the coordinates of the grid points
  real(8),      intent(inout):: PsiDotDefor(:,:,:)! Current time derivatives of the CRV of the nodes in the elements.
  real(8),      intent(out)  :: PosPsiTime(:,:)   ! Time-history of the position/rotation at selected nodes.
  real(8),      intent(out)  :: VelocTime (:,:)   ! Time-history of the time derivatives at selected nodes.
  real(8),      intent(out)  :: DynOut    (:,:)   ! Time-history of displacement of all nodes.
  logical,      intent(inout):: OutGrids  (:)     ! Output grids.
  type(xbopts),intent(in)    :: Options           ! Solver parameters.

! Local variables.
  real(8):: beta,gamma
  real(8):: dt                             ! Time step
  integer:: k                              ! Counters.
  integer:: iStep                          ! Current time step.
  integer:: NumE(1)                        ! Number of elements in the model.
  integer:: NumN                           ! Number of nodes in the model.

  real(8),allocatable:: X0(:),DX(:),DXDt(:),DXDDt(:)
  integer,allocatable:: ListIN     (:)     ! List of independent nodes.
  real(8),allocatable:: DeltaPos   (:,:)   ! Initial/final position vector of grid points
  real(8),allocatable:: DeltaPsi   (:,:,:) ! Initial/final CRV of the nodes in the elements.
  real(8),allocatable:: DeltaPosDot(:,:)   ! Current time derivatives of the coordinates of the grid points
  real(8),allocatable:: DeltaPsiDot(:,:,:) ! Current time derivatives of the CRV of the nodes in the elements.

  integer:: as,cs,ks,ms,fs
  type(sparse),allocatable  :: Asys   (:)     ! System matrix for implicit Newmark method.
  type(sparse),allocatable  :: Cglobal(:)     ! Sparse damping matrix.
  type(sparse),allocatable  :: Kglobal(:)     ! Global stiffness matrix in sparse storage.
  type(sparse),allocatable  :: Mglobal(:)     ! Global mass matrix in sparse storage.
  type(sparse),allocatable  :: Fglobal(:)
  real(8),allocatable       :: Qglobal(:)     ! Global vector of discrete generalize forces.
  real(8),allocatable       :: Mvel(:,:)      ! Mass and damping from the motions of reference system.
  real(8),allocatable       :: Cvel(:,:)      ! Mass and damping from the motions of reference system.

  real(8) :: Cao(3,3)           ! Rotation operator for body-fixed frame using quaternions
  real(8) :: Quat(4)
  real(8) :: Temp(4,4)

  character(len=80)  :: Text          ! Text with current time information to print out.
  type(outopts)      :: OutOptions    ! Output options.
  real(8),allocatable:: Displ(:,:),Veloc(:,:)

! Initialize.
  NumN=size(Node)
  NumE(1)=size(Elem)
  allocate (ListIN (NumN));
  do k=1,NumN
    ListIN(k)=Node(k)%Vdof
  end do
  gamma=1.d0/2.d0+Options%NewmarkDamp
  beta =1.d0/4.d0*(gamma+0.5d0)*(gamma+0.5d0)

! Allocate memory for solver (Use a conservative estimate of the size of the matrices).
  allocate (Asys   (DimMat*NumDof)); call sparse_zero (as,Asys)
  allocate (Mglobal(DimMat*NumDof)); call sparse_zero (ms,Mglobal)
  allocate (Cglobal(DimMat*NumDof)); call sparse_zero (cs,Cglobal)
  allocate (Kglobal(DimMat*NumDof)); call sparse_zero (ks,Kglobal)
  allocate (Fglobal(DimMat*NumDof)); call sparse_zero (fs,Fglobal)

  allocate (Qglobal(NumDof));   Qglobal= 0.d0
  allocate (Mvel   (NumDof,6)); Mvel   = 0.d0
  allocate (Cvel   (NumDof,6)); Cvel   = 0.d0

  allocate (X0     (NumDof));   X0     = 0.d0
  allocate (DX     (NumDof));   DX     = 0.d0
  allocate (DXDt   (NumDof));   DXDt   = 0.d0
  allocate (DXDDt  (NumDof));   DXDDt  = 0.d0

  allocate (Displ      (NumN,         6));      Displ   =    0.d0
  allocate (Veloc      (NumN,         6));      Veloc   =    0.d0
  allocate (DeltaPos   (NumN,         3));      DeltaPos=    0.d0
  allocate (DeltaPsi   (NumE(1),MaxElNod,3));   DeltaPsi=    0.d0
  allocate (DeltaPosDot(NumN,         3));      DeltaPosDot= 0.d0
  allocate (DeltaPsiDot(NumE(1),MaxElNod,3));   DeltaPsiDot= 0.d0

! Allocate quaternions and rotation operator for initially undeformed system
  Quat = (/1.d0,0.d0,0.d0,0.d0/); Cao = Unit; Temp = Unit4

!Find initial conditions
  call cbeam3_solv_disp2state (Node,Coords,Psi0,PosDotDefor,PsiDotDefor,X0,DXDt)
  call cbeam3_solv_disp2state (Node,PosDefor,PsiDefor,PosDotDefor,PsiDotDefor,DX,DXDt)

  DX=DX-X0

! Compute tangent matrices at initial time.
  call cbeam3_asbly_dynamic (Elem,Node,Coords,Psi0,Coords,Psi0,PosDotDefor,                                   &
&                            PsiDotDefor,0.d0*PosDefor,0.d0*PsiDefor,F0+Ftime(1)*Fa,Vrel(1,:),VrelDot(1,:),   &
&                            ms,Mglobal,Mvel,cs,Cglobal,Cvel,ks,Kglobal,fs,Fglobal,Qglobal,Options,Unit)

! Find initial acceleration.
  Qglobal= sparse_matvmul(fs,Fglobal,NumDof,fem_m2v(F0+Ftime(1)*Fa,NumDof,Filter=ListIN))   &
&          - matmul(Mvel,VrelDot(1,:)) - matmul(Cvel,Vrel(1,:))                             &
&          - sparse_matvmul(cs,Cglobal,NumDof,DXDt) - sparse_matvmul(ks,Kglobal,NumDof,DX)

  call sparse_addsparse(0,0,ms,Mglobal,as,Asys)
#ifdef NOLAPACK
  call lu_sparse(as,Asys,Qglobal,DXDDt)
#else
  call lapack_sparse (as,Asys,Qglobal,DXDDt)
#endif

! Loop in the time steps.
  do iStep=1,size(Time)-1
    dt     =Time(iStep+1)-Time(iStep)
    call out_time(iStep,Time(iStep+1),Text)
    write (*,'(5X,A)') trim(Text)

! Update transformation matrix for given angular velocity
    call lu_invers ((Unit4+0.25d0*xbeam_QuadSkew(Vrel(iStep+1,4:6))*(Time(iStep+1)-Time(iStep))),Temp)
    Quat=matmul(Temp,matmul((Unit4-0.25d0*xbeam_QuadSkew(Vrel(iStep,4:6))*(Time(iStep+1)-Time(iStep))),Quat))
    Cao = xbeam_Rot(Quat)

! Predictor step.
    DX= DX + dt*DXDt + (0.5d0-beta)*dt*dt*DXDDt
    DXDt=DXDt + (1-gamma)*dt*DXDDt

! Compute system functionals and matrices. Only updating Felast.
    call sparse_zero (ms,Mglobal); Mvel=0.d0
    call sparse_zero (cs,Cglobal); Cvel=0.d0
    call sparse_zero (ks,Kglobal);
    call sparse_zero (fs,Fglobal);

    call cbeam3_asbly_dynamic (Elem,Node,Coords,Psi0,PosDefor,PsiDefor,PosDotDefor,PsiDotDefor,                &
&                              0.d0*PosDefor,0.d0*PsiDefor,F0+Ftime(iStep+1)*Fa,Vrel(iStep+1,:),VrelDot(iStep+1,:),  &
&                              ms,Mglobal,Mvel,cs,Cglobal,Cvel,ks,Kglobal,fs,Fglobal,Qglobal,Options,Cao)

! Compute right-hand side of the equation.
    Qglobal= sparse_matvmul(fs,Fglobal,NumDof,fem_m2v(F0+Ftime(iStep+1)*Fa,NumDof,Filter=ListIN))
    Qglobal= Qglobal - matmul(Mvel,Vreldot(iStep+1,:)) &
&                    - matmul(Cvel,Vrel   (iStep+1,:))
    Qglobal= Qglobal - sparse_matvmul(cs,Cglobal,NumDof,DXDt) &
&                    - sparse_matvmul(ks,Kglobal,NumDof,DX)

! Compute left-hand side of the equation (if dt=constant then Asys=constant too, but this was not
! assumed here).
    call sparse_zero (as,Asys)
    call sparse_addsparse(0,0,ms,Mglobal,as,Asys,Factor=1.d0)
    call sparse_addsparse(0,0,cs,Cglobal,as,Asys,Factor=gamma*dt)
    call sparse_addsparse(0,0,ks,Kglobal,as,Asys,Factor=beta*dt*dt)

! Solve equation.
#ifdef NOLAPACK
    call lu_sparse(as,Asys,Qglobal,DXDDt)
#else
    call lapack_sparse (as,Asys,Qglobal,DXDDt)
#endif
    DX    = DX   + beta *dt*dt*DXDDt
    DXDt  = DXDt + gamma*dt   *DXDDt

! Update solution and store relevant information at current time step.
    call cbeam3_solv_update_lindyn (Elem,Node,PsiDefor,DX,DXDt,DeltaPos,DeltaPsi,DeltaPosDot,DeltaPsiDot)

! Store data in output variables.
    do k=1,NumN
      Veloc(k,1:3)= DeltaPosDot(k,:) + Vrel(iStep,1:3) &
&                 + rot_cross(Vrel(iStep,4:6),PosDefor(k,:)+DeltaPos(k,:))
    end do
    VelocTime (iStep+1, : )= Veloc   (1:NumN,3)

    PosPsiTime(iStep+1,1:3)= Coords(NumN,:) + DeltaPos(NumN,:)
    PosPsiTime(iStep+1,4:6)= Psi0(NumE(1),Elem(NumE(1))%NumNodes,:) &
&                          + DeltaPsi(NumE(1),Elem(NumE(1))%NumNodes,:)

!  Write output information in output file.
    if (any(OutGrids)) then
      do k=1,NumN
        if (OutGrids(k)) then
          Displ(k,1:3)= DeltaPos(k,:)
          Displ(k,4:6)= DeltaPsi(Node(k)%Master(1),Node(k)%Master(2),:)
        end if
      end do
      call out_title   (iOut,trim(Text))
      call out_outgrid (iOut,'NODE',OutOptions,1,NumE,6,OutGrids,DISPL=Displ,VELOC=Veloc)
    end if

    do k=1,NumN
        DynOut(iStep*NumN+k,:) = PosDefor(k,:) + DeltaPos(k,:)
    end do

  end do

! Write information at last time step.
  PosDefor= Coords + DeltaPos
  PsiDefor= Psi0 + DeltaPsi

  deallocate (ListIn,Mvel,Cvel)
  deallocate (Asys,Fglobal,Mglobal)
  deallocate (Kglobal,Cglobal,Qglobal)
  deallocate (DX,DXdt,DXDDt)
  deallocate (Displ,Veloc)
  deallocate (DeltaPos,DeltaPsi,DeltaPosDot,DeltaPsiDot)
  return
 end subroutine cbeam3_solv_lindyn


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-> Subroutine CBEAM3_SOLV_UPDATE_STATIC
!
!-> Description:
!
!    Update results from static solution.
!
!-> Remarks.-
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine cbeam3_solv_update_static (Elem,Node,Psi0,DeltaX,PosDefor,PsiDefor)
  use lib_fem
  use lib_rotvect
  use lib_bgeom
  use lib_cbeam3

! I/O Variables.
  type(xbelem),intent(in) :: Elem      (:)       ! Element information.
  type(xbnode),intent(in) :: Node      (:)       ! Nodal information.
  real(8),intent(in)      :: Psi0      (:,:,:)   ! Initial rotation vector at element grids.
  real(8),intent(in)      :: DeltaX    (:)       ! Incremental State Vector.
  real(8),intent(inout)   :: PosDefor  (:,:)     ! Current nodal position.
  real(8),intent(inout)   :: PsiDefor  (:,:,:)   ! Current rotation at element nodes.

! Local variables.
  integer :: i,j,k               ! Counters.
  integer :: ix                  ! Counter on the degrees of freedom.
  integer :: iElem               ! Counter on the elements.
  integer :: iNode               ! Counter on the nodes.

! Store current state in CBEAM3 elements.
  ix=0
  do iNode=1,size(PosDefor,DIM=1)
    iElem=Node(iNode)%Master(1)
    if ((Node(iNode)%Vdof.ne.0).and.(Elem(iElem)%MemNo.eq.0)) then
      k=Node(iNode)%Master(2)

! Nodal displacements.
      PosDefor(iNode,:)   = PosDefor(iNode,:)   + DeltaX(ix+1:ix+3)

! Cartesian rotation vector at master nodes.
      PsiDefor(iElem,k,:) = PsiDefor(iElem,k,:) + DeltaX(ix+4:ix+6)
      ix=ix+6
    end if
  end do

!!! Post-processing.
! Compute rotation vector at slave nodes of CBEAM3 elements.
  do i=1,size(Elem)
    if (Elem(i)%MemNo.eq.0) then

! Copy rotation from master node for each slave node.
      do j=1,Elem(i)%NumNodes
        if (Elem(i)%Master(j,1).ne.0) then
          PsiDefor(i,j,:)=PsiDefor(Elem(i)%Master(j,1),Elem(i)%Master(j,2),:)
        end if
      end do

! Include master-to-slave initial rotations from the undeformed configuration.
      call cbeam3_projm2s (Elem(i)%NumNodes,Elem(i)%Master,Psi0(i,:,:),Psi0,PsiDefor(i,:,:))
    end if
  end do

  return
 end subroutine cbeam3_solv_update_static



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-> Subroutine CBEAM3_SOLV_UPDATE_lindyn
!
!-> Description:
!
!    Update results from linear dynamic solution.
!
!-> Remarks.-
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine cbeam3_solv_update_lindyn (Elem,Node,Psi0,DX,DXDt,Pos,Psi,PosDot,PsiDot)
  use lib_fem
  use lib_rotvect
  use lib_bgeom
  use lib_cbeam3

! I/O Variables.
  type(xbelem),intent(in) :: Elem      (:)       ! Element information.
  type(xbnode),intent(in) :: Node      (:)       ! Nodal information.
  real(8),intent(in)      :: Psi0      (:,:,:)   ! Initial rotation vector at element grids.
  real(8),intent(in)      :: DX        (:)       ! Incremental State Vector.
  real(8),intent(in)      :: DXDt      (:)       ! Time derivative of DX.
  real(8),intent(inout)   :: Pos       (:,:)     ! Delta nodal position.
  real(8),intent(inout)   :: Psi       (:,:,:)   ! Delta rotation at element nodes.
  real(8),intent(inout)   :: PosDot    (:,:)     ! Delta nodal velocities.
  real(8),intent(inout)   :: PsiDot    (:,:,:)   ! Delta nodal angular velocities.

! Local variables.
  integer :: i,j,k               ! Counters.
  integer :: ix                  ! Counter on the degrees of freedom.
  integer :: iElem               ! Counter on the elements.
  integer :: iNode               ! Counter on the nodes.

! Store current displacement and its time derivative at all nodes and the
! rotations and its first derivatives at the master nodes of each element.
  ix=0
  do iNode=1,size(Pos,DIM=1)
    iElem=Node(iNode)%Master(1)
    k=Node(iNode)%Master(2)

    ! Constrained nodes.
    if (Node(iNode)%Vdof.eq.0) then
      Pos   (iNode,:)  = 0.d0
      Psi   (iElem,k,:)= 0.d0
      PosDot(iNode,:)  = 0.d0
      PsiDot(iElem,k,:)= 0.d0

    ! Unconstrained nodes.
    else
      Pos   (iNode,:)= DX(ix+1:ix+3)
      PosDot(iNode,:)= DXDt(ix+1:ix+3)

      Psi   (iElem,k,:)= DX(ix+4:ix+6)
      PsiDot(iElem,k,:)= DXDt(ix+4:ix+6)

      ix=ix+6
    end if
  end do

!!! Post-processing.
! Compute rotation vector at slave nodes of CBEAM3 elements.
  do i=1,size(Elem)

! Copy rotation and derivative from master node for each slave node.
    do j=1,Elem(i)%NumNodes
      if (Elem(i)%Master(j,1).ne.0) then
        Psi   (i,j,:)= Psi0 (Elem(i)%Master(j,1),Elem(i)%Master(j,2),:) &
&                    + Psi  (Elem(i)%Master(j,1),Elem(i)%Master(j,2),:)
        PsiDot(i,j,:)=PsiDot(Elem(i)%Master(j,1),Elem(i)%Master(j,2),:)
      end if
    end do

    ! Include master-to-slave initial rotations from the undeformed configuration.
    call cbeam3_projm2s (Elem(i)%NumNodes,Elem(i)%Master,Psi0(i,:,:),Psi0,Psi(i,:,:))
    
    ! Compute the delta value.
    do j=1,Elem(i)%NumNodes
      if (Elem(i)%Master(j,1).ne.0) then
        Psi   (i,j,:)= Psi (Elem(i)%Master(j,1),Elem(i)%Master(j,2),:) &
&                    - Psi0(Elem(i)%Master(j,1),Elem(i)%Master(j,2),:)
      end if
    end do
  end do


  return
 end subroutine cbeam3_solv_update_lindyn



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-> Subroutine CBEAM3_SOLV_STATE2DISP
!
!-> Description:
!
!    Extract current positiion, orientation and velocities from the
!    state vector.
!
!-> Remarks.-
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine cbeam3_solv_state2disp (Elem,Node,Coords,Psi0,X,dXdt,Pos,Psi,PosDot,PsiDot)
  use lib_fem
  use lib_rotvect
  use lib_bgeom
  use lib_cbeam3

! I/O Variables.
  type(xbelem),intent(in) :: Elem    (:)       ! Element information.
  type(xbnode),intent(in) :: Node    (:)       ! Nodal information.
  real(8),intent(in)      :: Coords  (:,:)     ! Initial coordinates of the grid points.
  real(8),intent(in)      :: Psi0    (:,:,:)   ! Initial rotation vector at element grids.
  real(8),intent(in)      :: X       (:)       ! Current generalized coordinates.
  real(8),intent(in)      :: dXdt    (:)       ! Time derivatives of X.
  real(8),intent(out)     :: Pos     (:,:)     ! Current nodal position.
  real(8),intent(out)     :: Psi     (:,:,:)   ! Current rotation at element nodes.
  real(8),intent(out)     :: PosDot  (:,:)     ! Current nodal position.
  real(8),intent(out)     :: PsiDot  (:,:,:)   ! Current rotation at element nodes.

! Local variables.
  integer :: i,j,k               ! Counters.
  integer :: ix                  ! Counter on the degrees of freedom.
  integer :: iElem               ! Counter on the elements.
  integer :: iNode               ! Counter on the nodes.

! Store current displacement and its time derivative at all nodes and the
! rotations and its first derivatives at the master nodes of each element.
  ix=0
  do iNode=1,size(Pos,DIM=1)
    iElem=Node(iNode)%Master(1)
    k=Node(iNode)%Master(2)

    ! Constrained nodes.
    if (Node(iNode)%Vdof.eq.0) then
      Pos   (iNode,:)= Coords(iNode,:)
      PosDot(iNode,:)= 0.d0

      Psi   (iElem,k,:)= Psi0(iElem,k,:)
      PsiDot(iElem,k,:)= 0.d0

    ! Unconstrained nodes.
    else
      Pos   (iNode,:)= X   (ix+1:ix+3)
      PosDot(iNode,:)= dXdt(ix+1:ix+3)

      Psi   (iElem,k,:)= X   (ix+4:ix+6)
      PsiDot(iElem,k,:)= dXdt(ix+4:ix+6)

      ix=ix+6
    end if
  end do

! Compute rotation vector and time derivative at slave nodes within elements.
  do i=1,size(Elem)
    do j=1,Elem(i)%NumNodes

! Copy rotation and derivative from master node for each slave node.
      if (Elem(i)%Master(j,1).ne.0) then
        Psi(i,j,:)   =Psi   (Elem(i)%Master(j,1),Elem(i)%Master(j,2),:)
        PsiDot(i,j,:)=PsiDot(Elem(i)%Master(j,1),Elem(i)%Master(j,2),:)
      end if
    end do

! Include master-to-slave initial rotations from the undeformed configuration.
    call cbeam3_projm2s (Elem(i)%NumNodes,Elem(i)%Master,Psi0(i,:,:),Psi0,Psi(i,:,:))
  end do

  return
 end subroutine cbeam3_solv_state2disp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-> Subroutine CBEAM3_SOLV_STATE2ACCEL
!
!-> Description:
!
!    Extract current acceleration from the state vector.
!
!-> Remarks:
!-> Author: Rob Simpson
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine cbeam3_solv_state2accel (Elem,Node,dXddt,PosDDot,PsiDDot)
  use lib_fem
  use lib_rotvect
  use lib_bgeom
  use lib_cbeam3

! I/O Variables.
  type(xbelem),intent(in) :: Elem    (:)       ! Element information.
  type(xbnode),intent(in) :: Node    (:)       ! Nodal information.
  real(8),intent(in)      :: dXddt    (:)       ! Time derivatives of X.
  real(8),intent(out)     :: PosDDot  (:,:)     ! Current nodal position.
  real(8),intent(out)     :: PsiDDot  (:,:,:)   ! Current rotation at element nodes.

! Local variables.
  integer :: i,j,k               ! Counters.
  integer :: ix                  ! Counter on the degrees of freedom.
  integer :: iElem               ! Counter on the elements.
  integer :: iNode               ! Counter on the nodes.

! Store current accelerations at all nodes and the
! rotational accel. at the master nodes of each element.
  ix=0
  do iNode=1,size(PosDDot,DIM=1)
    iElem=Node(iNode)%Master(1)
    k=Node(iNode)%Master(2)

    ! Constrained nodes.
    if (Node(iNode)%Vdof.eq.0) then

      PosDDot(iNode,:) = 0.d0
      PsiDDot(iElem,k,:) = 0.d0

    ! Unconstrained nodes.
    else
      PosDDot(iNode,:)=dXddt(ix+1:ix+3)
      PsiDDot(iElem,k,:) = dXddt(ix+4:ix+6)

      ix=ix+6
    end if
  end do

! Compute rotation vector and time derivative at slave nodes within elements.
  do i=1,size(Elem)
    do j=1,Elem(i)%NumNodes

! Copy rotational accel. from master node for each slave node.
      if (Elem(i)%Master(j,1).ne.0) then
        PsiDDot(i,j,:)=PsiDDot(Elem(i)%Master(j,1),Elem(i)%Master(j,2),:)
      end if
    end do
  end do

  return
 end subroutine cbeam3_solv_state2accel


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-> Subroutine CBEAM3_SOLV_DISP2STATE
!
!-> Description:
!
!    Write current state vector from current displacements and rotations.
!
!-> Remarks.-
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine cbeam3_solv_disp2state (Node,Pos,Psi,PosDot,PsiDot,X,dXdt)
  use lib_fem

! I/O Variables.
  type(xbnode),intent(in) :: Node      (:)       ! Nodal information.
  real(8),intent(in)      :: Pos       (:,:)     ! Current nodal position.
  real(8),intent(in)      :: Psi       (:,:,:)   ! Current rotation at element nodes.
  real(8),intent(in)      :: PosDot    (:,:)     ! Current nodal position.
  real(8),intent(in)      :: PsiDot    (:,:,:)   ! Current rotation at element nodes.
  real(8),intent(out)     :: X         (:)       ! Current displacements/rotations.
  real(8),intent(out)     :: dXdt      (:)       ! Current time derivative of X.


! Local variables.
  integer :: k                   ! Counters.
  integer :: ix                  ! Counter on the degrees of freedom.
  integer :: iElem               ! Counter on the elements.
  integer :: iNode               ! Counter on the nodes.

! Loop in all nodes in the model.
  ix=0
  do iNode=1,size(Pos,DIM=1)
    iElem=Node(iNode)%Master(1)
    if (Node(iNode)%Vdof.ne.0) then
      k=Node(iNode)%Master(2)

! Current nodal displacements and derivatives.
      X   (ix+1:ix+3)= Pos(iNode,:)
      dXdt(ix+1:ix+3)= PosDot(iNode,:)

! Cartesian rotation vector at master nodes.
      X   (ix+4:ix+6)= Psi(iElem,k,:)
      dXdt(ix+4:ix+6)= PsiDot(iElem,k,:)
      ix=ix+6
    end if
  end do

  return
 end subroutine cbeam3_solv_disp2state


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-> Subroutine CBEAM3_SOLV_ACCEL2STATE
!
!-> Description:
!
!    Extract 2nd time-derivative of state vector from accelerations.
!
!-> Remarks:
!-> Author: Rob Simpson
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine cbeam3_solv_accel2state (Node,PosDDot,PsiDDot,dXddt)
  use lib_fem

! I/O Variables.
  type(xbnode),intent(in) :: Node      (:)       ! Nodal information.
  real(8),intent(in)      :: PosDDot    (:,:)     ! Current nodal acceleration.
  real(8),intent(in)      :: PsiDDot    (:,:,:)   ! Current rotational accel. at element nodes.
  real(8),intent(out)     :: dXddt      (:)       ! Current second time-derivative of X.

! Local variables.
  integer :: k                   ! Counters.
  integer :: ix                  ! Counter on the degrees of freedom.
  integer :: iElem               ! Counter on the elements.
  integer :: iNode               ! Counter on the nodes.

! Loop in all nodes in the model.
  ix=0
  do iNode=1,size(PosDDot,DIM=1)
    iElem=Node(iNode)%Master(1)
    if (Node(iNode)%Vdof.ne.0) then
      k=Node(iNode)%Master(2)

! Current nodal displacements and derivatives.
      dXddt(ix+1:ix+3)= PosDDot(iNode,:)

! Cartesian rotation vector at master nodes.
      dXddt(ix+4:ix+6)= PsiDDot(iElem,k,:)
      ix=ix+6
    end if
  end do

  return
 end subroutine cbeam3_solv_accel2state

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end module cbeam3_solv
