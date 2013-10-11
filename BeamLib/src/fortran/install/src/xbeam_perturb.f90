!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-> Module.- XBEAM_PERTURB  Henrik Hesse. 17/05/2011 - Last Update 05/07/2011
!
!-> Language: FORTRAN90, Free format.
!
!-> Description.-
!
!   This module includes the high-level routines to solve flexibel-body dynamics problems for small geometrically-linear
!   deformations and large nonlinear rigid-body motion. This module includes the assembly routines for the system matrices and
!   the solver routine to solve the nonlinear system. The solution is based on the Newmark method using the full linearisation
!   at each subiteration in module xbeam_asbly.
!
!-> Subroutines.-
!
!   -xbeam_pert_asbly0:      Initial assembly of constant system matrices for linear elastic simulation.
!   -xbeam_perturb_asbly:    Assembly of dynamic perturbed system matrices.
!   -xbeam_perturb_solv:     Flexible-body dynamic solution for small deformations and nonlinear rigid-body motion.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module xbeam_perturb
 use xbeam_shared
 implicit none

! Shared variables within the module.
 integer,private,parameter:: MaxNodCB3=3               ! Max number of nodes per element is 3.
 integer,private,parameter:: DimMat=28                 ! Memory index for sparse matrices.

 real(8),private,parameter,dimension(4,4):: Unit4= &   ! 4x4 Unit matrix.
&         reshape((/1.d0,0.d0,0.d0,0.d0,0.d0,1.d0,0.d0,0.d0,0.d0,0.d0,1.d0,0.d0,0.d0,0.d0,0.d0,1.d0/),(/4,4/))

 contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-> Subroutine XBEAM_PERTURB_ASBLY0
!
!-> Description:
!
!   Initial assembly of constant system matrices for linear elastic simulation.
!
!-> Remarks:
!
!   - M and Kstiff is constant and only computed once in xbeam_perturb_asbly0.
!   - M is symmetric and hence MRS=transpose(MRS)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine xbeam_perturb_asbly0 (Elem,Node,NumDof,Coords,Psi0,iml,Mlin,ikstif,Kstif,Options,Cao)
  use lib_rotvect
  use lib_fem
  use lib_sparse
  use lib_cbeam3
  use lib_xbeam

! I/O variables.
  type(xbelem), intent(in) :: Elem(:)               ! Element information.
  type(xbnode), intent(in) :: Node(:)               ! List of independent nodes.
  integer,      intent(in) :: NumDof                ! Number of independent DoFs.
  real(8),      intent(in) :: Coords    (:,:)       ! Initial coordinates of the grid points.
  real(8),      intent(in) :: Psi0      (:,:,:)     ! Initial CRV of the nodes in the elements.

  integer,      intent(out):: iml               ! Size of the sparse mass matrix.
  type(sparse), intent(out):: Mlin  (:)         ! Sparse mass matrix.
  integer,      intent(out):: ikstif            ! Size of the sparse stiffness matrix.
  type(sparse), intent(out):: Kstif (:)         ! Sparse stiffness matrix.
  type(xbopts), intent(in) :: Options           ! Solver parameters.
  real(8),      intent(in) :: Cao   (:,:)    ! Rotation operator from reference to inertial frame

! Local variables.
  integer:: i,j,i1,j1                      ! Counters.
  integer:: iElem                          ! Counter on the finite elements.
  integer:: NumE                           ! Number of elements in the model.
  integer:: NumNE                          ! Number of nodes in an element.
  integer:: NumGaussMass                   ! Number of Gaussian points in the inertia terms.

  real(8):: MSSelem (6*MaxElNod,6*MaxElNod)  ! Element mass matrix.
  real(8):: MSRelem (6*MaxElNod,6)           ! Element reference-system mass matrix.
  real(8):: MRRelem (6,6)                    ! Mass from the motions of reference system.
  real(8):: KSSelem (6*MaxElNod,6*MaxElNod)  ! Element mass matrix.

  real(8):: rElem0(MaxElNod,6)             ! Initial Coordinates/CRV of nodes in the element.
  real(8):: SB2B1 (6*MaxElNod,6*MaxElNod)  ! Transformation from master to rigid node orientations.

! Allocate memory for solver (Use a conservative estimate of the size of the matrices).
  call sparse_zero (iml,Mlin);
  call sparse_zero (ikstif,Kstif);

! Loop in all elements in the model.
  NumE=size(Elem)

  do iElem=1,NumE
    MSSelem=0.d0; MSRelem=0.d0; MRRelem=0.d0;
    KSSelem=0.d0; SB2B1=0.d0

! Extract coords of elem nodes and determine if they are master (Flag=T) or slave.
    call fem_glob2loc_extract(Elem(iElem)%Conn,Coords,rElem0(:,1:3),NumNE)
    rElem0(:,4:6)= Psi0(iElem,:,:)

! Prevent singularities in mass matrices when using reduced integration.
	if (NumNE.gt.Options%NumGauss) then
		NumGaussMass=NumNE
	else
		NumGaussMass=Options%NumGauss
	end if

! Contributions of the structural mass to the linearised inertia matrices.
    call cbeam3_mass (NumNE,rElem0,rElem0,Elem(iElem)%Mass,MSSelem,NumGaussMass)
    call cbeam3_mvel (NumNE,rElem0,rElem0,Elem(iElem)%Mass,MSRelem,Options%NumGauss)
    call xbeam_mrr   (NumNE,rElem0,rElem0,Elem(iElem)%Mass,MRRelem,Options%NumGauss)
    
! Compute the element tangent stiffness matrix.
    call cbeam3_kmat  (NumNE,rElem0,rElem0,Elem(iElem)%Stiff,KSSelem,Options%NumGauss)
    call cbeam3_kgeom (NumNE,rElem0,rElem0,Elem(iElem)%Stiff,KSSelem,Options%NumGauss)

! Compute the element mass tangent stiffness matrix (can be neglected).
   !call xbeam_kmass  (NumNE,rElem0,rElem,rElemDDot,VrelDot,Elem(iElem)%Mass,KRSelem,Options%NumGauss)

! Project slave degrees of freedom to the orientation of the "master" ones.
    call cbeam3_projs2m (NumNE,Elem(iElem)%Master(:,:),Psi0(iElem,:,:),Psi0,SB2B1)
    MSSelem=matmul(transpose(SB2B1),matmul(MSSelem,SB2B1))
    MSRelem=matmul(transpose(SB2B1),MSRelem)
    KSSelem=matmul(transpose(SB2B1),matmul(KSSelem,SB2B1))

! Add to global matrix. Remove columns and rows at clamped points.
    call sparse_addmat (NumDof,NumDof,MRRelem,iml,Mlin)

    do i=1,NumNE
      i1=Node(Elem(iElem)%Conn(i))%Vdof
      if (i1.ne.0) then
        call sparse_addmat (6*(i1-1),NumDof,MSRelem(6*(i-1)+1:6*i,:),iml,Mlin)
        call sparse_addmat (NumDof,6*(i1-1),transpose(MSRelem(6*(i-1)+1:6*i,:)),iml,Mlin)

        do j=1,NumNE
          j1=Node(Elem(iElem)%Conn(j))%Vdof
          if (j1.ne.0) then
            call sparse_addmat (6*(i1-1),6*(j1-1),MSSelem(6*(i-1)+1:6*i,6*(j-1)+1:6*j),iml   ,Mlin)
            call sparse_addmat (6*(i1-1),6*(j1-1),KSSelem(6*(i-1)+1:6*i,6*(j-1)+1:6*j),ikstif,Kstif)
          end if
        end do

      end if
    end do

  end do

! Include quaternion matrix.
  call sparse_addmat (NumDof+6,NumDof+6,Unit4,iml,Mlin)

  return
 end subroutine xbeam_perturb_asbly0


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-> Subroutine XBEAM_PERTURB_ASBLY
!
!-> Description:
!
!   Assembly of dynamic perturbed system matrices.
!
!-> Remarks:
!
!   - Kstiff is constant and only computed once in xbeam_perturb_asbly0 and inputed here.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine xbeam_perturb_asbly (Elem,Node,NumDof,Coords,Psi0,PosDefor,PsiDefor,Vrel,VrelDot,Quat,                      &
&                               ikstif,Kstif,icp,Cpert,ikp,Kpert,icl,Clin,ikl,Klin,ifs,Felast,ifr,Frigid,Options,Cao)
  use lib_rotvect
  use lib_fem
  use lib_sparse
  use lib_cbeam3
  use lib_xbeam
  use lib_perturb

! I/O variables.
  type(xbelem), intent(in) :: Elem(:)               ! Element information.
  type(xbnode), intent(in) :: Node(:)               ! List of independent nodes.
  integer,      intent(in) :: NumDof                ! Number of independent DoFs.
  real(8),      intent(in) :: Coords    (:,:)       ! Initial coordinates of the grid points.
  real(8),      intent(in) :: Psi0      (:,:,:)     ! Initial CRV of the nodes in the elements.
  real(8),      intent(in) :: PosDefor  (:,:)       ! Current coordinates of the grid points
  real(8),      intent(in) :: PsiDefor  (:,:,:)     ! Current CRV of the nodes in the elements.
  real(8),      intent(in) :: Vrel(6), VrelDot(6)   ! Velocity of reference frame and derivative.
  real(8),      intent(in) :: Quat(4)               ! Quaternions.
  integer,      intent(in) :: ikstif                ! Size of the sparse stiffness matrix.
  type(sparse), intent(in) :: Kstif (:)             ! Sparse stiffness matrix.

  integer,      intent(out):: icp                 ! Size of the sparse perturbed damping matrix.
  type(sparse), intent(out):: Cpert(:)            ! Sparse perturbed damping matrix.
  integer,      intent(out):: ikp                 ! Size of the sparse perturbed stiffness matrix.
  type(sparse), intent(out):: Kpert(:)            ! Sparse perturbed stiffness matrix.
  integer,      intent(out):: icl                 ! Size of the sparse damping matrix.
  type(sparse), intent(out):: Clin (:)            ! Sparse damping matrix.
  integer,      intent(out):: ikl                 ! Size of the sparse stiffness matrix.
  type(sparse), intent(out):: Klin (:)            ! Sparse stiffness matrix.
  integer,      intent(out):: ifs                 ! Size of the sparse stiffness matrix.
  type(sparse), intent(out):: Felast (:)          ! Elastic influence coefficients matrix for applied forces.
  integer,      intent(out):: ifr                 ! Size of the sparse stiffness matrix.
  type(sparse), intent(out):: Frigid (:)          ! Rigid influence coefficients matrix for applied forces.

  type(xbopts), intent(in) :: Options             ! Solver parameters.
  real(8),      intent(in) :: Cao    (:,:)        ! Rotation operator from reference to inertial frame

! Local variables.
  logical:: Flags(MaxElNod)                  ! Auxiliary flags.
  integer:: i,j,i1,j1                        ! Counters.
  integer:: iElem                            ! Counter on the finite elements.
  integer:: NumE                             ! Number of elements in the model.
  integer:: NumNE                            ! Number of nodes in an element.
  
  ! Define element matrices
  real(8):: CSSelem (6*MaxElNod,6*MaxElNod)  ! Element damping matrix.
  real(8):: CSRelem (6*MaxElNod,6)           ! Element reference-system damping matrix.
  real(8):: CRSelem (6,6*MaxElNod)           ! Element damping matrix.
  real(8):: CRRelem (6,6)                    ! Reference-system damping matrix.
  real(8):: KSSelem (6*MaxElNod,6*MaxElNod)  ! Element tangent stiffness matrix.
  real(8):: KRSelem (6,6*MaxElNod)           ! Element tangent stiffness matrix.

  real(8):: FSelem (6*MaxElNod,6*MaxElNod)   ! Element force influence coefficients.
  real(8):: FRelem (6,6*MaxElNod)            ! Element force influence coefficients.
  
  real(8):: rElem0(MaxElNod,6)               ! Initial Coordinates/CRV of nodes in the element.
  real(8):: rElem (MaxElNod,6)               ! Current Coordinates/CRV of nodes in the element.
  real(8):: SB2B1 (6*MaxElNod,6*MaxElNod)    ! Transformation from master to rigid node orientations.

  real(8):: CQR(4,6),CQQ(4,4)                ! Tangent matrices from linearisation of quaternion equation.

! Allocate memory for solver (Use a conservative estimate of the size of the matrices).
  call sparse_zero (icl,Clin)
  call sparse_zero (ikl,Klin)
  call sparse_zero (icp,Cpert)
  call sparse_zero (ikp,Kpert)
  call sparse_zero (ifs,Felast)
  call sparse_zero (ifr,Frigid)
  CQR = 0.d0; CQQ = 0.d0

! Loop in all elements in the model.
  NumE=size(Elem)

  do iElem=1,NumE
    CSSelem=0.d0; CSRelem=0.d0; CRSelem=0.d0; CRRelem=0.d0; KSSelem=0.d0; KRSelem=0.d0;
    FSelem=0.d0; FRelem=0.d0; SB2B1=0.d0

! Extract coords of elem nodes and determine if they are master (Flag=T) or slave.
    Flags=.false.
    do i=1,Elem(iElem)%NumNodes
      if (Node(Elem(iElem)%Conn(i))%Master(1).eq.iElem) Flags(i)=.true.
    end do

    call fem_glob2loc_extract(Elem(iElem)%Conn,PosDefor,rElem (:,1:3),NumNE)
    call fem_glob2loc_extract(Elem(iElem)%Conn,Coords  ,rElem0(:,1:3),NumNE)
    rElem0(:,4:6)= Psi0    (iElem,:,:)
    rElem (:,4:6)= PsiDefor(iElem,:,:)

! Compute matrices for linearised system.
!=========================================
    call cbeam3_cgyr (NumNE,rElem0,rElem0,0*rElem0,Vrel,Elem(iElem)%Mass,CSSelem,Options%NumGauss)
    call cbeam3_cvel (NumNE,rElem0,rElem0,0*rElem0,Vrel,Elem(iElem)%Mass,CSRelem,Options%NumGauss)
    call xbeam_cgyr  (NumNE,rElem0,rElem0,0*rElem0,Vrel,Elem(iElem)%Mass,CRSelem,Options%NumGauss)
    call xbeam_crr   (NumNE,rElem0,rElem0,0*rElem0,Vrel,Elem(iElem)%Mass,CRRelem,Options%NumGauss)

    call cbeam3_kgyr (NumNE,rElem0,rElem0,0*rElem0,0*rElem0,Vrel,VrelDot,Elem(iElem)%Mass,KSSelem,Options%NumGauss)
    call xbeam_kgyr  (NumNE,rElem0,rElem0,0*rElem0,0*rElem0,Vrel,VrelDot,Elem(iElem)%Mass,KRSelem,Options%NumGauss)

    ! Project slave degrees of freedom to the orientation of the "master" ones.
    call cbeam3_projs2m (NumNE,Elem(iElem)%Master(:,:),Psi0(iElem,:,:),Psi0,SB2B1)
    CSSelem=matmul(transpose(SB2B1),matmul(CSSelem,SB2B1))
    CSRelem=matmul(transpose(SB2B1),CSRelem)
    CRSelem=matmul(CRSelem,SB2B1)
    KSSelem=matmul(transpose(SB2B1),matmul(KSSelem,SB2B1))
    KRSelem=matmul(KRSelem,SB2B1)

    ! Compute the influence coefficients multiplying the vector of external forces.
    call cbeam3_fext (NumNE,rElem,Flags(1:NumNE),FSelem,Options%FollowerForce,Options%FollowerForceRig,Cao)
    call xbeam_fext  (NumNE,rElem,Flags(1:NumNE),FRelem,Options%FollowerForce,Options%FollowerForceRig,Cao)

    ! Add to global linear matrices.
    call sparse_addmat (NumDof,NumDof,CRRelem,icl,Clin)

    do i=1,NumNE
      i1=Node(Elem(iElem)%Conn(i))%Vdof
      call sparse_addmat (0,6*(i1),FRelem(:,6*(i-1)+1:6*i),ifr,Frigid)
      if (i1.ne.0) then
        call sparse_addmat (6*(i1-1),NumDof,CSRelem(6*(i-1)+1:6*i,:),icl,Clin)
        call sparse_addmat (NumDof,6*(i1-1),CRSelem(:,6*(i-1)+1:6*i),icl,Clin)
        call sparse_addmat (NumDof,6*(i1-1),KRSelem(:,6*(i-1)+1:6*i),ikl,Klin)

        do j=1,NumNE
          j1=Node(Elem(iElem)%Conn(j))%Vdof
          if (j1.ne.0) then
            call sparse_addmat (6*(i1-1),6*(j1-1),CSSelem(6*(i-1)+1:6*i,6*(j-1)+1:6*j),icl,Clin)
            call sparse_addmat (6*(i1-1),6*(j1-1),KSSelem(6*(i-1)+1:6*i,6*(j-1)+1:6*j),ikl,Klin)
            call sparse_addmat (6*(i1-1),6*(j1-1),FSelem (6*(i-1)+1:6*i,6*(j-1)+1:6*j),ifs,Felast)
          end if
        end do
      end if
    end do 

! Compute matrices for perturbed system.
!=========================================
    CSSelem=0.d0; CSRelem=0.d0; CRSelem=0.d0; CRRelem=0.d0; KSSelem=0.d0; KRSelem=0.d0;

    call perturb_css (NumNE,rElem0,rElem0,Vrel,Elem(iElem)%Mass,CSSelem,Options%NumGauss)
    call perturb_csr (NumNE,rElem0,rElem0,Vrel,Elem(iElem)%Mass,CSRelem,Options%NumGauss)
    call perturb_crs (NumNE,rElem0,rElem0,Vrel,Elem(iElem)%Mass,CRSelem,Options%NumGauss)
    call perturb_crr (NumNE,rElem0,rElem0,Vrel,Elem(iElem)%Mass,CRRelem,Options%NumGauss)

    call perturb_kss (NumNE,rElem0,rElem0,Vrel,Elem(iElem)%Mass,KSSelem,Options%NumGauss)
    call perturb_krs (NumNE,rElem0,rElem0,Vrel,Elem(iElem)%Mass,KRSelem,Options%NumGauss)

    ! Project slave degrees of freedom to the orientation of the "master" ones.
    CSSelem=matmul(transpose(SB2B1),matmul(CSSelem,SB2B1))
    CSRelem=matmul(transpose(SB2B1),CSRelem)
    CRSelem=matmul(CRSelem,SB2B1)
    KSSelem=matmul(transpose(SB2B1),matmul(KSSelem,SB2B1))
    KRSelem=matmul(KRSelem,SB2B1)

    ! Add to global perturbed matrices.
    call sparse_addmat (NumDof,NumDof,CRRelem,icp,Cpert)

    do i=1,NumNE
      i1=Node(Elem(iElem)%Conn(i))%Vdof
      if (i1.ne.0) then
        call sparse_addmat (6*(i1-1),NumDof,CSRelem(6*(i-1)+1:6*i,:),icp,Cpert)
        call sparse_addmat (NumDof,6*(i1-1),CRSelem(:,6*(i-1)+1:6*i),icp,Cpert)
        call sparse_addmat (NumDof,6*(i1-1),KRSelem(:,6*(i-1)+1:6*i),ikp,Kpert)
        do j=1,NumNE
          j1=Node(Elem(iElem)%Conn(j))%Vdof
          if (j1.ne.0) then
            call sparse_addmat (6*(i1-1),6*(j1-1),CSSelem(6*(i-1)+1:6*i,6*(j-1)+1:6*j),icp,Cpert)
            call sparse_addmat (6*(i1-1),6*(j1-1),KSSelem(6*(i-1)+1:6*i,6*(j-1)+1:6*j),ikp,Kpert)
          end if
        end do
      end if
    end do

  end do

! Include initial elastic stiffness matrix
  call sparse_addsparse (0,0,ikstif,Kstif,ikl,Klin)
  call sparse_addsparse (0,0,ikstif,Kstif,ikp,Kpert)

! Compute tangent matrices for quaternion equations
  CQR(:,:)=0.d0
  CQR(1,4)= Quat(2); CQR(1,5)= Quat(3); CQR(1,6)= Quat(4)
  CQR(2,4)=-Quat(1); CQR(2,5)= Quat(4); CQR(2,6)=-Quat(3)
  CQR(3,4)=-Quat(4); CQR(3,5)=-Quat(1); CQR(3,6)= Quat(2)
  CQR(4,4)=-Quat(3); CQR(4,5)=-Quat(2); CQR(4,6)=-Quat(1)
  CQR=0.5d0*CQR

  CQQ=0.5d0*xbeam_QuadSkew(Vrel(4:6))

! Contribution of quaternions to linear damping matrix
  call sparse_addmat    (NumDof+6,NumDof  ,CQR,icl,Clin)
  call sparse_addmat    (NumDof+6,NumDof+6,CQQ,icl,Clin)

! Contribution of quaternions to perturbed damping matrix
  call sparse_addmat    (NumDof+6,NumDof+6,CQQ,icp,Cpert)

  return
 end subroutine xbeam_perturb_asbly


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-> Subroutine XBEAM_PERTURB_SOLV
!
!-> Description:
!
!    Flexible-body dynamic solution of multibeam problem under applied forces,
!    coupled with geometrically-linear deformations and nonlinear rigid-body motion
!
!-> Remarks.-
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine xbeam_perturb_solv (iOut,NumDof,Time,Elem,Node,F0,Fa,Ftime,             &
&                               Vrel,VrelDot,Quat,Coords,Psi0,PosDefor,PsiDefor,    &
&                               PosDotDefor,PsiDotDefor,DynOut,Options)
  use lib_fem
  use lib_rot
  use lib_rotvect
  use lib_out
  use lib_sparse
  use lib_xbeam
  use cbeam3_solv
  use lib_lu

! I/O Variables.
  integer,      intent(in)   :: iOut              ! Output file.
  integer,      intent(in)   :: NumDof            ! Number of independent DoFs.
  real(8),      intent(in)   :: Time      (:)     ! Time steps.
  type(xbelem), intent(in)   :: Elem      (:)     ! Element information.
  type(xbnode), intent(in)   :: Node      (:)     ! Nodal information.
  real(8),      intent(in)   :: F0        (:,:)   ! Applied static nodal forces.
  real(8),      intent(in)   :: Fa        (:,:)   ! Amplitude of the dynamic nodal forces.
  real(8),      intent(in)   :: Ftime     (:)     ! Time history of the applied forces.
  real(8),      intent(out)  :: Vrel      (:,:)   ! Time history of the velocities of the reference frame.
  real(8),      intent(out)  :: VrelDot   (:,:)   ! Time history of the accelerations of the reference frame.
  real(8),      intent(inout):: Quat      (:)     ! Quaternions to describes motion of reference system.
  real(8),      intent(in)   :: Coords    (:,:)   ! Initial coordinates of the grid points.
  real(8),      intent(in)   :: Psi0      (:,:,:) ! Initial CRV of the nodes in the elements.
  real(8),      intent(inout):: PosDefor  (:,:)   ! Current coordinates of the grid points
  real(8),      intent(inout):: PsiDefor  (:,:,:) ! Current CRV of the nodes in the elements.
  real(8),      intent(inout):: PosDotDefor(:,:)  ! Current time derivatives of the coordinates of the grid points
  real(8),      intent(inout):: PsiDotDefor(:,:,:)! Current time derivatives of the CRV of the nodes in the elements.
  real(8),      intent(out)  :: DynOut    (:,:)   ! Time-history of displacement of all nodes.
  type(xbopts), intent(in)   :: Options           ! Solver parameters.

! Local variables.
  real(8):: beta,gamma                     ! Newmark coefficients.
  real(8):: dt                             ! Time step
  integer:: k                              ! Counters.
  integer:: iStep,Iter                     ! Counters on time steps and subiterations.
  real(8):: MinDelta                       ! Value of Delta for convergence.
  integer:: NumE(1)                        ! Number of elements in the model.
  integer:: NumN                           ! Number of nodes in the model.

  real(8),allocatable:: dQdt(:),dQddt(:)   ! Generalized coordinates and derivatives for coupled system.
  real(8),allocatable:: Q(:), DQ(:)

  integer,allocatable:: ListIN     (:)     ! List of independent nodes.
  real(8),allocatable:: DeltaPos   (:,:)   ! Initial/final position vector of grid points
  real(8),allocatable:: DeltaPsi   (:,:,:) ! Initial/final CRV of the nodes in the elements.
  real(8),allocatable:: DeltaPosDot(:,:)   ! Current time derivatives of the coordinates of the grid points
  real(8),allocatable:: DeltaPsiDot(:,:,:) ! Current time derivatives of the CRV of the nodes in the elements.

  ! Define variables for rigid-body motion
  real(8) :: Cao(3,3),ACoa(6,6)            ! Rotation operator from reference to inertial frame

  ! Define variables for system matrices.
  integer:: ias,iml,ikstif,icl,ikl,icp,ikp,ifs,ifr
  type(sparse),allocatable:: Asys(:)      ! System matrix for implicit Newmark method.
  type(sparse),allocatable:: Mlin(:)      ! Total mass matrix in sparse storage.
  type(sparse),allocatable:: Clin(:)      ! Total Sparse damping matrix for linear system.
  type(sparse),allocatable:: Cpert(:)     ! Total Sparse damping matrix for perturbed system.
  type(sparse),allocatable:: Kstif(:)     ! Stiffness matrix in sparse storage for complete system.
  type(sparse),allocatable:: Klin(:)      ! Total stiffness matrix in sparse storage for linear system.
  type(sparse),allocatable:: Kpert(:)     ! Total stiffness matrix in sparse storage for perturbed system.
  real(8),allocatable     :: Qpert(:)    ! Total vector of discrete generalize forces.
  type(sparse),allocatable:: Felast(:)    ! Applied external forces on structure
  type(sparse),allocatable:: Frigid(:)    ! rigid matrix of applied forces in sparse format

  ! Define vaiables for output information.
  character(len=80)  ::  Text          ! Text with current time information to print out.
  real(8),allocatable::  Displ(:,:)    ! Current nodal displacements/rotations.
  real(8),allocatable::  Veloc(:,:)    ! Current nodal velocities.

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
  allocate (Asys(4*DimMat*NumDof)); call sparse_zero (ias,Asys)

  allocate (Mlin  (DimMat*NumDof));   call sparse_zero (iml,Mlin)
  allocate (Clin  (DimMat*NumDof));   call sparse_zero (icl,Clin)
  allocate (Klin  (DimMat*NumDof));   call sparse_zero (ikl,Klin)
  allocate (Kstif (DimMat*NumDof));   call sparse_zero (ikstif,Kstif)
  allocate (Felast(DimMat*NumDof));   call sparse_zero (ifs,Felast)
  allocate (Frigid(DimMat*NumDof));   call sparse_zero (ifr,Frigid)
  allocate (Cpert(DimMat*NumDof));    call sparse_zero (icp,Cpert)
  allocate (Kpert(DimMat*NumDof));    call sparse_zero (icp,Kpert)
  allocate (Qpert(NumDof+6+4));       Qpert= 0.d0

  ! Updated state vector with rigid body states and quaternions
  allocate (Q     (NumDof+6+4)); Q      = 0.d0
  allocate (dQdt  (NumDof+6+4)); dQdt   = 0.d0
  allocate (dQddt (NumDof+6+4)); dQddt  = 0.d0
  allocate (DQ    (NumDof+6+4)); DQ     = 0.d0

! Compute system information at initial condition.
  allocate (Displ      (NumN,            6)); Displ       = 0.d0
  allocate (Veloc      (NumN,            6)); Veloc       = 0.d0
  allocate (DeltaPos   (NumN,            3)); DeltaPos    = 0.d0
  allocate (DeltaPsi   (NumE(1),MaxElNod,3)); DeltaPsi    = 0.d0
  allocate (DeltaPosDot(NumN,            3)); DeltaPosDot = 0.d0
  allocate (DeltaPsiDot(NumE(1),MaxElNod,3)); DeltaPsiDot = 0.d0

  ! Extract initial displacements, velocities and orientation.
  dQdt(NumDof+1:NumDof+6) = Vrel(1,:)
  dQdt(NumDof+7:NumDof+10)= Quat

  ACoa= 0.d0
  Cao = xbeam_Rot(dQdt(NumDof+7:NumDof+10))

! Compute initial acceleration (we are neglecting qdotdot in Kmass).
  call xbeam_perturb_asbly0 (Elem,Node,NumDof,Coords,Psi0,iml,Mlin,ikstif,Kstif,Options,Cao)

  call xbeam_perturb_asbly (Elem,Node,NumDof,Coords,Psi0,PosDefor,PsiDefor,                                         &
&                           dQdt(NumDof+1:NumDof+6),0.d0*dQddt(NumDof+1:NumDof+6),dQdt(NumDof+7:NumDof+10),         &
&                           ikstif,Kstif,icp,Cpert,ikp,Kpert,icl,Clin,ikl,Klin,ifs,Felast,ifr,Frigid,Options,Cao)

! Compute initial accelerations
  Qpert(1:NumDof)          = sparse_matvmul(ifs,Felast,NumDof,fem_m2v(F0+Ftime(iStep+1)*Fa,NumDof,Filter=ListIN))
  Qpert(NumDof+1:NumDof+6) = sparse_matvmul(ifr,Frigid,6,fem_m2v(F0+Ftime(iStep+1)*Fa,NumDof+6))

  call lu_sparse (iml,Mlin,Qpert,dQddt)

! Loop in the time steps.
  do iStep=1,size(Time)-1
    dt= Time(iStep+1)-Time(iStep)
    call out_time(iStep,Time(iStep+1),Text)
    write (*,'(5X,A,$)') trim(Text)

! Predictor step.
    Q    = Q + dt*dQdt + (0.5d0-beta)*dt*dt*dQddt
    dQdt = dQdt + (1.d0-gamma)*dt*dQddt
    dQddt= 0.d0

! Iteration until convergence.
    do Iter=1,Options%MaxIterations+1
      if (Iter.gt.Options%MaxIterations) STOP 'Solution did not converge (18235)'

      ! Update nodal positions, velocities and orientation.
      call cbeam3_solv_update_lindyn (Elem,Node,PsiDefor,Q(1:NumDof),DQDt(1:NumDof),DeltaPos,DeltaPsi,DeltaPosDot,DeltaPsiDot)
      Cao = xbeam_Rot(dQdt(NumDof+7:NumDof+10))

      ! Compute system functionals and matrices.
      call sparse_zero (ifs,Felast)
      call sparse_zero (ifr,Frigid)
      call sparse_zero (icp,Cpert)
      call sparse_zero (icp,Kpert)
      call sparse_zero (icl,Clin)
      call sparse_zero (ikl,Klin)

      call xbeam_perturb_asbly (Elem,Node,NumDof,Coords,Psi0,PosDefor+DeltaPos,PsiDefor+DeltaPsi,                       &
&                               dQdt(NumDof+1:NumDof+6),0.d0*dQddt(NumDof+1:NumDof+6),dQdt(NumDof+7:NumDof+10),         &
&                               ikstif,Kstif,icp,Cpert,ikp,Kpert,icl,Clin,ikl,Klin,ifs,Felast,ifr,Frigid,Options,Cao)

      ! Compute the residual.
      Qpert= sparse_matvmul(iml,Mlin,NumDof+6+4,dQddt)
      Qpert= Qpert + sparse_matvmul(icp,Cpert,NumDof+6+4,dQdt)
      Qpert= Qpert + sparse_matvmul(ikp,Kpert,NumDof+6+4,Q)

      ! Compute admissible error.
      MinDelta=Options%MinDelta*max(1.d0,maxval(abs(Qpert(1:NumDof))))

      Qpert(1:NumDof) = Qpert(1:NumDof) - sparse_matvmul(ifs,Felast,NumDof,fem_m2v(F0+Ftime(iStep+1)*Fa,NumDof,Filter=ListIN))
      Qpert(NumDof+1:NumDof+6) = Qpert(NumDof+1:NumDof+6) - sparse_matvmul(ifr,Frigid,6,fem_m2v(F0+Ftime(iStep+1)*Fa,NumDof+6))

      ! Check convergence.
      if (maxval(abs(Qpert)).lt.MinDelta) then
        write (*,'(5X,A,I4,A,1PE12.3)') 'Subiteration',Iter, '  Delta=', maxval(abs(Qpert))
        exit
      end if

      ! Compute Jacobian
      call sparse_zero (ias,Asys)
      call sparse_addsparse(0,0,ikl,Klin,ias,Asys,Factor=1.d0)
      call sparse_addsparse(0,0,icl,Clin,ias,Asys,Factor=gamma/(beta*dt))
      call sparse_addsparse(0,0,iml,Mlin,ias,Asys,Factor=1.d0/(beta*dt*dt))

      ! Calculation of the correction.
      call lu_sparse (ias,Asys,-Qpert,DQ)

      Q     = Q     + DQ
      dQdt  = dQdt  + gamma/(beta*dt)*DQ
      dQddt = dQddt + 1.d0/(beta*dt*dt)*DQ
    end do

! Postprocessing
    Quat= dQdt(NumDof+7:NumDof+10)
    Cao = xbeam_Rot(Quat)
    ACoa(1:3,1:3) = transpose(Cao)
    ACoa(4:6,4:6) = transpose(Cao)

! Export velocities and accelerations in body frame
    if (Options%OutInaframe) then
        Vrel   (iStep+1,:) = dQdt (NumDof+1:NumDof+6)
        VrelDot(iStep+1,:) = dQddt(NumDof+1:NumDof+6)
! Export velocities and accelerations in inertial frame
    else
        Vrel   (iStep+1,:) = matmul(ACoa,dQdt (NumDof+1:NumDof+6))
        VrelDot(iStep+1,:) = matmul(ACoa,dQddt(NumDof+1:NumDof+6))
    end if

    do k=1,NumN
        DynOut((iStep-1)*NumN+k,:) = PosDefor(k,:) + DeltaPos(k,:)
    end do

  end do

! Write information at last time step.
  PosDefor= PosDefor + DeltaPos
  PsiDefor= PsiDefor + DeltaPsi

  deallocate (Asys,Mlin,Klin,Clin,Felast,Frigid)
  deallocate (Kpert,Cpert,Qpert)
  deallocate (Q,DQ,dQdt,dQddt)
  deallocate (ListIN,Displ,Veloc)
  return
 end subroutine xbeam_perturb_solv

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end module xbeam_perturb

