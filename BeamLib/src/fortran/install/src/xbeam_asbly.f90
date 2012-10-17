!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-> Module.- XBEAM_ASBLY Henrik Hesse. 07/01/2011 - Last Update 05/07/2011
!
!-> Language: FORTRAN90, Free format.
!
!-> Description.-
!
!  Assembly rigid-body components of beam equations.
!
!-> Subroutines.-
!
!    -xbeam_asbly_dynamic:        Assembly rigid-body matrices for the dynamic problem.
!
!-> Remarks.-
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module xbeam_asbly
 use xbeam_shared
 implicit none

 ! Shared variables within the module.
 integer,private,parameter:: MaxNodCB3=3               ! Max number of nodes per element is 3.

 contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-> Subroutine XBEAM_ASBLY_DYNAMIC
!
!-> Description:
!
!   Assembly rigid-body matrices for the dynamic problem.
!
!-> Remarks.-
!
!   - Check influence of mass stiffness matrix
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine xbeam_asbly_dynamic (Elem,Node,Coords,Psi0,PosDefor,PsiDefor,PosDeforDot,PsiDeforDot,PosDeforDDot,PsiDeforDDot,  &
&                               Vrel,VrelDot,Quat,ms,MRS,MRR,cs,CRS,CRR,CQR,CQQ,ks,KRS,fs,Frigid,Qrigid,Options,Cao)
  use lib_rotvect
  use lib_fem
  use lib_sparse
  use lib_cbeam3
  use lib_xbeam
!  use xbeam_fdiff

! I/O variables.
  type(xbelem), intent(in) :: Elem(:)               ! Element information.
  type(xbnode), intent(in) :: Node(:)               ! List of independent nodes.
  real(8),      intent(in) :: Coords    (:,:)       ! Initial coordinates of the grid points.
  real(8),      intent(in) :: Psi0      (:,:,:)     ! Initial CRV of the nodes in the elements.
  real(8),      intent(in) :: PosDefor  (:,:)       ! Current coordinates of the grid points
  real(8),      intent(in) :: PsiDefor  (:,:,:)     ! Current CRV of the nodes in the elements.
  real(8),      intent(in) :: PosDeforDot  (:,:)    ! Current coordinates of the grid points
  real(8),      intent(in) :: PsiDeforDot  (:,:,:)  ! Current CRV of the nodes in the elements.
  real(8),      intent(in) :: PosDeforDDot  (:,:)   ! Current coordinates of the grid points
  real(8),      intent(in) :: PsiDeforDDot  (:,:,:) ! Current CRV of the nodes in the elements.
  real(8),      intent(in) :: Vrel(6), VrelDot(6)   ! Velocity of reference frame and derivative.
  real(8),      intent(in) :: Quat(4)               ! Quaternions.

  integer,      intent(out):: ms                ! Size of the sparse mass matrix.
  type(sparse), intent(out):: MRS(:)            ! Sparse mass matrix.
  real(8),      intent(out):: MRR(:,:)          ! Reference system mass matrix.
  integer,      intent(out):: cs                ! Size of the sparse damping matrix.
  type(sparse), intent(out):: CRS(:)            ! Sparse damping matrix.
  real(8),      intent(out):: CRR(:,:)          ! Reference system damping matrix.
  real(8),      intent(out):: CQR(:,:),CQQ(:,:) ! Tangent matrices from linearisation of quaternion equation.
  integer,      intent(out):: ks                ! Size of the sparse stiffness matrix.
  type(sparse), intent(out):: KRS(:)            ! Sparse stiffness matrix.
  integer,      intent(out):: fs                ! Size of the sparse stiffness matrix.
  type(sparse), intent(out):: Frigid   (:)      ! Influence coefficients matrix for applied forces.
  real(8),      intent(out):: Qrigid   (:)      ! Stiffness and gyroscopic force vector.
  type(xbopts), intent(in) :: Options           ! Solver parameters.
  real(8),      intent(in) :: Cao      (:,:)    ! Rotation operator from reference to inertial frame

! Local variables.
  logical:: Flags(MaxElNod)                ! Auxiliary flags.
  integer:: i,i1                           ! Counters.
  integer:: iElem                          ! Counter on the finite elements.
  integer:: NumE                           ! Number of elements in the model.
  integer:: NumNE                          ! Number of nodes in an element.

  real(8):: CRSelem (6,6*MaxElNod)        ! Element damping matrix.
  real(8):: Felem (6,6*MaxElNod)          ! Element force influence coefficients.
  real(8):: KRSelem (6,6*MaxElNod)        ! Element tangent stiffness matrix.
  real(8):: MRSelem (6,6*MaxElNod)        ! Element mass matrix.
  real(8):: Qelem (6)                     ! Total generalized forces on the element.
  real(8):: MRRelem(6,6)                  ! Element reference-system mass matrix.
  real(8):: CRRelem(6,6)                  ! Element reference-system damping matrix.

  real(8):: rElem0(MaxElNod,6)             ! Initial Coordinates/CRV of nodes in the element.
  real(8):: rElem (MaxElNod,6)             ! Current Coordinates/CRV of nodes in the element.
  real(8):: rElemDot (MaxElNod,6)          ! Current Coordinates/CRV of nodes in the element.
  real(8):: rElemDDot (MaxElNod,6)         ! Current Coordinates/CRV of nodes in the element.
  real(8):: SB2B1 (6*MaxElNod,6*MaxElNod)  ! Transformation from master to rigid node orientations.

! Loop in all elements in the model.
  NumE=size(Elem)

  do iElem=1,NumE
    MRSelem=0.d0; CRSelem=0.d0; KRSelem=0.d0; Felem=0.d0; Qelem=0.d0
    MRRelem=0.d0; CRRelem=0.d0; SB2B1=0.d0

! Extract coords of elem nodes and determine if they are master (Flag=T) or slave.
    call fem_glob2loc_extract (Elem(iElem)%Conn,Coords,rElem0(:,1:3),NumNE)

    Flags=.false.
    do i=1,Elem(iElem)%NumNodes
      if (Node(Elem(iElem)%Conn(i))%Master(1).eq.iElem) Flags(i)=.true.
    end do

    call fem_glob2loc_extract (Elem(iElem)%Conn,PosDefor,    rElem    (:,1:3),NumNE)
    call fem_glob2loc_extract (Elem(iElem)%Conn,PosDeforDot, rElemDot (:,1:3),NumNE)
    call fem_glob2loc_extract (Elem(iElem)%Conn,PosDeforDDot,rElemDDot(:,1:3),NumNE)

    rElem0   (:,4:6)= Psi0        (iElem,:,:)
    rElem    (:,4:6)= PsiDefor    (iElem,:,:)
    rElemDot (:,4:6)= PsiDeforDot (iElem,:,:)
    rElemDDot(:,4:6)= PsiDeforDDot(iElem,:,:)

! CHECK LINEARISED FORM
!    call fdiff_check (NumNE,rElem0,rElem,rElemDot,rElemDDot,Vrel,VrelDot,Elem(iElem)%Mass,Elem(iElem)%Stiff,Options%NumGauss)

! Compute linearized inertia matrices.
    call xbeam_mrs  (NumNE,rElem0,rElem,                                Elem(iElem)%Mass,MRSelem,Options%NumGauss)
    call xbeam_cgyr (NumNE,rElem0,rElem,rElemDot,Vrel,                  Elem(iElem)%Mass,CRSelem,Options%NumGauss)
    call xbeam_kgyr (NumNE,rElem0,rElem,rElemDot,rElemDDot,Vrel,VrelDot,Elem(iElem)%Mass,KRSelem,Options%NumGauss)

! Compute the gyroscopic force vector.
    call xbeam_fgyr (NumNE,rElem0,rElem,rElemDot,Vrel,Elem(iElem)%Mass,Qelem,Options%NumGauss)

! Compute the element mass tangent stiffness matrix (can be neglected).
    call xbeam_kmass  (NumNE,rElem0,rElem,rElemDDot,VrelDot,Elem(iElem)%Mass,KRSelem,Options%NumGauss)

! Compute the element contribution to the mass and damping in the motion of the reference frame.
    call xbeam_mrr  (NumNE,rElem0,rElem              ,Elem(iElem)%Mass,MRRelem,Options%NumGauss)
    call xbeam_crr  (NumNE,rElem0,rElem,rElemDot,Vrel,Elem(iElem)%Mass,CRRelem,Options%NumGauss)

! Add contributions of non-structural (lumped) mass.
	if (any(Elem(iElem)%RBMass.ne.0.d0)) then
      call xbeam_rbmrs  (NumNE,rElem0,rElem,                                Elem(iElem)%RBMass,MRSelem)
      call xbeam_rbcgyr (NumNE,rElem0,rElem,rElemDot,          Vrel,        Elem(iElem)%RBMass,CRSelem)
      call xbeam_rbkgyr (NumNE,rElem0,rElem,rElemDot,rElemDDot,Vrel,VrelDot,Elem(iElem)%RBMass,KRSelem)
      call xbeam_rbfgyr (NumNE,rElem0,rElem,rElemDot,          Vrel,        Elem(iElem)%RBMass,Qelem)
      call xbeam_rbkmass(NumNE,rElem0,rElem,         rElemDDot,     VrelDot,Elem(iElem)%RBMass,KRSelem)
      call xbeam_rbmrr  (NumNE,rElem0,rElem,                                Elem(iElem)%RBMass,MRRelem)
      call xbeam_rbcrr  (NumNE,rElem0,rElem,rElemDot,          Vrel,        Elem(iElem)%RBMass,CRRelem)
    end if

! Project slave degrees of freedom to the orientation of the "master" ones.
    call cbeam3_projs2m (NumNE,Elem(iElem)%Master(:,:),Psi0(iElem,:,:),Psi0,SB2B1)
    MRSelem=matmul(MRSelem,SB2B1)
    CRSelem=matmul(CRSelem,SB2B1)
    KRSelem=matmul(KRSelem,SB2B1)

! Compute the influence coefficients multiplying the vector of external forces.
    call xbeam_fext  (NumNE,rElem,Flags(1:NumNE),Felem,Options%FollowerForce,Options%FollowerForceRig,Cao)

! Add to global matrix. Remove columns and rows at clamped points.
    Qrigid      = Qrigid + Qelem

    MRR  = MRR + MRRelem
    CRR  = CRR + CRRelem

    do i=1,NumNE
      i1=Node(Elem(iElem)%Conn(i))%Vdof
      call sparse_addmat (0,6*(i1),Felem(:,6*(i-1)+1:6*i),fs,Frigid)
      if (i1.ne.0) then
        call sparse_addmat (0,6*(i1-1),MRSelem(:,6*(i-1)+1:6*i),ms,MRS)
        call sparse_addmat (0,6*(i1-1),CRSelem(:,6*(i-1)+1:6*i),cs,CRS)
        call sparse_addmat (0,6*(i1-1),KRSelem(:,6*(i-1)+1:6*i),ks,KRS)
      end if
    end do

  end do

! Compute tangent matrices for quaternion equations
  CQR(:,:)=0.d0
  CQR(1,4)= Quat(2); CQR(1,5)= Quat(3); CQR(1,6)= Quat(4)
  CQR(2,4)=-Quat(1); CQR(2,5)= Quat(4); CQR(2,6)=-Quat(3)
  CQR(3,4)=-Quat(4); CQR(3,5)=-Quat(1); CQR(3,6)= Quat(2)
  CQR(4,4)=-Quat(3); CQR(4,5)=-Quat(2); CQR(4,6)=-Quat(1)
  CQR=0.5d0*CQR

  CQQ=0.5d0*xbeam_QuadSkew(Vrel(4:6))

  return
 end subroutine xbeam_asbly_dynamic


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-> Subroutine XBEAM_ASBLY_ORIENT
!
!-> Description:
!
!   Assembly rigid-body matrices to account for change in orientation. This is already done in xbeam_asbly_dynamic
!   but separated in this routine for the linear case, which still requires recomputation of force matrices and CQR and CQQ.
!
!-> Remarks.-
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine xbeam_asbly_orient (Elem,Node,PosDefor,PsiDefor,Vrel,Quat,CQR,CQQ,fs,Frigid,Options,Cao)
  use lib_rotvect
  use lib_fem
  use lib_sparse
  use lib_cbeam3
  use lib_xbeam

! I/O variables.
  type(xbelem), intent(in) :: Elem(:)               ! Element information.
  type(xbnode), intent(in) :: Node(:)               ! List of independent nodes.
  real(8),      intent(in) :: PosDefor  (:,:)       ! Current coordinates of the grid points
  real(8),      intent(in) :: PsiDefor  (:,:,:)     ! Current CRV of the nodes in the elements.
  real(8),      intent(in) :: Vrel(6)               ! Velocity of reference frame and derivative.
  real(8),      intent(in) :: Quat(4)               ! Quaternions.

  real(8),      intent(out):: CQR(:,:),CQQ(:,:) ! Tangent matrices from linearisation of quaternion equation.
  integer,      intent(out):: fs                ! Size of the sparse stiffness matrix.
  type(sparse), intent(out):: Frigid   (:)      ! Influence coefficients matrix for applied forces.
  type(xbopts), intent(in) :: Options           ! Solver parameters.
  real(8),      intent(in) :: Cao      (:,:)    ! Rotation operator from reference to inertial frame

! Local variables.
  logical:: Flags(MaxElNod)                ! Auxiliary flags.
  integer:: i,i1                           ! Counters.
  integer:: iElem                          ! Counter on the finite elements.
  integer:: NumE                           ! Number of elements in the model.
  integer:: NumNE                          ! Number of nodes in an element.

  real(8):: Felem (6,6*MaxElNod)           ! Element force influence coefficients.
  real(8):: rElem (MaxElNod,6)             ! Current Coordinates/CRV of nodes in the element.

! Loop in all elements in the model.
  NumE=size(Elem)

  do iElem=1,NumE
    Felem=0.d0;

    ! Extract coords of elem nodes and determine if they are master (Flag=T) or slave.
    Flags=.false.
    do i=1,Elem(iElem)%NumNodes
      if (Node(Elem(iElem)%Conn(i))%Master(1).eq.iElem) Flags(i)=.true.
    end do

    call fem_glob2loc_extract (Elem(iElem)%Conn,PosDefor,rElem(:,1:3),NumNE)
    rElem(:,4:6)= PsiDefor(iElem,:,:)

    ! Compute the influence coefficients multiplying the vector of external forces.
    call xbeam_fext  (NumNE,rElem,Flags(1:NumNE),Felem,Options%FollowerForce,Options%FollowerForceRig,Cao)

    ! Add to global matrix. Remove columns and rows at clamped points.
    do i=1,NumNE
      i1=Node(Elem(iElem)%Conn(i))%Vdof
      call sparse_addmat (0,6*(i1),Felem(:,6*(i-1)+1:6*i),fs,Frigid)
    end do

  end do

! Compute tangent matrices for quaternion equations
  CQR(:,:)=0.d0
  CQR(1,4)= Quat(2); CQR(1,5)= Quat(3); CQR(1,6)= Quat(4)
  CQR(2,4)=-Quat(1); CQR(2,5)= Quat(4); CQR(2,6)=-Quat(3)
  CQR(3,4)=-Quat(4); CQR(3,5)=-Quat(1); CQR(3,6)= Quat(2)
  CQR(4,4)=-Quat(3); CQR(4,5)=-Quat(2); CQR(4,6)=-Quat(1)
  CQR=0.5d0*CQR

  CQQ=0.5d0*xbeam_QuadSkew(Vrel(4:6))

  return
 end subroutine xbeam_asbly_orient


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end module xbeam_asbly
