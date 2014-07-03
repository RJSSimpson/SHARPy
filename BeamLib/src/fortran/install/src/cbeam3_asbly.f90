!-> Copyright by Imperial College London, 2009
!
!-> Module.- CBEAM3_ASBLY Rafa Palacios. 15Jul2009 - Last Update 20/09/2011
!
!-> Language: FORTRAN90, Free format.
!
!-> Description.-
!
!  Assembly beam equations.
!
!-> Subroutines.-
!
!    -cbeam3_asbly_static : Assembly matrices for nonlinear static problem.
!    -cbeam3_asbly_modal  : Assembly matrices for modal analysis.
!    -cbeam3_asbly_dynamic: Assembly matrices for nonlinear dynamic problems.
!
!-> Remarks.-
!  1) HH (20.09.2011) Changes made according to new version of nlabs r3.0,
!     which include derivatives of follower forces and new slave2master trans
!
!  2) HH (01.11.2013) Need to use full integration in assembly of mass matrix.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module cbeam3_asbly
  use xbeam_shared
  implicit none

 contains


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-> Subroutine CBEAM3_ASBLY_STATIC
!
!-> Description:
!
!   Assembly matrices for a nonlinear static problem.
!
!-> Remarks.-
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine cbeam3_asbly_static (Elem,Node,Coords,Psi0,PosDefor,PsiDefor,Force, &
&                                ks,Kglobal,fs,Fglobal,Qglobal,Options)
  use lib_rotvect
  use lib_fem
  use lib_sparse
  use lib_cbeam3

! I/O variables.
  type(xbelem), intent(in) :: Elem(:)           ! Element information.
  type(xbnode), intent(in) :: Node(:)           ! List of independent nodes.
  real(8),      intent(in) :: Coords    (:,:)   ! Initial coordinates of the grid points.
  real(8),      intent(in) :: Psi0      (:,:,:) ! Initial CRV of the nodes in the elements.
  real(8),      intent(in) :: PosDefor  (:,:)   ! Current coordinates of the grid points
  real(8),      intent(in) :: PsiDefor  (:,:,:) ! Current CRV of the nodes in the elements.
  real(8),      intent(in) :: Force     (:,:)   ! Force vector.
  integer,      intent(out):: ks                ! Size of the sparse stiffness matrix.
  type(sparse), intent(out):: Kglobal   (:)     ! Sparse stiffness matrix.
  real(8),      intent(out):: Qglobal   (:)     ! Discrete force vector.
  integer,      intent(out):: fs                ! Size of the sparse stiffness matrix.
  type(sparse), intent(out):: Fglobal   (:)     ! Influence coefficients matrix for applied forces.
  type(xbopts), intent(in) :: Options           ! Solver parameters.

! Local variables.
  logical:: Flags(MaxElNod)                ! Auxiliary flags.
  integer:: i,j,i1,j1                      ! Counters.
  integer:: iElem                          ! Counter on the finite elements.
  integer:: NumE                           ! Number of elements in the model.
  integer:: NumNE                          ! Number of nodes in an element.
  real(8):: Felem (6*MaxElNod,6*MaxElNod)  ! Element force influence coefficients.  
  real(8):: Kelem (6*MaxElNod,6*MaxElNod)  ! Element tangent stiffness matrix.  
  real(8):: Qelem (6*MaxElNod)             ! Total generalized forces on the element.
  real(8):: rElem0(MaxElNod,6)             ! Initial Coordinates/CRV of nodes in the element.
  real(8):: rElem (MaxElNod,6)             ! Current Coordinates/CRV of nodes in the element.
  real(8):: ForceElem (MaxElNod,6)         ! Current forces/moments of nodes in the element.
  real(8):: SB2B1 (6*MaxElNod,6*MaxElNod)  ! Transformation from master to global node orientations.

! Loop in all elements in the model.
  NumE=size(Elem)

  do iElem=1,NumE
    Kelem=0.d0; Felem=0.d0; Qelem=0.d0; SB2B1=0.d0

! Determine if the element nodes are master (Flag=T) or slave.
    Flags=.false.
    do i=1,Elem(iElem)%NumNodes
      if (Node(Elem(iElem)%Conn(i))%Master(1).eq.iElem) Flags(i)=.true.
    end do

! Extract components of the displacement and rotation vector at the element nodes
! and for the reference and current configurations.
    call fem_glob2loc_extract (Elem(iElem)%Conn,Coords, rElem0(:,1:3),NumNE)
    call fem_glob2loc_extract (Elem(iElem)%Conn,PosDefor,rElem(:,1:3),NumNE)

    rElem0(:,4:6)= Psi0    (iElem,:,:)
    rElem (:,4:6)= PsiDefor(iElem,:,:)
    call rotvect_boundscheck2(rElem(1,4:6),rElem(2,4:6))
    if (NumNE.eq.3) &
&   call rotvect_boundscheck2(rElem(3,4:6),rElem(2,4:6))

! Extract current applied forces/moments at the element nodes.
    call fem_glob2loc_extract (Elem(iElem)%Conn,Force,ForceElem,NumNE)

! Compute the element tangent stiffness matrix and force vectors.
    call cbeam3_kmat  (NumNE,rElem0,rElem,Elem(iElem)%Stiff,Kelem,Options%NumGauss)
    call cbeam3_kgeom (NumNE,rElem0,rElem,Elem(iElem)%Stiff,Kelem,Options%NumGauss)
    call cbeam3_fstif (NumNE,rElem0,rElem,Elem(iElem)%Stiff,Qelem,Options%NumGauss)

! Compute the influence coefficients multiplying the vector of external forces.
! (rotate if follower forces and filter out slave nodes).
    call cbeam3_fext  (NumNE,rElem,Flags(1:NumNE),Felem,Options%FollowerForce,Options%FollowerForceRig,Unit)
    call cbeam3_dqext (NumNE,rElem,ForceElem,Flags(1:NumNE),Kelem,Options%FollowerForce)

! Project equations to the orientation of the "master" degrees of freedom.
!    call cbeam3_projs2m (NumNE,Elem(iElem)%Master(:,:),Psi0(iElem,:,:),Psi0,SB2B1)
    call cbeam3_slave2master (NumNE,Elem(iElem)%Master(:,:),rElem0(:,4:6),Psi0,rElem(:,4:6),PsiDefor,SB2B1)
    Kelem=matmul(transpose(SB2B1),matmul(Kelem,SB2B1))
    Felem=matmul(transpose(SB2B1),Felem)
    Qelem=matmul(transpose(SB2B1),Qelem)

! Add to global matrix. Remove columns and rows at clamped points.
    do i=1,NumNE
      i1=Node(Elem(iElem)%Conn(i))%Vdof
      if (i1.ne.0) then

        Qglobal(6*(i1-1)+1:6*i1)= Qglobal(6*(i1-1)+1:6*i1)+Qelem(6*(i-1)+1:6*i)

        do j=1,NumNE
          j1=Node(Elem(iElem)%Conn(j))%Vdof
          if (j1.ne.0) then
            call sparse_addmat (6*(i1-1),6*(j1-1),Kelem(6*(i-1)+1:6*i,6*(j-1)+1:6*j),&
&                               ks,Kglobal)
            call sparse_addmat (6*(i1-1),6*(j1-1),Felem(6*(i-1)+1:6*i,6*(j-1)+1:6*j),&
&                               fs,Fglobal)
          end if
        end do
      end if
    end do

  end do

  return
 end subroutine cbeam3_asbly_static



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-> Subroutine CBEAM3_ASBLY_MODAL
!
!-> Description:
!
!   Assembly matrices for a modal analysis.
!
!-> Remarks.-
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine cbeam3_asbly_modal (Elem,Node,Coords,Psi0,PosDefor,PsiDefor,  &
&                               Vrel,ms,Mglobal,cs,Cglobal,ks,Kglobal,Options)
  use lib_rotvect
  use lib_fem
  use lib_sparse
  use lib_cbeam3

! I/O variables.
  type(xbelem),intent(in) :: Elem(:)               ! Element information.
  type(xbnode),intent(in) :: Node(:)               ! List of independent nodes.
  real(8),     intent(in) :: Coords    (:,:)       ! Initial coordinates of the grid points.
  real(8),     intent(in) :: Psi0      (:,:,:)     ! Initial CRV of the nodes in the elements.
  real(8),     intent(in) :: PosDefor  (:,:)       ! Current coordinates of the grid points
  real(8),     intent(in) :: PsiDefor  (:,:,:)     ! Current CRV of the nodes in the elements.
  real(8),     intent(in) :: Vrel(6)               ! Velocity of reference frame .

  integer,     intent(out):: ms                    ! Size of the sparse mass matrix.
  type(sparse),intent(out):: Mglobal   (:)         ! Sparse mass matrix.
  integer,     intent(out):: cs                    ! Size of the sparse damping matrix.
  type(sparse),intent(out):: Cglobal   (:)         ! Sparse damping matrix.
  integer,     intent(out):: ks                    ! Size of the sparse stiffness matrix.
  type(sparse),intent(out):: Kglobal   (:)         ! Sparse stiffness matrix.
  type(xbopts),intent(in) :: Options               ! Solver parameters.

! Local variables.
  integer:: i,j,i1,j1                      ! Counters.
  integer:: iElem                          ! Counter on the finite elements.
  integer:: NumE                           ! Number of elements in the model.
  integer:: NumNE                          ! Number of nodes in an element.
  integer:: NumGaussMass                   ! Number of Gaussian points in the inertia terms.
  real(8):: Celem (6*MaxElNod,6*MaxElNod)  ! Element damping matrix.
  real(8):: Kelem (6*MaxElNod,6*MaxElNod)  ! Element tangent stiffness matrix.
  real(8):: Melem (6*MaxElNod,6*MaxElNod)  ! Element mass matrix.
  real(8):: rElem0(MaxElNod,6)             ! Initial Coordinates/CRV of nodes in the element.
  real(8):: rElem (MaxElNod,6)             ! Current Coordinates/CRV of nodes in the element.
  real(8):: rElemDot (MaxElNod,6)          ! Current Coordinates/CRV of nodes in the element.
  real(8):: rElemDDot (MaxElNod,6)         ! Current Coordinates/CRV of nodes in the element.
  real(8):: VrelDot  (6)                   ! Time derivative of Vrel.
  real(8):: SB2B1 (6*MaxElNod,6*MaxElNod)  ! Transformation from master to global node orientations.

! Initialize.
  NumE=size(Elem)
  rElemDot = 0.d0
  rElemDDot= 0.d0
  VrelDot  = 0.d0

  do iElem=1,NumE
    Melem=0.d0; Celem=0.d0; Kelem=0.d0; SB2B1=0.d0

! Extract coords of elem nodes and determine if they are master (Flag=T) or slave.
    call fem_glob2loc_extract (Elem(iElem)%Conn,Coords,  rElem0(:,1:3),NumNE)
    call fem_glob2loc_extract (Elem(iElem)%Conn,PosDefor,rElem (:,1:3),NumNE)
    rElem0 (:,4:6)= Psi0     (iElem,:,:)
    rElem  (:,4:6)= PsiDefor (iElem,:,:)

! Use full integration for mass matrix.
    if (NumNE.eq.2) then
        NumGaussMass=NumNE
    elseif (NumNE.eq.3) then
        NumGaussMass=NumNE+1
    end if

! Compute linearized inertia matrices.
    call cbeam3_mass (NumNE,rElem0,rElem,                                Elem(iElem)%Mass,Melem,NumGaussMass)
    call cbeam3_cgyr (NumNE,rElem0,rElem,rElemDot,Vrel,                  Elem(iElem)%Mass,Celem,Options%NumGauss)
    call cbeam3_kgyr (NumNE,rElem0,rElem,rElemDot,rElemDDot,Vrel,VrelDot,Elem(iElem)%Mass,Kelem,Options%NumGauss)

    ! Add contributions of non-structural (lumped) mass.
    if (any(Elem(iElem)%RBMass.ne.0.d0)) then
      call cbeam3_rbmass (NumNE,rElem0,rElem,                                Elem(iElem)%RBMass,Melem)
      call cbeam3_rbcgyr (NumNE,rElem0,rElem,rElemDot,          Vrel,        Elem(iElem)%RBMass,Celem)
      call cbeam3_rbkgyr (NumNE,rElem0,rElem,rElemDot,rElemDDot,Vrel,VrelDot,Elem(iElem)%RBMass,Kelem)
    end if

! Compute the element tangent stiffness matrix and force vectors.
    call cbeam3_kmat  (NumNE,rElem0,rElem,Elem(iElem)%Stiff,Kelem,Options%NumGauss)
    call cbeam3_kgeom (NumNE,rElem0,rElem,Elem(iElem)%Stiff,Kelem,Options%NumGauss)

! Project slave degrees of freedom to the orientation of the "master" ones.
    call cbeam3_projs2m (NumNE,Elem(iElem)%Master(:,:),Psi0(iElem,:,:),Psi0,SB2B1)
    Melem=matmul(transpose(SB2B1),matmul(Melem,SB2B1))
    Celem=matmul(transpose(SB2B1),matmul(Celem,SB2B1))
    Kelem=matmul(transpose(SB2B1),matmul(Kelem,SB2B1))

! Add to global matrix. Remove columns and rows at clamped points.
    do i=1,NumNE
      i1=Node(Elem(iElem)%Conn(i))%Vdof
      if (i1.ne.0) then
        do j=1,NumNE
          j1=Node(Elem(iElem)%Conn(j))%Vdof
          if (j1.ne.0) then
            call sparse_addmat (6*(i1-1),6*(j1-1),Melem(6*(i-1)+1:6*i,6*(j-1)+1:6*j),ms,Mglobal)
            call sparse_addmat (6*(i1-1),6*(j1-1),Celem(6*(i-1)+1:6*i,6*(j-1)+1:6*j),cs,Cglobal)
            call sparse_addmat (6*(i1-1),6*(j1-1),Kelem(6*(i-1)+1:6*i,6*(j-1)+1:6*j),ks,Kglobal)
          end if
        end do
      end if
    end do
  end do

  return
 end subroutine cbeam3_asbly_modal


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-> Subroutine CBEAM3_ASBLY_DYNAMIC
!
!-> Description:
!
!   Assembly matrices for a dynamic problem.
!
!-> Remarks.-
!
!  1) Acceleration is used to compute the linearized stiffness matrix of the
!     inertia terms, although this can often be neglected.
!  2) RPN (20.04.2011) added lumped masses to the model.
!  3) RPN (20.04.2011) The number of Gauss points has been increased in all the
!     inertia matrices so that they are consistent among them.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine cbeam3_asbly_dynamic (Elem,Node,Coords,Psi0,PosDefor,PsiDefor,         &
&                                 PosDeforDot,PsiDeforDot,PosDeforDDot,            &
&                                 PsiDeforDDot,Force,Vrel,VrelDot,ms,Mglobal,Mvel, &
&                                 cs,Cglobal,Cvel,ks,Kglobal,fs,Fglobal,           &
&                                 Qglobal,Options,Cao)
  use lib_rotvect
  use lib_fem
  use lib_sparse
  use lib_cbeam3

! I/O variables.
  type(xbelem),intent(in) :: Elem(:)               ! Element information.
  type(xbnode),intent(in) :: Node(:)               ! List of independent nodes.
  real(8),     intent(in) :: Coords    (:,:)       ! Initial coordinates of the grid points.
  real(8),     intent(in) :: Psi0      (:,:,:)     ! Initial CRV of the nodes in the elements.
  real(8),     intent(in) :: PosDefor  (:,:)       ! Current coordinates of the grid points
  real(8),     intent(in) :: PsiDefor  (:,:,:)     ! Current CRV of the nodes in the elements.
  real(8),     intent(in) :: PosDeforDot  (:,:)    ! Current coordinates of the grid points
  real(8),     intent(in) :: PsiDeforDot  (:,:,:)  ! Current CRV of the nodes in the elements.
  real(8),     intent(in) :: PosDeforDDot (:,:)    ! Current coordinates of the grid points
  real(8),     intent(in) :: PsiDeforDDot (:,:,:)  ! Current CRV of the nodes in the elements.
  real(8),     intent(in) :: Force     (:,:)       ! Force vector.
  real(8),     intent(in) :: Vrel(6), VrelDot(6)   ! Velocity of reference frame and derivative.
  
  integer,     intent(out):: ms                ! Size of the sparse mass matrix.
  type(sparse),intent(out):: Mglobal   (:)     ! Sparse mass matrix.
  real(8),     intent(out):: Mvel      (:,:)   ! Reference system mass matrix.
  integer,     intent(out):: cs                ! Size of the sparse damping matrix.
  type(sparse),intent(out):: Cglobal   (:)     ! Sparse damping matrix.
  real(8),     intent(out):: Cvel      (:,:)   ! Reference system damping matrix.
  integer,     intent(out):: ks                ! Size of the sparse stiffness matrix.
  type(sparse),intent(out):: Kglobal   (:)     ! Sparse stiffness matrix.
  integer,     intent(out):: fs                ! Size of the sparse stiffness matrix.
  type(sparse),intent(out):: Fglobal   (:)     ! Influence coefficients matrix for applied forces.
  real(8),     intent(out):: Qglobal   (:)     ! Stiffness and gyroscopic force vector.
  type(xbopts),intent(in) :: Options           ! Solver parameters.
  real(8),     intent(in) :: Cao       (:,:)   ! Rotation operator

! Local variables.
  logical:: Flags(MaxElNod)                ! Auxiliary flags.
  integer:: i,j,i1,j1                      ! Counters.
  integer:: iElem                          ! Counter on the finite elements.
  integer:: NumE                           ! Number of elements in the model.
  integer:: NumNE                          ! Number of nodes in an element.
  integer:: NumGaussMass                   ! Number of Gaussian points in the inertia terms.
  real(8):: Celem (6*MaxElNod,6*MaxElNod)  ! Element damping matrix.
  real(8):: Felem (6*MaxElNod,6*MaxElNod)  ! Element force influence coefficients.
  real(8):: Kelem (6*MaxElNod,6*MaxElNod)  ! Element tangent stiffness matrix.
  real(8):: Melem (6*MaxElNod,6*MaxElNod)  ! Element mass matrix.
  real(8):: Qelem (6*MaxElNod)             ! Total generalized forces on the element.
  real(8):: Mvelelem(6*MaxElNod,6)         ! Element reference-system mass matrix.
  real(8):: Cvelelem(6*MaxElNod,6)         ! Element reference-system damping matrix.
  real(8):: rElem0(MaxElNod,6)             ! Initial Coordinates/CRV of nodes in the element.
  real(8):: rElem (MaxElNod,6)             ! Current Coordinates/CRV of nodes in the element.
  real(8):: rElemDot (MaxElNod,6)          ! Current Coordinates/CRV of nodes in the element.
  real(8):: rElemDDot (MaxElNod,6)         ! Current Coordinates/CRV of nodes in the element.
  real(8):: ForceElem (MaxElNod,6)         ! Current forces/moments of nodes in the element.
  real(8):: SB2B1 (6*MaxElNod,6*MaxElNod)  ! Transformation from master to global node orientations.

! Loop in all elements in the model.
  NumE=size(Elem)

  do iElem=1,NumE
    Melem=0.d0; Celem=0.d0; Kelem=0.d0; Felem=0.d0; Qelem=0.d0
    Mvelelem=0.d0; Cvelelem=0.d0; SB2B1=0.d0

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

    call rotvect_boundscheck2(rElem(1,4:6),rElem(2,4:6))
    if (NumNE.eq.3) &
&   call rotvect_boundscheck2(rElem(3,4:6),rElem(2,4:6))

! Extract current applied forces/moments at the element nodes.
    call fem_glob2loc_extract (Elem(iElem)%Conn,Force,ForceElem,NumNE)

! Use full integration for mass matrix.
    if (NumNE.eq.2) then
        NumGaussMass=NumNE
    elseif (NumNE.eq.3) then
        NumGaussMass=NumNE+1
    end if

! Compute the element contribution to the mass and damping in the motion of the reference frame.
    call cbeam3_mvel (NumNE,rElem0,rElem,              Elem(iElem)%Mass,Mvelelem,NumGaussMass)
    call cbeam3_cvel (NumNE,rElem0,rElem,rElemDot,Vrel,Elem(iElem)%Mass,Cvelelem,Options%NumGauss)

! Contributions of the structural mass to the linearized inertia matrices.
    call cbeam3_mass (NumNE,rElem0,rElem,                                Elem(iElem)%Mass,Melem,NumGaussMass)
    call cbeam3_cgyr (NumNE,rElem0,rElem,rElemDot,          Vrel,        Elem(iElem)%Mass,Celem,Options%NumGauss)
    call cbeam3_kgyr (NumNE,rElem0,rElem,rElemDot,rElemDDot,Vrel,VrelDot,Elem(iElem)%Mass,Kelem,Options%NumGauss)

! Compute the gyroscopic force vector.
    call cbeam3_fgyr (NumNE,rElem0,rElem,rElemDot,Vrel,Elem(iElem)%Mass,Qelem,Options%NumGauss)

! Add contributions of non-structural (lumped) mass.
    if (any(Elem(iElem)%RBMass.ne.0.d0)) then
      call cbeam3_rbmvel (NumNE,rElem0,rElem,                                Elem(iElem)%RBMass,Mvelelem)
      call cbeam3_rbcvel (NumNE,rElem0,rElem,rElemDot,          Vrel,        Elem(iElem)%RBMass,Cvelelem)
      call cbeam3_rbmass (NumNE,rElem0,rElem,                                Elem(iElem)%RBMass,Melem)
      call cbeam3_rbcgyr (NumNE,rElem0,rElem,rElemDot,          Vrel,        Elem(iElem)%RBMass,Celem)
      call cbeam3_rbkgyr (NumNE,rElem0,rElem,rElemDot,rElemDDot,Vrel,VrelDot,Elem(iElem)%RBMass,Kelem)
      call cbeam3_rbfgyr (NumNE,rElem0,rElem,rElemDot,          Vrel,        Elem(iElem)%RBMass,Qelem)
    end if

! Compute the element tangent stiffness matrix and force vectors.
    call cbeam3_kmat  (NumNE,rElem0,rElem,Elem(iElem)%Stiff,Kelem,Options%NumGauss)
    call cbeam3_kgeom (NumNE,rElem0,rElem,Elem(iElem)%Stiff,Kelem,Options%NumGauss)
    call cbeam3_fstif (NumNE,rElem0,rElem,Elem(iElem)%Stiff,Qelem,Options%NumGauss)

! Add external forces to the residual, and their derivatives with respect to the nodal
! degrees of freedom to the stiffness.
    if (any(ForceElem.ne.0.d0)) then
      call cbeam3_fext (NumNE,rElem,Flags(1:NumNE),Felem,Options%FollowerForce,Options%FollowerForceRig,Cao)
      call cbeam3_dqext(NumNE,rElem,ForceElem,Flags(1:NumNE),Kelem,Options%FollowerForce)
    end if

! Project slave degrees of freedom to the orientation of the "master" ones.
!    call cbeam3_projs2m (NumNE,Elem(iElem)%Master(:,:),Psi0(iElem,:,:),Psi0,SB2B1)
    call cbeam3_slave2master (NumNE,Elem(iElem)%Master(:,:),rElem0(:,4:6),Psi0,rElem(:,4:6),PsiDefor,SB2B1)
    Melem=matmul(transpose(SB2B1),matmul(Melem,SB2B1))
    Celem=matmul(transpose(SB2B1),matmul(Celem,SB2B1))
    Kelem=matmul(transpose(SB2B1),matmul(Kelem,SB2B1))
    Qelem=matmul(transpose(SB2B1),Qelem)
    Felem=matmul(transpose(SB2B1),Felem)
    Mvelelem=matmul(transpose(SB2B1),Mvelelem)
    Cvelelem=matmul(transpose(SB2B1),Cvelelem)

! Add to global matrix. Remove columns and rows at clamped points.
    do i=1,NumNE
      i1=Node(Elem(iElem)%Conn(i))%Vdof
      if (i1.ne.0) then

        Qglobal(6*(i1-1)+1:6*i1)   = Qglobal(6*(i1-1)+1:6*i1)    + Qelem   (6*(i-1)+1:6*i)
        Mvel   (6*(i1-1)+1:6*i1,:) = Mvel   (6*(i1-1)+1:6*i1,:)  + Mvelelem(6*(i-1)+1:6*i,:)
        Cvel   (6*(i1-1)+1:6*i1,:) = Cvel   (6*(i1-1)+1:6*i1,:)  + Cvelelem(6*(i-1)+1:6*i,:)

        do j=1,NumNE
          j1=Node(Elem(iElem)%Conn(j))%Vdof
          if (j1.ne.0) then
            call sparse_addmat (6*(i1-1),6*(j1-1),Melem(6*(i-1)+1:6*i,6*(j-1)+1:6*j),ms,Mglobal)
            call sparse_addmat (6*(i1-1),6*(j1-1),Celem(6*(i-1)+1:6*i,6*(j-1)+1:6*j),cs,Cglobal)
            call sparse_addmat (6*(i1-1),6*(j1-1),Kelem(6*(i-1)+1:6*i,6*(j-1)+1:6*j),ks,Kglobal)
            call sparse_addmat (6*(i1-1),6*(j1-1),Felem(6*(i-1)+1:6*i,6*(j-1)+1:6*j),fs,Fglobal)
          end if
        end do

      end if
    end do
  end do

  return
 end subroutine cbeam3_asbly_dynamic
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-> Subroutine CBEAM3_ASBLY_FGLOBAL
!
!-> Description:
!
!   Separate assembly of influence coefficients matrix for 
!   applied follower and dead loads
!
!-> Remarks.-
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine cbeam3_asbly_Fglobal (Elem,Node,Coords,Psi0,PosDefor,PsiDefor,Force,           &
&                                 ksf,Kglobal_foll,fsf,Fglobal_foll,fsd,Fglobal_dead,CAG)
  use lib_rotvect
  use lib_fem
  use lib_sparse
  use lib_cbeam3

! I/O variables.
  type(xbelem), intent(in) :: Elem(:)           ! Element information.
  type(xbnode), intent(in) :: Node(:)           ! List of independent nodes.
  real(8),      intent(in) :: Coords    (:,:)   ! Initial coordinates of the grid points.
  real(8),      intent(in) :: Psi0      (:,:,:) ! Initial CRV of the nodes in the elements.
  real(8),      intent(in) :: PosDefor  (:,:)   ! Current coordinates of the grid points
  real(8),      intent(in) :: PsiDefor  (:,:,:) ! Current CRV of the nodes in the elements.
  real(8),      intent(in) :: Force     (:,:)   ! Force vector.
  integer,      intent(out):: ksf               ! Size of the sparse stiffness matrix.
  type(sparse), intent(out):: Kglobal_foll (:)  ! Sparse stiffness matrix.
  integer,      intent(out):: fsf               ! Size of the sparse force matrix, Fglobal_foll.
  type(sparse), intent(out):: Fglobal_foll (:)  ! Influence coefficients matrix for follower forces.
  integer,      intent(out):: fsd               ! Size of the sparse force matrix, Fglobal_dead.
  type(sparse), intent(out):: Fglobal_dead (:)  ! Influence coefficients matrix for dead forces.
  real(8),      intent(in) :: CAG       (:,:)   ! Rotation operator

! Local variables.
  logical:: Flags(MaxElNod)                     ! Auxiliary flags.
  integer:: i,j,i1,j1                           ! Counters.
  integer:: iElem                               ! Counter on the finite elements.
  integer:: NumE                                ! Number of elements in the model.
  integer:: NumNE                               ! Number of nodes in an element.
  real(8):: Kelem_foll (6*MaxElNod,6*MaxElNod)  ! Element tangent stiffness matrix.
  real(8):: Felem_foll (6*MaxElNod,6*MaxElNod)  ! Element force influence coefficients.
  real(8):: Felem_dead (6*MaxElNod,6*MaxElNod)  ! Element force influence coefficients.
  real(8):: rElem0(MaxElNod,6)                  ! Initial Coordinates/CRV of nodes in the element.
  real(8):: rElem (MaxElNod,6)                  ! Current Coordinates/CRV of nodes in the element.
  real(8):: ForceElem (MaxElNod,6)              ! Current forces/moments of nodes in the element.
  real(8):: SB2B1 (6*MaxElNod,6*MaxElNod)       ! Transformation from master to global node orientations.

! Initialise
  call sparse_zero(ksf,Kglobal_foll)
  call sparse_zero(fsf,Fglobal_foll)
  call sparse_zero(fsd,Fglobal_dead)
  
! Loop in all elements in the model.
  NumE=size(Elem)

  do iElem=1,NumE
    Kelem_foll=0.d0; Felem_foll=0.d0; Felem_dead=0.d0;

    ! Determine if the element nodes are master (Flag=T) or slave.
    Flags=.false.
    do i=1,Elem(iElem)%NumNodes
      if (Node(Elem(iElem)%Conn(i))%Master(1).eq.iElem) Flags(i)=.true.
    end do

    ! Extract components of the displacement and rotation vector at the element nodes
    ! and for the reference and current configurations.
    call fem_glob2loc_extract (Elem(iElem)%Conn,Coords, rElem0(:,1:3),NumNE)
    call fem_glob2loc_extract (Elem(iElem)%Conn,PosDefor,rElem(:,1:3),NumNE)

    rElem0(:,4:6)= Psi0    (iElem,:,:)
    rElem (:,4:6)= PsiDefor(iElem,:,:)
    call rotvect_boundscheck2(rElem(1,4:6),rElem(2,4:6))
    if (NumNE.eq.3) &
&   call rotvect_boundscheck2(rElem(3,4:6),rElem(2,4:6))

    ! Compute the influence coefficients multiplying the vector of external forces.
    call cbeam3_fext (NumNE,rElem,Flags(1:NumNE),Felem_foll,.true.,.true.,CAG)
    call cbeam3_fext (NumNE,rElem,Flags(1:NumNE),Felem_dead,.false.,.false.,CAG)
    
    ! Compute contribution to Kglobal from follower forces
    call fem_glob2loc_extract (Elem(iElem)%Conn,Force,ForceElem,NumNE)
    if (any(ForceElem.ne.0.d0)) call cbeam3_dqext (NumNE,rElem,ForceElem,Flags(1:NumNE),Kelem_foll,.true.)
    
    ! Project equations to the orientation of the "master" degrees of freedom.
    call cbeam3_slave2master (NumNE,Elem(iElem)%Master(:,:),rElem0(:,4:6),Psi0,rElem(:,4:6),PsiDefor,SB2B1)
    Kelem_foll=matmul(transpose(SB2B1),matmul(Kelem_foll,SB2B1))
    Felem_foll=matmul(transpose(SB2B1),Felem_foll)
    Felem_dead=matmul(transpose(SB2B1),Felem_dead)

    ! Add to global matrix. Remove columns and rows at clamped points.
    do i=1,NumNE
      i1=Node(Elem(iElem)%Conn(i))%Vdof
      if (i1.ne.0) then
        do j=1,NumNE
          j1=Node(Elem(iElem)%Conn(j))%Vdof
          if (j1.ne.0) then
            call sparse_addmat (6*(i1-1),6*(j1-1),Kelem_foll(6*(i-1)+1:6*i,6*(j-1)+1:6*j),&
&                               ksf,Kglobal_foll)
            call sparse_addmat (6*(i1-1),6*(j1-1),Felem_foll(6*(i-1)+1:6*i,6*(j-1)+1:6*j),&
&                               fsf,Fglobal_foll)
            call sparse_addmat (6*(i1-1),6*(j1-1),Felem_dead(6*(i-1)+1:6*i,6*(j-1)+1:6*j),&
&                               fsd,Fglobal_dead)
          end if
        end do
      end if
    end do

  end do

  return
 end subroutine cbeam3_asbly_Fglobal
 
end module cbeam3_asbly
