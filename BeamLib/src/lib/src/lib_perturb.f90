!-> Copyright by Imperial College London, 2011
!
!-> Module.- LIB_CBEAM3 Rafa Palacios. 15Jul2008 - Last Update 20Apr2011
!
!-> Language: FORTRAN90, Free format.
!
!-> Description.-
!
!  Compute perturbed matrices for the 3-noded geometrically-exact beam element. This is the contribution coming from the perturbed
!  gyroscopic forces to the elasically linearised system with large nonlinear rigid-body motion.
!
!-> Subroutines.-
!
!    -perturb_css:  Compute tangent perturbed gyroscopic damping matrix for structural dof.
!    -perturb_csr:  Compute tangent perturbed gyroscopic damping matrix for rigid-body contribution to structural dynamics.
!    -perturb_crs:  Compute tangent perturbed gyroscopic damping matrix for structural contribution to rigid-body motion.
!    -perturb_crr:  Compute tangent perturbed gyroscopic damping matrix for rigid-body dof.
!    -perturb_kss:  Compute tangent gyroscopic stiffness matrix for structural dof.
!    -perturb_krs:  Compute tangent gyroscopic stiffness matrix for structural contribution to rigid-body dynamics.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module lib_perturb
  implicit none

! Shared variables within the module.
 integer,private,parameter:: MaxNodCB3=3               ! Max number of nodes per element is 3.
 real(8),private,parameter,dimension(3,3):: Unit= &    ! Unit matrix.
&         reshape((/1.d0,0.d0,0.d0,0.d0,1.d0,0.d0,0.d0,0.d0,1.d0/),(/3,3/))

 contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-> Subroutine PERTURB_CSS
!
!-> Description:
!
!    Compute tangent perturbed gyroscopic damping matrix for structural dof
!
!-> Remarks.-
!
!   1) New values are added to those already in the matrix.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine perturb_css (NumNodesElem,r0,Ri,Vrel,ElemMass,CSSpert,NumGauss)
  use lib_fem
  use lib_rot
  use lib_rotvect

! I/O Variables.
  integer,intent(in)    :: NumNodesElem     ! Number of nodes in the element.
  real(8),intent(in)    :: r0       (:,:)   ! Initial position/orientation of grid points.
  real(8),intent(in)    :: Ri       (:,:)   ! Current position/orientation of grid points.
  real(8),intent(in)    :: Vrel     (6)     ! Linear/angular velocity of reference frame.
  real(8),intent(in)    :: ElemMass (6,6)   ! Inertia properties in the element.
  real(8),intent(inout) :: CSSpert  (:,:)   ! Tangent gyroscopic damping matrix.
  integer,intent(in)    :: NumGauss         ! Number of Gauss points in the element.

! Local variables
  integer :: i,j                    ! Counters.
  integer :: iGauss                 ! Counter on the Gaussian points.
  real(8) :: Jacobian               ! ds/deta, with eta the nondimensional arclength.
  real(8) :: ShapeFun(MaxNodCB3)    ! Shape functions in a gauss point.
  real(8) :: ShapeDer(MaxNodCB3)    ! Derivatives of ShapeFun in a gauss point.
  real(8),allocatable :: CoordGauss (:) ! Coords of Gauss points.
  real(8),allocatable :: WeightGauss(:) ! Coords of Gauss points.

  real(8) :: dr0_g (6)              ! Spacial derivative of r0 at Gauss point.
  real(8) :: Ra(3)                  ! Displacement at Gauss point.
  real(8) :: Psi(3)                 ! Rotation vector at Gauss point.
  real(8) :: VRB(3),WRB(3)          ! Linear and angular rigid-body velocities in the material frame.
  real(8) :: PRB(3),HRB(3)          ! Linear and angular rigid-body momenta in the material frame.

  real(8) :: CBa(3,3)               ! Coordinate transformation matrix.
  real(8) :: Rot(3,3)               ! Tangential operator.
  real(8) :: Lambda(6,6)            ! Lambda operator.
  real(8) :: ARVW (6,6)             ! A^R_{Omega Omega} matrix.
  real(8) :: ARdeltaWMass (6,6)     ! A^R_{Delta Omega}*ElemMass matrix.
  real(8) :: N (6,6*MaxNodCB3)      ! Element shape function matrix.

! Define Gauss points and loop on them.
  allocate (CoordGauss(NumGauss))
  allocate (WeightGauss(NumGauss))
  call fem_1d_gauss_val (NumGauss,CoordGauss,WeightGauss)

  do iGauss=1,NumGauss

! Obtain the shape functions and their derivatives.
    call fem_1d_shapefun (NumNodesElem,CoordGauss(iGauss),ShapeFun,ShapeDer)

! Global coordinates and Jacobian at current Gauss point.
    do i=1,6
      dr0_g(i)= dot_product (ShapeDer(1:NumNodesElem),r0(1:NumNodesElem,i))
    end do
    Jacobian =sqrt(dot_product (dr0_g(1:3),dr0_g(1:3)))

! Compute the current position vector and rotations, and their derivatives.
    do i=1,3
      Ra    (i)=dot_product (ShapeFun(1:NumNodesElem),Ri   (1:NumNodesElem,i))
      Psi   (i)=dot_product (ShapeFun(1:NumNodesElem),Ri   (1:NumNodesElem,i+3))
    end do

! Compute element shape functions.
    N =0.d0
    do i=1,6
      do j=1,NumNodesElem
        N (i,i+(j-1)*6)= ShapeFun(j)
      end do
    end do

! Compute the current coordinate transformation matrix and rotational operator.
    CBa= rotvect_psi2mat  (Psi)
    Rot= rotvect_psi2rot  (Psi)

! Compute the current linear and angular velocity.
    VRB= matmul(CBa,rot_cross(Vrel(4:6),Ra)+Vrel(1:3))
    WRB= matmul(CBa,Vrel(4:6))

    PRB= matmul(ElemMass(1:3,1:3),VRB)+matmul(ElemMass(1:3,4:6),WRB)
    HRB= matmul(ElemMass(4:6,1:3),VRB)+matmul(ElemMass(4:6,4:6),WRB)

!Operators for tangent matrix.
    Lambda=0.0d0
    Lambda(1:3,1:3)=CBa
    Lambda(4:6,4:6)=Rot

    ARVW(1:3,1:3)=rot_skew(WRB)
    ARVW(1:3,4:6)=0.d0
    ARVW(4:6,1:3)=rot_skew(VRB)
    ARVW(4:6,4:6)=rot_skew(WRB)

    ARdeltaWMass=matmul(ARVW,ElemMass)

    ARdeltaWMass(1:3,4:6)= ARdeltaWMass(1:3,4:6) - rot_skew(PRB)
    ARdeltaWMass(4:6,1:3)= ARdeltaWMass(4:6,1:3) - rot_skew(PRB)
    ARdeltaWMass(4:6,4:6)= ARdeltaWMass(4:6,4:6) - rot_skew(HRB)

! Compute perturbed elastic gyroscopic damping matrix.
    CSSpert=CSSpert + WeightGauss(iGauss)*Jacobian*    &
&                     matmul(transpose(matmul(Lambda,N)),matmul(ARdeltaWMass-matmul(ElemMass,transpose(ARVW)),matmul(Lambda,N)))

  end do
  deallocate (CoordGauss,WeightGauss)
  return
 end subroutine perturb_css


subroutine perturb_rbcss (NumNodesElem,r0,Ri,Vrel,NodalMass,CSSpert)
  use lib_fem
  use lib_rot
  use lib_rotvect

! I/O Variables.
  integer,intent(in)    :: NumNodesElem     ! Number of nodes in the element.
  real(8),intent(in)    :: r0       (:,:)   ! Initial position/orientation of grid points.
  real(8),intent(in)    :: Ri       (:,:)   ! Current position/orientation of grid points.
  real(8),intent(in)    :: Vrel     (6)     ! Linear/angular velocity of reference frame.
  real(8),intent(in)    :: NodalMass(:,:,:) ! Inertia tensor of lumped masses at element nodes.
  real(8),intent(inout) :: CSSpert  (:,:)   ! Tangent gyroscopic damping matrix.

! Local variables
  integer :: i                      ! Counters.
  integer :: iNode                  ! Counter on the element nodes.

  real(8) :: Ra(3)                  ! Displacement at Gauss point.
  real(8) :: Psi(3)                 ! Rotation vector at Gauss point.
  real(8) :: VRB(3),WRB(3)          ! Linear and angular rigid-body velocities in the material frame.
  real(8) :: PRB(3),HRB(3)          ! Linear and angular rigid-body momenta in the material frame.

  real(8) :: CBa(3,3)               ! Coordinate transformation matrix.
  real(8) :: Rot(3,3)               ! Tangential operator.
  real(8) :: Lambda(6,6)            ! Lambda operator.
  real(8) :: ARVW (6,6)             ! A^R_{Omega Omega} matrix.
  real(8) :: ARdeltaWMass (6,6)     ! A^R_{Delta Omega}*ElemMass matrix.
  real(8) :: N (6,6*MaxNodCB3)      ! Element shape function matrix.

! Loop in the element nodes.
  do iNode=1,NumNodesElem

! Compute the current position vector and rotations, and their derivatives.
	Ra    = Ri   (iNode,1:3)
	Psi   = Ri   (iNode,4:6)

! Compute element shape functions at the corrent node.
    N =0.d0
    do i=1,6
      N (i,i+(iNode-1)*6)= 1.d0
    end do

! Compute the current coordinate transformation matrix and rotational operator.
    CBa= rotvect_psi2mat  (Psi)
    Rot= rotvect_psi2rot  (Psi)

! Compute the current linear and angular velocity.
    VRB= matmul(CBa,rot_cross(Vrel(4:6),Ra)+Vrel(1:3))
    WRB= matmul(CBa,Vrel(4:6))

    PRB= matmul(NodalMass(iNode,1:3,1:3),VRB)+matmul(NodalMass(iNode,1:3,4:6),WRB)
    HRB= matmul(NodalMass(iNode,4:6,1:3),VRB)+matmul(NodalMass(iNode,4:6,4:6),WRB)

!Operators for tangent matrix.
    Lambda=0.0d0
    Lambda(1:3,1:3)=CBa
    Lambda(4:6,4:6)=Rot

    ARVW(1:3,1:3)=rot_skew(WRB)
    ARVW(1:3,4:6)=0.d0
    ARVW(4:6,1:3)=rot_skew(VRB)
    ARVW(4:6,4:6)=rot_skew(WRB)

    ARdeltaWMass=matmul(ARVW,NodalMass(iNode,:,:))

    ARdeltaWMass(1:3,4:6)= ARdeltaWMass(1:3,4:6) - rot_skew(PRB)
    ARdeltaWMass(4:6,1:3)= ARdeltaWMass(4:6,1:3) - rot_skew(PRB)
    ARdeltaWMass(4:6,4:6)= ARdeltaWMass(4:6,4:6) - rot_skew(HRB)

! Compute perturbed elastic gyroscopic damping matrix.
    CSSpert=CSSpert + matmul(transpose(matmul(Lambda,N)), &
&                     matmul(ARdeltaWMass-matmul(NodalMass(iNode,:,:),transpose(ARVW)),matmul(Lambda,N)))

  end do
  return
 end subroutine perturb_rbcss

 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-> Subroutine PERTURB_CSR
!
!-> Description:
!
!    Compute element tangent perturbed gyroscopic damping matrix for rigid-body contribution to structural dynamics
!
!-> Remarks.-
!
!   1) New values are added to those already in the matrix.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine perturb_csr (NumNodesElem,r0,Ri,Vrel,ElemMass,CSRpert,NumGauss)
  use lib_fem
  use lib_rot
  use lib_rotvect

! I/O Variables.
  integer,intent(in)    :: NumNodesElem     ! Number of nodes in the element.
  real(8),intent(in)    :: r0       (:,:)   ! Initial position/orientation of grid points.
  real(8),intent(in)    :: Ri       (:,:)   ! Current position/orientation of grid points.
  real(8),intent(in)    :: Vrel     (6)     ! Linear/angular velocity of reference frame.
  real(8),intent(in)    :: ElemMass (6,6)   ! Inertia properties in the element.
  real(8),intent(inout) :: CSRpert  (:,:)   ! Tangent gyroscopic damping matrix.
  integer,intent(in)    :: NumGauss         ! Number of Gauss points in the element.

! Local variables
  integer :: i,j                    ! Counters.
  integer :: iGauss                 ! Counter on the Gaussian points.
  real(8) :: Jacobian               ! ds/deta, with eta the nondimensional arclength.
  real(8) :: ShapeFun(MaxNodCB3)    ! Shape functions in a gauss point.
  real(8) :: ShapeDer(MaxNodCB3)    ! Derivatives of ShapeFun in a gauss point.
  real(8),allocatable :: CoordGauss (:) ! Coords of Gauss points.
  real(8),allocatable :: WeightGauss(:) ! Coords of Gauss points.

  real(8) :: dr0_g (6)              ! Spacial derivative of r0 at Gauss point.
  real(8) :: Ra(3)                  ! Displacement at Gauss point.
  real(8) :: Psi(3)                 ! Rotation vector at Gauss point.
  real(8) :: VRB(3),WRB(3)          ! Linear and angular rigid-body velocities in the material frame.

  real(8) :: CBa(3,3)               ! Coordinate transformation matrix.
  real(8) :: Rot(3,3)               ! Tangential operator.
  real(8) :: Lambda(6,6)            ! Lambda operator.
  real(8) :: ARC(6,6),ARVW(6,6)     ! A_{RC} and A^R_{Omega Omega} matrices.
  real(8) :: N (6,6*MaxNodCB3)      ! Element shape function matrix.

! Define Gauss points and loop on them.
  allocate (CoordGauss(NumGauss))
  allocate (WeightGauss(NumGauss))
  call fem_1d_gauss_val (NumGauss,CoordGauss,WeightGauss)

  do iGauss=1,NumGauss

! Obtain the shape functions and their derivatives.
    call fem_1d_shapefun (NumNodesElem,CoordGauss(iGauss),ShapeFun,ShapeDer)

! Global coordinates and Jacobian at current Gauss point.
    do i=1,6
      dr0_g(i)= dot_product (ShapeDer(1:NumNodesElem),r0(1:NumNodesElem,i))
    end do
    Jacobian =sqrt(dot_product (dr0_g(1:3),dr0_g(1:3)))

! Compute the current position vector and rotations, and their derivatives.
    do i=1,3
      Ra    (i)=dot_product (ShapeFun(1:NumNodesElem),Ri   (1:NumNodesElem,i))
      Psi   (i)=dot_product (ShapeFun(1:NumNodesElem),Ri   (1:NumNodesElem,i+3))
    end do

! Compute element shape functions.
    N =0.d0
    do i=1,6
      do j=1,NumNodesElem
        N (i,i+(j-1)*6)= ShapeFun(j)
      end do
    end do

! Compute the current coordinate transformation matrix and rotational operator.
    CBa= rotvect_psi2mat  (Psi)
    Rot= rotvect_psi2rot  (Psi)

! Compute the current linear and angular velocity.
    VRB= matmul(CBa,rot_cross(Vrel(4:6),Ra)+Vrel(1:3))
    WRB= matmul(CBa,Vrel(4:6))

!Operators for tangent matrix.
    Lambda=0.0d0
    Lambda(1:3,1:3)=CBa
    Lambda(4:6,4:6)=Rot

    ARC=0.d0
    ARC(1:3,1:3)=transpose(CBa)
    ARC(4:6,1:3)=matmul(rot_skew(Ra),transpose(CBa))
    ARC(4:6,4:6)=transpose(CBa)

    ARVW(1:3,1:3)=rot_skew(WRB)
    ARVW(1:3,4:6)=0.d0
    ARVW(4:6,1:3)=rot_skew(VRB)
    ARVW(4:6,4:6)=rot_skew(WRB)

! Compute perturbed elastic gyroscopic damping matrix.
    CSRpert=CSRpert + WeightGauss(iGauss)*Jacobian*    &
&                     matmul(transpose(matmul(Lambda,N)),matmul(matmul(ARVW,ElemMass),transpose(ARC)))

  end do
  deallocate (CoordGauss,WeightGauss)
  return
 end subroutine perturb_csr


 subroutine perturb_rbcsr (NumNodesElem,r0,Ri,Vrel,NodalMass,CSRpert)
  use lib_fem
  use lib_rot
  use lib_rotvect

! I/O Variables.
  integer,intent(in)    :: NumNodesElem     ! Number of nodes in the element.
  real(8),intent(in)    :: r0       (:,:)   ! Initial position/orientation of grid points.
  real(8),intent(in)    :: Ri       (:,:)   ! Current position/orientation of grid points.
  real(8),intent(in)    :: Vrel     (6)     ! Linear/angular velocity of reference frame.
  real(8),intent(in)    :: NodalMass(:,:,:) ! Inertia tensor of lumped masses at element nodes.
  real(8),intent(inout) :: CSRpert  (:,:)   ! Tangent gyroscopic damping matrix.

! Local variables
  integer :: i                      ! Counters.
  integer :: iNode                  ! Counter on the element nodes.

  real(8) :: Ra(3)                  ! Displacement at Gauss point.
  real(8) :: Psi(3)                 ! Rotation vector at Gauss point.
  real(8) :: VRB(3),WRB(3)          ! Linear and angular rigid-body velocities in the material frame.

  real(8) :: CBa(3,3)               ! Coordinate transformation matrix.
  real(8) :: Rot(3,3)               ! Tangential operator.
  real(8) :: Lambda(6,6)            ! Lambda operator.
  real(8) :: ARC(6,6),ARVW(6,6)     ! A_{RC} and A^R_{Omega Omega} matrices.
  real(8) :: N (6,6*MaxNodCB3)      ! Element shape function matrix.

! Loop in the element nodes.
  do iNode=1,NumNodesElem

! Compute the current position vector and rotations, and their derivatives.
	Ra    = Ri   (iNode,1:3)
	Psi   = Ri   (iNode,4:6)

! Compute element shape functions at the corrent node.
    N =0.d0
    do i=1,6
      N (i,i+(iNode-1)*6)= 1.d0
    end do

! Compute the current coordinate transformation matrix and rotational operator.
    CBa= rotvect_psi2mat  (Psi)
    Rot= rotvect_psi2rot  (Psi)

! Compute the current linear and angular velocity.
    VRB= matmul(CBa,rot_cross(Vrel(4:6),Ra)+Vrel(1:3))
    WRB= matmul(CBa,Vrel(4:6))

!Operators for tangent matrix.
    Lambda=0.0d0
    Lambda(1:3,1:3)=CBa
    Lambda(4:6,4:6)=Rot

    ARC=0.d0
    ARC(1:3,1:3)=transpose(CBa)
    ARC(4:6,1:3)=matmul(rot_skew(Ra),transpose(CBa))
    ARC(4:6,4:6)=transpose(CBa)

    ARVW(1:3,1:3)=rot_skew(WRB)
    ARVW(1:3,4:6)=0.d0
    ARVW(4:6,1:3)=rot_skew(VRB)
    ARVW(4:6,4:6)=rot_skew(WRB)

! Compute perturbed elastic gyroscopic damping matrix.
    CSRpert=CSRpert + matmul(transpose(matmul(Lambda,N)),matmul(matmul(ARVW,NodalMass(iNode,:,:)),transpose(ARC)))

  end do
  return
 end subroutine perturb_rbcsr


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-> Subroutine PERTURB_CRS
!
!-> Description:
!
!    Compute tangent perturbed gyroscopic damping matrix for structural contribution to rigid-body motion
!
!-> Remarks.-
!
!   1) New values are added to those already in the matrix.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine perturb_crs (NumNodesElem,r0,Ri,Vrel,ElemMass,CRSpert,NumGauss)
  use lib_fem
  use lib_rot
  use lib_rotvect

! I/O Variables.
  integer,intent(in)    :: NumNodesElem     ! Number of nodes in the element.
  real(8),intent(in)    :: r0       (:,:)   ! Initial position/orientation of grid points.
  real(8),intent(in)    :: Ri       (:,:)   ! Current position/orientation of grid points.
  real(8),intent(in)    :: Vrel     (6)     ! Linear/angular velocity of reference frame.
  real(8),intent(in)    :: ElemMass (6,6)   ! Inertia properties in the element.
  real(8),intent(inout) :: CRSpert  (:,:)   ! Tangent gyroscopic damping matrix.
  integer,intent(in)    :: NumGauss         ! Number of Gauss points in the element.

! Local variables
  integer :: i,j                    ! Counters.
  integer :: iGauss                 ! Counter on the Gaussian points.
  real(8) :: Jacobian               ! ds/deta, with eta the nondimensional arclength.
  real(8) :: ShapeFun(MaxNodCB3)    ! Shape functions in a gauss point.
  real(8) :: ShapeDer(MaxNodCB3)    ! Derivatives of ShapeFun in a gauss point.
  real(8),allocatable :: CoordGauss (:) ! Coords of Gauss points.
  real(8),allocatable :: WeightGauss(:) ! Coords of Gauss points.

  real(8) :: dr0_g (6)              ! Spacial derivative of r0 at Gauss point.
  real(8) :: Ra(3)                  ! Displacement at Gauss point.
  real(8) :: Psi(3)                 ! Rotation vector at Gauss point.
  real(8) :: VRB(3),WRB(3)          ! Linear and angular rigid-body velocities in the material frame.
  real(8) :: PRB(3),HRB(3)          ! Linear and angular rigid-body momenta in the material frame.

  real(8) :: CBa(3,3)               ! Coordinate transformation matrix.
  real(8) :: Rot(3,3)               ! Tangential operator.
  real(8) :: Lambda(6,6)            ! Lambda operator.
  real(8) :: ARC(6,6),ARVW(6,6)     ! A_{RC} and A^R_{Omega Omega} matrices.
  real(8) :: ARdeltaWMass (6,6)     ! A^R_{Delta Omega}*ElemMass matrix.
  real(8) :: N (6,6*MaxNodCB3)      ! Element shape function matrix.

! Define Gauss points and loop on them.
  allocate (CoordGauss(NumGauss))
  allocate (WeightGauss(NumGauss))
  call fem_1d_gauss_val (NumGauss,CoordGauss,WeightGauss)

  do iGauss=1,NumGauss

! Obtain the shape functions and their derivatives.
    call fem_1d_shapefun (NumNodesElem,CoordGauss(iGauss),ShapeFun,ShapeDer)

! Global coordinates and Jacobian at current Gauss point.
    do i=1,6
      dr0_g(i)= dot_product (ShapeDer(1:NumNodesElem),r0(1:NumNodesElem,i))
    end do
    Jacobian =sqrt(dot_product (dr0_g(1:3),dr0_g(1:3)))

! Compute the current position vector and rotations, and their derivatives.
    do i=1,3
      Ra    (i)=dot_product (ShapeFun(1:NumNodesElem),Ri   (1:NumNodesElem,i))
      Psi   (i)=dot_product (ShapeFun(1:NumNodesElem),Ri   (1:NumNodesElem,i+3))
    end do

! Compute element shape functions.
    N =0.d0
    do i=1,6
      do j=1,NumNodesElem
        N (i,i+(j-1)*6)= ShapeFun(j)
      end do
    end do

! Compute the current coordinate transformation matrix and rotational operator.
    CBa= rotvect_psi2mat  (Psi)
    Rot= rotvect_psi2rot  (Psi)

! Compute the current linear and angular velocity.
    VRB= matmul(CBa,rot_cross(Vrel(4:6),Ra)+Vrel(1:3))
    WRB= matmul(CBa,Vrel(4:6))

    PRB= matmul(ElemMass(1:3,1:3),VRB)+matmul(ElemMass(1:3,4:6),WRB)
    HRB= matmul(ElemMass(4:6,1:3),VRB)+matmul(ElemMass(4:6,4:6),WRB)

!Operators for tangent matrix.
    Lambda=0.0d0
    Lambda(1:3,1:3)=CBa
    Lambda(4:6,4:6)=Rot

    ARC=0.d0
    ARC(1:3,1:3)=transpose(CBa)
    ARC(4:6,1:3)=matmul(rot_skew(Ra),transpose(CBa))
    ARC(4:6,4:6)=transpose(CBa)

    ARVW(1:3,1:3)=rot_skew(WRB)
    ARVW(1:3,4:6)=0.d0
    ARVW(4:6,1:3)=rot_skew(VRB)
    ARVW(4:6,4:6)=rot_skew(WRB)

    ARdeltaWMass=matmul(ARVW,ElemMass)

    ARdeltaWMass(1:3,4:6)= ARdeltaWMass(1:3,4:6) - rot_skew(PRB)
    ARdeltaWMass(4:6,1:3)= ARdeltaWMass(4:6,1:3) - rot_skew(PRB)
    ARdeltaWMass(4:6,4:6)= ARdeltaWMass(4:6,4:6) - rot_skew(HRB)

! Compute perturbed elastic gyroscopic damping matrix.
    CRSpert=CRSpert + WeightGauss(iGauss)*Jacobian*    &
&                     matmul(ARC,matmul(ARdeltaWMass-matmul(ElemMass,transpose(ARVW)),matmul(Lambda,N)))

  end do
  deallocate (CoordGauss,WeightGauss)
  return
 end subroutine perturb_crs


 subroutine perturb_rbcrs (NumNodesElem,r0,Ri,Vrel,NodalMass,CRSpert)
  use lib_fem
  use lib_rot
  use lib_rotvect

! I/O Variables.
  integer,intent(in)    :: NumNodesElem     ! Number of nodes in the element.
  real(8),intent(in)    :: r0       (:,:)   ! Initial position/orientation of grid points.
  real(8),intent(in)    :: Ri       (:,:)   ! Current position/orientation of grid points.
  real(8),intent(in)    :: Vrel     (6)     ! Linear/angular velocity of reference frame.
  real(8),intent(in)    :: NodalMass(:,:,:) ! Inertia tensor of lumped masses at element nodes.
  real(8),intent(inout) :: CRSpert  (:,:)   ! Tangent gyroscopic damping matrix.

! Local variables
  integer :: i                      ! Counters.
  integer :: iNode                  ! Counter on the element nodes.

  real(8) :: Ra(3)                  ! Displacement at Gauss point.
  real(8) :: Psi(3)                 ! Rotation vector at Gauss point.
  real(8) :: VRB(3),WRB(3)          ! Linear and angular rigid-body velocities in the material frame.
  real(8) :: PRB(3),HRB(3)          ! Linear and angular rigid-body momenta in the material frame.

  real(8) :: CBa(3,3)               ! Coordinate transformation matrix.
  real(8) :: Rot(3,3)               ! Tangential operator.
  real(8) :: Lambda(6,6)            ! Lambda operator.
  real(8) :: ARC(6,6),ARVW(6,6)     ! A_{RC} and A^R_{Omega Omega} matrices.
  real(8) :: ARdeltaWMass (6,6)     ! A^R_{Delta Omega}*ElemMass matrix.
  real(8) :: N (6,6*MaxNodCB3)      ! Element shape function matrix.

! Loop in the element nodes.
  do iNode=1,NumNodesElem

! Compute the current position vector and rotations, and their derivatives.
	Ra    = Ri   (iNode,1:3)
	Psi   = Ri   (iNode,4:6)

! Compute element shape functions at the corrent node.
    N =0.d0
    do i=1,6
      N (i,i+(iNode-1)*6)= 1.d0
    end do

! Compute the current coordinate transformation matrix and rotational operator.
    CBa= rotvect_psi2mat  (Psi)
    Rot= rotvect_psi2rot  (Psi)

! Compute the current linear and angular velocity.
    VRB= matmul(CBa,rot_cross(Vrel(4:6),Ra)+Vrel(1:3))
    WRB= matmul(CBa,Vrel(4:6))

    PRB= matmul(NodalMass(iNode,1:3,1:3),VRB)+matmul(NodalMass(iNode,1:3,4:6),WRB)
    HRB= matmul(NodalMass(iNode,4:6,1:3),VRB)+matmul(NodalMass(iNode,4:6,4:6),WRB)

!Operators for tangent matrix.
    Lambda=0.0d0
    Lambda(1:3,1:3)=CBa
    Lambda(4:6,4:6)=Rot

    ARC=0.d0
    ARC(1:3,1:3)=transpose(CBa)
    ARC(4:6,1:3)=matmul(rot_skew(Ra),transpose(CBa))
    ARC(4:6,4:6)=transpose(CBa)

    ARVW(1:3,1:3)=rot_skew(WRB)
    ARVW(1:3,4:6)=0.d0
    ARVW(4:6,1:3)=rot_skew(VRB)
    ARVW(4:6,4:6)=rot_skew(WRB)

    ARdeltaWMass=matmul(ARVW,NodalMass(iNode,:,:))

    ARdeltaWMass(1:3,4:6)= ARdeltaWMass(1:3,4:6) - rot_skew(PRB)
    ARdeltaWMass(4:6,1:3)= ARdeltaWMass(4:6,1:3) - rot_skew(PRB)
    ARdeltaWMass(4:6,4:6)= ARdeltaWMass(4:6,4:6) - rot_skew(HRB)

! Compute perturbed elastic gyroscopic damping matrix.
    CRSpert=CRSpert + matmul(ARC,matmul(ARdeltaWMass-matmul(NodalMass(iNode,:,:),transpose(ARVW)),matmul(Lambda,N)))

  end do
  return
 end subroutine perturb_rbcrs


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-> Subroutine PERTURB_CRR
!
!-> Description:
!
!    Compute tangent perturbed gyroscopic damping matrix for rigid-body dof
!
!-> Remarks.-
!
!   1) New values are added to those already in the matrix.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine perturb_crr (NumNodesElem,r0,Ri,Vrel,ElemMass,CRRpert,NumGauss)
  use lib_fem
  use lib_rot
  use lib_rotvect

! I/O Variables.
  integer,intent(in)    :: NumNodesElem     ! Number of nodes in the element.
  real(8),intent(in)    :: r0       (:,:)   ! Initial position/orientation of grid points.
  real(8),intent(in)    :: Ri       (:,:)   ! Current position/orientation of grid points.
  real(8),intent(in)    :: Vrel     (6)     ! Linear/angular velocity of reference frame.
  real(8),intent(in)    :: ElemMass (6,6)   ! Inertia properties in the element.
  real(8),intent(inout) :: CRRpert  (:,:)   ! Tangent gyroscopic damping matrix.
  integer,intent(in)    :: NumGauss         ! Number of Gauss points in the element.

! Local variables
  integer :: i                      ! Counters.
  integer :: iGauss                 ! Counter on the Gaussian points.
  real(8) :: Jacobian               ! ds/deta, with eta the nondimensional arclength.
  real(8) :: ShapeFun(MaxNodCB3)    ! Shape functions in a gauss point.
  real(8) :: ShapeDer(MaxNodCB3)    ! Derivatives of ShapeFun in a gauss point.
  real(8),allocatable :: CoordGauss (:) ! Coords of Gauss points.
  real(8),allocatable :: WeightGauss(:) ! Coords of Gauss points.

  real(8) :: dr0_g (6)              ! Spacial derivative of r0 at Gauss point.
  real(8) :: Ra(3)                  ! Displacement at Gauss point.
  real(8) :: Psi(3)                 ! Rotation vector at Gauss point.
  real(8) :: VRB(3),WRB(3)          ! Linear and angular rigid-body velocities in the material frame.

  real(8) :: CBa(3,3)               ! Coordinate transformation matrix.
  real(8) :: ARC(6,6),ARVW(6,6)     ! A_{RC} and A^R_{Omega Omega} matrices.

! Define Gauss points and loop on them.
  allocate (CoordGauss(NumGauss))
  allocate (WeightGauss(NumGauss))
  call fem_1d_gauss_val (NumGauss,CoordGauss,WeightGauss)

  do iGauss=1,NumGauss

! Obtain the shape functions and their derivatives.
    call fem_1d_shapefun (NumNodesElem,CoordGauss(iGauss),ShapeFun,ShapeDer)

! Global coordinates and Jacobian at current Gauss point.
    do i=1,6
      dr0_g(i)= dot_product (ShapeDer(1:NumNodesElem),r0(1:NumNodesElem,i))
    end do
    Jacobian =sqrt(dot_product (dr0_g(1:3),dr0_g(1:3)))

! Compute the current position vector and rotations, and their derivatives.
    do i=1,3
      Ra    (i)=dot_product (ShapeFun(1:NumNodesElem),Ri   (1:NumNodesElem,i))
      Psi   (i)=dot_product (ShapeFun(1:NumNodesElem),Ri   (1:NumNodesElem,i+3))
    end do

! Compute the current coordinate transformation matrix.
    CBa= rotvect_psi2mat  (Psi)

! Compute the current linear and angular velocity.
    VRB= matmul(CBa,rot_cross(Vrel(4:6),Ra)+Vrel(1:3))
    WRB= matmul(CBa,Vrel(4:6))

!Operators for tangent matrix.
    ARC=0.d0
    ARC(1:3,1:3)=transpose(CBa)
    ARC(4:6,1:3)=matmul(rot_skew(Ra),transpose(CBa))
    ARC(4:6,4:6)=transpose(CBa)

    ARVW(1:3,1:3)=rot_skew(WRB)
    ARVW(1:3,4:6)=0.d0
    ARVW(4:6,1:3)=rot_skew(VRB)
    ARVW(4:6,4:6)=rot_skew(WRB)

! Compute perturbed elastic gyroscopic damping matrix.
    CRRpert=CRRpert + WeightGauss(iGauss)*Jacobian*matmul(matmul(matmul(ARC,ARVW),ElemMass),transpose(ARC))

  end do
  deallocate (CoordGauss,WeightGauss)
  return
 end subroutine perturb_crr


 subroutine perturb_rbcrr (NumNodesElem,r0,Ri,Vrel,NodalMass,CRRpert)
  use lib_fem
  use lib_rot
  use lib_rotvect

! I/O Variables.
  integer,intent(in)    :: NumNodesElem     ! Number of nodes in the element.
  real(8),intent(in)    :: r0       (:,:)   ! Initial position/orientation of grid points.
  real(8),intent(in)    :: Ri       (:,:)   ! Current position/orientation of grid points.
  real(8),intent(in)    :: Vrel     (6)     ! Linear/angular velocity of reference frame.
  real(8),intent(in)    :: NodalMass(:,:,:) ! Inertia tensor of lumped masses at element nodes.
  real(8),intent(inout) :: CRRpert  (:,:)   ! Tangent gyroscopic damping matrix.

! Local variables
  integer :: iNode                  ! Counter on the element nodes.

  real(8) :: Ra(3)                  ! Displacement at Gauss point.
  real(8) :: Psi(3)                 ! Rotation vector at Gauss point.
  real(8) :: VRB(3),WRB(3)          ! Linear and angular rigid-body velocities in the material frame.

  real(8) :: CBa(3,3)               ! Coordinate transformation matrix.
  real(8) :: ARC(6,6),ARVW(6,6)     ! A_{RC} and A^R_{Omega Omega} matrices.

! Loop in the element nodes.
  do iNode=1,NumNodesElem

! Compute the current position vector and rotations, and their derivatives.
	Ra    = Ri   (iNode,1:3)
	Psi   = Ri   (iNode,4:6)

! Compute the current coordinate transformation matrix.
    CBa= rotvect_psi2mat  (Psi)

! Compute the current linear and angular velocity.
    VRB= matmul(CBa,rot_cross(Vrel(4:6),Ra)+Vrel(1:3))
    WRB= matmul(CBa,Vrel(4:6))

!Operators for tangent matrix.
    ARC=0.d0
    ARC(1:3,1:3)=transpose(CBa)
    ARC(4:6,1:3)=matmul(rot_skew(Ra),transpose(CBa))
    ARC(4:6,4:6)=transpose(CBa)

    ARVW(1:3,1:3)=rot_skew(WRB)
    ARVW(1:3,4:6)=0.d0
    ARVW(4:6,1:3)=rot_skew(VRB)
    ARVW(4:6,4:6)=rot_skew(WRB)

! Compute perturbed elastic gyroscopic damping matrix.
    CRRpert=CRRpert + matmul(matmul(matmul(ARC,ARVW),NodalMass(iNode,:,:)),transpose(ARC))

  end do
  return
 end subroutine perturb_rbcrr


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-> Subroutine PERTURB_KSS
!
!-> Description:
!
!    Compute tangent gyroscopic stiffness matrix for structural dof
!
!-> Remarks.-
!
!   1) New values are added to those already in the matrix.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine perturb_kss (NumNodesElem,r0,Ri,Vrel,ElemMass,KSSpert,NumGauss)
  use lib_fem
  use lib_rot
  use lib_rotvect

! I/O Variables.
  integer,intent(in)    :: NumNodesElem     ! Number of nodes in the element.
  real(8),intent(in)    :: r0       (:,:)   ! Initial position/orientation of grid points.
  real(8),intent(in)    :: Ri       (:,:)   ! Current position/orientation of grid points.
  real(8),intent(in)    :: Vrel     (6)     ! Linear/angular velocity of reference frame.
  real(8),intent(in)    :: ElemMass (6,6)   ! Inertia properties in the element.
  real(8),intent(inout) :: KSSpert  (:,:)   ! Tangent gyroscopic damping matrix.
  integer,intent(in)    :: NumGauss         ! Number of Gauss points in the element.

! Local variables
  integer :: i,j                    ! Counters.
  integer :: iGauss                 ! Counter on the Gaussian points.
  real(8) :: Jacobian               ! ds/deta, with eta the nondimensional arclength.
  real(8) :: ShapeFun(MaxNodCB3)    ! Shape functions in a gauss point.
  real(8) :: ShapeDer(MaxNodCB3)    ! Derivatives of ShapeFun in a gauss point.
  real(8),allocatable :: CoordGauss (:) ! Coords of Gauss points.
  real(8),allocatable :: WeightGauss(:) ! Coords of Gauss points.

  real(8) :: dr0_g (6)              ! Spacial derivative of r0 at Gauss point.
  real(8) :: Ra(3)                  ! Displacement at Gauss point.
  real(8) :: Psi(3)                 ! Rotation vector at Gauss point.
  real(8) :: VRB(3),WRB(3)          ! Linear and angular rigid-body velocities in the material frame.
  real(8) :: PRB(3),HRB(3)          ! Linear and angular rigid-body momenta in the material frame.

  real(8) :: CBa(3,3)               ! Coordinate transformation matrix.
  real(8) :: Rot(3,3)               ! Tangential operator.
  real(8) :: Lambda(6,6)            ! Lambda operator.
  real(8) :: LambdaBar(6,6)         ! Lambda operator.
  real(8) :: ARVW (6,6)             ! A^R_{Omega Omega} matrix.
  real(8) :: ARdeltaWMass (6,6)     ! A^R_{Delta Omega}*ElemMass matrix.
  real(8) :: N (6,6*MaxNodCB3)      ! Element shape function matrix.

! Define Gauss points and loop on them.
  allocate (CoordGauss(NumGauss))
  allocate (WeightGauss(NumGauss))
  call fem_1d_gauss_val (NumGauss,CoordGauss,WeightGauss)

  do iGauss=1,NumGauss

! Obtain the shape functions and their derivatives.
    call fem_1d_shapefun (NumNodesElem,CoordGauss(iGauss),ShapeFun,ShapeDer)

! Global coordinates and Jacobian at current Gauss point.
    do i=1,6
      dr0_g(i)= dot_product (ShapeDer(1:NumNodesElem),r0(1:NumNodesElem,i))
    end do
    Jacobian =sqrt(dot_product (dr0_g(1:3),dr0_g(1:3)))

! Compute the current position vector and rotations, and their derivatives.
    do i=1,3
      Ra    (i)=dot_product (ShapeFun(1:NumNodesElem),Ri   (1:NumNodesElem,i))
      Psi   (i)=dot_product (ShapeFun(1:NumNodesElem),Ri   (1:NumNodesElem,i+3))
    end do

! Compute element shape functions.
    N =0.d0
    do i=1,6
      do j=1,NumNodesElem
        N (i,i+(j-1)*6)= ShapeFun(j)
      end do
    end do

! Compute the current coordinate transformation matrix and rotational operator.
    CBa= rotvect_psi2mat  (Psi)
    Rot= rotvect_psi2rot  (Psi)

! Compute the current linear and angular velocity.
    VRB= matmul(CBa,rot_cross(Vrel(4:6),Ra)+Vrel(1:3))
    WRB= matmul(CBa,Vrel(4:6))

    PRB= matmul(ElemMass(1:3,1:3),VRB)+matmul(ElemMass(1:3,4:6),WRB)
    HRB= matmul(ElemMass(4:6,1:3),VRB)+matmul(ElemMass(4:6,4:6),WRB)

!Operators for tangent matrix.
    Lambda=0.0d0
    Lambda(1:3,1:3)=CBa
    Lambda(4:6,4:6)=Rot

    LambdaBar=0.0d0
    LambdaBar(1:3,4:6)=-matmul(matmul(transpose(CBa),rot_skew(matmul(rot_skew(WRB),PRB))),Rot)
    LambdaBar(4:6,4:6)=-rotvect_a1(-Psi,matmul(rot_skew(VRB),PRB)+matmul(rot_skew(WRB),HRB))

    ARVW=0.d0
    ARVW(1:3,1:3)=rot_skew(WRB)
    ARVW(1:3,4:6)=0.d0
    ARVW(4:6,1:3)=rot_skew(VRB)
    ARVW(4:6,4:6)=rot_skew(WRB)

    ARdeltaWMass=matmul(ARVW,ElemMass)

    ARdeltaWMass(1:3,4:6)= ARdeltaWMass(1:3,4:6) - rot_skew(PRB)
    ARdeltaWMass(4:6,1:3)= ARdeltaWMass(4:6,1:3) - rot_skew(PRB)
    ARdeltaWMass(4:6,4:6)= ARdeltaWMass(4:6,4:6) - rot_skew(HRB)

! Compute perturbed elastic gyroscopic damping matrix.
    KSSpert=KSSpert + WeightGauss(iGauss)*Jacobian*matmul(transpose(N),matmul(LambdaBar,N))

    KSSpert=KSSpert - WeightGauss(iGauss)*Jacobian*    &
&                     matmul(transpose(matmul(Lambda,N)),matmul(matmul(ARdeltaWMass,transpose(ARVW)),matmul(Lambda,N)))

  end do
  deallocate (CoordGauss,WeightGauss)
  return
 end subroutine perturb_kss

 
 subroutine perturb_rbkss (NumNodesElem,r0,Ri,Vrel,NodalMass,KSSpert)
  use lib_fem
  use lib_rot
  use lib_rotvect

! I/O Variables.
  integer,intent(in)    :: NumNodesElem     ! Number of nodes in the element.
  real(8),intent(in)    :: r0       (:,:)   ! Initial position/orientation of grid points.
  real(8),intent(in)    :: Ri       (:,:)   ! Current position/orientation of grid points.
  real(8),intent(in)    :: Vrel     (6)     ! Linear/angular velocity of reference frame.
  real(8),intent(in)    :: NodalMass(:,:,:) ! Inertia tensor of lumped masses at element nodes.
  real(8),intent(inout) :: KSSpert  (:,:)   ! Tangent gyroscopic damping matrix.

! Local variables
  integer :: i                      ! Counters.
  integer :: iNode                  ! Counter on the element nodes.

  real(8) :: Ra(3)                  ! Displacement at Gauss point.
  real(8) :: Psi(3)                 ! Rotation vector at Gauss point.
  real(8) :: VRB(3),WRB(3)          ! Linear and angular rigid-body velocities in the material frame.
  real(8) :: PRB(3),HRB(3)          ! Linear and angular rigid-body momenta in the material frame.

  real(8) :: CBa(3,3)               ! Coordinate transformation matrix.
  real(8) :: Rot(3,3)               ! Tangential operator.
  real(8) :: Lambda(6,6)            ! Lambda operator.
  real(8) :: LambdaBar(6,6)         ! Lambda operator.
  real(8) :: ARVW (6,6)             ! A^R_{Omega Omega} matrix.
  real(8) :: ARdeltaWMass (6,6)     ! A^R_{Delta Omega}*ElemMass matrix.
  real(8) :: N (6,6*MaxNodCB3)      ! Element shape function matrix.

! Loop in the element nodes.
  do iNode=1,NumNodesElem

! Compute the current position vector and rotations, and their derivatives.
	Ra    = Ri   (iNode,1:3)
	Psi   = Ri   (iNode,4:6)

! Compute element shape functions at the corrent node.
    N =0.d0
    do i=1,6
      N (i,i+(iNode-1)*6)= 1.d0
    end do

! Compute the current coordinate transformation matrix and rotational operator.
    CBa= rotvect_psi2mat  (Psi)
    Rot= rotvect_psi2rot  (Psi)

! Compute the current linear and angular velocity.
    VRB= matmul(CBa,rot_cross(Vrel(4:6),Ra)+Vrel(1:3))
    WRB= matmul(CBa,Vrel(4:6))

    PRB= matmul(NodalMass(iNode,1:3,1:3),VRB)+matmul(NodalMass(iNode,1:3,4:6),WRB)
    HRB= matmul(NodalMass(iNode,4:6,1:3),VRB)+matmul(NodalMass(iNode,4:6,4:6),WRB)

!Operators for tangent matrix.
    Lambda=0.0d0
    Lambda(1:3,1:3)=CBa
    Lambda(4:6,4:6)=Rot

    LambdaBar=0.0d0
    LambdaBar(1:3,4:6)=-matmul(matmul(transpose(CBa),rot_skew(matmul(rot_skew(WRB),PRB))),Rot)
    LambdaBar(4:6,4:6)=-rotvect_a1(-Psi,matmul(rot_skew(VRB),PRB)+matmul(rot_skew(WRB),HRB))

    ARVW=0.d0
    ARVW(1:3,1:3)=rot_skew(WRB)
    ARVW(1:3,4:6)=0.d0
    ARVW(4:6,1:3)=rot_skew(VRB)
    ARVW(4:6,4:6)=rot_skew(WRB)

    ARdeltaWMass=matmul(ARVW,NodalMass(iNode,:,:))

    ARdeltaWMass(1:3,4:6)= ARdeltaWMass(1:3,4:6) - rot_skew(PRB)
    ARdeltaWMass(4:6,1:3)= ARdeltaWMass(4:6,1:3) - rot_skew(PRB)
    ARdeltaWMass(4:6,4:6)= ARdeltaWMass(4:6,4:6) - rot_skew(HRB)

! Compute perturbed elastic gyroscopic damping matrix.
    KSSpert=KSSpert + matmul(transpose(N),matmul(LambdaBar,N))

    KSSpert=KSSpert - matmul(transpose(matmul(Lambda,N)),matmul(matmul(ARdeltaWMass,transpose(ARVW)),matmul(Lambda,N)))

  end do
  return
 end subroutine perturb_rbkss

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-> Subroutine PERTURB_KRS
!
!-> Description:
!
!    Compute tangent gyroscopic stiffness matrix for structural contribution to rigid-body dynamics
!
!-> Remarks.-
!
!   1) New values are added to those already in the matrix.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine perturb_krs (NumNodesElem,r0,Ri,Vrel,ElemMass,KRSpert,NumGauss)
  use lib_fem
  use lib_rot
  use lib_rotvect

! I/O Variables.
  integer,intent(in)    :: NumNodesElem     ! Number of nodes in the element.
  real(8),intent(in)    :: r0       (:,:)   ! Initial position/orientation of grid points.
  real(8),intent(in)    :: Ri       (:,:)   ! Current position/orientation of grid points.
  real(8),intent(in)    :: Vrel     (6)     ! Linear/angular velocity of reference frame.
  real(8),intent(in)    :: ElemMass (6,6)   ! Inertia properties in the element.
  real(8),intent(inout) :: KRSpert  (:,:)   ! Tangent gyroscopic damping matrix.
  integer,intent(in)    :: NumGauss         ! Number of Gauss points in the element.

! Local variables
  integer :: i,j                    ! Counters.
  integer :: iGauss                 ! Counter on the Gaussian points.
  real(8) :: Jacobian               ! ds/deta, with eta the nondimensional arclength.
  real(8) :: ShapeFun(MaxNodCB3)    ! Shape functions in a gauss point.
  real(8) :: ShapeDer(MaxNodCB3)    ! Derivatives of ShapeFun in a gauss point.
  real(8),allocatable :: CoordGauss (:) ! Coords of Gauss points.
  real(8),allocatable :: WeightGauss(:) ! Coords of Gauss points.

  real(8) :: dr0_g (6)              ! Spacial derivative of r0 at Gauss point.
  real(8) :: Ra(3)                  ! Displacement at Gauss point.
  real(8) :: Psi(3)                 ! Rotation vector at Gauss point.
  real(8) :: VRB(3),WRB(3)          ! Linear and angular rigid-body velocities in the material frame.
  real(8) :: PRB(3),HRB(3)          ! Linear and angular rigid-body momenta in the material frame.

  real(8) :: CBa(3,3)                   ! Coordinate transformation matrix.
  real(8) :: Rot(3,3)                   ! Tangential operator.
  real(8) :: Lambda(6,6)                ! Lambda operator.
  real(8) :: ARC(6,6),ARVW(6,6)         ! A_{RC} and A^R_{Omega Omega} matrices.
  real(8) :: MassBRdeltaW(6,6)          ! A^R_{Delta Omega}*ElemMass matrix.
  real(8) :: N (6,6*MaxNodCB3)          ! Element shape function matrix.

! Define Gauss points and loop on them.
  allocate (CoordGauss(NumGauss))
  allocate (WeightGauss(NumGauss))
  call fem_1d_gauss_val (NumGauss,CoordGauss,WeightGauss)

  do iGauss=1,NumGauss

! Obtain the shape functions and their derivatives.
    call fem_1d_shapefun (NumNodesElem,CoordGauss(iGauss),ShapeFun,ShapeDer)

! Global coordinates and Jacobian at current Gauss point.
    do i=1,6
      dr0_g(i)= dot_product (ShapeDer(1:NumNodesElem),r0(1:NumNodesElem,i))
    end do
    Jacobian =sqrt(dot_product (dr0_g(1:3),dr0_g(1:3)))

! Compute the current position vector and rotations, and their derivatives.
    do i=1,3
      Ra    (i)=dot_product (ShapeFun(1:NumNodesElem),Ri   (1:NumNodesElem,i))
      Psi   (i)=dot_product (ShapeFun(1:NumNodesElem),Ri   (1:NumNodesElem,i+3))
    end do

! Compute element shape functions.
    N =0.d0
    do i=1,6
      do j=1,NumNodesElem
        N (i,i+(j-1)*6)= ShapeFun(j)
      end do
    end do

! Compute the current coordinate transformation matrix and rotational operator.
    CBa= rotvect_psi2mat  (Psi)
    Rot= rotvect_psi2rot  (Psi)

! Compute the current linear and angular velocity.
    VRB= matmul(CBa,rot_cross(Vrel(4:6),Ra)+Vrel(1:3))
    WRB= matmul(CBa,Vrel(4:6))

    PRB= matmul(ElemMass(1:3,1:3),VRB)+matmul(ElemMass(1:3,4:6),WRB)
    HRB= matmul(ElemMass(4:6,1:3),VRB)+matmul(ElemMass(4:6,4:6),WRB)

!Operators for tangent matrix.
    Lambda=0.0d0
    Lambda(1:3,1:3)=CBa
    Lambda(4:6,4:6)=Rot

    ARC=0.d0
    ARC(1:3,1:3)=transpose(CBa)
    ARC(4:6,1:3)=matmul(rot_skew(Ra),transpose(CBa))
    ARC(4:6,4:6)=transpose(CBa)

    ARVW=0.d0
    ARVW(1:3,1:3)=rot_skew(WRB)
    ARVW(1:3,4:6)=0.d0
    ARVW(4:6,1:3)=rot_skew(VRB)
    ARVW(4:6,4:6)=rot_skew(WRB)

    MassBRdeltaW=matmul(ElemMass,-transpose(ARVW))

    MassBRdeltaW(1:3,4:6)= MassBRdeltaW(1:3,4:6) - rot_skew(PRB)
    MassBRdeltaW(4:6,1:3)= MassBRdeltaW(4:6,1:3) - rot_skew(PRB)
    MassBRdeltaW(4:6,4:6)= MassBRdeltaW(4:6,4:6) - rot_skew(HRB)

! Compute perturbed elastic gyroscopic damping matrix.
    KRSpert=KRSpert + WeightGauss(iGauss)*Jacobian*matmul(matmul(matmul(ARC,ARVW),MassBRdeltaW),matmul(Lambda,N))

  end do
  deallocate (CoordGauss,WeightGauss)
  return
 end subroutine perturb_krs
 
 
 subroutine perturb_rbkrs (NumNodesElem,r0,Ri,Vrel,NodalMass,KRSpert)
  use lib_fem
  use lib_rot
  use lib_rotvect

! I/O Variables.
  integer,intent(in)    :: NumNodesElem     ! Number of nodes in the element.
  real(8),intent(in)    :: r0       (:,:)   ! Initial position/orientation of grid points.
  real(8),intent(in)    :: Ri       (:,:)   ! Current position/orientation of grid points.
  real(8),intent(in)    :: Vrel     (6)     ! Linear/angular velocity of reference frame.
  real(8),intent(in)    :: NodalMass(:,:,:) ! Inertia tensor of lumped masses at element nodes.
  real(8),intent(inout) :: KRSpert  (:,:)   ! Tangent gyroscopic damping matrix.

! Local variables
  integer :: i                      ! Counters.
  integer :: iNode                  ! Counter on the element nodes.

  real(8) :: Ra(3)                  ! Displacement at Gauss point.
  real(8) :: Psi(3)                 ! Rotation vector at Gauss point.
  real(8) :: VRB(3),WRB(3)          ! Linear and angular rigid-body velocities in the material frame.
  real(8) :: PRB(3),HRB(3)          ! Linear and angular rigid-body momenta in the material frame.

  real(8) :: CBa(3,3)                   ! Coordinate transformation matrix.
  real(8) :: Rot(3,3)                   ! Tangential operator.
  real(8) :: Lambda(6,6)                ! Lambda operator.
  real(8) :: ARC(6,6),ARVW(6,6)         ! A_{RC} and A^R_{Omega Omega} matrices.
  real(8) :: MassBRdeltaW(6,6)          ! A^R_{Delta Omega}*ElemMass matrix.
  real(8) :: N (6,6*MaxNodCB3)          ! Element shape function matrix.

! Loop in the element nodes.
  do iNode=1,NumNodesElem

! Compute the current position vector and rotations, and their derivatives.
	Ra    = Ri   (iNode,1:3)
	Psi   = Ri   (iNode,4:6)

! Compute element shape functions at the corrent node.
    N =0.d0
    do i=1,6
      N (i,i+(iNode-1)*6)= 1.d0
    end do

! Compute the current coordinate transformation matrix and rotational operator.
    CBa= rotvect_psi2mat  (Psi)
    Rot= rotvect_psi2rot  (Psi)

! Compute the current linear and angular velocity.
    VRB= matmul(CBa,rot_cross(Vrel(4:6),Ra)+Vrel(1:3))
    WRB= matmul(CBa,Vrel(4:6))

    PRB= matmul(NodalMass(iNode,1:3,1:3),VRB)+matmul(NodalMass(iNode,1:3,4:6),WRB)
    HRB= matmul(NodalMass(iNode,4:6,1:3),VRB)+matmul(NodalMass(iNode,4:6,4:6),WRB)

!Operators for tangent matrix.
    Lambda=0.0d0
    Lambda(1:3,1:3)=CBa
    Lambda(4:6,4:6)=Rot

    ARC=0.d0
    ARC(1:3,1:3)=transpose(CBa)
    ARC(4:6,1:3)=matmul(rot_skew(Ra),transpose(CBa))
    ARC(4:6,4:6)=transpose(CBa)

    ARVW=0.d0
    ARVW(1:3,1:3)=rot_skew(WRB)
    ARVW(1:3,4:6)=0.d0
    ARVW(4:6,1:3)=rot_skew(VRB)
    ARVW(4:6,4:6)=rot_skew(WRB)

    MassBRdeltaW=matmul(NodalMass(iNode,:,:),-transpose(ARVW))

    MassBRdeltaW(1:3,4:6)= MassBRdeltaW(1:3,4:6) - rot_skew(PRB)
    MassBRdeltaW(4:6,1:3)= MassBRdeltaW(4:6,1:3) - rot_skew(PRB)
    MassBRdeltaW(4:6,4:6)= MassBRdeltaW(4:6,4:6) - rot_skew(HRB)

! Compute perturbed elastic gyroscopic damping matrix.
    KRSpert=KRSpert + matmul(matmul(matmul(ARC,ARVW),MassBRdeltaW),matmul(Lambda,N))

  end do
  return
 end subroutine perturb_rbkrs


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end module lib_perturb

