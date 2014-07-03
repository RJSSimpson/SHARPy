!-> Copyright by Imperial College London, 2011
!
!-> Module.- LIB_CBEAM3 Rafa Palacios. 15Jul2008 - Last Update 20Apr2011
!
!-> Language: FORTRAN90, Free format.
!
!-> Description.-
!
!  Compute element matrices for the 3-noded geometrically-exact beam element.
!
!-> Subroutines.-
!
!    -cbeam3_kmat:      Compute material tangent stiffness matrix.
!    -cbeam3_kgeom:     Compute geometrical tangent stiffness matrix.
!    -cbeam3_fstif:     Compute discrete stiffness forces.
!    -cbeam3_kfindiff:  Compute tangent stiffness matrix by finite differences.
!
!    -cbeam3_mass:      Compute tangent mass matrix.
!    -cbeam3_cgyr:      Compute tangent gyroscopic damping matrix.
!    -cbeam3_kgyr:      Compute tangent gyroscopic stiffness matrix.
!    -cbeam3_fgyr:      Compute discrete gyroscopic forces.
!
!    -cbeam3_mvel:      Compute reference system mass matrix.
!    -cbeam3_cvel:      Compute reference system damping matrix.
!
!    -cbeam3_rbmass:    Compute lumped mass contribution to tangent mass matrix.
!    -cbeam3_rbcgyr:    Compute lumped mass contribution to tangent gyroscopic damping matrix.
!    -cbeam3_rbkgyr:    Compute lumped mass contribution to tangent gyroscopic stiffness matrix.
!    -cbeam3_rbfgyr:    Compute lumped mass contribution to discrete gyroscopic forces.
!    -cbeam3_rbmvel:    Compute lumped mass contribution to reference system mass matrix.
!    -cbeam3_rbcvel:    Compute lumped mass contribution to reference system damping matrix.
!
!    -cbeam3_fext:      Compute influence coefficients for force vector.
!    -cbeam3_projs2m:   Compute transformation from slave to master equation.
!    -cbeam3_projm2s:   Convert rotations from master to slave nodes.
!    -cbeam3_proj2a:    Project vector from local to element frame.
!    -cbeam3_glob2loc:  Compute transformation from global to local frame.
!   
!-> Remarks.-
!
!  1) The degrees of freedom in the element are the displacements and the
!     Cartesian Rotation Vector (CRV) at each of the nodes)
!
!  2) The element can have 2 or 3 nodes.
!
!  3) Update 20.04.2011: Added subroutines for lumped mass contributions.
!
!-> Modifications.-
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module lib_cbeam3
  implicit none

! Shared variables within the module.
 integer,private,parameter:: MaxNodCB3=3               ! Max number of nodes per element is 3.
 real(8),private,parameter,dimension(3,3):: Unit= &    ! Unit matrix.
&         reshape((/1.d0,0.d0,0.d0,0.d0,1.d0,0.d0,0.d0,0.d0,1.d0/),(/3,3/))

 contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-> Subroutine CBEAM3_KMAT
!
!-> Description:
!
!    Compute the material tangent stiffness matrix.
!
!-> Remarks.-
!
!   1) New values are added to those already in the matrix.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine cbeam3_kmat (NumNodesElem,r0,Ri,ElemStiff,Kmat,NumGauss)
  use lib_fem
  use lib_rot
  use lib_rotvect

! I/O Variables.
  integer,intent(in)   :: NumNodesElem      ! Number of nodes in the element.
  real(8),intent(in)   :: r0       (:,:)    ! Initial position/orientation of grid points.
  real(8),intent(in)   :: Ri       (:,:)    ! Current position/orientation of grid points.
  real(8),intent(in)   :: ElemStiff(:,:)    ! Stiffness properties in the element.
  real(8),intent(inout):: Kmat     (:,:)   ! Material tangent stiffness matrix.
  integer,intent(in)   :: NumGauss          ! Number of Gauss points in the element.

! Local variables.
  integer :: i,j                   ! Counters.
  integer :: iGauss                ! Counter on the Gaussian points.
  real(8) :: Jacobian              ! ds/deta, with eta the nondimensional arclength.
  real(8),allocatable :: CoordGauss (:) ! Coords of Gauss points.
  real(8),allocatable :: WeightGauss(:) ! Coords of Gauss points.
  real(8) :: ShapeFun(MaxNodCB3)    ! Shape functions in a gauss point.
  real(8) :: ShapeDer(MaxNodCB3)    ! Derivatives of ShapeFun in a gauss point.

  real(8) :: Ri_g (6)              ! Displacement/rotation at Gauss point.
  real(8) :: r0_g (6)              ! Position vector at the Gauss point.
  real(8) :: dRi_g(6)              ! Derivative of Ri_g.
  real(8) :: dr0_g(6)              ! Derivative of r0_g.

  real(8) :: CBa(3,3)              ! Coordinate transformation matrix.
  real(8) :: KB(3)                 ! Curvature in the deformed configuration at Gauss point.
  real(8) :: gamma(3)              ! Strain at Gauss point.
  real(8) :: Rot(3,3),dRot(3,3)    ! Tangential operator and spatial derivative.

  real(8) :: Ak(6,6),D(6,6)        ! A_K and D.
  real(8) :: B(6,6*MaxNodCB3)       ! Strain matrix operator.
  real(8) :: Yp(6,6),dYp(6,6)      ! Ypsilon and its derivative.
  real(8) :: N (6,6*MaxNodCB3)      ! Element shape function matrix
  real(8) :: dN(6,6*MaxNodCB3)      ! Derivatives of shape function.


! Define Gauss points and loop on them.
  allocate (CoordGauss(NumGauss))
  allocate (WeightGauss(NumGauss))
  call fem_1d_gauss_val (NumGauss,CoordGauss,WeightGauss)

  do iGauss=1,NumGauss

! Obtain the shape functions and their derivatives.
    call fem_1d_shapefun (NumNodesElem,CoordGauss(iGauss),ShapeFun,ShapeDer)

! Global coordinates and their derivatives at the gauss points.
    do i=1,6
      r0_g(i) = dot_product (ShapeFun(1:NumNodesElem),r0(1:NumNodesElem,i))
      dr0_g(i)= dot_product (ShapeDer(1:NumNodesElem),r0(1:NumNodesElem,i))
    end do
    Jacobian =sqrt(dot_product (dr0_g(1:3),dr0_g(1:3)))

! Rescale the derivatives to be in the physical coordinates.
    ShapeDer=ShapeDer/Jacobian
    dr0_g=dr0_g/Jacobian

! Compute the current position vector and rotations, and their derivatives.
    do i=1,6
      Ri_g(i) = dot_product (ShapeFun(1:NumNodesElem),Ri(1:NumNodesElem,i))
      dRi_g(i)= dot_product (ShapeDer(1:NumNodesElem),Ri(1:NumNodesElem,i))
    end do

! Compute element shape function.
    N =0.d0
    dN=0.d0
    do i=1,6
      do j=1,NumNodesElem
        N (i,i+(j-1)*6)= ShapeFun(j)
        dN(i,i+(j-1)*6)= ShapeDer(j)
      end do
    end do

! Compute the current coordinate transformation matrix and rotational operator.
    CBa= rotvect_psi2mat  (Ri_g(4:6))
    Rot= rotvect_psi2rot  (Ri_g(4:6))
    dRot=rotvect_drotdpsib(Ri_g(4:6),dRi_g(4:6))

! Compute the current curvature and strain.
    gamma=matmul(CBa,dRi_g(1:3))-matmul(rotvect_psi2mat(r0_g(4:6)),dr0_g(1:3))
    KB=matmul(Rot,dRi_g(4:6))

!Operators for the strain matrix.
    D=0.d0
    D(1:3,1:3)=transpose(CBa)
    D(4:6,4:6)=Unit

    Yp=0.d0
    Yp(1:3,1:3)=Unit
    Yp(4:6,4:6)=Rot

    dYp=0.d0
    dYp(4:6,4:6)=dRot

    AK=0.d0
    AK(4:6,1:3)=-rot_skew(gamma+Unit(1,1:3))
    AK(4:6,4:6)=-rot_skew(KB)

! Compute strain matrix operator and material tangent stiffness.
    B= matmul(transpose(D), matmul(Yp,dN)) &
&    + matmul(transpose(AK),matmul(Yp, N)) &
&    + matmul(dYp,N)

    Kmat=Kmat + WeightGauss(iGauss)*Jacobian &
&             * matmul(transpose(B),matmul(ElemStiff,B))

  end do
  deallocate (CoordGauss,WeightGauss)

  return
 end subroutine cbeam3_kmat



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-> Subroutine CBEAM3_KGEOM
!
!-> Description:
!
!    Compute geometric tangent stiffness matrix
!
!-> Remarks.-
!
!   1) New values are added to those already in the matrix.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine cbeam3_kgeom (NumNodesElem,r0,Ri,ElemStiff,Kgeom,NumGauss)
  use lib_fem
  use lib_rot
  use lib_rotvect

! I/O Variables.
  integer,intent(in)      :: NumNodesElem     ! Number of nodes in the element.
  real(8),intent(in)      :: r0       (:,:)   ! Initial position/orientation of grid points.
  real(8),intent(in)      :: Ri       (:,:)   ! Current position/orientation of grid points.
  real(8),intent(in)      :: ElemStiff(:,:)   ! Stiffness properties in the element.
  real(8),intent(inout)   :: Kgeom    (:,:)   ! Geometric tangent stiffness matrix.
  integer,intent(in)      :: NumGauss         ! Number of Gauss points in the element.

! Local variables.
  integer :: i,j                    ! Counters.
  integer :: iGauss                 ! Counter on the Gaussian points.
  real(8) :: Jacobian               ! ds/deta, with eta the nondimensional arclength.
  real(8) :: ShapeFun(MaxNodCB3)    ! Shape functions in a gauss point.
  real(8) :: ShapeDer(MaxNodCB3)    ! Derivatives of ShapeFun in a gauss point.
  real(8),allocatable :: CoordGauss (:) ! Coords of Gauss points.
  real(8),allocatable :: WeightGauss(:) ! Coords of Gauss points.

  real(8) :: r0_g (6)              ! Position vector at the Gauss point.
  real(8) :: dr0_g(6)              ! Derivative of r0_g.

  real(8) :: KB(3)                 ! Curvature in the deformed configuration at Gauss point.
  real(8) :: gamma(3),kappa(3)     ! Strains at Gauss point.
  real(8) :: CBa(3,3)              ! Coordinate transformation matrix.
  real(8) :: Rot(3,3)              ! Tangential operator and spatial derivative.
  real(8) :: FB(3),MB(3)           ! Force and moment.
  real(8) :: Psi(3),dPSi(3)        ! Current orientation and derivative.
  real(8) :: dRa(3)                ! Derivative of the position vector.

  real(8) :: Aq0q0(6,6)            ! A_qq.
  real(8) :: Aq0q1(6,6)            ! A_qq'.
  real(8) :: Aq1q0(6,6)            ! A_q'q.
  real(8) :: N (6,6*MaxNodCB3)     ! Element shape function matrix
  real(8) :: dN(6,6*MaxNodCB3)     ! Derivatives of shape function.


! Define Gauss points and loop on them.
  allocate (CoordGauss(NumGauss))
  allocate (WeightGauss(NumGauss))
  call fem_1d_gauss_val (NumGauss,CoordGauss,WeightGauss)

  do iGauss=1,NumGauss

! Obtain the shape functions and their derivatives.
    call fem_1d_shapefun (NumNodesElem,CoordGauss(iGauss),ShapeFun,ShapeDer)

! Global coordinates and their derivatives at the gauss points.
    do i=1,6
      r0_g(i) = dot_product (ShapeFun(1:NumNodesElem),r0(1:NumNodesElem,i))
      dr0_g(i)= dot_product (ShapeDer(1:NumNodesElem),r0(1:NumNodesElem,i))
    end do
    Jacobian =sqrt(dot_product (dr0_g(1:3),dr0_g(1:3)))

! Rescale the derivatives to be in the physical coordinates.
    ShapeDer=ShapeDer/Jacobian
    dr0_g=dr0_g/Jacobian

! Compute the current position vector and rotations, and their derivatives.
    do i=1,3
      dRa (i)=dot_product (ShapeDer(1:NumNodesElem),Ri(1:NumNodesElem,i))
      Psi (i)=dot_product (ShapeFun(1:NumNodesElem),Ri(1:NumNodesElem,i+3))
      dPsi(i)=dot_product (ShapeDer(1:NumNodesElem),Ri(1:NumNodesElem,i+3))
    end do

! Compute element shape function.
    N =0.d0
    dN=0.d0
    do i=1,6
      do j=1,NumNodesElem
        N (i,i+(j-1)*6)= ShapeFun(j)
        dN(i,i+(j-1)*6)= ShapeDer(j)
      end do
    end do

! Compute the current coordinate transformation matrix and rotational operator.
    CBa= rotvect_psi2mat  (Psi)
    Rot= rotvect_psi2rot  (Psi)

! Compute the current curvature and strain.
    gamma= matmul(CBa,dRa) - matmul(rotvect_psi2mat(r0_g(4:6)),dr0_g(1:3))
    KB=    matmul(Rot,dPsi)
    kappa= KB - matmul(rotvect_psi2rot(r0_g(4:6)),dr0_g(4:6))

    FB=matmul(ElemStiff(1:3,1:3),gamma)+matmul(ElemStiff(1:3,4:6),kappa)
    MB=matmul(ElemStiff(4:6,1:3),gamma)+matmul(ElemStiff(4:6,4:6),kappa)

! Operators for the tangent matrix.
    Aq0q0= 0.d0
    Aq0q0(4:6,4:6)= rotvect_A2(Psi,rot_cross(FB,gamma+Unit(1,1:3))+rot_cross(MB,KB)) &
&                 + rotvect_B2(Psi,dPsi,MB)                                          &
&                 + matmul(matmul(transpose(Rot),rot_skew(FB)),                      &
&                          matmul(rot_skew(gamma+Unit(1,1:3)),Rot))                  &
&                 + matmul(matmul(transpose(Rot),rot_skew(MB)),rotvect_A1(Psi,dPsi))

    Aq0q1= 0.d0
    Aq0q1(4:6,1:3)= matmul(transpose(Rot),matmul(rot_skew(FB),CBa))
    Aq0q1(4:6,4:6)= matmul(transpose(Rot),matmul(rot_skew(MB),Rot)) + rotvect_A2(Psi,MB)

    Aq1q0= 0.d0
    Aq1q0(1:3,4:6)= -matmul(transpose(CBa),matmul(rot_skew(FB),Rot))
    Aq1q0(4:6,4:6)= rotvect_A2(Psi,MB)

! Compute strain matrix operator and material tangent stiffness.
    Kgeom=Kgeom + WeightGauss(iGauss)*Jacobian                                   &
&               *(matmul(transpose(N), matmul((Aq0q0+transpose(Aq0q0))/2.d0,N))  &
&                +matmul(transpose(N), matmul((Aq0q1+transpose(Aq1q0))/2.d0,dN)) &
&                +matmul(transpose(dN),matmul((Aq1q0+transpose(Aq0q1))/2.d0,N)))

  end do
  deallocate (CoordGauss,WeightGauss)
  return
 end subroutine cbeam3_kgeom


 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-> Subroutine CBEAM3_FSTIF
!
!-> Description:
!
!    Compute element stiffness forces.
!
!-> Remarks.-
!
!   1) New values are added to those already in the vector.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine cbeam3_fstif (NumNodesElem,r0,Ri,ElemStiff,Qstiff,NumGauss)
  use lib_fem
  use lib_rot
  use lib_rotvect

! I/O Variables.
  integer,intent(in)   :: NumNodesElem       ! Number of nodes in the element.
  real(8),intent(in)   :: r0       (:,:)     ! Initial position/orientation of grid points.
  real(8),intent(in)   :: Ri       (:,:)     ! Current position/orientation of grid points.
  real(8),intent(in)   :: ElemStiff(:,:)     ! Stiffness properties in the element.
  real(8),intent(inout):: Qstiff   (:)       ! Discrete generalized stiffness forces.
  integer,intent(in)   :: NumGauss           ! Number of Gauss points in the element.

! Local variables.
  integer :: i,j                   ! Counters.
  integer :: iGauss                ! Counter on the Gaussian points.
  real(8) :: Jacobian              ! ds/deta, with eta the nondimensional arclength.
  real(8),allocatable :: CoordGauss (:) ! Coords of Gauss points.
  real(8),allocatable :: WeightGauss(:) ! Coords of Gauss points.
  real(8) :: ShapeFun(MaxNodCB3)    ! Shape functions in a gauss point.
  real(8) :: ShapeDer(MaxNodCB3)    ! Derivatives of ShapeFun in a gauss point.

  real(8) :: Ri_g (6)              ! Displacement/rotation at Gauss point.
  real(8) :: r0_g (6)              ! Position vector at the Gauss point.
  real(8) :: dRi_g(6)              ! Derivative of Ri_g.
  real(8) :: dr0_g(6)              ! Derivative of r0_g.

  real(8) :: CBa(3,3)              ! Coordinate transformation matrix.
  real(8) :: KB(3)                 ! Curvature in the deformed configuration at Gauss point.
  real(8) :: gamma(3)              ! Strain at Gauss point.
  real(8) :: Forces(6)             ! Forces and moments.
  real(8) :: Rot(3,3),dRot(3,3)    ! Tangential operator and spatial derivative.
  real(8) :: Strain(6)             ! Beam strains at the Gauss point.

  real(8) :: Ak(6,6),D(6,6)        ! A_K and D.
  real(8) :: B(6,6*MaxNodCB3)       ! Strain matrix operator.
  real(8) :: Yp(6,6),dYp(6,6)      ! Ypsilon and its derivative.
  real(8) :: N (6,6*MaxNodCB3)      ! Element shape function matrix
  real(8) :: dN(6,6*MaxNodCB3)      ! Derivatives of shape function.


! Define Gauss points and loop on them.
  allocate (CoordGauss(NumGauss))
  allocate (WeightGauss(NumGauss))
  call fem_1d_gauss_val (NumGauss,CoordGauss,WeightGauss)

  do iGauss=1,NumGauss

! Obtain the shape functions and their derivatives.
    call fem_1d_shapefun (NumNodesElem,CoordGauss(iGauss),ShapeFun,ShapeDer)

! Global coordinates and their derivatives at the gauss points.
    do i=1,6
      r0_g(i) = dot_product (ShapeFun(1:NumNodesElem),r0(1:NumNodesElem,i))
      dr0_g(i)= dot_product (ShapeDer(1:NumNodesElem),r0(1:NumNodesElem,i))
    end do
    Jacobian =sqrt(dot_product (dr0_g(1:3),dr0_g(1:3)))

! Rescale the derivatives to be in the physical coordinates.
    ShapeDer=ShapeDer/Jacobian
    dr0_g=dr0_g/Jacobian

! Compute the current position vector and rotations, and their derivatives.
    do i=1,6
      Ri_g(i) = dot_product (ShapeFun(1:NumNodesElem),Ri(1:NumNodesElem,i))
      dRi_g(i)= dot_product (ShapeDer(1:NumNodesElem),Ri(1:NumNodesElem,i))
    end do

! Compute element shape function.
    N =0.d0
    dN=0.d0
    do i=1,6
      do j=1,NumNodesElem
        N (i,i+(j-1)*6)= ShapeFun(j)
        dN(i,i+(j-1)*6)= ShapeDer(j)
      end do
    end do

! Compute the current coordinate transformation matrix, the rotational operator, and
! its spatial derivative at the Gauss point.
    CBa= rotvect_psi2mat  (Ri_g(4:6))
    Rot= rotvect_psi2rot  (Ri_g(4:6))
    dRot=rotvect_drotdpsib(Ri_g(4:6),dRi_g(4:6))

! Compute the current curvature and strain at the Gaussian point.
    gamma=matmul(CBa,dRi_g(1:3)) - matmul(rotvect_psi2mat(r0_g(4:6)),dr0_g(1:3))
    KB=matmul(Rot,dRi_g(4:6))

!Operators for the strain matrix.
    D=0.d0
    D(1:3,1:3)=transpose(CBa)
    D(4:6,4:6)=Unit

    Yp=0.d0
    Yp(1:3,1:3)=Unit
    Yp(4:6,4:6)=Rot

    dYp=0.d0
    dYp(4:6,4:6)=dRot

    AK=0.d0
    AK(4:6,1:3)=-rot_skew(gamma+Unit(1,1:3))
    AK(4:6,4:6)=-rot_skew(KB)

! Compute strain matrix operator and material tangent stiffness.
    B= matmul(transpose(D), matmul(Yp,dN)) &
&    + matmul(transpose(AK),matmul(Yp, N)) &
&    + matmul(dYp,N)

! Compute force/moment strains.
    Strain(1:3)= gamma
    Strain(4:6)= KB - matmul(rotvect_psi2rot(r0_g(4:6)),dr0_g(4:6))
    Forces=matmul(ElemStiff,Strain)

! Compute stiffness forces.
    Qstiff= Qstiff + WeightGauss(iGauss)*Jacobian &
&                  * matmul(transpose(B),Forces)

  end do

  deallocate (CoordGauss,WeightGauss)
  return
 end subroutine cbeam3_fstif



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-> Subroutine CBEAM3_KFINDIFF
!
!-> Description:
!
!    Compute tangent stiffness matrix using finite differences.
!
!-> Remarks.-
!
!   1) This routine is only used to test the implementation of the
!      analytical stiffness matrices.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine cbeam3_kfindiff (NumNodesElem,r0,Ri,ElemStiff,Kfindiff,NumGauss)
  use lib_fem
  use lib_rot
  use lib_rotvect

! I/O Variables.
  integer,intent(in) :: NumNodesElem    ! Number of nodes in the element.
  real(8),intent(in) :: r0    (:,:)     ! Initial position/orientation of grid points.
  real(8),intent(in) :: Ri    (:,:)     ! Current position/orientation of grid points.
  real(8),intent(in) :: ElemStiff(:,:)  ! Stiffness properties in the element.
  real(8),intent(out):: Kfindiff (:,:)  ! Geometric tangent stiffness matrix.
  integer,intent(in) :: NumGauss        ! Number of Gauss points in the element..

! Local variables.
  integer            :: i,j,k           ! Counters.
  real(8),parameter  :: Delta=1.d-7     ! Delta for finite-differences.
  real(8)            :: dRi(MaxNodCB3,6) ! Delta R.
  real(8),allocatable:: Q0(:),Q1(:)     ! Values of the element stiffness functional.

! Initialize.
  allocate (Q0(6*MaxNodCB3)); Q0=0.d0
  allocate (Q1(6*MaxNodCB3)); Q1=0.d0

! Compute force vector at reference point.
  call cbeam3_fstif (NumNodesElem,r0,Ri,ElemStiff,Q0,NumGauss)

! Loop in the element degrees of freedom.
  k=0
  do i=1,NumNodesElem
    do j=1,6
      k=k+1
      dRi=0.d0
      dRi(i,j)=Delta
      call cbeam3_fstif (NumNodesElem,r0,Ri+dRi,ElemStiff,Q1,NumGauss)
      Kfindiff(k,:)=(Q1-Q0)/Delta
    end do
  end do

  return
 end subroutine cbeam3_kfindiff



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-> Subroutine CBEAM3_MASS
!
!-> Description:
!
!    Compute the element tangent mass matrix.
!
!-> Remarks.-
!
!   1) New values are added to those already in the matrix.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine cbeam3_mass (NumNodesElem,r0,Ri,ElemMass,Mass,NumGauss)
  use lib_fem
  use lib_rot
  use lib_rotvect

! I/O Variables.
  integer,intent(in)   :: NumNodesElem      ! Number of nodes in the element.
  real(8),intent(in)   :: r0       (:,:)    ! Initial position/orientation of grid points.
  real(8),intent(in)   :: Ri       (:,:)    ! Current position/orientation of grid points.
  real(8),intent(in)   :: ElemMass (:,:)    ! Inertial properties in the element.
  real(8),intent(inout):: Mass     (:,:)    ! Tangent mass matrix.
  integer,intent(in)   :: NumGauss          ! Number of Gauss points in the element.

! Local variables.
  integer :: i,j                   ! Counters.
  integer :: iGauss                ! Counter on the Gaussian points.
  real(8) :: Jacobian              ! ds/deta, with eta the nondimensional arclength.
  real(8),allocatable :: CoordGauss (:) ! Coords of Gauss points.
  real(8),allocatable :: WeightGauss(:) ! Coords of Gauss points.
  real(8) :: ShapeFun(MaxNodCB3)    ! Shape functions in a gauss point.
  real(8) :: ShapeDer(MaxNodCB3)    ! Derivatives of ShapeFun in a gauss point.

  real(8) :: Ri_g (6)              ! Displacement/rotation at Gauss point.
  real(8) :: dr0_g(6)              ! Derivative of r0_g.
  real(8) :: CBa(3,3)              ! Coordinate transformation matrix.
  real(8) :: Rot(3,3)              ! Tangential operator.
  real(8) :: Yp(6,6),D(6,6)        ! Ypsilon and D.
  real(8) :: DTYpN(6,6*MaxNodCB3)  ! Strain matrix operator.
  real(8) :: N (6,6*MaxNodCB3)     ! Element shape function matrix

! Define Gauss points and loop on them.
  allocate (CoordGauss(NumGauss))
  allocate (WeightGauss(NumGauss))
  call fem_1d_gauss_val (NumGauss,CoordGauss,WeightGauss)

  do iGauss=1,NumGauss

! Obtain the shape functions and their derivatives.
    call fem_1d_shapefun (NumNodesElem,CoordGauss(iGauss),ShapeFun,ShapeDer)

! Global coordinates at the gauss points and element Jacobian.
    do i=1,6
      dr0_g(i)= dot_product (ShapeDer(1:NumNodesElem),r0(1:NumNodesElem,i))
    end do
    Jacobian =sqrt(dot_product (dr0_g(1:3),dr0_g(1:3)))

! Compute the current position vector and rotations.
    do i=1,6
      Ri_g(i) = dot_product (ShapeFun(1:NumNodesElem),Ri(1:NumNodesElem,i))
    end do

! Compute element shape function.
    N =0.d0
    do i=1,6
      do j=1,NumNodesElem
        N (i,i+(j-1)*6)= ShapeFun(j)
      end do
    end do

! Compute the current coordinate transformation matrix and rotational operator.
    CBa= rotvect_psi2mat  (Ri_g(4:6))
    Rot= rotvect_psi2rot  (Ri_g(4:6))

!Operators for the tangent mass matrix.
    D=0.d0
    D(1:3,1:3)=transpose(CBa)
    D(4:6,4:6)=Unit

    Yp=0.d0
    Yp(1:3,1:3)=Unit
    Yp(4:6,4:6)=Rot

    DTYpN=matmul(matmul(transpose(D),Yp),N)

! Compute mass tangent stiffness.
    Mass= Mass + WeightGauss(iGauss)*Jacobian &
&              * matmul(transpose(DTYpN),matmul(ElemMass,DTYpN))
  end do
  deallocate (CoordGauss,WeightGauss)

  return
 end subroutine cbeam3_mass


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-> Subroutine CBEAM3_RBMASS
!
!-> Description:
!
!    Compute the contribution of lumped masses to the element tangent mass matrix.
!
!-> Remarks.-
!
!   1) New values are added to those already in the matrix.
!   2) The lumped masses are defined by the inertia tensor of the nonstructural masses
!      in the section of the node, and expressed with respect to the local B frame.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine cbeam3_rbmass (NumNodesElem,r0,Ri,NodalMass,Mass)
  use lib_rot
  use lib_rotvect

! I/O Variables.
  integer,intent(in)   :: NumNodesElem      ! Number of nodes in the element.
  real(8),intent(in)   :: r0       (:,:)    ! Initial position/orientation of grid points.
  real(8),intent(in)   :: Ri       (:,:)    ! Current position/orientation of grid points.
  real(8),intent(in)   :: NodalMass(:,:,:)  ! Inertia tensor of lumped masses at element nodes.
  real(8),intent(inout):: Mass     (:,:)    ! Tangent mass matrix.

! Local variables.
  integer :: i,iNode               ! Counters.
  real(8) :: CBa(3,3)              ! Coordinate transformation matrix.
  real(8) :: Rot(3,3)              ! Tangential operator.
  real(8) :: Yp(6,6),D(6,6)        ! Ypsilon and D.
  real(8) :: DTYpN(6,6*MaxNodCB3)  ! Strain matrix operator.
  real(8) :: N (6,6*MaxNodCB3)     ! Element shape function matrix

! Loop in the nodes.
  do iNode=1,NumNodesElem

! Compute the nodal coordinate transformation matrix and rotational operator.
    CBa= rotvect_psi2mat  (Ri(iNode,4:6))
    Rot= rotvect_psi2rot  (Ri(iNode,4:6))

!Operators for the tangent mass matrix.
    D=0.d0
    D(1:3,1:3)=transpose(CBa)
    D(4:6,4:6)=Unit

    Yp=0.d0
    Yp(1:3,1:3)=Unit
    Yp(4:6,4:6)=Rot

! Compute element shape function at the current node.
    N =0.d0
    do i=1,6
        N (i,i+(iNode-1)*6)= 1.d0
    end do

! Compute mass tangent stiffness.
    DTYpN=matmul(matmul(transpose(D),Yp),N)
    Mass= Mass + matmul(transpose(DTYpN),matmul(NodalMass(iNode,:,:),DTYpN))
  end do

  return
 end subroutine cbeam3_rbmass



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-> Subroutine CBEAM3_CGYR
!
!-> Description:
!
!    Compute element tangent gyroscopic damping matrix
!
!-> Remarks.-
!
!   1) New values are added to those already in the matrix.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine cbeam3_cgyr (NumNodesElem,r0,Ri,RiDot,Vrel,ElemMass,Cgyr,NumGauss)
  use lib_fem
  use lib_rot
  use lib_rotvect

! I/O Variables.
  integer,intent(in)    :: NumNodesElem     ! Number of nodes in the element.
  real(8),intent(in)    :: r0       (:,:)   ! Initial position/orientation of grid points.
  real(8),intent(in)    :: Ri       (:,:)   ! Current position/orientation of grid points.
  real(8),intent(in)    :: RiDot    (:,:)   ! Current time derivative of position/CRV of points.
  real(8),intent(in)    :: Vrel     (6)     ! Linear/angular velocity of reference frame.
  real(8),intent(in)    :: ElemMass (6,6)   ! Inertia properties in the element.
  real(8),intent(inout) :: Cgyr     (:,:)   ! Tangent gyroscopic damping matrix.
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
  real(8) :: Ra(3),RaDot(3)         ! Displacement and its derivative at Gauss point.
  real(8) :: Psi(3),PsiDot(3)       ! Rotation vector and its derivative at Gauss point.
  real(8) :: VB(3),WB(3)            ! Linear and angular velocity in the material frame..
  real(8) :: PB(3),HB(3)            ! Linear and angular momenta in the material frame..

  real(8) :: CBa(3,3)               ! Coordinate transformation matrix.
  real(8) :: Rot(3,3)               ! Tangential operator.
  real(8) :: Yp(6,6),D(6,6)         ! Ypsilon and D.
  real(8) :: ADeltaW   (6,6)        ! A_delta_Omega matrix.
  real(8) :: Cbar      (6,6)        ! Cbar at current Gauss point.
  real(8) :: dVgyrdqdot(6,6)        ! d(Vdot)/d(qdot) at Gauss point.
  real(8) :: dVdqdot   (6,6)        ! dV/d(qdot) at Gauss point.
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
      RaDot (i)=dot_product (ShapeFun(1:NumNodesElem),RiDot(1:NumNodesElem,i))
      Psi   (i)=dot_product (ShapeFun(1:NumNodesElem),Ri   (1:NumNodesElem,i+3))
      PsiDot(i)=dot_product (ShapeFun(1:NumNodesElem),RiDot(1:NumNodesElem,i+3))
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
    VB= matmul(CBa,RaDot+rot_cross(Vrel(4:6),Ra)+Vrel(1:3))
    WB= matmul(Rot,PsiDot) + matmul(CBa,Vrel(4:6))

    PB= matmul(ElemMass(1:3,1:3),VB)+matmul(ElemMass(1:3,4:6),WB)
    HB= matmul(ElemMass(4:6,1:3),VB)+matmul(ElemMass(4:6,4:6),WB)

!Operators in Cbar.
    D=0.0d0
    D(1:3,1:3)=transpose(CBa)
    D(4:6,4:6)=Unit

    Yp=0.0d0
    Yp(1:3,1:3)=Unit
    Yp(4:6,4:6)=Rot

    dVdqdot=0.0d0
    dVdqdot(1:3,1:3)= CBa
    dVdqdot(4:6,4:6)= Rot

    dVgyrdqdot(1:3,1:3)= 2.d0*matmul(CBa,rot_skew(Vrel(4:6)))-matmul(rot_skew(WB),CBa)
    dVgyrdqdot(1:3,4:6)= matmul(rot_skew(VB),Rot)
    dVgyrdqdot(4:6,1:3)= 0.d0
    dVgyrdqdot(4:6,4:6)= 2.d0*rotvect_drotdpsib(Psi,PsiDot)+matmul(rot_skew(WB),Rot)

    ADeltaW(1:3,1:3)= rot_skew(WB)
    ADeltaW(1:3,4:6)= 0.d0
    ADeltaW(4:6,1:3)= rot_skew(VB)
    ADeltaW(4:6,4:6)= rot_skew(WB)

    ADeltaW=matmul(ADeltaW,ElemMass)

    AdeltaW(1:3,4:6)= AdeltaW(1:3,4:6) - rot_skew(PB)
    AdeltaW(4:6,1:3)= AdeltaW(4:6,1:3) - rot_skew(PB)
    AdeltaW(4:6,4:6)= AdeltaW(4:6,4:6) - rot_skew(HB)

! Operators for the tangent matrix.
    Cbar= matmul(matmul(transpose(Yp),D),      &
&                (matmul(ElemMass,dVgyrdqdot)  &
&                +matmul(ADeltaW,dVdqdot)))

! Compute strain matrix operator and material tangent stiffness.
    Cgyr= Cgyr+WeightGauss(iGauss)*Jacobian*matmul(transpose(N),matmul(Cbar,N))

  end do
  deallocate (CoordGauss,WeightGauss)
  return
 end subroutine cbeam3_cgyr

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-> Subroutine CBEAM3_RBCGYR
!
!-> Description:
!
!    Compute the contribution of lumped masses to the element tangent
!    gyroscopic matrix.
!
!-> Remarks.-
!
!   1) New values are added to those already in the matrix.
!   2) The lumped masses are defined by the inertia tensor of the nonstructural masses
!      in the section of the node, and expressed with respect to the local B frame.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine cbeam3_rbcgyr (NumNodesElem,r0,Ri,RiDot,Vrel,NodalMass,Cgyr)
  use lib_rot
  use lib_rotvect

! I/O Variables.
  integer,intent(in)    :: NumNodesElem     ! Number of nodes in the element.
  real(8),intent(in)    :: r0       (:,:)   ! Initial position/orientation of grid points.
  real(8),intent(in)    :: Ri       (:,:)   ! Current position/orientation of grid points.
  real(8),intent(in)    :: RiDot    (:,:)   ! Current time derivative of position/CRV of points.
  real(8),intent(in)    :: Vrel     (6)     ! Linear/angular velocity of reference frame.
  real(8),intent(in)    :: NodalMass(:,:,:) ! Inertia tensor of lumped masses at element nodes.
  real(8),intent(inout) :: Cgyr     (:,:)   ! Tangent gyroscopic damping matrix.

! Local variables
  integer :: i                      ! Counters.
  integer :: iNode                  ! Counter on the element nodes.

  real(8) :: Ra(3),RaDot(3)         ! Displacement and its derivative at Gauss point.
  real(8) :: Psi(3),PsiDot(3)       ! Rotation vector and its derivative at Gauss point.
  real(8) :: VB(3),WB(3)            ! Linear and angular velocity in the material frame.
  real(8) :: PB(3),HB(3)            ! Linear and angular momenta in the material frame.

  real(8) :: CBa(3,3)               ! Coordinate transformation matrix.
  real(8) :: Rot(3,3)               ! Tangential operator.
  real(8) :: Yp(6,6),D(6,6)         ! Ypsilon and D.
  real(8) :: ADeltaW   (6,6)        ! A_delta_Omega matrix.
  real(8) :: Cbar      (6,6)        ! Cbar at current Gauss point.
  real(8) :: dVgyrdqdot(6,6)        ! d(Vdot)/d(qdot) at Gauss point.
  real(8) :: dVdqdot   (6,6)        ! dV/d(qdot) at Gauss point.
  real(8) :: N (6,6*MaxNodCB3)      ! Element shape function matrix.

! Loop in the element ndoes.
  do iNode=1,NumNodesElem

! Compute the current position vector and rotations, and their derivatives.
		Ra    = Ri   (iNode,1:3)
		RaDot = RiDot(iNode,1:3)
		Psi   = Ri   (iNode,4:6)
		PsiDot= RiDot(iNode,4:6)

! Compute element shape functions at the corrent node.
    N =0.d0
    do i=1,6
      N (i,i+(iNode-1)*6)= 1.d0
    end do

! Compute the current coordinate transformation matrix and rotational operator.
    CBa= rotvect_psi2mat  (Psi)
    Rot= rotvect_psi2rot  (Psi)

! Compute the current linear and angular velocity.
    VB= matmul(CBa,RaDot+rot_cross(Vrel(4:6),Ra)+Vrel(1:3))
    WB= matmul(Rot,PsiDot) + matmul(CBa,Vrel(4:6))

    PB= matmul(NodalMass(iNode,1:3,1:3),VB)+matmul(NodalMass(iNode,1:3,4:6),WB)
    HB= matmul(NodalMass(iNode,4:6,1:3),VB)+matmul(NodalMass(iNode,4:6,4:6),WB)

!Operators in Cbar.
    D=0.0d0
    D(1:3,1:3)=transpose(CBa)
    D(4:6,4:6)=Unit

    Yp=0.0d0
    Yp(1:3,1:3)=Unit
    Yp(4:6,4:6)=Rot

    dVdqdot=0.0d0
    dVdqdot(1:3,1:3)= CBa
    dVdqdot(4:6,4:6)= Rot

    dVgyrdqdot(1:3,1:3)= 2.d0*matmul(CBa,rot_skew(Vrel(4:6)))-matmul(rot_skew(WB),CBa)
    dVgyrdqdot(1:3,4:6)= matmul(rot_skew(VB),Rot)
    dVgyrdqdot(4:6,1:3)= 0.d0
    dVgyrdqdot(4:6,4:6)= 2.d0*rotvect_drotdpsib(Psi,PsiDot)+matmul(rot_skew(WB),Rot)

    ADeltaW(1:3,1:3)= rot_skew(WB)
    ADeltaW(1:3,4:6)= 0.d0
    ADeltaW(4:6,1:3)= rot_skew(VB)
    ADeltaW(4:6,4:6)= rot_skew(WB)

    ADeltaW=matmul(ADeltaW,NodalMass(iNode,:,:))

    AdeltaW(1:3,4:6)= AdeltaW(1:3,4:6) - rot_skew(PB)
    AdeltaW(4:6,1:3)= AdeltaW(4:6,1:3) - rot_skew(PB)
    AdeltaW(4:6,4:6)= AdeltaW(4:6,4:6) - rot_skew(HB)

! Operators for the tangent matrix.
    Cbar= matmul(matmul(transpose(Yp),D),      &
&                (matmul(NodalMass(iNode,:,:),dVgyrdqdot)  &
&                +matmul(ADeltaW,dVdqdot)))

! Compute strain matrix operator and material tangent stiffness.
    Cgyr= Cgyr + matmul(transpose(N),matmul(Cbar,N))

  end do
 return
 end subroutine cbeam3_rbcgyr

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-> Subroutine CBEAM3_KGYR
!
!-> Description:
!
!    Compute tangent gyroscopic stiffness matrix
!
!-> Remarks.-
!
!   1) New values are added to those already in the matrix.
!
!   2) The contribution due to the derivative of the mass matrix with the
!      degrees of freedom are included (thus the dependency with the
!      accelerations).
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine cbeam3_kgyr (NumNodesElem,r0,Ri,RiDot,RiDDot,Vrel,VrelDot, &
&                                  ElemMass,Kgyr,NumGauss)
  use lib_fem
  use lib_rot
  use lib_rotvect

! I/O Variables.
  integer,intent(in)    :: NumNodesElem     ! Number of nodes in the element.
  real(8),intent(in)    :: r0       (:,:)   ! Initial position/orientation of grid points.
  real(8),intent(in)    :: Ri       (:,:)   ! Current position/orientation of grid points.
  real(8),intent(in)    :: RiDot    (:,:)   ! Current time derivative of position/CRV of points.
  real(8),intent(in)    :: RiDDot   (:,:)   ! Second time derivative of Ri.
  real(8),intent(in)    :: Vrel     (6)     ! Linear/angular velocity of reference frame.
  real(8),intent(in)    :: VrelDot  (6)     ! Currnet time derivative of Vrel.
  real(8),intent(in)    :: ElemMass (6,6)   ! Inertia properties in the element.
  real(8),intent(inout) :: Kgyr     (:,:)   ! Tangent gyroscopic stiffness matrix.
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
  real(8) :: Ra(3),Psi(3)           ! Displacement and Rotation vector at Gauss point.
  real(8) :: RaDot(3),PsiDot(3)     ! Time derivatives of Ra and Psi.
  real(8) :: RaDDot(3),PsiDDot(3)   ! Time derivatives of RaDot and PsiDot.
  real(8) :: RotPsiDot(3)           ! Rot*PsiDot.
  real(8) :: VB(3),WB(3)            ! Linear and angular velocity in the material frame.
  real(8) :: VBdot(3),WBdot(3)      ! Time derivatives of VB and WB.
  real(8) :: Vhat(6),Phat(6)        ! Vectors with linear/angular velocities/momenta.

  real(8) :: CBa(3,3)               ! Coordinate transformation matrix.
  real(8) :: Rot(3,3),RotDot(3,3)   ! Tangential operator and time derivative.
  real(8) :: Yp(6,6),D(6,6)         ! Ypsilon and D.
  real(8) :: ADeltaW (6,6)          ! A_delta_Omega matrix.
  real(8) :: AWW      (6,6)         ! A_Omega_Omega matrix.
  real(8) :: Kbar     (6,6)         ! Kbar at current Gauss point.
  real(8) :: dVdotdq  (6,6)         ! d(Vdot)/dq at Gauss point.
  real(8) :: dVdq     (6,6)         ! dV/dq at Gauss point.
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
      Ra     (i)=dot_product (ShapeFun(1:NumNodesElem),Ri    (1:NumNodesElem,i))
      RaDot  (i)=dot_product (ShapeFun(1:NumNodesElem),RiDot (1:NumNodesElem,i))
      RaDDot (i)=dot_product (ShapeFun(1:NumNodesElem),RiDDot(1:NumNodesElem,i))
      Psi    (i)=dot_product (ShapeFun(1:NumNodesElem),Ri    (1:NumNodesElem,i+3))
      PsiDot (i)=dot_product (ShapeFun(1:NumNodesElem),RiDot (1:NumNodesElem,i+3))
      PsiDDot(i)=dot_product (ShapeFun(1:NumNodesElem),RiDDot(1:NumNodesElem,i+3))
    end do

! Compute element shape functions.
    N =0.d0
    do i=1,6
      do j=1,NumNodesElem
        N (i,i+(j-1)*6)= ShapeFun(j)
      end do
    end do

! Compute the current coordinate transformation matrix and rotational operator.
    CBa= rotvect_psi2mat(Psi)
    Rot= rotvect_psi2rot(Psi)
    RotDot=rotvect_drotdpsib(Psi,PsiDot)
    RotPsiDot=matmul(Rot,PsiDot)

! Compute the current linear and angular velocity.
    VB= matmul(CBa,RaDot+rot_cross(Vrel(4:6),Ra)+Vrel(1:3))
    WB= RotPsiDot + matmul(CBa,Vrel(4:6))

    Phat= matmul(ElemMass(:,1:3),VB)+matmul(ElemMass(:,4:6),WB)

    VBdot= matmul(CBa,RaDDot+rot_cross(Vrel(4:6),RaDot)+rot_cross(VrelDot(4:6),Ra)+VrelDot(1:3)) &
&        + rot_cross(VB,RotPsiDot)
    WBdot= matmul(Rot,PsiDDot) + matmul(RotDot,PsiDot) + rot_cross(WB,RotPsiDot) &
&        + matmul(CBa,VrelDot(4:6))

! Operators in Kbar from the derivatives of velocities.
    D=0.0d0
    D(1:3,1:3)=transpose(CBa)
    D(4:6,4:6)=Unit

    Yp=0.0d0
    Yp(1:3,1:3)=Unit
    Yp(4:6,4:6)=Rot

    AWW(1:3,1:3)=rot_skew(WB)
    AWW(1:3,4:6)=0.d0
    AWW(4:6,1:3)=rot_skew(VB)
    AWW(4:6,4:6)=rot_skew(WB)
    
    dVdq(1:3,1:3)= matmul(CBa,rot_skew(Vrel(4:6)))
    dVdq(1:3,4:6)= matmul(rot_skew(VB),Rot)
    dVdq(4:6,1:3)= 0.0d0
    dVdq(4:6,4:6)= RotDot + matmul(rot_skew(WB),Rot)

    dVdotdq(1:3,1:3)= matmul(CBa,rot_skew(VrelDot(4:6))) &
&                   - matmul(rot_skew(RotPsiDot),matmul(CBa,rot_skew(Vrel(4:6))))
    dVdotdq(1:3,4:6)= matmul(rot_skew(VBdot),Rot) + matmul(rot_skew(VB),RotDot)
    dVdotdq(4:6,1:3)= 0.d0
    dVdotdq(4:6,4:6)= rotvect_ddrotdpsib(Psi,PsiDot,PsiDot) &
&                   + rotvect_drotdpsib (Psi,PsiDDot)       &
&                   + matmul(rot_skew(WB),RotDot)           &
&                   + matmul(rot_skew(WBdot),Rot)

    ADeltaW=matmul(AWW,ElemMass)

    AdeltaW(1:3,4:6)= AdeltaW(1:3,4:6) - rot_skew(Phat(1:3))
    AdeltaW(4:6,1:3)= AdeltaW(4:6,1:3) - rot_skew(Phat(1:3))
    AdeltaW(4:6,4:6)= AdeltaW(4:6,4:6) - rot_skew(Phat(4:6))

! Operators in Kbar from the derivatives of matrices D and T.
    
    Vhat(1:3)=VB
    Vhat(4:6)=WB
    Phat=matmul(AWW,matmul(ElemMass,Vhat))

    Vhat(1:3)=VBdot
    Vhat(4:6)=WBdot
    Phat= Phat + matmul(ElemMass,Vhat)
    
! Compute Kbar for the tangent inertial stiffness at the current Gauss point.
    Kbar= 0.d0
    Kbar(1:3,4:6)=-matmul(transpose(CBa),matmul(rot_skew(Phat(1:3)),Rot))
    Kbar(4:6,4:6)= rotvect_a2(Psi,Phat(4:6))
    Kbar= Kbar + matmul(matmul(transpose(Yp),D),(matmul(ElemMass,dVdotdq)   &
&                                               +matmul( ADeltaW,dVdq   )))

! Add the contribution of the current Gauss point to inertial stiffness matrix.
    Kgyr= Kgyr + WeightGauss(iGauss) * Jacobian &
&       * matmul(transpose(N),matmul((Kbar+transpose(Kbar))/2.d0,N))

  end do
  deallocate (CoordGauss,WeightGauss)
  return
 end subroutine cbeam3_kgyr


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-> Subroutine CBEAM3_RBKGYR
!
!-> Description:
!
!    Compute lumped masses contribution to the tangent gyroscopic stiffness
!    matrix
!
!-> Remarks.-
!
!   1) New values are added to those already in the matrix.
!   2) The contribution due to the derivative of the mass matrix with the
!      degrees of freedom are included (thus the dependency with the
!      accelerations).
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine cbeam3_rbkgyr (NumNodesElem,r0,Ri,RiDot,RiDDot,Vrel,VrelDot, &
&                          NodalMass,Kgyr)
  use lib_rot
  use lib_rotvect

! I/O Variables.
  integer,intent(in)    :: NumNodesElem     ! Number of nodes in the element.
  real(8),intent(in)    :: r0       (:,:)   ! Initial position/orientation of grid points.
  real(8),intent(in)    :: Ri       (:,:)   ! Current position/orientation of grid points.
  real(8),intent(in)    :: RiDot    (:,:)   ! Current time derivative of position/CRV of points.
  real(8),intent(in)    :: RiDDot   (:,:)   ! Second time derivative of Ri.
  real(8),intent(in)    :: Vrel     (6)     ! Linear/angular velocity of reference frame.
  real(8),intent(in)    :: VrelDot  (6)     ! Currnet time derivative of Vrel.
  real(8),intent(in)    :: NodalMass(:,:,:) ! Inertia tensor of lumped masses at element nodes.
  real(8),intent(inout) :: Kgyr     (:,:)   ! Tangent gyroscopic stiffness matrix.

! Local variables
  integer :: i,iNode                ! Counters.
  real(8) :: Ra(3),Psi(3)           ! Displacement and Rotation vector at Gauss point.
  real(8) :: RaDot(3),PsiDot(3)     ! Time derivatives of Ra and Psi.
  real(8) :: RaDDot(3),PsiDDot(3)   ! Time derivatives of RaDot and PsiDot.
  real(8) :: RotPsiDot(3)           ! Rot*PsiDot.
  real(8) :: VB(3),WB(3)            ! Linear and angular velocity in the material frame.
  real(8) :: VBdot(3),WBdot(3)      ! Time derivatives of VB and WB.
  real(8) :: Vhat(6),Phat(6)        ! Vectors with linear/angular velocities/momenta.

  real(8) :: CBa(3,3)               ! Coordinate transformation matrix.
  real(8) :: Rot(3,3),RotDot(3,3)   ! Tangential operator and time derivative.
  real(8) :: Yp(6,6),D(6,6)         ! Ypsilon and D.
  real(8) :: ADeltaW (6,6)          ! A_delta_Omega matrix.
  real(8) :: AWW      (6,6)         ! A_Omega_Omega matrix.
  real(8) :: Kbar     (6,6)         ! Kbar at current Gauss point.
  real(8) :: dVdotdq  (6,6)         ! d(Vdot)/dq at Gauss point.
  real(8) :: dVdq     (6,6)         ! dV/dq at Gauss point.
  real(8) :: N (6,6*MaxNodCB3)      ! Element shape function matrix.

! Loop in the element ndoes.
  do iNode=1,NumNodesElem

! Compute the current position vector and rotations, and their derivatives.
		Ra     = Ri    (iNode,1:3)
		RaDot  = RiDot (iNode,1:3)
		RaDDot = RiDDot(iNode,1:3)
		Psi    = Ri    (iNode,4:6)
		PsiDot = RiDot (iNode,4:6)
		PsiDDot= RiDDot(iNode,4:6)

! Compute element shape functions at the corrent node.
    N =0.d0
    do i=1,6
      N (i,i+(iNode-1)*6)= 1.d0
    end do

! Compute the current coordinate transformation matrix and rotational operator.
    CBa= rotvect_psi2mat(Psi)
    Rot= rotvect_psi2rot(Psi)
    RotDot=rotvect_drotdpsib(Psi,PsiDot)
    RotPsiDot=matmul(Rot,PsiDot)

! Compute the current linear and angular velocity.
    VB= matmul(CBa,RaDot+rot_cross(Vrel(4:6),Ra)+Vrel(1:3))
    WB= RotPsiDot + matmul(CBa,Vrel(4:6))

    Phat= matmul(NodalMass(iNode,:,1:3),VB)+matmul(NodalMass(iNode,:,4:6),WB)

    VBdot= matmul(CBa,RaDDot+rot_cross(Vrel(4:6),RaDot)+rot_cross(VrelDot(4:6),Ra)+VrelDot(1:3)) &
&        + rot_cross(VB,RotPsiDot)
    WBdot= matmul(Rot,PsiDDot) + matmul(RotDot,PsiDot) + rot_cross(WB,RotPsiDot) &
&        + matmul(CBa,VrelDot(4:6))

! Operators in Kbar from the derivatives of velocities.
    D=0.0d0
    D(1:3,1:3)=transpose(CBa)
    D(4:6,4:6)=Unit

    Yp=0.0d0
    Yp(1:3,1:3)=Unit
    Yp(4:6,4:6)=Rot

    AWW(1:3,1:3)=rot_skew(WB)
    AWW(1:3,4:6)=0.d0
    AWW(4:6,1:3)=rot_skew(VB)
    AWW(4:6,4:6)=rot_skew(WB)

    dVdq(1:3,1:3)= matmul(CBa,rot_skew(Vrel(4:6)))
    dVdq(1:3,4:6)= matmul(rot_skew(VB),Rot)
    dVdq(4:6,1:3)= 0.0d0
    dVdq(4:6,4:6)= RotDot + matmul(rot_skew(WB),Rot)

    dVdotdq(1:3,1:3)= matmul(CBa,rot_skew(VrelDot(4:6))) &
&                   - matmul(rot_skew(RotPsiDot),matmul(CBa,rot_skew(Vrel(4:6))))
    dVdotdq(1:3,4:6)= matmul(rot_skew(VBdot),Rot) + matmul(rot_skew(VB),RotDot)
    dVdotdq(4:6,1:3)= 0.d0
    dVdotdq(4:6,4:6)= rotvect_ddrotdpsib(Psi,PsiDot,PsiDot) &
&                   + rotvect_drotdpsib (Psi,PsiDDot)       &
&                   + matmul(rot_skew(WB),RotDot)           &
&                   + matmul(rot_skew(WBdot),Rot)

    ADeltaW=matmul(AWW,NodalMass(iNode,:,:))

    AdeltaW(1:3,4:6)= AdeltaW(1:3,4:6) - rot_skew(Phat(1:3))
    AdeltaW(4:6,1:3)= AdeltaW(4:6,1:3) - rot_skew(Phat(1:3))
    AdeltaW(4:6,4:6)= AdeltaW(4:6,4:6) - rot_skew(Phat(4:6))

! Operators in Kbar from the derivatives of matrices D and T.

    Vhat(1:3)=VB
    Vhat(4:6)=WB
    Phat=matmul(AWW,matmul(NodalMass(iNode,:,:),Vhat))

    Vhat(1:3)=VBdot
    Vhat(4:6)=WBdot
    Phat= Phat + matmul(NodalMass(iNode,:,:),Vhat)

! Compute Kbar for the tangent inertial stiffness at the current Gauss point.
    Kbar= 0.d0
    Kbar(1:3,4:6)=-matmul(transpose(CBa),matmul(rot_skew(Phat(1:3)),Rot))
    Kbar(4:6,4:6)= rotvect_a2(Psi,Phat(4:6))
    Kbar= Kbar + matmul(matmul(transpose(Yp),D),(matmul(NodalMass(iNode,:,:),dVdotdq)   &
&                                               +matmul(ADeltaW,dVdq )))

! Add the contribution of the current Gauss point to inertial stiffness matrix.
    Kgyr= Kgyr + matmul(transpose(N),matmul((Kbar+transpose(Kbar))/2.d0,N))
  end do

	return
 end subroutine cbeam3_rbkgyr


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-> Subroutine CBEAM3_FGYR
!
!-> Description:
!
!    Compute generalized gyroscopic forces.
!
!-> Remarks.-
!
!   1) New values are added to those already in the matrix.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine cbeam3_fgyr (NumNodesElem,r0,Ri,RiDot,Vrel,ElemMass,Qgyr,NumGauss)
  use lib_fem
  use lib_rot
  use lib_rotvect

! I/O Variables.
  integer,intent(in)    :: NumNodesElem     ! Number of nodes in the element.
  real(8),intent(in)    :: r0       (:,:)   ! Initial position/orientation of grid points.
  real(8),intent(in)    :: Ri       (:,:)   ! Current position/orientation of grid points.
  real(8),intent(in)    :: RiDot    (:,:)   ! Current time derivative of position/CRV of points.
  real(8),intent(in)    :: Vrel     (6)     ! Linear/angular velocity of reference frame.
  real(8),intent(in)    :: ElemMass (6,6)   ! Inertia properties in the element.
  real(8),intent(inout) :: Qgyr     (:)     ! Vector of discrete gyroscopic forces.
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
  real(8) :: Ra(3),Psi(3)           ! Displacement and Rotation vector at Gauss point.
  real(8) :: RaDot(3),PsiDot(3)     ! Time derivatives of Ra and Psi.
  real(8) :: VB(3),WB(3)            ! Linear and angular velocity in the material frame.
  real(8)  ::VBgyr(3),WBgyr(3)      ! Time derivatives of VB and WB.
  real(8) :: Vhat(6),Phat(6)        ! Vectors with linear/angular velocities/momenta.

  real(8) :: CBa(3,3)               ! Coordinate transformation matrix.
  real(8) :: Rot(3,3),RotDot(3,3)   ! Tangential operator and time derivative.
  real(8) :: Yp(6,6),D(6,6)         ! Ypsilon and D.
  real(8) :: AWW      (6,6)         ! A_Omega_Omega matrix.
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
      Ra     (i)=dot_product (ShapeFun(1:NumNodesElem),Ri    (1:NumNodesElem,i))
      RaDot  (i)=dot_product (ShapeFun(1:NumNodesElem),RiDot (1:NumNodesElem,i))
      Psi    (i)=dot_product (ShapeFun(1:NumNodesElem),Ri    (1:NumNodesElem,i+3))
      PsiDot (i)=dot_product (ShapeFun(1:NumNodesElem),RiDot (1:NumNodesElem,i+3))
    end do

! Compute element shape functions.
    N =0.d0
    do i=1,6
      do j=1,NumNodesElem
        N (i,i+(j-1)*6)= ShapeFun(j)
      end do
    end do

! Compute the current coordinate transformation matrix and rotational operator.
    CBa= rotvect_psi2mat(Psi)
    Rot= rotvect_psi2rot(Psi)
    RotDot=rotvect_drotdpsib(Psi,PsiDot)

! Compute the current linear and angular velocity.
    VB= matmul(CBa,RaDot+rot_cross(Vrel(4:6),Ra)+Vrel(1:3))
    WB= matmul(Rot,PsiDot) + matmul(CBa,Vrel(4:6))

    VBgyr= matmul(CBa,rot_cross(Vrel(4:6),RaDot))+ rot_cross(VB,WB-matmul(CBa,Vrel(4:6)))
    WBgyr= matmul(RotDot,PsiDot) - rot_cross(WB,matmul(CBa,Vrel(4:6)))

! Operators in Kbar from the derivatives of velocities.
    D=0.0d0
    D(1:3,1:3)=transpose(CBa)
    D(4:6,4:6)=Unit

    Yp=0.0d0
    Yp(1:3,1:3)=Unit
    Yp(4:6,4:6)=Rot

    AWW(1:3,1:3)=rot_skew(WB)
    AWW(1:3,4:6)=0.d0
    AWW(4:6,1:3)=rot_skew(VB)
    AWW(4:6,4:6)=rot_skew(WB)

! Operators in Kbar from the derivatives of matrices D and T.
    Vhat(1:3)=VB
    Vhat(4:6)=WB
    Phat=matmul(AWW,matmul(ElemMass,Vhat))

    Vhat(1:3)=VBgyr
    Vhat(4:6)=WBgyr
    Phat= Phat + matmul(ElemMass,Vhat)

! Add the contribution of the current Gauss point to inertial force.
    Qgyr= Qgyr + WeightGauss(iGauss) * Jacobian &
&        * matmul(transpose(N),matmul(matmul(transpose(Yp),D),Phat))

  end do
  deallocate (CoordGauss,WeightGauss)
  return
 end subroutine cbeam3_fgyr


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-> Subroutine CBEAM3_RBFGYR
!
!-> Description:
!
!    Compute generalized gyroscopic forces.
!
!-> Remarks.-
!
!   1) New values are added to those already in the matrix.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine cbeam3_rbfgyr (NumNodesElem,r0,Ri,RiDot,Vrel,NodalMass,Qgyr)
  use lib_rot
  use lib_rotvect

! I/O Variables.
  integer,intent(in)   :: NumNodesElem     ! Number of nodes in the element.
  real(8),intent(in)   :: r0       (:,:)   ! Initial position/orientation of grid points.
  real(8),intent(in)   :: Ri       (:,:)   ! Current position/orientation of grid points.
  real(8),intent(in)   :: RiDot    (:,:)   ! Current time derivative of position/CRV of points.
  real(8),intent(in)   :: Vrel     (6)     ! Linear/angular velocity of reference frame.
  real(8),intent(in)   :: NodalMass(:,:,:) ! Inertia tensor of lumped masses at element nodes.
  real(8),intent(inout):: Qgyr     (:)     ! Vector of discrete gyroscopic forces.

! Local variables
  integer :: i                      ! Counters.
  integer :: iNode                  ! Counter on the element nodes.
  real(8) :: Ra(3),Psi(3)           ! Displacement and Rotation vector at Gauss point.
  real(8) :: RaDot(3),PsiDot(3)     ! Time derivatives of Ra and Psi.
  real(8) :: VB(3),WB(3)            ! Linear and angular velocity in the material frame.
  real(8)  ::VBgyr(3),WBgyr(3)      ! Time derivatives of VB and WB.
  real(8) :: Vhat(6),Phat(6)        ! Vectors with linear/angular velocities/momenta.

  real(8) :: CBa(3,3)               ! Coordinate transformation matrix.
  real(8) :: Rot(3,3),RotDot(3,3)   ! Tangential operator and time derivative.
  real(8) :: Yp(6,6),D(6,6)         ! Ypsilon and D.
  real(8) :: AWW      (6,6)         ! A_Omega_Omega matrix.
  real(8) :: N (6,6*MaxNodCB3)      ! Element shape function matrix.

! Loop on the nodes in the current elemetn.
  do iNode=1,NumNodesElem

! Compute the current position vector and rotations, and their derivatives.
    Ra    = Ri   (iNode,1:3)
    Psi   = Ri   (iNode,4:6)
    RaDot = RiDot(iNode,1:3)
    PsiDot= RiDot(iNode,4:6)

! Compute element shape functions.
    N =0.d0
    do i=1,6
      N (i,i+(iNode-1)*6)= 1.d0
    end do

! Compute the current coordinate transformation matrix, the rotational operator and
! its derivative.
    CBa= rotvect_psi2mat(Psi)
    Rot= rotvect_psi2rot(Psi)
    RotDot=rotvect_drotdpsib(Psi,PsiDot)

! Compute the current linear and angular velocity (note that the relations between
! Wb and Rot habe been used in Wgyr).
    VB= matmul(CBa,RaDot+rot_cross(Vrel(4:6),Ra)+Vrel(1:3))
    WB= matmul(Rot,PsiDot) + matmul(CBa,Vrel(4:6))

    VBgyr= matmul(CBa,rot_cross(Vrel(4:6),RaDot))+ rot_cross(VB,WB-matmul(CBa,Vrel(4:6)))
    WBgyr= matmul(RotDot,PsiDot) - rot_cross(WB,matmul(CBa,Vrel(4:6)))

! Operators in Kbar from the derivatives of velocities.
    D=0.0d0
    D(1:3,1:3)=transpose(CBa)
    D(4:6,4:6)=Unit

    Yp=0.0d0
    Yp(1:3,1:3)=Unit
    Yp(4:6,4:6)=Rot

    AWW(1:3,1:3)=rot_skew(WB)
    AWW(1:3,4:6)=0.d0
    AWW(4:6,1:3)=rot_skew(VB)
    AWW(4:6,4:6)=rot_skew(WB)

! Contribution from the time derivatives of matrices D and T.
    Vhat(1:3)=VB
    Vhat(4:6)=WB
    Phat=matmul(AWW,matmul(NodalMass(iNode,:,:),Vhat))

    Vhat(1:3)=VBgyr
    Vhat(4:6)=WBgyr
    Phat= Phat + matmul(NodalMass(iNode,:,:),Vhat)

! Add the contribution of the current node to inertial force.
    Qgyr= Qgyr + matmul(transpose(N),matmul(matmul(transpose(Yp),D),Phat))

  end do

  return
 end subroutine cbeam3_rbfgyr


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-> Subroutine CBEAM3_MVEL
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
 subroutine cbeam3_mvel (NumNodesElem,r0,Ri,ElemMass,Mvel,NumGauss)
  use lib_fem
  use lib_rot
  use lib_rotvect

! I/O Variables.
  integer,intent(in)   :: NumNodesElem      ! Number of nodes in the element.
  real(8),intent(in)   :: r0       (:,:)    ! Initial position/orientation of grid points.
  real(8),intent(in)   :: Ri       (:,:)    ! Current position/orientation of grid points.
  real(8),intent(in)   :: ElemMass (:,:)    ! Inertial properties in the element.
  real(8),intent(inout):: Mvel     (:,:)    ! Tangent mass matrix.
  integer,intent(in)   :: NumGauss          ! Number of Gauss points in the element.

! Local variables.
  integer :: i,j                        ! Counters.
  integer :: iGauss                     ! Counter on the Gaussian points.
  real(8) :: Jacobian                   ! ds/deta, with eta the nondimensional arclength.
  real(8),allocatable :: CoordGauss (:) ! Coords of Gauss points.
  real(8),allocatable :: WeightGauss(:) ! Coords of Gauss points.
  real(8) :: ShapeFun(MaxNodCB3)        ! Shape functions in a gauss point.
  real(8) :: ShapeDer(MaxNodCB3)        ! Derivatives of ShapeFun in a gauss point.

  real(8) :: Ri_g (6)              ! Displacement/rotation at Gauss point.
  real(8) :: dr0_g(6)              ! Derivative of r0_g.
  real(8) :: CBa(3,3)              ! Coordinate transformation matrix.
  real(8) :: Rot(3,3)              ! Tangential operator.
  real(8) :: Yp(6,6),D(6,6)        ! Ypsilon and D.
  real(8) :: dVdotdvrel(6,6)       ! d(Vdot)/d(vreldot)
  real(8) :: DTYpN(6,6*MaxNodCB3)  ! Strain matrix operator.
  real(8) :: N (6,6*MaxNodCB3)     ! Element shape function matrix

! Define Gauss points and loop on them.
  allocate (CoordGauss(NumGauss))
  allocate (WeightGauss(NumGauss))
  call fem_1d_gauss_val (NumGauss,CoordGauss,WeightGauss)

  do iGauss=1,NumGauss

! Obtain the shape functions and their derivatives.
    call fem_1d_shapefun (NumNodesElem,CoordGauss(iGauss),ShapeFun,ShapeDer)

! Global coordinates at the gauss points and element Jacobian.
    do i=1,6
      dr0_g(i)= dot_product (ShapeDer(1:NumNodesElem),r0(1:NumNodesElem,i))
    end do
    Jacobian =sqrt(dot_product (dr0_g(1:3),dr0_g(1:3)))

! Compute the current position vector and rotations.
    do i=1,6
      Ri_g(i) = dot_product (ShapeFun(1:NumNodesElem),Ri(1:NumNodesElem,i))
    end do

! Compute element shape function.
    N =0.d0
    do i=1,6
      do j=1,NumNodesElem
        N (i,i+(j-1)*6)= ShapeFun(j)
      end do
    end do

! Compute the current coordinate transformation matrix and rotational operator.
    CBa= rotvect_psi2mat  (Ri_g(4:6))
    Rot= rotvect_psi2rot  (Ri_g(4:6))

!Operators for the tangent mass matrix.
    D=0.d0
    D(1:3,1:3)=transpose(CBa)
    D(4:6,4:6)=Unit

    Yp=0.d0
    Yp(1:3,1:3)=Unit
    Yp(4:6,4:6)=Rot

    DTYpN=matmul(matmul(transpose(D),Yp),N)

    dVdotdvrel=0.d0
    dVdotdvrel(1:3,1:3)= CBa
    dVdotdvrel(1:3,4:6)=-matmul(CBa,rot_skew(Ri_g(1:3)))
    dVdotdvrel(4:6,4:6)= CBa

! Compute mass tangent stiffness.
    Mvel= Mvel + WeightGauss(iGauss)*Jacobian &
&              * matmul(transpose(DTYpN),matmul(ElemMass,dVdotdvrel))
  end do
  deallocate (CoordGauss,WeightGauss)

  return
 end subroutine cbeam3_mvel

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-> Subroutine CBEAM3_RBMVEL
!
!-> Description:
!
!    Compute the contribution of lumped masses to the mass matrix due to
!    the reference system accelerations.
!
!-> Remarks.-
!
!   1) New values are added to those already in the matrix.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine cbeam3_rbmvel (NumNodesElem,r0,Ri,NodalMass,Mvel)
  use lib_rot
  use lib_rotvect

! I/O Variables.
  integer,intent(in)   :: NumNodesElem      ! Number of nodes in the element.
  real(8),intent(in)   :: r0       (:,:)    ! Initial position/orientation of grid points.
  real(8),intent(in)   :: Ri       (:,:)    ! Current position/orientation of grid points.
  real(8),intent(in)   :: NodalMass(:,:,:)  ! Inertia tensor of lumped masses at element nodes.
  real(8),intent(inout):: Mvel     (:,:)    ! Tangent mass matrix.

! Local variables.
  integer :: i                     ! Counters.
  integer :: iNode                 ! Counter on the element nodes.
  real(8) :: CBa(3,3)              ! Coordinate transformation matrix.
  real(8) :: Rot(3,3)              ! Tangential operator.
  real(8) :: Yp(6,6),D(6,6)        ! Ypsilon and D.
  real(8) :: dVdotdvrel(6,6)       ! d(Vdot)/d(vreldot)
  real(8) :: DTYpN(6,6*MaxNodCB3)  ! Strain matrix operator.
  real(8) :: N (6,6*MaxNodCB3)     ! Element shape function matrix

! Loop on the element nodes.
	do iNode=1,NumNodesElem

! Compute element shape function.
    N =0.d0
    do i=1,6
      N(i,i+(iNode-1)*6)= 1.d0
    end do

! Compute the current coordinate transformation matrix and rotational operator.
    CBa= rotvect_psi2mat  (Ri(iNode,4:6))
    Rot= rotvect_psi2rot  (Ri(iNode,4:6))

!Operators for the tangent mass matrix.
    D=0.d0
    D(1:3,1:3)=transpose(CBa)
    D(4:6,4:6)=Unit

    Yp=0.d0
    Yp(1:3,1:3)=Unit
    Yp(4:6,4:6)=Rot

    DTYpN=matmul(matmul(transpose(D),Yp),N)

    dVdotdvrel=0.d0
    dVdotdvrel(1:3,1:3)= CBa
    dVdotdvrel(1:3,4:6)=-matmul(CBa,rot_skew(Ri(iNode,1:3)))
    dVdotdvrel(4:6,4:6)= CBa

! Compute mass tangent stiffness.
    Mvel= Mvel + matmul(transpose(DTYpN),matmul(NodalMass(iNode,:,:),dVdotdvrel))
  end do
  return
 end subroutine cbeam3_rbmvel


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-> Subroutine CBEAM3_CVEL
!
!-> Description:
!
!    Compute tangent damping matrix for motions of reference system.
!
!-> Remarks.-
!
!   1) New values are added to those already in the matrix.
!
!   2) Corrections made to dVgyrdvrel
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine cbeam3_cvel (NumNodesElem,r0,Ri,RiDot,Vrel,ElemMass,Cvel,NumGauss)
  use lib_fem
  use lib_rot
  use lib_rotvect

! I/O Variables.
  integer,intent(in)    :: NumNodesElem     ! Number of nodes in the element.
  real(8),intent(in)    :: r0       (:,:)   ! Initial position/orientation of grid points.
  real(8),intent(in)    :: Ri       (:,:)   ! Current position/orientation of grid points.
  real(8),intent(in)    :: RiDot    (:,:)   ! Current time derivative of position/CRV of points.
  real(8),intent(in)    :: Vrel     (6)     ! Linear/angular velocity of reference frame.
  real(8),intent(in)    :: ElemMass (6,6)   ! Inertia properties in the element.
  real(8),intent(inout) :: Cvel     (:,:)   ! Tangent gyroscopic damping matrix.
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
  real(8) :: Ra(3),RaDot(3)         ! Displacement and its derivative at Gauss point.
  real(8) :: Psi(3),PsiDot(3)       ! Rotation vector and its derivative at Gauss point.
  real(8) :: VB(3),WB(3)            ! Linear and angular velocity in the material frame..
  real(8) :: PB(3),HB(3)            ! Linear and angular momenta in the material frame..

  real(8) :: CBa(3,3)               ! Coordinate transformation matrix.
  real(8) :: Rot(3,3)               ! Tangential operator.
  real(8) :: Yp(6,6),D(6,6)         ! Ypsilon and D.
  real(8) :: ADeltaW   (6,6)        ! A_delta_Omega matrix.
  real(8) :: dVgyrdvrel(6,6)        ! d(Vdot)/d(vrel) at Gauss point.
  real(8) :: dVdvrel   (6,6)        ! dV/d(vrel) at Gauss point.
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
      RaDot (i)=dot_product (ShapeFun(1:NumNodesElem),RiDot(1:NumNodesElem,i))
      Psi   (i)=dot_product (ShapeFun(1:NumNodesElem),Ri   (1:NumNodesElem,i+3))
      PsiDot(i)=dot_product (ShapeFun(1:NumNodesElem),RiDot(1:NumNodesElem,i+3))
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
    VB= matmul(CBa,RaDot+rot_cross(Vrel(4:6),Ra)+Vrel(1:3))
    WB= matmul(Rot,PsiDot) + matmul(CBa,Vrel(4:6))

    PB= matmul(ElemMass(1:3,1:3),VB)+matmul(ElemMass(1:3,4:6),WB)
    HB= matmul(ElemMass(4:6,1:3),VB)+matmul(ElemMass(4:6,4:6),WB)

!Operators in Cbar.
    D=0.0d0
    D(1:3,1:3)=transpose(CBa)
    D(4:6,4:6)=Unit

    Yp=0.0d0
    Yp(1:3,1:3)=Unit
    Yp(4:6,4:6)=Rot

    dVdvrel=0.0d0
    dVdvrel(1:3,1:3)= CBa
    dVdvrel(1:3,4:6)= -matmul(CBa,rot_skew(Ra))
    dVdvrel(4:6,4:6)= CBa

    dVgyrdvrel(1:3,1:3)= -matmul(rot_skew(WB),CBa) + matmul(CBa,rot_skew(Vrel(4:6)))
    dVgyrdvrel(1:3,4:6)= -matmul(CBa,rot_skew(Radot))+matmul(rot_skew(WB),matmul(CBa,rot_skew(Ra))) &
&                        -matmul(CBa,matmul(rot_skew(Vrel(4:6)),rot_skew(Ra)))
    dVgyrdvrel(4:6,1:3)=  0.d0
    dVgyrdvrel(4:6,4:6)= -matmul(rot_skew(WB),CBa) + matmul(CBa,rot_skew(Vrel(4:6)))

    ADeltaW(1:3,1:3)= rot_skew(WB)
    ADeltaW(1:3,4:6)= 0.d0
    ADeltaW(4:6,1:3)= rot_skew(VB)
    ADeltaW(4:6,4:6)= rot_skew(WB)

    ADeltaW=matmul(ADeltaW,ElemMass)

    AdeltaW(1:3,4:6)= AdeltaW(1:3,4:6) - rot_skew(PB)
    AdeltaW(4:6,1:3)= AdeltaW(4:6,1:3) - rot_skew(PB)
    AdeltaW(4:6,4:6)= AdeltaW(4:6,4:6) - rot_skew(HB)

! Compute strain matrix operator and material tangent stiffness.
    Cvel= Cvel + WeightGauss(iGauss) * Jacobian                                 &
&              * matmul(transpose(N),matmul(matmul(transpose(Yp),D),     &
&                                                 (matmul(ElemMass,dVgyrdvrel)  &
&                                                 +matmul(ADeltaW,dVdvrel))))

  end do
  deallocate (CoordGauss,WeightGauss)
  return
 end subroutine cbeam3_cvel

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-> Subroutine CBEAM3_RBCVEL
!
!-> Description:
!
!    Compute lumped masses contributions to the tangent damping matrix due to
!    accelerations of reference system.
!
!-> Remarks.-
!
!   1) New values are added to those already in the matrix.
!
!   2) Corrections made to dVgyrdvrel
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine cbeam3_rbcvel (NumNodesElem,r0,Ri,RiDot,Vrel,NodalMass,Cvel)
  use lib_rot
  use lib_rotvect

! I/O Variables.
  integer,intent(in)   :: NumNodesElem     ! Number of nodes in the element.
  real(8),intent(in)   :: r0       (:,:)   ! Initial position/orientation of grid points.
  real(8),intent(in)   :: Ri       (:,:)   ! Current position/orientation of grid points.
  real(8),intent(in)   :: RiDot    (:,:)   ! Current time derivative of position/CRV of points.
  real(8),intent(in)   :: Vrel     (6)     ! Linear/angular velocity of reference frame.
  real(8),intent(in)   :: NodalMass(:,:,:) ! Inertia tensor of lumped masses at element nodes.
  real(8),intent(inout):: Cvel     (:,:)   ! Tangent gyroscopic damping matrix.

! Local variables
  integer :: i                      ! Counters.
  integer :: iNode                  ! Counter on the nodes in the current element.
  real(8) :: Ra(3),RaDot(3)         ! Displacement and its derivative at Gauss point.
  real(8) :: Psi(3),PsiDot(3)       ! Rotation vector and its derivative at Gauss point.
  real(8) :: VB(3),WB(3)            ! Linear and angular velocity in the material frame..
  real(8) :: PB(3),HB(3)            ! Linear and angular momenta in the material frame..

  real(8) :: CBa(3,3)               ! Coordinate transformation matrix.
  real(8) :: Rot(3,3)               ! Tangential operator.
  real(8) :: Yp(6,6),D(6,6)         ! Ypsilon and D.
  real(8) :: ADeltaW   (6,6)        ! A_delta_Omega matrix.
  real(8) :: dVgyrdvrel(6,6)        ! d(Vdot)/d(vrel) at Gauss point.
  real(8) :: dVdvrel   (6,6)        ! dV/d(vrel) at Gauss point.
  real(8) :: N (6,6*MaxNodCB3)      ! Element shape function matrix.

! Loop in element nodes.
	do iNode=1,NumNodesElem

! Compute the current position vector and rotations, and their derivatives.
		Ra    = Ri   (iNode,1:3)
		RaDot = RiDot(iNode,1:3)
		Psi   = Ri   (iNode,4:6)
		PsiDot= RiDot(iNode,4:6)

! Compute element shape functions.
    N =0.d0
    do i=1,6
      N (i,i+(iNode-1)*6)= 1.d0
    end do

! Compute the current coordinate transformation matrix and rotational operator.
    CBa= rotvect_psi2mat(Psi)
    Rot= rotvect_psi2rot(Psi)

! Compute the current linear and angular velocity.
    VB= matmul(CBa,RaDot + rot_cross(Vrel(4:6),Ra) + Vrel(1:3))
    WB= matmul(Rot,PsiDot) + matmul(CBa,Vrel(4:6))

    PB= matmul(NodalMass(iNode,1:3,1:3),VB) + matmul(NodalMass(iNode,1:3,4:6),WB)
    HB= matmul(NodalMass(iNode,4:6,1:3),VB) + matmul(NodalMass(iNode,4:6,4:6),WB)

!Operators in Cbar.
    D=0.0d0
    D(1:3,1:3)=transpose(CBa)
    D(4:6,4:6)=Unit

    Yp=0.0d0
    Yp(1:3,1:3)=Unit
    Yp(4:6,4:6)=Rot

    dVdvrel=0.0d0
    dVdvrel(1:3,1:3)= CBa
    dVdvrel(1:3,4:6)=-matmul(CBa,rot_skew(Ra))
    dVdvrel(4:6,4:6)= CBa

    dVgyrdvrel(1:3,1:3)= -matmul(rot_skew(WB),CBa) + matmul(CBa,rot_skew(Vrel(4:6)))
    dVgyrdvrel(1:3,4:6)= -matmul(CBa,rot_skew(Radot))+matmul(rot_skew(WB),matmul(CBa,rot_skew(Ra))) &
&                        -matmul(CBa,matmul(rot_skew(Vrel(4:6)),rot_skew(Ra)))
    dVgyrdvrel(4:6,1:3)=  0.d0
    dVgyrdvrel(4:6,4:6)= -matmul(rot_skew(WB),CBa) + matmul(CBa,rot_skew(Vrel(4:6)))

    ADeltaW(1:3,1:3)= rot_skew(WB)
    ADeltaW(1:3,4:6)= 0.d0
    ADeltaW(4:6,1:3)= rot_skew(VB)
    ADeltaW(4:6,4:6)= rot_skew(WB)

    ADeltaW=matmul(ADeltaW,NodalMass(iNode,:,:))

    AdeltaW(1:3,4:6)= AdeltaW(1:3,4:6) - rot_skew(PB)
    AdeltaW(4:6,1:3)= AdeltaW(4:6,1:3) - rot_skew(PB)
    AdeltaW(4:6,4:6)= AdeltaW(4:6,4:6) - rot_skew(HB)

! Compute strain matrix operator and material tangent stiffness.
    Cvel= Cvel + matmul(transpose(N),matmul( matmul(transpose(Yp),D),                &
&		                                        (matmul(NodalMass(iNode,:,:),dVgyrdvrel) &
&                                           +matmul(ADeltaW,dVdvrel))))
  end do
	
  return
 end subroutine cbeam3_rbcvel



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-> Subroutine CBEAM3_FEXT
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
 subroutine cbeam3_fext (NumNodesElem,Ri,Flags,Fmat,FollowerForce,FollowerForceRig,Cao)
  use lib_fem
  use lib_rot
  use lib_rotvect

! I/O Variables.
  integer,intent(in)   :: NumNodesElem     ! Number of nodes in the element.
  real(8),intent(in)   :: Ri    (:,:)      ! Current position/orientation of grid points.
  logical,intent(in)   :: Flags (:)        ! Identify master nodes.
  real(8),intent(inout):: Fmat  (:,:)      ! Forces/moments on the element nodes.
  logical,intent(in)   :: FollowerForce    ! =T if follower forces.
  logical,intent(in)   :: FollowerForceRig ! =T if follower forces.
  real(8),intent(in)   :: Cao    (:,:)     ! Rotation operator

! Local variables.
  integer :: i                     ! Counters.
  integer :: iNode                 ! Counter on the element nodes.
  real(8) :: D(6,6)                ! D matrix.
  real(8) :: Yp(6,6)               ! Ypsilon.
  real(8) :: N (6,6*MaxNodCB3)     ! Element shape function matrix


!--> Case of follower forces.
  if (FollowerForceRig) then
    if (FollowerForce) then
     ! Loop in the element nodes.
     do iNode=1,NumNodesElem
      if (Flags(iNode)) then
        ! Compute shape function at the current node.
        N =0.d0
        do i=1,6
          N (i,i+(iNode-1)*6)= 1.d0
        end do

        ! Operators to project the forces on the dofs.
        D=0.d0
        D(1:3,1:3)=transpose(rotvect_psi2mat(Ri(iNode,4:6)))
        D(4:6,4:6)=Unit

        Yp=0.d0
        Yp(1:3,1:3)=Unit
        Yp(4:6,4:6)=rotvect_psi2rot(Ri(iNode,4:6))

        ! Compute influence coefficients for follower forces.
        Fmat(:,6*(iNode-1)+1:6*iNode)= Fmat(:,6*(iNode-1)+1:6*iNode) &
 &                                   + matmul(transpose(matmul(Yp,N)),D)
      end if
     end do

!--> Case of dead forces but following the a frame.
    else
     do iNode=1,NumNodesElem
      if (Flags(iNode)) then
        ! Compute shape function at the current node.
        N =0.d0
        do i=1,6
          N (i,i+(iNode-1)*6)= 1.d0
        end do

        Yp=0.d0
        Yp(1:3,1:3)=Unit
        Yp(4:6,4:6)=rotvect_psi2rot(Ri(iNode,4:6))

        ! Compute influence coefficients for follower forces.
        Fmat(:,6*(iNode-1)+1:6*iNode)= Fmat(:,6*(iNode-1)+1:6*iNode) + transpose(matmul(transpose(Yp),N))
      end if
    end do
  end if

!--> Case of dead forces. Not possible to have structural follower forces and rigid-body dead loads.
  else
! Loop in the element nodes.
     do iNode=1,NumNodesElem

      if (Flags(iNode)) then
! Compute shape function at the current node.
        N =0.d0
        do i=1,6
          N (i,i+(iNode-1)*6)= 1.d0
        end do

! Operators to project the forces on the dofs.
        D=0.d0
        D(1:3,1:3)=Cao
        D(4:6,4:6)=Cao

        Yp=0.d0
        Yp(1:3,1:3)=Unit
        Yp(4:6,4:6)=rotvect_psi2rot(Ri(iNode,4:6))

        ! Compute influence coefficients for follower forces.
        Fmat(:,6*(iNode-1)+1:6*iNode)= Fmat(:,6*(iNode-1)+1:6*iNode) + matmul(transpose(matmul(transpose(Yp),N)),D)

      end if
     end do
  end if

  return
 end subroutine cbeam3_fext


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
 subroutine cbeam3_dqext (NumNodesElem,Ri,Fi,Flags,Kmat,FollowerForce)
  use lib_fem
  use lib_rot
  use lib_rotvect

! I/O Variables.
  integer,intent(in)   :: NumNodesElem    ! Number of nodes in the element.
  real(8),intent(in)   :: Ri    (:,:)     ! Current position/orientation of grid points.
  real(8),intent(in)   :: Fi    (:,:)     ! Current forces/moments at grid points.
  logical,intent(in)   :: Flags (:)       ! Identify master nodes.
  real(8),intent(inout):: Kmat  (:,:)     ! Derivative of the external forces.
  logical,intent(in)   :: FollowerForce   ! =T if follower forces.

! Local variables.
  integer :: iNode                 ! Counter on the element nodes.
  real(8) :: CBa (3,3)             ! Coordinate transformation matrix.
  real(8) :: Rot (3,3)             ! Tangential operator.
  real(8) :: Kbar(6,6)             ! Matrix at current node.

! This routine gives an output only for follower forces.
  if (FollowerForce) then

! Loop in the element nodes.
  do iNode=1,NumNodesElem
  if (Flags(iNode)) then

! Compute the current coordinate transformation matrix and rotational operator.
    CBa= rotvect_psi2mat(Ri(iNode,4:6))
    Rot= rotvect_psi2rot(Ri(iNode,4:6))

! Compute Kbar for the tangent matrix at the current Gauss point.
    Kbar= 0.d0
    Kbar(1:3,4:6)= matmul(transpose(CBa),matmul(rot_skew(Fi(iNode,1:3)),Rot))
    Kbar(4:6,4:6)=-rotvect_a2(Ri(iNode,4:6),Fi(iNode,4:6))

! Compute influence coefficients for follower forces.
    Kmat(6*(iNode-1)+1:6*iNode,6*(iNode-1)+1:6*iNode)= Kmat(6*(iNode-1)+1:6*iNode,6*(iNode-1)+1:6*iNode) + Kbar

  end if
  end do
  end if

  return
 end subroutine cbeam3_dqext


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-> Subroutine CBEAM3_SLAVE2MASTER
!
!-> Description:
!
!     Compute transformation from slave to master nodal equations. This is
!     done to include discontinuities between adjacent elements
!
!-> Remarks.-
!
!   1) For a slave frame, 1 refers to the node on the master element and 2
!      is the same node on the slave element.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine cbeam3_slave2master (NNE,TreeConn,Psi02_0,AllPsi_0,Psi02_t, &
&                                AllPsi_t,S21)
  use lib_rotvect
  use lib_lu

! I/O Variables.
  integer,intent(in) :: NNE                ! Number of nodes in the element.
  integer,intent(in) :: TreeConn  (:,:)    ! Connectivity tree.
  real(8),intent(in) :: Psi02_0   (:,:)    ! Initial CRV from global to nodes in element.
  real(8),intent(in) :: AllPsi_0  (:,:,:)  ! Initial CRV from global to nodes in all elements.
  real(8),intent(in) :: Psi02_t   (:,:)    ! Current CRV from global to nodes in element.
  real(8),intent(in) :: AllPsi_t  (:,:,:)  ! Current CRV from global to nodes in all elements.
  real(8),intent(out):: S21       (:,:)    ! Element transformation matrix master2slave.

! Local variables.
  integer :: i1         ! Counter.
  integer :: Info       ! Error code
  integer :: iNode      ! Counter in the nodes in the element.
  real(8) :: C_12(3,3)  ! Coordinate transformation matrix from slave (2) to master (1) node.
  real(8) :: Psi21(3)   ! CRV from slave to master.
  real(8) :: Psi01(3)   ! CRV from global to slave.
  real(8) :: T01(3,3)   ! Tangent matrix at the master node.
  real(8) :: T02(3,3)   ! Tangent matrix at the slave node.
  real(8) :: invT02(3,3)! Inverse of T02.
  real(8) :: S21n(3,3)  ! Transformation matrix master-to-slave for the current node.

! Initialise
  S21(:,:)=0.d0

! Loop in the elements nodes.
  i1=0
  do iNode=1,NNE
    S21n= Unit

! For each slave node, determine the constant rotation with respect to its master.
    if (TreeConn(iNode,1).ne.0) then
      Psi01= AllPsi_0(TreeConn(iNode,1),TreeConn(iNode,2),:)
      Psi21= rotvect_addpsi(-Psi02_0(iNode,:),Psi01)

! If the rotation is non-zero, slave-node equations need to be transformed.
      if (dot_product(Psi21,Psi21).ne.0.d0) then
        C_12=  rotvect_psi2mat(Psi21)

! Transformation operator from local to global rotations.
        Psi01= AllPsi_t(TreeConn(iNode,1),TreeConn(iNode,2),:)
        T02= rotvect_psi2rot(Psi02_t(iNode,:))
        T01= rotvect_psi2rot(Psi01)

        call lu_invers(T02,invT02); Info=0        
        if (Info.ne.0) STOP 'Error: Could not compute the inverse of the tangent operator (77921)'

        S21n= matmul(invT02,matmul(transpose(C_12),T01))
      end if
    end if

! Only rotation equations on slave nodes are transformed.
    S21   (i1+1:i1+3,i1+1:i1+3)= Unit
    S21   (i1+4:i1+6,i1+4:i1+6)= S21n
    i1= i1+6
  end do

  return
 end subroutine cbeam3_slave2master



 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-> Subroutine CBEAM3_PROJS2M
!
!-> Description:
!
!     Compute transformation from slave to master nodal equations. This is done to
!     include discontinuities between adjacent elements
!
!-> Remarks.-
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine cbeam3_projs2m (NumNodesElem,TreeConn,Psi02,AllPsi,Telem)
  use lib_rotvect
  use lib_lu

! I/O Variables.
  integer,intent(in) :: NumNodesElem       ! Number of nodes in the element.
  integer,intent(in) :: TreeConn  (:,:)    ! Connectivity tree.
  real(8),intent(in) :: Psi02     (:,:)    ! Initial CRV from global to nodes in element.
  real(8),intent(in) :: AllPsi    (:,:,:)  ! Initial CRV from global to nodes in all elements.
  real(8),intent(out):: Telem     (:,:)    ! Element transformation matrix.

! Local variables.
  integer :: i1         ! Counter.
  integer :: Info       ! Error code
  integer :: iNode      ! Counter in the nodes in the element.
  real(8) :: C_12(3,3)  ! Coordinate transformation matrix from slave (2) to master (1) node.
  real(8) :: Psi21(3)   ! CRV from slave to master.
  real(8) :: Psi01(3)   ! CRV from global to slave.
  real(8) :: T01(3,3)   ! Tangent matrix at the master node.
  real(8) :: T02(3,3)   ! Tangent matrix at the slave node.
  real(8) :: invT02(3,3)! Inverse of T02.
  real(8) :: Tnode(3,3) ! Transformation matrix master-to-slave for the current node.

! Loop in the elements nodes.
  i1=0
  do iNode=1,NumNodesElem

! Slave-node equations need to be transformed.
    if (TreeConn(iNode,1).ne.0) then
      Psi01= AllPsi(TreeConn(iNode,1),TreeConn(iNode,2),:)
      Psi21= rotvect_addpsi(-Psi02(iNode,:),Psi01)
      C_12=  rotvect_psi2mat(Psi21)

! Transformation operator from local to global rotations.
      T02= rotvect_psi2rot(Psi02(iNode,:))
      T01= rotvect_psi2rot(Psi01)

      call lu_invers(T02,invT02); Info=0      
      if (Info.ne.0) STOP 'Error: Could not compute the inverse of the tangent operator (7921)'

      Tnode=matmul(invT02,matmul(transpose(C_12),T01))

! Master node equations are not transformed.
    else
      Tnode= Unit
    end if

! Only rotation equations are transformed.
    Telem(i1+1:i1+3,i1+1:i1+3)= Unit
    Telem(i1+4:i1+6,i1+4:i1+6)= Tnode
    i1= i1+6
  end do

  return
 end subroutine cbeam3_projs2m



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-> Subroutine CBEAM3_PROJM2S
!
!-> Description:
!
!     Compute transformation for common element equations. This is done to
!     include discontinuities between adjacent elements
!
!     Project rotation from master to slave grid node.
!
!-> Remarks.-
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine cbeam3_projm2s (NumNodesElem,TreeConn,Psi0,AllPsi,Psi)
  use lib_rotvect

! I/O Variables.
  integer,intent(in) :: NumNodesElem       ! Number of nodes in the element.
  integer,intent(in) :: TreeConn  (:,:)    ! Connectivity tree.
  real(8),intent(in) :: Psi0      (:,:)    ! Initial CRV from global to nodes in element.
  real(8),intent(in) :: AllPsi    (:,:,:)  ! Initial CRV from global to nodes in all elements.
  real(8),intent(out):: Psi       (:,:)    ! Current CRV in element.

! Local variables.
  integer :: iNode      ! Counter in the nodes in the element.
  real(8) :: Psi21(3)   ! CRV from slave to master.
  real(8) :: Psi01(3)   ! CRV of the master node.

! Loop in the elements nodes and identify slave ones.
  do iNode=1,NumNodesElem
    if (TreeConn(iNode,1).ne.0) then

      ! Find slave-to-master transformation in the initial configuration.
      Psi21= rotvect_addpsi(-Psi0(iNode,:),AllPsi(TreeConn(iNode,1),TreeConn(iNode,2),:))
      Psi01= Psi(iNode,:)

      ! Transform current rotation in the element.
      if (maxval(abs(Psi21)).gt.0.d0) Psi(iNode,:) = rotvect_addpsi(Psi01,-Psi21)

    end if
  end do

  return
 end subroutine cbeam3_projm2s


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-> Subroutine CBEAM3_PROJ2A
!
!-> Description:
!
!    Project position vector and rotations at element nodes into the
!    element frame.
!
!-> Remarks.-
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine cbeam3_proj2a (NumNodesElem,Psi,R)
  use lib_rot
  use lib_rotvect

! I/O Variables.
  integer,intent(in) :: NumNodesElem      ! Number of nodes in the element.
  real(8),intent(in) :: Psi    (:)        ! Element orientation.
  real(8),intent(out):: R      (:,:)      ! Current position vector/rotation at grid points.

! Local variables
  real(8):: Aux(3)   ! Auxiliary vector (against nasty compilers).
  integer:: i        ! Counter.

! Loop in the nodes of the element.
  do i=1,NumNodesElem

! Rotate the position vector of each grid.
    Aux=matmul(rotvect_psi2mat(Psi),R(i,1:3))
    R(i,1:3)=Aux

! Compute the rotation vector from the displaced element frame to the local deformed frame.
    Aux=rotvect_addpsi(-Psi,R(i,4:6))
    R(i,4:6)=Aux
  end do

  return
 end subroutine cbeam3_proj2a


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-> Subroutine CBEAM3_GLOB2LOC
!
!-> Description:
!
!     Compute transformation from local to global variations.
!
!-> Remarks.-
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine cbeam3_glob2loc (NumNodesElem,Psi_ga,Psi_aB,Sga)
  use lib_lu
  use lib_rotvect

! I/O Variables.
  integer,intent(in) :: NumNodesElem    ! Number of nodes in the element.
  real(8),intent(in) :: Psi_ga (:)      ! Rotation vector from global to element.
  real(8),intent(in) :: Psi_aB (:,:)    ! Rotation vector from element to local.
  real(8),intent(out):: Sga    (:,:)    ! Current position vector/rotation at grid points.

! Local variables.
  integer :: i, i1      ! Counters.
  integer :: Info       ! Error code.
  real(8) :: InvTaB(3,3)! Inverse of TaB.
  real(8) :: TaB(3,3)   ! Tangent operator in the rotation from element to local frame.
  real(8) :: TgB(3,3)   ! Tangent operator in the rotation from global to the local frame.

! Loop in the elements nodes.
  Sga= 0.d0
  i1 = 0
  do i=1,NumNodesElem

! Coordinate transformation matrix for the position vector at current node.
    Sga(i1+1:i1+3,i1+1:i1+3)= transpose(rotvect_psi2mat(Psi_ga))
    i1= i1+3

! Transformation operator from local to global rotations.
    TgB= rotvect_psi2rot(rotvect_addpsi(Psi_ga,Psi_aB(i,:)))
    TaB= rotvect_psi2rot(Psi_aB(i,:))

    call lu_invers(TaB,invTaB); Info=0
    if (Info.ne.0) STOP 'Error: Could not compute the inverse of the tangent operator (7921)'

    Sga(i1+1:i1+3,i1+1:i1+3)= transpose(matmul(invTaB,TgB))
    i1= i1+3

  end do

  return
 end subroutine cbeam3_glob2loc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end module lib_cbeam3
