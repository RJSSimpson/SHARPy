!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-> Module.- LIB_XBEAM Henrik Hesse. 07/01/2011 - Last Update 07/01/2011
!
!-> Language: FORTRAN90, Free format.
!
!-> Description.-
!
!   Compute element matrices for rigid-body components of the 3-noded geometrically-exact beam element.
!
!-> Subroutines.-
!
!    -xbeam_mrs:        Compute tangent mass matrix.
!    -xbeam_kmass:      Compute tangent mass stiffness matrix
!    -xbeam_cgyr:       Compute tangent gyroscopic damping matrix.
!    -xbeam_kgyr:       Compute tangent gyroscopic stiffness matrix.
!    -xbeam_fgyr:       Compute discrete gyroscopic forces.
!
!    -xbeam_mrr:        Compute reference system mass matrix.
!    -xbeam_crr:        Compute reference system damping matrix.
!
!    -xbeam_fext:       Compute influence coefficients for force vector.
!
!    -xbeam_rot:        Definition of rotation operator in Euler quaternions
!    -xbeam_quadskew:   See Shearer and Cesnik (2007)
!
!    -xbeam_lapack_inv: Computes the inverse of a matrix
!    -xbeam_2norm:      Computes the 2norm of a vector
!
!-> Remarks.-
!
!  1) The degrees of freedom in the element are the displacements and the
!     Cartesian Rotation Vector (CRV) at each of the nodes)
!
!  2) The element can have 2 or 3 nodes.
!
!-> Modifications.-
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module lib_xbeam
 implicit none

! Shared variables within the module.
 integer,private,parameter:: MaxNodCB3=3               ! Max number of nodes per element is 3.
 real(8),private,parameter,dimension(3,3):: Unit= &    ! Unit matrix.
&         reshape((/1.d0,0.d0,0.d0,0.d0,1.d0,0.d0,0.d0,0.d0,1.d0/),(/3,3/))

 contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-> Subroutine XBEAM_MRS
!
!-> Description:
!
!    Compute the tangent mass matrix for rigid equation.
!
!-> Remarks.-
!
!   1) New values are added to those already in the matrix.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine xbeam_mrs (NumNodesElem,r0,Ri,ElemMass,Mass,NumGauss)
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
  real(8) :: CaB(3,3)              ! Coordinate transformation matrix.
  real(8) :: Rot(3,3)              ! Tangential operator.
  real(8) :: Yp(6,6),D(6,6)        ! Ypsilon and D.
  real(8) :: ARC(6,6)              ! A_{RC}
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
    CaB= transpose(rotvect_psi2mat(Ri_g(4:6)))
    Rot= rotvect_psi2rot  (Ri_g(4:6))

!Operators for the tangent mass matrix.
    ARC=0.d0
    ARC(1:3,1:3)=CaB
    ARC(4:6,1:3)=matmul(rot_skew(Ri_g(1:3)),CaB)
    ARC(4:6,4:6)=CaB

    D=0.d0
    D(1:3,1:3)=CaB
    D(4:6,4:6)=Unit

    Yp=0.d0
    Yp(1:3,1:3)=Unit
    Yp(4:6,4:6)=Rot

    DTYpN=matmul(matmul(transpose(D),Yp),N)

! Compute mass tangent stiffness.
    Mass= Mass + WeightGauss(iGauss)*Jacobian*matmul(ARC,matmul(ElemMass,DTYpN))
  end do
  deallocate (CoordGauss,WeightGauss)

  return
 end subroutine xbeam_mrs

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-> Subroutine XBEAM_RBMRS
!
!-> Description:
!
!    Compute the lumped mass contributions to the tangent mass matrix for
!    rigid equation.
!
!-> Remarks.-
!
!   1) New values are added to those already in the matrix.
!   2) The lumped masses are defined by the inertia tensor of the nonstructural masses
!      in the section of the node, and expressed with respect to the local B frame.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine xbeam_rbmrs (NumNodesElem,r0,Ri,NodalMass,Mass)
  use lib_fem
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
  real(8) :: Ra(3)                 ! Current position vector
  real(8) :: CaB(3,3)              ! Coordinate transformation matrix.
  real(8) :: Rot(3,3)              ! Tangential operator.
  real(8) :: ARC(6,6)              ! A_{RC}
  real(8) :: Yp(6,6),D(6,6)        ! Ypsilon and D.
  real(8) :: DTYpN(6,6*MaxNodCB3)  ! Strain matrix operator.
  real(8) :: N (6,6*MaxNodCB3)     ! Element shape function matrix

! Loop in the nodes.
  do iNode=1,NumNodesElem

! Compute the nodal coordinate transformation matrix and rotational operator.
    CaB= transpose(rotvect_psi2mat(Ri(iNode,4:6)))
    Rot= rotvect_psi2rot(Ri(iNode,4:6))

! Compute element shape function at the current node.
    N =0.d0
    do i=1,6
        N (i,i+(iNode-1)*6)= 1.d0
    end do

!Operators for the tangent mass matrix.
    ARC=0.d0
    ARC(1:3,1:3)=CaB
    ARC(4:6,1:3)=matmul(rot_skew(Ra),CaB)
    ARC(4:6,4:6)=CaB

    D=0.d0
    D(1:3,1:3)=CaB
    D(4:6,4:6)=Unit

    Yp=0.d0
    Yp(1:3,1:3)=Unit
    Yp(4:6,4:6)=Rot

    DTYpN=matmul(matmul(transpose(D),Yp),N)

! Compute mass tangent stiffness.
    Mass= Mass + matmul(ARC,matmul(NodalMass(iNode,:,:),DTYpN))
  end do
  return
 end subroutine xbeam_rbmrs

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-> Subroutine XBEAM_KMASS
!
!-> Description:
!
!    Compute tangent damping matrix for motions of reference system.
!
!-> Remarks.-
!
!   1) New values are added to those already in the matrix.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine xbeam_kmass (NumNodesElem,r0,Ri,RiDDot,VrelDot,ElemMass,Kmass,NumGauss)
  use lib_fem
  use lib_rot
  use lib_rotvect

! I/O Variables.
  integer,intent(in)    :: NumNodesElem     ! Number of nodes in the element.
  real(8),intent(in)    :: r0       (:,:)   ! Initial position/orientation of grid points.
  real(8),intent(in)    :: Ri       (:,:)   ! Current position/orientation of grid points.
  real(8),intent(in)    :: RiDDot   (:,:)   ! Current time derivative of position/CRV of points.
  real(8),intent(in)    :: VrelDot  (6)     ! Linear/angular acceleration of reference frame.
  real(8),intent(in)    :: ElemMass (6,6)   ! Inertia properties in the element.
  real(8),intent(inout) :: Kmass     (:,:)   ! Tangent gyroscopic damping matrix.
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
  real(8) :: Ra(3),RaDDot(3)        ! Displacement and its derivative at Gauss point.
  real(8) :: Psi(3),PsiDDot(3)      ! Rotation vector and its derivative at Gauss point.
  real(8) :: ddq(6)                 ! state of system

  real(8) :: CBa(3,3),CaB(3,3)      ! Coordinate transformation matrix.
  real(8) :: Rot(3,3)               ! Tangential operator.
  real(8) :: Yp(6,6),D(6,6)         ! Ypsilon and D.
  real(8) :: ARC(6,6)               ! A_{RC} matrix
  real(8) :: N (6,6*MaxNodCB3)      ! Element shape function matrix.

  real(8) :: MDtYddq(6),dDtYddq(6,6),dARCdqMDtYddq(6,6)                  ! Derivatives d^2(.)/d(q)^2 computed at Gauss points.
  real(8) :: MARCtVrelDot(6),dARCdqMARCtVrelDot(6,6),dARCtVrelDot(6,6)   ! Derivatives d(.)/d(Vrel) computed at Gauss points.

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
      Ra     (i)=dot_product (ShapeFun(1:NumNodesElem),Ri   (1:NumNodesElem,i))
      RaDDot (i)=dot_product (ShapeFun(1:NumNodesElem),RiDDot(1:NumNodesElem,i))
      Psi    (i)=dot_product (ShapeFun(1:NumNodesElem),Ri   (1:NumNodesElem,i+3))
      PsiDDot(i)=dot_product (ShapeFun(1:NumNodesElem),RiDDot(1:NumNodesElem,i+3))

    end do

    ddq(1:3)= RaDDot
    ddq(4:6)= PsiDDot

! Compute element shape functions.
    N =0.d0
    do i=1,6
      do j=1,NumNodesElem
        N (i,i+(j-1)*6)= ShapeFun(j)
      end do
    end do

! Compute the current coordinate transformation matrix and rotational operator.
    CBa= rotvect_psi2mat  (Psi)
    CaB= transpose(CBa)
    Rot= rotvect_psi2rot  (Psi)

!Operators in d(Mtddq)/dq.
    D=0.0d0
    D(1:3,1:3)=CaB
    D(4:6,4:6)=Unit

    Yp=0.0d0
    Yp(1:3,1:3)=Unit
    Yp(4:6,4:6)=Rot

    ARC=0.d0
    ARC(1:3,1:3)=CaB
    ARC(4:6,1:3)=matmul(rot_skew(Ra),CaB)
    ARC(4:6,4:6)=CaB

    MDtYddq= matmul(ElemMass,matmul(transpose(D),matmul(Yp,ddq)))

    dARCdqMDtYddq(1:3,1:3)= 0.d0
    dARCdqMDtYddq(1:3,4:6)= -matmul(CaB,matmul(rot_skew(MDtYddq(1:3)),Rot))
    dARCdqMDtYddq(4:6,1:3)= -rot_skew(matmul(CaB,MDtYddq(1:3)))
    dARCdqMDtYddq(4:6,4:6)= -matmul(CaB,matmul(rot_skew(MDtYddq(4:6)),Rot))    &
&                           -matmul(rot_skew(Ra),matmul(CaB,matmul(rot_skew(MDtYddq(1:3)),Rot)))

    dDtYddq= 0.d0
    dDtYddq(1:3,4:6)= matmul(CBa,matmul(rot_skew(ddq(1:3)),transpose(Rot)))
    dDtYddq(4:6,4:6)= rotvect_a1 (Psi,ddq(4:6))

! Compute d(Mtddq)/dq part in Kmass
    Kmass= Kmass + WeightGauss(iGauss) * Jacobian                                 &
&                  * matmul((dARCdqMDtYddq + matmul(ARC,matmul(ElemMass,dDtYddq))),N)

! Additional operators in d(madVrel)/dq
    MARCtVrelDot=matmul(ElemMass,matmul(transpose(ARC),VrelDot))

    dARCdqMARCtVrelDot(1:3,1:3)= 0.d0
    dARCdqMARCtVrelDot(1:3,4:6)= -matmul(CaB,matmul(rot_skew(MARCtVrelDot(1:3)),Rot))
    dARCdqMARCtVrelDot(4:6,1:3)= -rot_skew(matmul(CaB,MARCtVrelDot(1:3)))
    dARCdqMARCtVrelDot(4:6,4:6)= -matmul(CaB,matmul(rot_skew(MARCtVrelDot(4:6)),Rot))    &
&                                -matmul(rot_skew(Ra),matmul(CaB,matmul(rot_skew(MARCtVrelDot(1:3)),Rot)))

    dARCtVrelDot(1:3,1:3)= matmul(CBa,rot_skew(VrelDot(4:6)))
    dARCtVrelDot(1:3,4:6)= matmul(CBa,matmul((rot_skew(VrelDot(1:3))+rot_skew(matmul(rot_skew(VrelDot(4:6)),Ra))),transpose(Rot)))
    dARCtVrelDot(4:6,1:3)= 0.d0
    dARCtVrelDot(4:6,4:6)= matmul(CBa,matmul(rot_skew(ddq(4:6)),transpose(Rot)))

! Compute d(madVrel)/dq and add to Kmass
    Kmass= Kmass + WeightGauss(iGauss) * Jacobian                                 &
&                  * matmul((dARCdqMARCtVrelDot + matmul(ARC,matmul(ElemMass,dARCtVrelDot))),N)

  end do
  deallocate (CoordGauss,WeightGauss)
  return
 end subroutine xbeam_kmass


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-> Subroutine XBEAM_RBKMASS
!
!-> Description:
!
!    Compute the lumped mass contributions to the tangent damping matrix
!    for motions of reference system.
!
!-> Remarks.-
!
!   1) New values are added to those already in the matrix.
!   2) The lumped masses are defined by the inertia tensor of the nonstructural masses
!      in the section of the node, and expressed with respect to the local B frame.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine xbeam_rbkmass (NumNodesElem,r0,Ri,RiDDot,VrelDot,NodalMass,Kmass)
  use lib_fem
  use lib_rot
  use lib_rotvect

! I/O Variables.
  integer,intent(in)   :: NumNodesElem     ! Number of nodes in the element.
  real(8),intent(in)   :: r0       (:,:)   ! Initial position/orientation of grid points.
  real(8),intent(in)   :: Ri       (:,:)   ! Current position/orientation of grid points.
  real(8),intent(in)   :: RiDDot   (:,:)   ! Current time derivative of position/CRV of points.
  real(8),intent(in)   :: VrelDot  (6)     ! Linear/angular acceleration of reference frame.
  real(8),intent(in)   :: NodalMass(:,:,:) ! Inertia tensor of lumped masses at element nodes.
  real(8),intent(inout):: Kmass     (:,:)  ! Tangent gyroscopic damping matrix.

! Local variables
  integer :: i,iNode               ! Counters.
  real(8) :: Ra(3),RaDDot(3)        ! Displacement and its derivative at Gauss point.
  real(8) :: Psi(3),PsiDDot(3)      ! Rotation vector and its derivative at Gauss point.
  real(8) :: ddq(6)                 ! state of system

  real(8) :: CBa(3,3),CaB(3,3)      ! Coordinate transformation matrix.
  real(8) :: Rot(3,3)               ! Tangential operator.
  real(8) :: Yp(6,6),D(6,6)         ! Ypsilon and D.
  real(8) :: ARC(6,6)               ! A_{RC} matrix
  real(8) :: N (6,6*MaxNodCB3)      ! Element shape function matrix.

  real(8) :: MDtYddq(6),dDtYddq(6,6),dARCdqMDtYddq(6,6)                 ! Derivatives d^2(.)/d(q)^2 computed at Gauss points.
  real(8) :: MARCtVrelDot(6),dARCdqMARCtVrelDot(6,6),dARCtVrelDot(6,6)  ! Derivatives d(.)/d(Vrel) computed at Gauss points.

! Loop in the nodes.
  do iNode=1,NumNodesElem

 ! Compute the current position vector and rotations, and their derivatives.
		Ra     = Ri    (iNode,1:3)
		RaDDot = RiDDot(iNode,1:3)
		Psi    = Ri    (iNode,4:6)
		PsiDDot= RiDDot(iNode,4:6)

    ddq(1:3)= RaDDot
    ddq(4:6)= PsiDDot

! Compute element shape functions at the corrent node.
    N =0.d0
    do i=1,6
      N (i,i+(iNode-1)*6)= 1.d0
    end do

! Compute the current coordinate transformation matrix and rotational operator.
    CBa= rotvect_psi2mat  (Psi)
    CaB= transpose(CBa)
    Rot= rotvect_psi2rot  (Psi)

!Operators in d(Mtddq)/dq.
    D=0.0d0
    D(1:3,1:3)=CaB
    D(4:6,4:6)=Unit

    Yp=0.0d0
    Yp(1:3,1:3)=Unit
    Yp(4:6,4:6)=Rot

    ARC=0.d0
    ARC(1:3,1:3)=CaB
    ARC(4:6,1:3)=matmul(rot_skew(Ra),CaB)
    ARC(4:6,4:6)=CaB

    MDtYddq= matmul(NodalMass(iNode,:,:),matmul(transpose(D),matmul(Yp,ddq)))

    dARCdqMDtYddq(1:3,1:3)= 0.d0
    dARCdqMDtYddq(1:3,4:6)= -matmul(CaB,matmul(rot_skew(MDtYddq(1:3)),Rot))
    dARCdqMDtYddq(4:6,1:3)= -rot_skew(matmul(CaB,MDtYddq(1:3)))
    dARCdqMDtYddq(4:6,4:6)= -matmul(CaB,matmul(rot_skew(MDtYddq(4:6)),Rot))    &
&                           -matmul(rot_skew(Ra),matmul(CaB,matmul(rot_skew(MDtYddq(1:3)),Rot)))

    dDtYddq= 0.d0
    dDtYddq(1:3,4:6)= matmul(CBa,matmul(rot_skew(ddq(1:3)),transpose(Rot)))
    dDtYddq(4:6,4:6)= rotvect_a1 (Psi,ddq(4:6))

! Compute d(Mtddq)/dq part in Kmass
    Kmass= Kmass + matmul((dARCdqMDtYddq + matmul(ARC,matmul(NodalMass(iNode,:,:),dDtYddq))),N)

! Additional operators in d(madVrel)/dq
    MARCtVrelDot=matmul(NodalMass(iNode,:,:),matmul(transpose(ARC),VrelDot))

    dARCdqMARCtVrelDot(1:3,1:3)= 0.d0
    dARCdqMARCtVrelDot(1:3,4:6)= -matmul(CaB,matmul(rot_skew(MARCtVrelDot(1:3)),Rot))
    dARCdqMARCtVrelDot(4:6,1:3)= -rot_skew(matmul(CaB,MARCtVrelDot(1:3)))
    dARCdqMARCtVrelDot(4:6,4:6)= -matmul(CaB,matmul(rot_skew(MARCtVrelDot(4:6)),Rot))    &
&                                -matmul(rot_skew(Ra),matmul(CaB,matmul(rot_skew(MARCtVrelDot(1:3)),Rot)))

    dARCtVrelDot(1:3,1:3)= matmul(CBa,rot_skew(VrelDot(4:6)))
    dARCtVrelDot(1:3,4:6)= matmul(CBa,matmul((rot_skew(VrelDot(1:3))+rot_skew(matmul(rot_skew(VrelDot(4:6)),Ra))),transpose(Rot)))
    dARCtVrelDot(4:6,1:3)= 0.d0
    dARCtVrelDot(4:6,4:6)= matmul(CBa,matmul(rot_skew(ddq(4:6)),transpose(Rot)))

! Compute d(madVrel)/dq and add to Kmass
    Kmass= Kmass + matmul((dARCdqMARCtVrelDot + matmul(ARC,matmul(NodalMass(iNode,:,:),dARCtVrelDot))),N)

  end do
  return
 end subroutine xbeam_rbkmass


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-> Subroutine XBEAM_KGYR
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
 subroutine xbeam_kgyr (NumNodesElem,r0,Ri,RiDot,RiDDot,Vrel,VrelDot, &
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
  real(8) :: RaDDot(3)              ! Time derivatives of RaDot.
  real(8) :: RotPsiDot(3)           ! Rot*PsiDot.
  real(8) :: VB(3),WB(3)            ! Linear and angular velocity in the material frame.
  real(8) :: VBdot(3)               ! Time derivatives of VB.
  real(8) :: Vgyr(6)                ! Vectors with angular velocities.
  real(8) :: Phat(6),Pgyr(6)        ! Vectors with linear/angular momenta

  real(8) :: CBa(3,3),CaB(3,3)      ! Coordinate transformation matrix.
  real(8) :: Rot(3,3),RotDot(3,3)   ! Tangential operator and time derivative.
  real(8) :: ARC(6,6),AWRC(6,6)     ! A_{RC} and A_{wRC}.
  real(8) :: dARCdqPgyr(6,6),dAWRCdqPhat(6,6)
  real(8) :: dVgyrdq  (6,6)         ! d(Vdot)/dq at Gauss point.
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
    CaB= transpose(CBa)
    Rot= rotvect_psi2rot(Psi)
    RotDot=rotvect_drotdpsib(Psi,PsiDot)
    RotPsiDot=matmul(Rot,PsiDot)

! Compute the current linear and angular velocity.
    VB= matmul(CBa,RaDot+matmul(rot_skew(Vrel(4:6)),Ra)+Vrel(1:3))
    WB= RotPsiDot + matmul(CBa,Vrel(4:6))

    Phat= matmul(ElemMass(:,1:3),VB)+matmul(ElemMass(:,4:6),WB)

    VBdot= matmul(CBa,RaDDot+rot_cross(Vrel(4:6),RaDot)+rot_cross(VrelDot(4:6),Ra)+VrelDot(1:3)) &
&          + rot_cross(VB,RotPsiDot)

    Vgyr(1:3)= matmul(CBa,matmul(rot_skew(Vrel(4:6)),RaDot)) + matmul(rot_skew(VB),matmul(Rot,PsiDot))
    Vgyr(4:6)= matmul(RotDot+matmul(rot_skew(WB),Rot),PsiDot)

    Pgyr= matmul(ElemMass,Vgyr)

! Operators in Kgyr from the derivatives of velocities.
    dVdq(1:3,1:3)= matmul(CBa,rot_skew(Vrel(4:6)))
    dVdq(1:3,4:6)= matmul(rot_skew(VB),Rot)
    dVdq(4:6,1:3)= 0.0d0
    dVdq(4:6,4:6)= RotDot + matmul(rot_skew(WB),Rot)

    dVgyrdq(1:3,1:3)= matmul(CBa,rot_skew(VrelDot(4:6))) &
&                   - matmul(rot_skew(RotPsiDot),matmul(CBa,rot_skew(Vrel(4:6))))
    dVgyrdq(1:3,4:6)= matmul(rot_skew(VBdot),Rot) + matmul(rot_skew(VB),RotDot)
    dVgyrdq(4:6,1:3)= 0.d0
    dVgyrdq(4:6,4:6)= rotvect_b1(Psi,PsiDot,PsiDot) + matmul(rot_skew(WB),rotvect_a1(Psi,PsiDot))    &
&                     - matmul(rot_skew(RotPsiDot),RotDot+matmul(rot_skew(WB),Rot))

! Operators in Kgyr from the derivatives of matrices ARC and AWRC.
    ARC=0.d0
    ARC(1:3,1:3)=CaB
    ARC(4:6,1:3)=matmul(rot_skew(Ra),CaB)
    ARC(4:6,4:6)=CaB

    dARCdqPgyr= 0.d0
    dARCdqPgyr(1:3,4:6)= -matmul(CaB,matmul(rot_skew(Pgyr(1:3)),Rot))
    dARCdqPgyr(4:6,1:3)= -matmul(CaB,matmul(rot_skew(Pgyr(1:3)),CBa))
    dARCdqPgyr(4:6,4:6)= -matmul(rot_skew(Ra),matmul(CaB,matmul(rot_skew(Pgyr(1:3)),Rot)))   &
&                        -matmul(CaB,matmul(rot_skew(Pgyr(4:6)),Rot))

    AWRC=0.d0
    AWRC(1:3,1:3) = matmul(CaB,rot_skew(WB))
    AWRC(4:6,1:3) = matmul(rot_skew(Ra),matmul(CaB,rot_skew(WB))) + matmul(CaB,rot_skew(VB))
    AWRC(4:6,4:6) = matmul(CaB,rot_skew(WB))

    dAWRCdqPhat= 0.d0
    dAWRCdqPhat(1:3,4:6)= -matmul(CaB,matmul(rot_skew(Phat(1:3)),RotDot)+matmul(rot_skew(WB),matmul(rot_skew(Phat(1:3)),Rot)))
    dAWRCdqPhat(4:6,1:3)= -matmul(CaB,matmul(matmul(rot_skew(WB),rot_skew(Phat(1:3))),CBa))          &
&                         +matmul(CaB,matmul(rot_skew(Phat(1:3)),matmul(rot_skew(WB),CBa)-matmul(CBa,rot_skew(Vrel(4:6)))))

    dAWRCdqPhat(4:6,4:6)= -matmul(rot_skew(Ra),matmul(CaB,matmul(rot_skew(Phat(1:3)),RotDot)                                   &
&                                                         +matmul(rot_skew(WB),matmul(rot_skew(Phat(1:3)),Rot))))              &
&                         -matmul(CaB,matmul(matmul(rot_skew(VB),rot_skew(Phat(1:3))),Rot))                                    &
&                         -matmul(CaB,matmul(rot_skew(Phat(4:6)),RotDot)+matmul(rot_skew(WB),matmul(rot_skew(Phat(4:6)),Rot)))

! Add the contribution of the current Gauss point to inertial stiffness matrix.
    Kgyr= Kgyr + WeightGauss(iGauss) * Jacobian &
&         * matmul((dARCdqPgyr + matmul(ARC,matmul(ElemMass,dVgyrdq)) + dAWRCdqPhat + matmul(AWRC,matmul(ElemMass,dVdq))),N)

  end do
  deallocate (CoordGauss,WeightGauss)
  return
 end subroutine xbeam_kgyr


 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-> Subroutine XBEAM_RBKGYR
!
!-> Description:
!
!    Compute the lumped mass contributions to the tangent gyroscopic
!    stiffness matrix
!
!-> Remarks.-
!
!   1) New values are added to those already in the matrix.
!   2) The contribution due to the derivative of the mass matrix with the
!      degrees of freedom are included (thus the dependency with the
!      accelerations).
!   3) The contribution due to the derivative of the mass matrix with the
!      degrees of freedom are included (thus the dependency with the
!      accelerations).
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine xbeam_rbkgyr (NumNodesElem,r0,Ri,RiDot,RiDDot,Vrel,VrelDot, &
&                         NodalMass,Kgyr)
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
  real(8) :: RaDDot(3)              ! Time derivatives of RaDot.
  real(8) :: RotPsiDot(3)           ! Rot*PsiDot.
  real(8) :: VB(3),WB(3)            ! Linear and angular velocity in the material frame.
  real(8) :: VBdot(3)               ! Time derivatives of VB.
  real(8) :: Vgyr(6)                ! Vectors with linear/angular velocities.
  real(8) :: Phat(6),Pgyr(6)        ! Vectors with linear/angular momenta

  real(8) :: CBa(3,3),CaB(3,3)      ! Coordinate transformation matrix.
  real(8) :: Rot(3,3),RotDot(3,3)   ! Tangential operator and time derivative.
  real(8) :: ARC(6,6),AWRC(6,6)     ! A_{RC} and A_{wRC}.
  real(8) :: dARCdqPgyr(6,6),dAWRCdqPhat(6,6)
  real(8) :: dVgyrdq  (6,6)         ! d(Vdot)/dq at Gauss point.
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

! Compute element shape functions at the corrent node.
    N =0.d0
    do i=1,6
      N (i,i+(iNode-1)*6)= 1.d0
    end do

! Compute the current coordinate transformation matrix and rotational operator.
    CBa= rotvect_psi2mat(Psi)
    CaB= transpose(CBa)
    Rot= rotvect_psi2rot(Psi)
    RotDot=rotvect_drotdpsib(Psi,PsiDot)
    RotPsiDot=matmul(Rot,PsiDot)

! Compute the current linear and angular velocity.
    VB= matmul(CBa,RaDot+matmul(rot_skew(Vrel(4:6)),Ra)+Vrel(1:3))
    WB= RotPsiDot + matmul(CBa,Vrel(4:6))

    Phat= matmul(NodalMass(iNode,:,1:3),VB)+matmul(NodalMass(iNode,:,4:6),WB)

    VBdot= matmul(CBa,RaDDot+rot_cross(Vrel(4:6),RaDot)+rot_cross(VrelDot(4:6),Ra)+VrelDot(1:3)) &
&          + rot_cross(VB,RotPsiDot)

    Vgyr(1:3)= matmul(CBa,matmul(rot_skew(Vrel(4:6)),RaDot)) + matmul(rot_skew(VB),matmul(Rot,PsiDot))
    Vgyr(4:6)= matmul(RotDot+matmul(rot_skew(WB),Rot),PsiDot)

    Pgyr= matmul(NodalMass(iNode,:,:),Vgyr)

! Operators in Kgyr from the derivatives of velocities.
    dVdq(1:3,1:3)= matmul(CBa,rot_skew(Vrel(4:6)))
    dVdq(1:3,4:6)= matmul(rot_skew(VB),Rot)
    dVdq(4:6,1:3)= 0.0d0
    dVdq(4:6,4:6)= RotDot + matmul(rot_skew(WB),Rot)

    dVgyrdq(1:3,1:3)= matmul(CBa,rot_skew(VrelDot(4:6))) &
&                   - matmul(rot_skew(RotPsiDot),matmul(CBa,rot_skew(Vrel(4:6))))
    dVgyrdq(1:3,4:6)= matmul(rot_skew(VBdot),Rot) + matmul(rot_skew(VB),RotDot)
    dVgyrdq(4:6,1:3)= 0.d0
    dVgyrdq(4:6,4:6)= rotvect_b1(Psi,PsiDot,PsiDot) + matmul(rot_skew(WB),rotvect_a1(Psi,PsiDot))    &
&                     - matmul(rot_skew(RotPsiDot),RotDot+matmul(rot_skew(WB),Rot))

! Operators in Kgyr from the derivatives of matrices ARC and AWRC.
    ARC=0.d0
    ARC(1:3,1:3)=CaB
    ARC(4:6,1:3)=matmul(rot_skew(Ra),CaB)
    ARC(4:6,4:6)=CaB

    dARCdqPgyr= 0.d0
    dARCdqPgyr(1:3,4:6)= -matmul(CaB,matmul(rot_skew(Pgyr(1:3)),Rot))
    dARCdqPgyr(4:6,1:3)= -matmul(CaB,matmul(rot_skew(Pgyr(1:3)),CBa))
    dARCdqPgyr(4:6,4:6)= -matmul(rot_skew(Ra),matmul(CaB,matmul(rot_skew(Pgyr(1:3)),Rot)))   &
&                        -matmul(CaB,matmul(rot_skew(Pgyr(4:6)),Rot))

    AWRC=0.d0
    AWRC(1:3,1:3) = matmul(CaB,rot_skew(WB))
    AWRC(4:6,1:3) = matmul(rot_skew(Ra),matmul(CaB,rot_skew(WB))) + matmul(CaB,rot_skew(VB))
    AWRC(4:6,4:6) = matmul(CaB,rot_skew(WB))

    dAWRCdqPhat= 0.d0
    dAWRCdqPhat(1:3,4:6)= -matmul(CaB,matmul(rot_skew(Phat(1:3)),RotDot)+matmul(rot_skew(WB),matmul(rot_skew(Phat(1:3)),Rot)))
    dAWRCdqPhat(4:6,1:3)= -matmul(CaB,matmul(matmul(rot_skew(WB),rot_skew(Phat(1:3))),CBa))          &
&                         +matmul(CaB,matmul(rot_skew(Phat(1:3)),matmul(rot_skew(WB),CBa)-matmul(CBa,rot_skew(Vrel(4:6)))))

    dAWRCdqPhat(4:6,4:6)= -matmul(rot_skew(Ra),matmul(CaB,matmul(rot_skew(Phat(1:3)),RotDot)                                   &
&                                                         +matmul(rot_skew(WB),matmul(rot_skew(Phat(1:3)),Rot))))              &
&                         -matmul(CaB,matmul(matmul(rot_skew(VB),rot_skew(Phat(1:3))),Rot))                                    &
&                         -matmul(CaB,matmul(rot_skew(Phat(4:6)),RotDot)+matmul(rot_skew(WB),matmul(rot_skew(Phat(4:6)),Rot)))

! Add the contribution of the current Gauss point to inertial stiffness matrix.
    Kgyr= Kgyr + matmul((dARCdqPgyr + matmul(ARC, matmul(NodalMass(iNode,:,:),dVgyrdq)) &
&                     + dAWRCdqPhat + matmul(AWRC,matmul(NodalMass(iNode,:,:),dVdq))),N)

  end do
  return
 end subroutine xbeam_rbkgyr


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-> Subroutine XBEAM_CGYR
!
!-> Description:
!
!    Compute tangent gyroscopic damping matrix
!
!-> Remarks.-
!
!   1) New values are added to those already in the matrix.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine xbeam_cgyr (NumNodesElem,r0,Ri,RiDot,Vrel,ElemMass,Cgyr,NumGauss)
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

  real(8) :: CBa(3,3),CaB(3,3)      ! Coordinate transformation matrix.
  real(8) :: Rot(3,3)               ! Tangential operator.
  real(8) :: dVgyrdqdot(6,6)        ! d(Vdot)/d(qdot) at Gauss point.
  real(8) :: dVdqdot   (6,6)        ! dV/d(qdot) at Gauss point.
  real(8) :: ARC(6,6),AWRC(6,6),dAWRCdqdotPhat(6,6)
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
    CaB= transpose(CBa)
    Rot= rotvect_psi2rot  (Psi)

! Compute the current linear and angular velocity.
    VB= matmul(CBa,RaDot+matmul(rot_skew(Vrel(4:6)),Ra)+Vrel(1:3))
    WB= matmul(Rot,PsiDot) + matmul(CBa,Vrel(4:6))

    PB= matmul(ElemMass(1:3,1:3),VB)+matmul(ElemMass(1:3,4:6),WB)
    HB= matmul(ElemMass(4:6,1:3),VB)+matmul(ElemMass(4:6,4:6),WB)

!Operators in Cgyr.
    dVdqdot=0.0d0
    dVdqdot(1:3,1:3)= CBa
    dVdqdot(4:6,4:6)= Rot

    dVgyrdqdot(1:3,1:3)= 2.d0*matmul(CBa,rot_skew(Vrel(4:6)))-matmul(rot_skew(WB),CBa)
    dVgyrdqdot(1:3,4:6)= matmul(rot_skew(VB),Rot)
    dVgyrdqdot(4:6,1:3)= 0.d0
    dVgyrdqdot(4:6,4:6)= 2.d0*rotvect_drotdpsib(Psi,PsiDot)+matmul(rot_skew(WB),Rot)

    ARC=0.d0
    ARC(1:3,1:3)=CaB
    ARC(4:6,1:3)=matmul(rot_skew(Ra),CaB)
    ARC(4:6,4:6)=CaB

    AWRC=0.d0
    AWRC(1:3,1:3) = matmul(CaB,rot_skew(WB))
    AWRC(4:6,1:3) = matmul(rot_skew(Ra),matmul(CaB,rot_skew(WB))) + matmul(CaB,rot_skew(VB))
    AWRC(4:6,4:6) = matmul(CaB,rot_skew(WB))

    dAWRCdqdotPhat(1:3,1:3)= 0.d0
    dAWRCdqdotPhat(1:3,4:6)= -matmul(CaB,matmul(rot_skew(PB),Rot))
    dAWRCdqdotPhat(4:6,1:3)= -rot_skew(matmul(CaB,PB))
    dAWRCdqdotPhat(4:6,4:6)= -matmul(rot_skew(Ra),matmul(CaB,matmul(rot_skew(PB),Rot)))   &
&                            -matmul(CaB,matmul(rot_skew(HB),Rot))

! Compute strain matrix operator and material tangent stiffness.
    Cgyr= Cgyr + WeightGauss(iGauss) * Jacobian   &
&                * matmul((matmul(ARC,matmul(ElemMass,dVgyrdqdot)) + dAWRCdqdotPhat + matmul(AWRC,matmul(ElemMass,dVdqdot))),N)

  end do
  deallocate (CoordGauss,WeightGauss)
  return
 end subroutine xbeam_cgyr


 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-> Subroutine XBEAM_RBCGYR
!
!-> Description:
!
!    Compute the lumped mass contributions to the tangent gyroscopic damping
!    matrix
!
!-> Remarks.-
!
!   1) New values are added to those already in the matrix.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine xbeam_rbcgyr (NumNodesElem,r0,Ri,RiDot,Vrel,NodalMass,Cgyr)
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
  integer :: i,iNode                ! Counters.
  real(8) :: Ra(3),RaDot(3)         ! Displacement and its derivative at Gauss point.
  real(8) :: Psi(3),PsiDot(3)       ! Rotation vector and its derivative at Gauss point.
  real(8) :: VB(3),WB(3)            ! Linear and angular velocity in the material frame..
  real(8) :: PB(3),HB(3)            ! Linear and angular momenta in the material frame..

  real(8) :: CBa(3,3),CaB(3,3)      ! Coordinate transformation matrix.
  real(8) :: Rot(3,3)               ! Tangential operator.
  real(8) :: dVgyrdqdot(6,6)        ! d(Vdot)/d(qdot) at Gauss point.
  real(8) :: dVdqdot   (6,6)        ! dV/d(qdot) at Gauss point.
  real(8) :: ARC(6,6),AWRC(6,6),dAWRCdqdotPhat(6,6)
  real(8) :: N (6,6*MaxNodCB3)      ! Element shape function matrix.

! Loop in the element ndoes.
  do iNode=1,NumNodesElem

! Compute the current position vector and rotations, and their derivatives.
		Ra     = Ri    (iNode,1:3)
		RaDot  = RiDot (iNode,1:3)
		Psi    = Ri    (iNode,4:6)
		PsiDot = RiDot (iNode,4:6)

! Compute element shape functions at the corrent node.
    N =0.d0
    do i=1,6
      N (i,i+(iNode-1)*6)= 1.d0
    end do

! Compute the current coordinate transformation matrix and rotational operator.
    CBa= rotvect_psi2mat  (Psi)
    CaB= transpose(CBa)
    Rot= rotvect_psi2rot  (Psi)

! Compute the current linear and angular velocity.
    VB= matmul(CBa,RaDot+matmul(rot_skew(Vrel(4:6)),Ra)+Vrel(1:3))
    WB= matmul(Rot,PsiDot) + matmul(CBa,Vrel(4:6))

    PB= matmul(NodalMass(iNode,1:3,1:3),VB) + matmul(NodalMass(iNode,1:3,4:6),WB)
    HB= matmul(NodalMass(iNode,4:6,1:3),VB) + matmul(NodalMass(iNode,4:6,4:6),WB)

!Operators in Cgyr.
    dVdqdot=0.0d0
    dVdqdot(1:3,1:3)= CBa
    dVdqdot(4:6,4:6)= Rot

    dVgyrdqdot(1:3,1:3)= 2.d0*matmul(CBa,rot_skew(Vrel(4:6)))-matmul(rot_skew(WB),CBa)
    dVgyrdqdot(1:3,4:6)= matmul(rot_skew(VB),Rot)
    dVgyrdqdot(4:6,1:3)= 0.d0
    dVgyrdqdot(4:6,4:6)= 2.d0*rotvect_drotdpsib(Psi,PsiDot)+matmul(rot_skew(WB),Rot)

    ARC=0.d0
    ARC(1:3,1:3)=CaB
    ARC(4:6,1:3)=matmul(rot_skew(Ra),CaB)
    ARC(4:6,4:6)=CaB

    AWRC=0.d0
    AWRC(1:3,1:3) = matmul(CaB,rot_skew(WB))
    AWRC(4:6,1:3) = matmul(rot_skew(Ra),matmul(CaB,rot_skew(WB))) + matmul(CaB,rot_skew(VB))
    AWRC(4:6,4:6) = matmul(CaB,rot_skew(WB))

    dAWRCdqdotPhat(1:3,1:3)= 0.d0
    dAWRCdqdotPhat(1:3,4:6)= -matmul(CaB,matmul(rot_skew(PB),Rot))
    dAWRCdqdotPhat(4:6,1:3)= -rot_skew(matmul(CaB,PB))
    dAWRCdqdotPhat(4:6,4:6)= -matmul(rot_skew(Ra),matmul(CaB,matmul(rot_skew(PB),Rot)))   &
&                            -matmul(CaB,matmul(rot_skew(HB),Rot))

! Compute strain matrix operator and material tangent stiffness.
    Cgyr= Cgyr + matmul((matmul(ARC,matmul(NodalMass(iNode,:,:),dVgyrdqdot)) &
&              + dAWRCdqdotPhat + matmul(AWRC,matmul(NodalMass(iNode,:,:),dVdqdot))),N)

  end do

  return
 end subroutine xbeam_rbcgyr


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-> Subroutine xbeam_FGYR
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
 subroutine xbeam_fgyr (NumNodesElem,r0,Ri,RiDot,Vrel,ElemMass,Qgyr,NumGauss)
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
  integer :: i                      ! Counters.
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
  real(8) :: Vhat(6),Vhatgyr(6)     ! Vectors with linear/angular velocities.

  real(8) :: CBa(3,3),CaB(3,3)      ! Coordinate transformation matrix.
  real(8) :: Rot(3,3),RotDot(3,3)   ! Tangential operator and time derivative.
  real(8) :: ARC(6,6),AWRC(6,6)     ! A_{RC} and A_{wRC} matrices.

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

! Compute the current coordinate transformation matrix and rotational operator.
    CBa    = rotvect_psi2mat(Psi)
    CaB    = transpose(CBa)
    Rot    = rotvect_psi2rot(Psi)
    RotDot = rotvect_drotdpsib(Psi,PsiDot)

! Compute the current linear and angular velocity.
    VB= matmul(CBa,RaDot+matmul(rot_skew(Vrel(4:6)),Ra)+Vrel(1:3))
    WB= matmul(Rot,PsiDot) + matmul(CBa,Vrel(4:6))

    VBgyr= matmul(CBa,matmul(rot_skew(Vrel(4:6)),RaDot)) + matmul(rot_skew(VB),matmul(Rot,PsiDot))
    WBgyr= matmul((RotDot+matmul(rot_skew(WB),Rot)),PsiDot)

! Operators in Qgyr
    ARC=0.d0
    ARC(1:3,1:3) = CaB
    ARC(4:6,1:3) = matmul(rot_skew(Ra),CaB)
    ARC(4:6,4:6) = CaB

    AWRC=0.d0
    AWRC(1:3,1:3) = matmul(CaB,rot_skew(WB))
    AWRC(4:6,1:3) = matmul(rot_skew(Ra),matmul(CaB,rot_skew(WB))) + matmul(CaB,rot_skew(VB))
    AWRC(4:6,4:6) = matmul(CaB,rot_skew(WB))

    Vhat(1:3)=VB
    Vhat(4:6)=WB

    Vhatgyr(1:3)=VBgyr
    Vhatgyr(4:6)=WBgyr

! Add the contribution of the current Gauss point to inertial force.
    Qgyr= Qgyr + WeightGauss(iGauss) * Jacobian * (matmul(ARC,matmul(ElemMass,Vhatgyr)) + matmul(AWRC,matmul(ElemMass,Vhat)))

  end do

  deallocate (CoordGauss,WeightGauss)
  return
 end subroutine xbeam_fgyr


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-> Subroutine xbeam_RBFGYR
!
!-> Description:
!
!    Lumped mass contribution to the generalized gyroscopic forces.
!
!-> Remarks.-
!
!   1) New values are added to those already in the matrix.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine xbeam_rbfgyr (NumNodesElem,r0,Ri,RiDot,Vrel,NodalMass,Qgyr)
  use lib_rot
  use lib_rotvect

! I/O Variables.
  integer,intent(in)    :: NumNodesElem     ! Number of nodes in the element.
  real(8),intent(in)    :: r0       (:,:)   ! Initial position/orientation of grid points.
  real(8),intent(in)    :: Ri       (:,:)   ! Current position/orientation of grid points.
  real(8),intent(in)    :: RiDot    (:,:)   ! Current time derivative of position/CRV of points.
  real(8),intent(in)    :: Vrel     (6)     ! Linear/angular velocity of reference frame.
  real(8),intent(in)    :: NodalMass(:,:,:) ! Inertia tensor of lumped masses at element nodes.
  real(8),intent(inout) :: Qgyr     (:)     ! Vector of discrete gyroscopic forces.

! Local variables
  integer :: iNode                  ! Counters.
  real(8) :: Ra(3),Psi(3)           ! Displacement and Rotation vector at Gauss point.
  real(8) :: RaDot(3),PsiDot(3)     ! Time derivatives of Ra and Psi.
  real(8) :: VB(3),WB(3)            ! Linear and angular velocity in the material frame.
  real(8)  ::VBgyr(3),WBgyr(3)      ! Time derivatives of VB and WB.
  real(8) :: Vhat(6),Vhatgyr(6)     ! Vectors with linear/angular velocities.

  real(8) :: CBa(3,3),CaB(3,3)      ! Coordinate transformation matrix.
  real(8) :: Rot(3,3),RotDot(3,3)   ! Tangential operator and time derivative.
  real(8) :: ARC(6,6),AWRC(6,6)     ! A_{RC} and A_{wRC} matrices.

! Loop in the element ndoes.
  do iNode=1,NumNodesElem

! Compute the current position vector and rotations, and their derivatives.
	Ra     = Ri    (iNode,1:3)
	RaDot  = RiDot (iNode,1:3)
	Psi    = Ri    (iNode,4:6)
	PsiDot = RiDot (iNode,4:6)

! Compute the current coordinate transformation matrix and rotational operator.
    CBa    = rotvect_psi2mat(Psi)
    CaB    = transpose(CBa)
    Rot    = rotvect_psi2rot(Psi)
    RotDot = rotvect_drotdpsib(Psi,PsiDot)

! Compute the current linear and angular velocity.
    VB= matmul(CBa,RaDot+matmul(rot_skew(Vrel(4:6)),Ra)+Vrel(1:3))
    WB= matmul(Rot,PsiDot) + matmul(CBa,Vrel(4:6))

    VBgyr= matmul(CBa,matmul(rot_skew(Vrel(4:6)),RaDot)) + matmul(rot_skew(VB),matmul(Rot,PsiDot))
    WBgyr= matmul((RotDot+matmul(rot_skew(WB),Rot)),PsiDot)

! Operators in Qgyr
    ARC=0.d0
    ARC(1:3,1:3) = CaB
    ARC(4:6,1:3) = matmul(rot_skew(Ra),CaB)
    ARC(4:6,4:6) = CaB

    AWRC=0.d0
    AWRC(1:3,1:3) = matmul(CaB,rot_skew(WB))
    AWRC(4:6,1:3) = matmul(rot_skew(Ra),matmul(CaB,rot_skew(WB))) + matmul(CaB,rot_skew(VB))
    AWRC(4:6,4:6) = matmul(CaB,rot_skew(WB))

    Vhat(1:3)=VB
    Vhat(4:6)=WB

    Vhatgyr(1:3)=VBgyr
    Vhatgyr(4:6)=WBgyr

! Add the contribution of the current Gauss point to inertial force.
    Qgyr= Qgyr +(matmul(ARC, matmul(NodalMass(iNode,:,:),Vhatgyr)) &
&              + matmul(AWRC,matmul(NodalMass(iNode,:,:),Vhat)))
  end do

  return
 end subroutine xbeam_rbfgyr


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
 subroutine xbeam_mrr (NumNodesElem,r0,Ri,ElemMass,MRR,NumGauss)
  use lib_fem
  use lib_rot
  use lib_rotvect

! I/O Variables.
  integer,intent(in)   :: NumNodesElem      ! Number of nodes in the element.
  real(8),intent(in)   :: r0       (:,:)    ! Initial position/orientation of grid points.
  real(8),intent(in)   :: Ri       (:,:)    ! Current position/orientation of grid points.
  real(8),intent(in)   :: ElemMass (:,:)    ! Inertial properties in the element.
  real(8),intent(inout):: MRR     (:,:)     ! Tangent mass matrix.
  integer,intent(in)   :: NumGauss          ! Number of Gauss points in the element.

! Local variables.
  integer :: i                          ! Counters.
  integer :: iGauss                     ! Counter on the Gaussian points.
  real(8) :: Jacobian                   ! ds/deta, with eta the nondimensional arclength.
  real(8),allocatable :: CoordGauss (:) ! Coords of Gauss points.
  real(8),allocatable :: WeightGauss(:) ! Coords of Gauss points.
  real(8) :: ShapeFun(MaxNodCB3)        ! Shape functions in a gauss point.
  real(8) :: ShapeDer(MaxNodCB3)        ! Derivatives of ShapeFun in a gauss point.

  real(8) :: Ri_g (6)              ! Displacement/rotation at Gauss point.
  real(8) :: dr0_g(6)              ! Derivative of r0_g.
  real(8) :: CaB(3,3)              ! Coordinate transformation matrix.
  real(8) :: ARC(6,6)              ! A_{RC} matrix

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

! Compute the current coordinate transformation matrix.
    CaB= transpose(rotvect_psi2mat(Ri_g(4:6)))

!Operators for the reference mass matrix.
    ARC=0.d0
    ARC(1:3,1:3)=CaB
    ARC(4:6,1:3)=matmul(rot_skew(Ri_g(1:3)),CaB)
    ARC(4:6,4:6)=CaB

! Compute reference mass matrix.
    MRR= MRR + WeightGauss(iGauss)*Jacobian*matmul(ARC,matmul(ElemMass,transpose(ARC)))
  end do

  deallocate (CoordGauss,WeightGauss)

  return
 end subroutine xbeam_mrr


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-> Subroutine XBEAM_RBMRR
!
!-> Description:
!
!    Compute the lumped mass contribution to the mass matrix corresponding
!    to the reference system accelerations.
!
!-> Remarks.-
!
!   1) New values are added to those already in the matrix.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine xbeam_rbmrr (NumNodesElem,r0,Ri,NodalMass,MRR)
  use lib_rot
  use lib_rotvect

! I/O Variables.
  integer,intent(in)   :: NumNodesElem     ! Number of nodes in the element.
  real(8),intent(in)   :: r0       (:,:)   ! Initial position/orientation of grid points.
  real(8),intent(in)   :: Ri       (:,:)   ! Current position/orientation of grid points.
  real(8),intent(in)   :: NodalMass(:,:,:) ! Inertia tensor of lumped masses at element nodes.
  real(8),intent(inout):: MRR     (:,:)    ! Tangent mass matrix.

! Local variables.
  integer :: iNode                 ! Counters.
  real(8) :: CaB(3,3)              ! Coordinate transformation matrix.
  real(8) :: ARC(6,6)              ! A_{RC} matrix

! Loop in the element ndoes.
  do iNode=1,NumNodesElem

! Compute the current coordinate transformation matrix.
    CaB= transpose(rotvect_psi2mat(Ri(iNode,4:6)))

!Operators for the reference mass matrix.
    ARC=0.d0
    ARC(1:3,1:3)=CaB
    ARC(4:6,1:3)=matmul(rot_skew(Ri(iNode,1:3)),CaB)
    ARC(4:6,4:6)=CaB

! Compute reference mass matrix.
    MRR= MRR + matmul(ARC,matmul(NodalMass(iNode,:,:),transpose(ARC)))
  end do
  return
 end subroutine xbeam_rbmrr


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-> Subroutine XBEAM_CRR
!
!-> Description:
!
!    Compute tangent damping matrix CRR for motions of reference system.
!
!-> Remarks.-
!
!   1) New values are added to those already in the matrix.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine xbeam_crr (NumNodesElem,r0,Ri,RiDot,Vrel,ElemMass,CRR,NumGauss)
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
  real(8),intent(inout) :: CRR     (:,:)    ! Tangent gyroscopic damping matrix.
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
  real(8) :: Ra(3),RaDot(3)         ! Displacement and its derivative at Gauss point.
  real(8) :: Psi(3),PsiDot(3)       ! Rotation vector and its derivative at Gauss point.
  real(8) :: VB(3),WB(3)            ! Linear and angular velocity in the material frame..
  real(8) :: PB(3),HB(3)            ! Linear and angular momenta in the material frame..

  real(8) :: CBa(3,3),CaB(3,3)      ! Coordinate transformation matrix.
  real(8) :: Rot(3,3)               ! Tangential operator.
  real(8) :: ARC(6,6),AWRC(6,6)
  real(8) :: dAWRCdvrelPhat(6,6)
  real(8) :: dVgyrdvrel(6,6)        ! d(Vdot)/d(vrel) at Gauss point.
  real(8) :: dVdvrel   (6,6)        ! dV/d(vrel) at Gauss point.

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

! Compute the current coordinate transformation matrix and rotational operator.
    CBa= rotvect_psi2mat  (Psi)
    CaB= transpose(CBa)
    Rot= rotvect_psi2rot  (Psi)

! Compute the current linear and angular velocity.
    VB= matmul(CBa,RaDot+matmul(rot_skew(Vrel(4:6)),Ra)+Vrel(1:3))
    WB= matmul(Rot,PsiDot) + matmul(CBa,Vrel(4:6))

    PB= matmul(ElemMass(1:3,1:3),VB)+matmul(ElemMass(1:3,4:6),WB)
    HB= matmul(ElemMass(4:6,1:3),VB)+matmul(ElemMass(4:6,4:6),WB)

!Operators in CRR.
    ARC=0.d0
    ARC(1:3,1:3)=CaB
    ARC(4:6,1:3)=matmul(rot_skew(Ra),CaB)
    ARC(4:6,4:6)=CaB

    AWRC=0.d0
    AWRC(1:3,1:3) = matmul(CaB,rot_skew(WB))
    AWRC(4:6,1:3) = matmul(rot_skew(Ra),matmul(CaB,rot_skew(WB))) + matmul(CaB,rot_skew(VB))
    AWRC(4:6,4:6) = matmul(CaB,rot_skew(WB))

    dVdvrel=transpose(ARC)

    dVgyrdvrel(1:3,1:3)= -matmul(rot_skew(WB),CBa) + matmul(CBa,rot_skew(Vrel(4:6)))
    dVgyrdvrel(1:3,4:6)= -matmul(CBa,rot_skew(Radot))+matmul(rot_skew(WB),matmul(CBa,rot_skew(Ra))) &
&                        -matmul(CBa,matmul(rot_skew(Vrel(4:6)),rot_skew(Ra)))
    dVgyrdvrel(4:6,1:3)=  0.d0
    dVgyrdvrel(4:6,4:6)= -matmul(rot_skew(WB),CBa) + matmul(CBa,rot_skew(Vrel(4:6)))

    dAWRCdvrelPhat=0.d0
    dAWRCdvrelPhat(1:3,4:6)= - matmul(CaB,matmul(rot_skew(PB),CBa))
    dAWRCdvrelPhat(4:6,1:3)= - matmul(CaB,matmul(rot_skew(PB),CBa))
    dAWRCdvrelPhat(4:6,4:6)= - matmul(CaB,matmul(rot_skew(HB),CBa)) &
&                            + matmul(matmul(CaB,matmul(rot_skew(PB),CBa)),rot_skew(Ra)) &
&                            - matmul(rot_skew(Ra),matmul(CaB,matmul(rot_skew(PB),CBa)))

! Compute strain matrix operator and material tangent stiffness.
    CRR= CRR + WeightGauss(iGauss) * Jacobian                                   &
&              * (matmul(ARC,matmul(ElemMass,dVgyrdvrel)) + dAWRCdvrelPhat + matmul(AWRC,matmul(ElemMass,dVdvrel)))

  end do
  deallocate (CoordGauss,WeightGauss)
  return
 end subroutine xbeam_crr


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-> Subroutine XBEAM_RBCRR
!
!-> Description:
!
!    Compute the lumped mass contribution to the tangent damping matrix
!    CRR for motions of reference system.
!
!-> Remarks.-
!
!   1) New values are added to those already in the matrix.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine xbeam_rbcrr (NumNodesElem,r0,Ri,RiDot,Vrel,NodalMass,CRR)
  use lib_rot
  use lib_rotvect

! I/O Variables.
  integer,intent(in)    :: NumNodesElem     ! Number of nodes in the element.
  real(8),intent(in)    :: r0       (:,:)   ! Initial position/orientation of grid points.
  real(8),intent(in)    :: Ri       (:,:)   ! Current position/orientation of grid points.
  real(8),intent(in)    :: RiDot    (:,:)   ! Current time derivative of position/CRV of points.
  real(8),intent(in)    :: Vrel     (6)     ! Linear/angular velocity of reference frame.
  real(8),intent(in)    :: NodalMass(:,:,:) ! Inertia tensor of lumped masses at element nodes.
  real(8),intent(inout) :: CRR     (:,:)    ! Tangent gyroscopic damping matrix.

! Local variables
  integer :: iNode                  ! Counters.
  real(8) :: Ra(3),RaDot(3)         ! Displacement and its derivative at Gauss point.
  real(8) :: Psi(3),PsiDot(3)       ! Rotation vector and its derivative at Gauss point.
  real(8) :: VB(3),WB(3)            ! Linear and angular velocity in the material frame..
  real(8) :: PB(3),HB(3)            ! Linear and angular momenta in the material frame..

  real(8) :: CBa(3,3),CaB(3,3)      ! Coordinate transformation matrix.
  real(8) :: Rot(3,3)               ! Tangential operator.
  real(8) :: ARC(6,6),AWRC(6,6)
  real(8) :: dAWRCdvrelPhat(6,6)
  real(8) :: dVgyrdvrel(6,6)        ! d(Vdot)/d(vrel) at Gauss point.
  real(8) :: dVdvrel   (6,6)        ! dV/d(vrel) at Gauss point.

! Loop in the element ndoes.
  do iNode=1,NumNodesElem

! Compute the current position vector and rotations, and their derivatives.
	Ra     = Ri    (iNode,1:3)
	RaDot  = RiDot (iNode,1:3)
	Psi    = Ri    (iNode,4:6)
	PsiDot = RiDot (iNode,4:6)

! Compute the current coordinate transformation matrix and rotational operator.
    CBa= rotvect_psi2mat  (Psi)
    CaB= transpose(CBa)
    Rot= rotvect_psi2rot  (Psi)

! Compute the current linear and angular velocity.
    VB= matmul(CBa,RaDot+matmul(rot_skew(Vrel(4:6)),Ra)+Vrel(1:3))
    WB= matmul(Rot,PsiDot) + matmul(CBa,Vrel(4:6))

    PB= matmul(NodalMass(iNode,1:3,1:3),VB) + matmul(NodalMass(iNode,1:3,4:6),WB)
    HB= matmul(NodalMass(iNode,4:6,1:3),VB) + matmul(NodalMass(iNode,4:6,4:6),WB)

!Operators in CRR.
    ARC=0.d0
    ARC(1:3,1:3)=CaB
    ARC(4:6,1:3)=matmul(rot_skew(Ra),CaB)
    ARC(4:6,4:6)=CaB

    AWRC=0.d0
    AWRC(1:3,1:3) = matmul(CaB,rot_skew(WB))
    AWRC(4:6,1:3) = matmul(rot_skew(Ra),matmul(CaB,rot_skew(WB))) + matmul(CaB,rot_skew(VB))
    AWRC(4:6,4:6) = matmul(CaB,rot_skew(WB))

    dVdvrel=transpose(ARC)

    dVgyrdvrel(1:3,1:3)= -matmul(rot_skew(WB),CBa) + matmul(CBa,rot_skew(Vrel(4:6)))
    dVgyrdvrel(1:3,4:6)= -matmul(CBa,rot_skew(Radot))+matmul(rot_skew(WB),matmul(CBa,rot_skew(Ra))) &
&                        -matmul(CBa,matmul(rot_skew(Vrel(4:6)),rot_skew(Ra)))
    dVgyrdvrel(4:6,1:3)=  0.d0
    dVgyrdvrel(4:6,4:6)= -matmul(rot_skew(WB),CBa) + matmul(CBa,rot_skew(Vrel(4:6)))

    dAWRCdvrelPhat=0.d0
    dAWRCdvrelPhat(1:3,4:6)= - matmul(CaB,matmul(rot_skew(PB),CBa))
    dAWRCdvrelPhat(4:6,1:3)= - matmul(CaB,matmul(rot_skew(PB),CBa))
    dAWRCdvrelPhat(4:6,4:6)= - matmul(CaB,matmul(rot_skew(HB),CBa)) &
&                            + matmul(matmul(CaB,matmul(rot_skew(PB),CBa)),rot_skew(Ra)) &
&                            - matmul(rot_skew(Ra),matmul(CaB,matmul(rot_skew(PB),CBa)))

! Compute strain matrix operator and material tangent stiffness.
    CRR= CRR +(matmul(ARC, matmul(NodalMass(iNode,:,:),dVgyrdvrel)) + dAWRCdvrelPhat      &
&            + matmul(AWRC,matmul(NodalMass(iNode,:,:),dVdvrel)))

  end do
  return
 end subroutine xbeam_rbcrr


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
 subroutine xbeam_fext (NumNodesElem,Ri,Flags,Fmat,FollowerForce,FollowerForceRig,Cao)
  use lib_fem
  use lib_rot
  use lib_rotvect

! I/O Variables.
  integer,intent(in)        :: NumNodesElem     ! Number of nodes in the element.
  real(8),intent(in)        :: Ri    (:,:)      ! Current position/orientation of grid points.
  logical,intent(in)        :: Flags (:)        ! Identify master nodes.
  real(8),intent(inout)     :: Fmat  (:,:)      ! Forces/moments on the element nodes.
  logical,intent(in)        :: FollowerForce    ! =T if follower forces.
  logical,intent(in)        :: FollowerForceRig ! =T if follower force in the body-fixed frame.
  real(8),intent(in)        :: Cao   (:,:)      ! Rotation operator from inertial frame to reference frame

! Local variables.
  integer :: iNode                 ! Counter on the element nodes.
  real(8) :: CaB(3,3)              ! Transformation matrix from B to a
  real(8) :: ARC(6,6)              ! Operators to project the forces on the dofs.

!--> Case of follower forces in the body-fixed frame.
  if (FollowerForceRig) then
    !--> Case of follower forces on the structure.
    if (FollowerForce) then
      ! Loop in the element nodes.
      do iNode=1,NumNodesElem
         if (Flags(iNode)) then
           ! Operators to project the forces on the dofs.
           CaB= transpose(rotvect_psi2mat(Ri(iNode,4:6)))
           ARC = 0.d0
           ARC(1:3,1:3)= CaB
           ARC(4:6,1:3)= matmul(rot_skew(Ri(iNode,1:3)),CaB)
           ARC(4:6,4:6)= CaB                                            !rotvect_psi2rot(Ri(iNode,4:6))

           ! Compute influence coefficients for follower forces.
           Fmat(:,6*(iNode-1)+1:6*iNode)= Fmat(:,6*(iNode-1)+1:6*iNode) + ARC
         end if
      end do

    !--> Case of dead forces on the structure but following the (a) frame.
    else
      do iNode=1,NumNodesElem
         if (Flags(iNode)) then
           ! Operators to project the forces on the dofs.
           ARC = 0.d0
           ARC(1:3,1:3)= Unit
           ARC(4:6,1:3)= rot_skew(Ri(iNode,1:3))
           ARC(4:6,4:6)= Unit

           ! Compute influence coefficients for follower forces.
           Fmat(:,6*(iNode-1)+1:6*iNode)= Fmat(:,6*(iNode-1)+1:6*iNode) + ARC
         end if
      end do
    end if

  !--> Case of dead forces in the inertial frame. Not possible to have structural follower forces and rigid-body dead loads.
  else
    do iNode=1,NumNodesElem
       if (Flags(iNode)) then
         ! Operators to project the forces on the dofs.
         ARC = 0.d0
         ARC(1:3,1:3)= Cao
         ARC(4:6,1:3)= matmul(rot_skew(Ri(iNode,1:3)),Cao)
         ARC(4:6,4:6)= Cao

         ! Compute influence coefficients for follower forces.
         Fmat(:,6*(iNode-1)+1:6*iNode)= Fmat(:,6*(iNode-1)+1:6*iNode) + ARC
       end if
    end do
  end if

  return
 end subroutine xbeam_fext

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!->Function XBEAM_ROT
!
!->Description.-
!
!   Definition of rotation operator in Euler quaternions
!   See Aircraft Control and Simulation by Stevens, Lewis
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 function xbeam_Rot (q)
! I/O Variables.
  real(8),  intent(inout):: q(:)           ! Quaternion
  real(8)             :: xbeam_Rot(3,3)

  xbeam_Rot = 0.d0

  xbeam_Rot(1,1) = q(1)*q(1) + q(2)*q(2) - q(3)*q(3) - q(4)*q(4)
  xbeam_Rot(2,2) = q(1)*q(1) - q(2)*q(2) + q(3)*q(3) - q(4)*q(4)
  xbeam_Rot(3,3) = q(1)*q(1) - q(2)*q(2) - q(3)*q(3) + q(4)*q(4)

  xbeam_Rot(1,2) = 2*(q(2)*q(3) + q(1)*q(4))
  xbeam_Rot(2,1) = 2*(q(2)*q(3) - q(1)*q(4))

  xbeam_Rot(1,3) = 2*(q(2)*q(4) - q(1)*q(3))
  xbeam_Rot(3,1) = 2*(q(2)*q(4) + q(1)*q(3))

  xbeam_Rot(2,3) = 2*(q(3)*q(4) + q(1)*q(2))
  xbeam_Rot(3,2) = 2*(q(3)*q(4) - q(1)*q(2))

  return
 end function xbeam_Rot


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!->Function XBEAM_QUADSKEW
!
!->Description.-
!
!  Quaternion ODE to compute orientation of body-fixed frame a.
!  See Shearer and Cesnik (2007) for definition
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 function xbeam_QuadSkew (w)

! I/O Variables.
  real(8),  intent(in):: w(:)           ! Angular velocity
  real(8)             :: xbeam_QuadSkew(4,4)

  xbeam_QuadSkew = reshape((/0.d0,w(1),w(2),w(3),-w(1),0.d0,-w(3),w(2),-w(2),w(3),0.d0,-w(1),-w(3),-w(2),w(1),0.d0/),   &
&                          (/4,4/), ORDER = (/ 2, 1 /))
  return
 end function xbeam_QuadSkew

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!->Subroutines XBEAM_2NORM
!
!->Description.-
!
!   computes the 2norm of a vector
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 function xbeam_2norm (V)
! I/O Variables.
  real(8),  intent(in):: V(:)           ! Vector
  real(8)             :: xbeam_2norm

! Local variables.
  integer            :: i

  xbeam_2norm = 0.d0
  do i=1,size(V)
      xbeam_2norm = xbeam_2norm + V(i)*V(i)
  end do
  xbeam_2norm = sqrt(xbeam_2norm)

  return
 end function xbeam_2norm


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end module lib_xbeam
