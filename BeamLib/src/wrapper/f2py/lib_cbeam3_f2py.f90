!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-> Subroutine CBEAM3_FSTIFZ
!
!-> Description:
!
!    Compute stiffness forces at local element ordinate z.
!
!-> Remarks.-
!
!   1) For three noded elements r0(i,:) corresponds to LHS elemnt for i=1, RHS 
!      for i=2, and center for i=3.
!   2) z = 0 is the center of the element with z in [-1,1].
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine cbeam3_fstifz (NumNodesElem,r0,Ri,ElemStiff,Forces,z)
  use lib_fem
  use lib_rot
  use lib_rotvect

! I/O Variables.
  integer,intent(in)   :: NumNodesElem       ! Number of nodes in the element.
  real(8),intent(in)   :: r0       (NumNodesElem,6)! Initial position/orient
  real(8),intent(in)   :: Ri       (NumNodesElem,6)! Current position/orient
  real(8),intent(in)   :: ElemStiff(6,6)     ! Stiffness properties in the element.
  real(8),intent(out):: Forces   (6)       ! Stiffness forces.
  real(8),intent(in)   :: z                  ! Local element coordinate

! Local variables.
  integer :: i                     ! Counter.
  real(8) :: Jacobian              ! ds/deta, with eta the nondimensional arclength.
  integer, parameter :: MaxNodCB3 = 3
  real(8) :: ShapeFun(MaxNodCB3)    ! Shape functions in a gauss point.
  real(8) :: ShapeDer(MaxNodCB3)    ! Derivatives of ShapeFun in a gauss point.

  real(8) :: Ri_z (6)              ! Displacement/rotation at Gauss point.
  real(8) :: r0_z (6)              ! Position vector at the Gauss point.
  real(8) :: dRi_z(6)              ! Derivative of Ri_z.
  real(8) :: dr0_z(6)              ! Derivative of r0_z.

  real(8) :: CBa(3,3)              ! Coordinate transformation matrix.
  real(8) :: KB(3)                 ! Curvature in the deformed configuration at Gauss point.
  real(8) :: gamma(3)              ! Strain at Gauss point.
  real(8) :: Rot(3,3),dRot(3,3)    ! Tangential operator and spatial derivative.
  real(8) :: Strain(6)             ! Beam strains at the Gauss point.

! Computation.

! Obtain the shape functions and their derivatives.
    call fem_1d_shapefun (NumNodesElem,z,ShapeFun,ShapeDer)

! Global coordinates and their derivatives at the gauss points.
    do i=1,6
      r0_z(i) = dot_product (ShapeFun(1:NumNodesElem),r0(1:NumNodesElem,i))
      dr0_z(i)= dot_product (ShapeDer(1:NumNodesElem),r0(1:NumNodesElem,i))
    end do
    Jacobian =sqrt(dot_product (dr0_z(1:3),dr0_z(1:3)))

! Rescale the derivatives to be in the physical coordinates.
    ShapeDer=ShapeDer/Jacobian
    dr0_z=dr0_z/Jacobian

! Compute the current position vector and rotations, and their derivatives.
    do i=1,6
      Ri_z(i) = dot_product (ShapeFun(1:NumNodesElem),Ri(1:NumNodesElem,i))
      dRi_z(i)= dot_product (ShapeDer(1:NumNodesElem),Ri(1:NumNodesElem,i))
    end do

! Compute the current coordinate transformation matrix, the rotational operator, and
! its spatial derivative at the Gauss point.
    CBa= rotvect_psi2mat  (Ri_z(4:6))
    Rot= rotvect_psi2rot  (Ri_z(4:6))
    dRot=rotvect_drotdpsib(Ri_z(4:6),dRi_z(4:6))

! Compute the current curvature and strain at the Gaussian point.
    gamma=matmul(CBa,dRi_z(1:3)) - matmul(rotvect_psi2mat(r0_z(4:6)),dr0_z(1:3))
    KB=matmul(Rot,dRi_z(4:6))

! Compute force/moment strains.
    Strain(1:3)= gamma
    Strain(4:6)= KB - matmul(rotvect_psi2rot(r0_z(4:6)),dr0_z(4:6))
    Forces=matmul(ElemStiff,Strain)

  return
 end subroutine cbeam3_fstifz
