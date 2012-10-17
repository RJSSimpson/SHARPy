!->Module LIB_BGEOM. Rafa Palacios. Last Update 6Aug2008.
!
!->Description.-
!
!    Subroutines for computing geometric magnitudes in beam elements.
!
!->Subroutines.-
!
!    -bgeom_elemlength  : Compute length of the element.
!    -bgeom_elemcurv2   : Compute curvature vector of a 2-noded element.
!    -bgeom_elemframe   : Compute element frame orientation.
!    -bgeom_nodeframe   : Compute frame orientation for nodes in an element.
!    -bgeom_strain2disp : Integrate strains along reference line.
!
!->Remarks.-
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module lib_bgeom
 implicit none

! (This module does not contain Public Variables).

!-> Private variables.
  real(8),private,parameter,dimension(3,3):: Unit= &         ! Unit matrix.
&         reshape((/1.d0,0.d0,0.d0,0.d0,1.d0,0.d0,0.d0,0.d0,1.d0/),(/3,3/))

 contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-> Subroutine BGEOM_ELEMLENGTH
!
!-> Description:
!
!    Compute length of an element using Gaussian integration.
!
!-> Remarks.-
!
!   1) The element can have an arbitrary number of nodes. Isoparametric 
!      interpolation is used.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine bgeom_elemlength (NumGauss,NumNodesElem,Coords,Arclength)
  use lib_fem

! I/O Variables.
  integer,intent(in) :: NumGauss          ! Number of Gauss integration points.
  integer,intent(in) :: NumNodesElem      ! Number of nodes in the element.
  real(8),intent(in) :: Coords    (:,:)   ! Coordinates of the grid points.
  real(8),intent(out):: ArcLength         ! Element total arclength.

! Local variables.
  integer :: i                     ! Counters.
  integer :: iGauss                ! Counter on the Gaussian points.
 real(8) :: dXgauss(3)            ! Coordinate at the Gauss point.
  real(8),allocatable :: ShapeFun(:)    ! Shape functions for a gauss CoordGauss.
  real(8),allocatable :: ShapeDer(:)    ! Derivatives of ShapeFun in a gauss CoordGauss.
  real(8),allocatable :: CoordGauss (:) ! Coords of Gauss points.
  real(8),allocatable :: WeightGauss(:) ! Coords of Gauss points.

! Initialize.
  Arclength=0.d0
  allocate (CoordGauss (NumGauss))
  allocate (WeightGauss(NumGauss))
  allocate (ShapeFun   (NumNodesElem))
  allocate (ShapeDer   (NumNodesElem))

! Define Gauss points and loop on them.
  call fem_1d_gauss_val (NumGauss,CoordGauss,WeightGauss)

  do iGauss=1,NumGauss

! Obtain the shape functions and their derivatives.
    call fem_1d_shapefun (NumNodesElem,CoordGauss(iGauss),ShapeFun,ShapeDer)

! Global coordinates and their derivatives at the gauss points.
    do i=1,3
      dXgauss(i)= dot_product (ShapeDer(1:NumNodesElem),Coords(1:NumNodesElem,i))
    end do

! Jacobian of the arclength function at the Gauss point.
    ArcLength=ArcLength+WeightGauss(iGauss)*(sqrt(dot_product(dXgauss,dXgauss)))

  end do

  deallocate (CoordGauss,WeightGauss)
  deallocate (ShapeFun,ShapeDer)
  return
 end subroutine bgeom_elemlength



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-> Subroutine BGEOM_ELEMCURV2
!
!-> Description:
!
!    Compute curvature at the center of a 2-noded beam element for given 
!    cartesian rotation vector at its nodes.
!
!-> Remarks.-
!
!  1) This is vector kb in the undeformed frame, b.
!
!  2) Element goes from node A to node B.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine bgeom_elemcurv2 (PosGrids,PsiGrids,kb)
  use lib_fem
  use lib_rotvect

! I/O Variables.
  real(8),intent(in) :: PosGrids (:,:) ! Position vector at nodes.
  real(8),intent(in) :: PsiGrids (:,:) ! CRV at nodes
  real(8),intent(out):: kb     (3)     ! Curvature of the element.

! Local variables.
  real(8) :: dx(3)                     ! Derivatives at each point.
  real(8) :: Psi0(3),dPsi0(3)          ! Interpolated rotation vector and its derivative.
  real(8) :: Jacobian                  ! ds/deta.
  integer :: NumNodesElem              ! Number of nodes in the element.
  real(8),allocatable :: ShapeFun(:)   ! Shape functions for a gauss CoordGauss.
  real(8),allocatable :: ShapeDer(:)   ! Derivatives of ShapeFun in a gauss CoordGauss.

! Initialize
  NumNodesElem=2
  allocate (ShapeFun(NumNodesElem)); ShapeFun=0.d0
  allocate (ShapeDer(NumNodesElem)); ShapeDer=0.d0

! Obtain the shape functions derivatives at the element center point.
  call fem_1d_shapefun (NumNodesElem,0.d0,ShapeFun,ShapeDer)

! Compute Jacobian and the differential operator in physical coordinates.
  dx= matmul(ShapeDer,PosGrids(1:NumNodesElem,:))
  Jacobian=dsqrt(dot_product(dx,dx))
  ShapeDer=ShapeDer/Jacobian

! Compute value of the rotation vector and its derivative at the center point.
  Psi0 = matmul(ShapeFun,PsiGrids(1:NumNodesElem,:))
  dPsi0= matmul(ShapeDer,PsiGrids(1:NumNodesElem,:))
  
  kb= matmul(rotvect_psi2rot (Psi0),dPsi0)

  deallocate (ShapeFun,ShapeDer)
  return
 end subroutine bgeom_elemcurv2



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-> Subroutine BGEOM_ELEMFRAME
!
!-> Description:
!
!    Compute element frame as defined by MSC conventions.
!
!-> Remarks.-
!
!  1) Element frame is the 'a' frame in the co-rotational approach.
!
!  2) X axis is defined by the two end nodes of the elements.  The X-Z plane
!     is orthogonal to the given V vector.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine bgeom_elemframe (r,V,Psi,Delta)
  use lib_rot
  use lib_rotvect

! I/O Variables.
  real(8),intent(in) :: r(2,3)         ! Position vector of initial and final node.
  real(8),intent(in) :: V(3)           ! Cross-sectional orientation vector.
  real(8),intent(out):: Psi(3)         ! CRV of the element orientation.
  real(8),intent(in) :: Delta          ! Minimum angle for two unit vectors to be parallel.

! Local variables.
  real(8)::  Cag (3,3)                 ! Coordinate transformation matrix from global to element frame.
  real(8)::  d(3)                      ! Direction of base vector.

! X-axis.
  d=r(2,:)-r(1,:)
  Cag(1,:)= d/sqrt(dot_product(d,d))

! Z-axis.
  d=rot_cross(Cag(1,:),V)

  if (sqrt(dot_product(d,d)).le.Delta) &
&    STOP 'Element Vector V is not orthogonal to reference line (51105)'

  Cag(3,:)= d/sqrt(dot_product(d,d))

! Y-axis.
  Cag(2,:)= rot_cross(Cag(3,:),Cag(1,:))

! Convert to rotation vector with respect to the reference (global) frame.
  Psi= rotvect_mat2psi(Cag)

  return
 end subroutine bgeom_elemframe


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-> Subroutine BGEOM_NODEFRAME
!
!-> Description:
!
!    For a given node, compute undeformed coordinate system ('b' frame) at 
!    each of the element nodes.
!
!-> Remarks.-
!
!  1) V is only used for collinear points in an element. This is defined
!     by the convergence parameter (Delta).
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine bgeom_nodeframe (NumNodesElem,Coords,Phi,V,Psi,Delta)
  use lib_fem
  use lib_rot
  use lib_rotvect

! I/O Variables.
  integer,intent(in) :: NumNodesElem   ! Number of nodes in the element.
  real(8),intent(in) :: Coords  (:,:)  ! Coordinates of the grid points.
  real(8),intent(in) :: Phi     (:)    ! Pretwist at grid points.
  real(8),intent(in) :: V       (:)    ! Element cross-sectional orientation.
  real(8),intent(out):: Psi     (:,:)  ! CRV at each node.
  real(8),intent(in) :: Delta          ! Minimum angle for two unit vectors to be parallel.

! Local variables.
  real(8) :: AuxV(3)                   ! Auxiliary vector.
  real(8) :: Cgb(3,3)                  ! Coordinate transformation matrix from global to local.
  integer :: iNode                     ! Counter on the nodes in the element.
  integer :: j                         ! Counters.
  real(8) :: dx(3)                     ! Derivatives of the position vector at each point.
  real(8) :: e_t(3),e_n(3),e_b(3)      ! Tangent, normal, and binormal vectors.
  real(8) :: Phi0                      ! Defines direction of zero pretwist.
  real(8) :: v12(3),v13(3),v23(3)      ! Auxiliary vectors.

  real(8),allocatable:: Eta     (:)    ! Nodal position in the 1-D nondimensional coordinate.
  real(8),allocatable:: ShapeFun(:)    ! Shape functions for a gauss CoordGauss.
  real(8),allocatable:: ShapeDer(:)    ! Derivatives of ShapeFun in a gauss CoordGauss.


! Initialize
  allocate (Eta     (NumNodesElem)); Eta     = 0.d0
  allocate (ShapeFun(NumNodesElem)); ShapeFun= 0.d0
  allocate (ShapeDer(NumNodesElem)); ShapeFun= 0.d0

  v12=Coords(2,:)-Coords(1,:)
  v12=v12/sqrt(dot_product(v12,v12))

  select case (NumNodesElem)
  case (3)
! If there are three nodes, compute the vector in the normal plane.
    Eta=(/-1.d0,1.d0,0.d0/)
    v13=Coords(3,:)-Coords(1,:)
    v13=v13/sqrt(dot_product(v13,v13))
    v23=rot_cross(v13,v12)

  case (2)
    Eta=(/-1.d0,1.d0/)
    v23=0.d0

  case default
    STOP 'Wrong number of nodes in the element (51248)'
  end select

! Loop through the element nodes.

  do iNode=1,NumNodesElem

! Obtain the shape functions derivatives.
    call fem_1d_shapefun (NumNodesElem,Eta(iNode),ShapeFun,ShapeDer)

! Compute Jacobian.
    do j=1,3
      dx(j)= dot_product (ShapeDer,Coords(1:NumNodesElem,j))
    end do

! Tangent vector.
    e_t=dx/sqrt(dot_product(dx,dx))

! Binormal vector for collinear and non-collinear cases.
    if (sqrt(dot_product(v23,v23)).le.Delta) then
      AuxV=rot_cross(v12,V)
      e_b=AuxV/dsqrt(dot_product(AuxV,AuxV))
    else
      e_b=v23/dsqrt(dot_product(v23,v23))
    end if

! Normal vector (normalization included to avoid numerical problems).
    e_n=rot_cross(e_b,e_t)
    e_n=e_n/dsqrt(dot_product(e_n,e_n))

! Grid A is used to establish the reference for the pretwist angle, which is defined
! as the angle that vector e_n needs to be rotated around e_t so that it lies in the
! element X-Y plane.
    if (iNode.eq.1) then
      AuxV=rot_cross(rot_cross(v12,V),e_t)
      AuxV=AuxV/sqrt(dot_product(AuxV,AuxV))
      if (sqrt(dot_product(AuxV,AuxV)).gt.Delta) then
        Phi0= sign(acos(max(-1.d0,min(dot_product(e_n,AuxV),1.d0))), &
&                  dot_product(e_t,rot_cross(v12,V)))
      else
        AuxV=rot_cross(v12,rot_cross(v12,V))
        AuxV=AuxV/sqrt(dot_product(AuxV,AuxV))
        Phi0= sign(acos(max(-1.d0,min(dot_product(e_n,AuxV),1.d0))), &
&                  dot_product(e_n,V))
      end if
    end if

! Rotate vectors with the pretwist angle.    
    Cgb(1,:)= e_t
    Cgb(2,:)= e_n*cos(Phi0+Phi(iNode))+e_b*sin(Phi0+Phi(iNode))
    Cgb(3,:)=-e_n*sin(Phi0+Phi(iNode))+e_b*cos(Phi0+Phi(iNode))

! Compute the cartesian rotation vector that defines 
    Psi(iNode,:)= rotvect_mat2psi(Cgb)
  end do

  return
 end subroutine bgeom_nodeframe


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-> Subroutine BGEOM_STRAIN2DISP
!
!-> Description:
!
!    Integrate strain/curvature along a continuous reference line and
!    compute the nodal values of the position and rotation vectors.
!
!-> Remarks.-
!
!  1) The strain is assumed constant within each element.
!
!  2) This routine is only valid for 2-noded elements along a single line.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine bgeom_strain2disp (NumElems,PsiIni,Lengths,Gamma,Curvature,PosFin,PsiFin)                                 
  use lib_fem
  use lib_rotvect

! I/O Variables.
  integer,intent(in) :: NumElems         ! Number of elements in the model.
  real(8),intent(in) :: PsiIni  (:,:,:)  ! Undeformed orientation vector (CRV).
  real(8),intent(in) :: Lengths (:)      ! Element lengths.
  real(8),intent(in) :: Gamma   (:,:)    ! Extensional strain within each element.
  real(8),intent(in) :: Curvature (:,:)  ! Total curvature in the element (KB).
  real(8),intent(out):: PosFin  (:,:)    ! Position vector at end nodes.
  real(8),intent(out):: PsiFin  (:,:,:)  ! CRV at the nodes.

! Local variables.
  integer:: iElem              ! Counter on the finite elements.
  real(8):: C_aB(3,3)          ! Coordinate transformation matrix at current node.
  real(8):: C_aB1(3,3)         ! Trafo matrix at end of previous element.
  real(8):: C_aB2(3,3)         ! Trafo matrix at beginning of current element.
  real(8):: DeltaPsi(3)        ! KB*DeltaS.
  real(8):: H0(3,3),H1(3,3)    ! Element rotation operators.
  real(8):: Ra(3)              ! Position vector at current node.


! Initialize the kinematic chain.
  C_aB= Unit
  Ra = 0.d0

! Loop in all elements in the model.
  do iElem=1,NumElems

! Compute orientation change at the first node of the element
! from the previous element.
    if (iElem.gt.1) then
      C_aB1= rotvect_psi2mat(-PsiIni(iElem-1,2,:))
    else
      C_aB1= rotvect_psi2mat(-PsiIni(iElem,1,:))
    end if
    C_aB2= rotvect_psi2mat(-PsiIni(iElem,1,:))
    C_aB=matmul(C_aB,matmul(transpose(C_aB1),C_aB2))
    PsiFin(iElem,1,:)= rotvect_mat2psi(transpose(C_aB))

! Integrate displacements.
    DeltaPsi= Lengths(iElem) * Curvature(iElem,:)
    H1= Lengths(iElem)*rotvect_psi2intmat(-DeltaPsi)
    Ra= Ra + matmul(C_aB,matmul(H1,Unit(1:3,1)+Gamma(iElem,:)))

! Integrate rotation.
    H0  = rotvect_psi2mat(-DeltaPsi)
    C_aB= matmul(C_aB,H0)

! Store values.
    PosFin(iElem,:)  = Ra
    PsiFin(iElem,2,:)= rotvect_mat2psi(transpose(C_aB))
  end do

  return
 end subroutine bgeom_strain2disp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end module lib_bgeom
