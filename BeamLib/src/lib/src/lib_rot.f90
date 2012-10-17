!->Module LIB_ROT. Rafa Palacios. 07May2008
!
!->Description.-
!
!  Routines with basic operations for finite rotations.
!
!->Subroutines:
!
!   rot_beta2mat:    Convert cosines with b1 to Rotation Matrix.
!   rot_phi2mat:     Convert rotation angles to rotation matrix.
!   rot_points2mat:  Convert two points to Rotation Matrix.
!   rot_cross:       Cross product of two vectors.
!   rot_skew:        Convert a vector to its dual skew-symmetric matrix.
!   rot_outprod:     Outer product of two vectors.
!
!-> Remarks.
!
!  1) Parameter rot_epsilon defines linearization of the finite-rotation
!     equations throughout rot and rotvect libraries. Its value affects the
!     computational costs in the code, but can also affect accuracy if it is
!     set too large.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module lib_rot
 implicit none

! Public variables.
  real(8),public,parameter:: rot_Epsilon   =1.d-3  ! Rotations below this are
                                                   ! linearized
  real(8),public,parameter:: rot_EpsilonMin=1.d-9  ! Limit used for Rodrigues'
                                                   ! singularity at 180 degrees.
! Private Variables.
  real(8), private, parameter,dimension(3,3) :: Unit= &
&     reshape((/1.d0,0.d0,0.d0,0.d0,1.d0,0.d0,0.d0,0.d0,1.d0/),(/3,3/))

 contains


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!-> Subrutine ROT_BETA2MAT
! 
!-> Description.- 
! 
!  Find the transformation matrix from the director cosines for obliqueness.
!
!-> Remarks:
!
!  1) Definition of obliqueness follows the formulation of Popescu, Hodges
!     and Cesnik (2000): Beta is the vector of director cosines between the
!     the oblique axes and the beam longitudinal axis (normal axis to the
!     normal cross section, b1). They satisfy:
!
!               b1= Beta1*c1+Beta2*c2+Beta3*c3
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine rot_beta2mat (Beta,RotMatrix)
!
!-> I/O Variables.
!
  real(8),intent(in) :: Beta(:)        ! Obliqueness Angles.
  real(8),intent(out):: RotMatrix(:,:) ! Rot matrix from b to c.
!
!-> Local variables.
!
  real(8):: Aux      ! Auxiliar Variable.
!
  Aux= 1.d0/(1.d0+Beta(1))
!
  RotMatrix(1,1)= Beta(1)
  RotMatrix(1,2)=-Beta(2)
  RotMatrix(1,3)=-Beta(3)
!
  RotMatrix(2,1)= Beta(2)
  RotMatrix(2,2)= 1.d0 - Beta(2)*Beta(2)*Aux
  RotMatrix(2,3)=-Beta(2)*Beta(3)*Aux
!
  RotMatrix(3,1)= Beta(3)
  RotMatrix(3,2)= RotMatrix(2,3)  
  RotMatrix(3,3)= 1.d0 - Beta(3)*Beta(3)*Aux
!
  return
 end subroutine rot_beta2mat


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!-> Function ROT_PHI2MAT.
!
!-> Description.-
!
!  Convert simple rotations along the coordinate axis to rotation matrices.
!
!-> Remarks.
!
!  1) The elemental rotations are applied in the following order: Z,Y,X.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 function rot_phi2mat (Phi)
!
! I/O Variables.
!
  real(8),intent(in):: Phi(3)     ! Elemental rotations.
  real(8)::   rot_phi2mat(3,3)   ! Rotation matrix.
!
! Local variables
!
  real(8):: RotX(3,3)      ! Simple rotations.
  real(8):: RotY(3,3)
  real(8):: RotZ(3,3)

  RotX(1,1)= 1.d0
  RotX(2,2)= dcos(Phi(1))
  RotX(3,3)= dcos(Phi(1))
  RotX(2,3)=-dsin(Phi(1))
  RotX(3,2)= dsin(Phi(1))
!
  RotY(2,2)= 1.d0
  RotY(3,3)= dcos(Phi(2))
  RotY(1,1)= dcos(Phi(2))
  RotY(3,1)=-dsin(Phi(2))
  RotY(1,3)= dsin(Phi(2))
!
  RotZ(3,3)= 1.d0
  RotZ(1,1)= dcos(Phi(3))
  RotZ(2,2)= dcos(Phi(3))
  RotZ(1,2)=-dsin(Phi(3))
  RotZ(2,1)= dsin(Phi(3))
!
  rot_phi2mat=matmul(matmul(RotX,RotY),RotZ)
!
  return
 end function rot_phi2mat


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!->Subroutine ROT_POINTS2MAT.
!
!->Description.-
!
!   Convert the location of two points in frame A to a Rotation Matrix from A to B.
!
!   Point 1 is in the Z axis of B.
!   Point 2 is in the X-Z plane of B.
!
!->Remarks.-
!
!   1) Both coordinate reference systems have the same origin, 0.
!
!   2) The rotation matrix from A to B is defined as:
!
!          Bi=(C_BA)ij�Aj
!
!          Bi and Aj are the base vectors of frames B and A, respectively.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine rot_points2mat (Point1, Point2, RotMatrix)

! I/O Variables.
  real(8),intent(in) :: Point1(3)       ! Point 1 coordinates.
  real(8),intent(in) :: Point2(3)       ! Point 2 coordinates.
  real(8),intent(out):: RotMatrix(3,3)  ! Rotation Matrix.

! Local Variables.
  real(8) :: b1(3)    ! Unit vector of B along x1.
  real(8) :: b2(3)    ! Unit vector of B along x2.
  real(8) :: b3(3)    ! Unit vector of B along x3.
  real(8) :: Unit1(3) ! Unit vector in the direction of point 1.
  real(8) :: Unit2(3) ! Unit vector in the direction of point 2.
!
! Get unit vectors in the direction of point 1 and point2.
!
  Unit1=Point1/dsqrt(dot_product(Point1,Point1))
  Unit2=Point2/dsqrt(dot_product(Point2,Point2))
!
! Vector 0-1 gives b3.
!
  b3= Unit1
!
! b2 is normal to 0-1 and 0-2.
!
  b2= matmul(rot_skew(Unit1),Unit2)
  b2= b2/dot_product(b2,b2)
!
! b1 is normal to b2 and b3.
!
  b1= matmul(rot_skew(b2),b3)
!
! The Rotation matrix from A to B is defined by the components of this
! orthonormal base.
!
  RotMatrix(1,:)=b1
  RotMatrix(2,:)=b2
  RotMatrix(3,:)=b3
!  
  return
 end subroutine rot_points2mat


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!->Subroutine ROT_CROSS.
!
!->Description.-
!
!   Cross Product of two vectors.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 function rot_cross (Vec1, Vec2)
!
! I/O Variables.
!
  real(8),intent(in)  :: Vec1 (3)
  real(8),intent(in)  :: Vec2 (3)
  real(8),dimension(3):: rot_cross
!
! Local variables.
!
  rot_cross(1)= Vec1(2)*Vec2(3) - Vec1(3)*Vec2(2)
  rot_cross(2)= Vec1(3)*Vec2(1) - Vec1(1)*Vec2(3)
  rot_cross(3)= Vec1(1)*Vec2(2) - Vec1(2)*Vec2(1)
!
  return
 end function rot_cross


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!->Subroutine ROT_SKEW.
!
!->Description.-
!
!   Compute the Skew-Symmetric Matrix given a vector.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 function rot_skew (Vector)

! I/O Variables.

  real(8), intent(in)    :: Vector(3)
  real(8), dimension(3,3):: rot_skew

  rot_skew=0.0d0
  rot_skew(1,2)=-Vector(3)
  rot_skew(1,3)= Vector(2)
  rot_skew(2,1)= Vector(3)
  rot_skew(2,3)=-Vector(1)
  rot_skew(3,1)=-Vector(2)
  rot_skew(3,2)= Vector(1)

  return
 end function rot_skew
 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!->Subroutine ROT_VECT.
!
!->Description.-
!
!   Extract the vector of a skew-symmetric matrix.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 function rot_vect (SkewMatrix)

! I/O Variables.
  real(8), intent(in)  :: SkewMatrix(3,3)
  real(8), dimension(3):: rot_vect

  rot_vect(1)= SkewMatrix(3,2)
  rot_vect(2)= SkewMatrix(1,3)
  rot_vect(3)= SkewMatrix(2,1)

  return
 end function rot_vect


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!->Subroutine ROT_OUTPROD.
!
!->Description.-
!
!   Compute the outer product of two vectors: a*transpose(b)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 function rot_outprod (a,b)

! I/O Variables.
  real(8),intent(in)    :: a(3)
  real(8),intent(in)    :: b(3)
  real(8),dimension(3,3):: rot_outprod

! Local variables.
  integer :: i,j 
  
  do i=1,3
    do j=1,3
      rot_outprod(i,j)= a(i)*b(j)
    end do
  end do

  return
 end function rot_outprod
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end module lib_rot
