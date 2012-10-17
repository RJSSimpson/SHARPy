!->Module LIB_ROTVECT. Rafa Palacios. 07May2008-25May2011
!
!->Description.-
!
!  This module defines functions to operate with Cartesian Rotation
!  Vector (CRV) on finite rotations.
!
!-> List of functions:
!
!   rotvect_bounds:     Redefine a CRV for its norm to be within [-pi,pi].
!   rotvect_addpsi:     Compose rotations given by CRV.
!   rotvect_addquat:    Compose rotations given by quaternions.
!   rotvect_deltapsi:   Substract two Cartesian rotation vectors.
!   rotvect_lintpsi:    Interpolate CRV through linint of quaternions.
!   rotvect_lintquat:   Linear interpolation of quaternions.
!
!   rotvect_psi2theta:  Transform CRV to Rodrigues parameters.
!   rotvect_theta2psi:  Transform Rodrigues parameters to CRV.
!   rotvect_psi2quat:   Transform CRV to quaternions.
!   rotvect_quat2psi:   Transform quaternions to CRV.
!   rotvect_mat2quat:   Transform coordinate transf matrix to quaternions.
!   rotvect_psi2mat:    Transform CRV to the coordinate transf matrix.
!   rotvect_mat2psi:    Transform coordinate transf matrix to CRV.
!   rotvect_psi2rot:    Compute rotational operator for a CRV.
!   rotvect_psi2intmat: Transform CRV to the integral of the rotation matrix.
!
!   rotvect_drotdpsib:  Derivative of the rotational operator with the CRV.
!   rotvect_ddrotdpsib: Derivative of the time/space derivative of rot with the CRV.
!   rotvect_a1:         A1 operator (delta(T)*b=a1(Psi,b)*Delta(Psi))
!   rotvect_a2:         A2 operator (A2(Psi,b)=-A1(-Psi,b))
!   rotvect_b1:         B1 operator.
!   rotvect_b2:         B2 operator.
!
!-> Remarks.
!
!   1) The rotation vector is defined between -180 and +180 degrees. For
!      angles beyond this limits, it is mapped into this range. Note that
!      the discontinuity needs to be considered in later operations with
!      rotations.
!
!-> History.
!   25.05.2011: Review to set bounds on the CRV. 
!               Added bounds deltapsi lintpsi lintquat
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module lib_rotvect
 use lib_rot
 implicit none
!
!   (There are no public variables in this module).
!
!-> Private Variables.
!
  real(8),private,parameter,dimension(3,3) :: Unit= &
&       reshape((/1.d0,0.d0,0.d0,0.d0,1.d0,0.d0,0.d0,0.d0,1.d0/),(/3,3/))

 contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!->Function ROTVECT_BOUNDS.
!
!->Description.-
!
!   Obtain the norm of the the Cartesian Rotation Vector and force the vector
!   to be in the [-pi,pi] interval.
!
!->Remarks.-
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine rotvect_bounds (Psi)

! I/O Variables.
  real(8),intent(inout)   :: Psi(3)

! Local variables.
  real(8):: OldNormPsi   ! Original norm of the CRV.
  real(8):: NormPsi      ! New norm within bounds.
  real(8):: pi

! Initialize
  pi=4.d0*datan(1.d0)

! Compute original norm.
  OldNormPsi=dsqrt(dot_product(Psi,Psi))

! Force the norm to be in [-pi,pi].
  NormPsi=OldNormPsi-2.d0*pi*int(OldNormPsi/(2.d0*pi))
  if (NormPsi.eq.0.d0) then
    Psi=0.d0
  else
    if (NormPsi.gt.pi) then
      NormPsi=NormPsi - 2.d0*pi
    elseif (NormPsi.le.-pi) then
      NormPsi=NormPsi + 2.d0*pi
    end if
    Psi=Psi*(NormPsi/OldNormPsi)
  end if
  
  return
 end subroutine rotvect_bounds


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!->Function ROTVECT_BOUNDSCHECK2.
!
!->Description.-
!
!   Check if two rotations are on both sides of the -pi / pi discontinuity and
!   add 2pi to one of them.
!
!->Remarks.-
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine rotvect_boundscheck2 (PsiA,PsiB)

! I/O Variables.
  real(8),intent(inout):: PsiA (3)
  real(8),intent(in)   :: PsiB (3)

! Local variables.
  real(8)              :: Norm
  real(8)              :: pi

! Initialize
  pi=4.d0*datan(1.d0)
  Norm=dsqrt(dot_product(PsiA,PsiA))

  if (Norm.ne.0.d0) then
! Check if the closest "distance" from the first to the second vector is obtained
! by adding or substracting 2*pi. This should remove the discontinuity at +/- pi,
! but needs some extensive checking.
    if (dot_product(PsiB-PsiA*(Norm+2.d0*pi)/Norm,PsiB-PsiA*(Norm+2.d0*pi)/Norm).lt.pi*pi) then
      PsiA = PsiA*(Norm+2.d0*pi)/Norm
    else if (dot_product(PsiB-PsiA*(Norm-2.d0*pi)/Norm,PsiB-PsiA*(Norm-2.d0*pi)/Norm).lt.pi*pi) then
      PsiA = PsiA*(Norm-2.d0*pi)/Norm
    end if
  end if

  return
 end subroutine rotvect_boundscheck2


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!->Function ROTVECT_ADDPSI.
!
!->Description.-
!
!   Compose two rotations given by Cartesian Rotation vectors (right rotation).
!
!->Remarks.-
!  1) The second rotation is given with reference on the first one (a right
!     rotation).
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 function rotvect_addpsi (Psi01, Psi12)

! I/O Variables.
  real(8),intent(in)   :: Psi01 (3)
  real(8),intent(in)   :: Psi12 (3)   
  real(8),dimension(3) :: rotvect_addpsi ! Composed rotation (Psi02)

! Local variables.
  real(8),dimension(4) :: Quat01,Quat12,Quat02  ! Quaternions.
  real(8),dimension(3) :: Psi                   ! CRV.
  
! Transform to quaternions, add them, and convert the result back.
  Quat01=rotvect_psi2quat(Psi01)
  Quat12=rotvect_psi2quat(Psi12)

  Quat02=rotvect_addquat(Quat01,Quat12)

  Psi=rotvect_quat2psi(Quat02)
  call rotvect_bounds (Psi)

  rotvect_addpsi=Psi
  
  return
 end function rotvect_addpsi


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!->Function ROTVECT_ADDQUAT.
!
!->Description.-
!
!   Compose two rotations defined by quaternions.
!
!->Remarks.-
!  1) The second rotation is given with reference on the first one, which
!     corresponds to a right rotation.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 function rotvect_addquat (Quat01, Quat12)

! I/O Variables.
  real(8),intent(in)   :: Quat01 (4)
  real(8),intent(in)   :: Quat12 (4)   
  real(8),dimension(4) :: rotvect_addquat ! Composed rotation (Quat02)

! Local variables.
  real(8),dimension(4,4):: B              ! B matrix.
  
! Define the B matrix based on the right rotation (Quat12)
  B(1,1)     = Quat12(1)
  B(1,2:4)   =-Quat12(2:4)
  B(2:4,1)   = Quat12(2:4)
  B(2:4,2:4) = Quat12(1)*Unit - rot_skew(Quat12(2:4))
  
! Compute composite rotation.
  rotvect_addquat=matmul(B,Quat01)

  return
 end function rotvect_addquat


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!->Function ROTVECT_DELTAPSI.
!
!->Description.-
!
!  Substract two Cartesian rotation vectors (PsiB-PsiA).
!
!->Remarks.-
!
!   1) This routine simple checks that the two vectors don't lie on both sides
!      of the +/- pi discontinuity.
!
!   2) It needs further work to ensure "objectivity" (invariance in the evaluation
!      of derivatives of the CRV with rigid-body rotations of the global system.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 function rotvect_deltapsi (PsiA, PsiB)

! I/O Variables.
  real(8),intent(in)   :: PsiA (3)
  real(8),intent(in)   :: PsiB (3)
  real(8),dimension(3) :: rotvect_deltapsi

! Local variables.
  real(8),dimension(3) :: DeltaPsi       ! Increment between both CRVs.
  real(8)              :: Norm
  real(8)              :: pi

! Initialize
  pi=4.d0*datan(1.d0)
  Norm=dsqrt(dot_product(PsiB,PsiB))

! Transform to quaternions. Check if the closest "distance" from the first to the second
! vector is obtained by adding or substracting the 2*pi. In FEM applications this should
! remove the discontinuity at +/- pi, but needs some extensive checking.
  if (dot_product(PsiB*(Norm+2.d0*pi)/Norm-PsiA,PsiB*(Norm+2.d0*pi)/Norm-PsiA).lt.Pi*Pi) then
    DeltaPsi = PsiB*(Norm+2.d0*pi)/Norm - PsiA
  else if (dot_product(PsiB*(Norm-2.d0*pi)/Norm-PsiA,PsiB*(Norm-2.d0*pi)/Norm-PsiA).lt.Pi*pi) then
    DeltaPsi = PsiB*(Norm-2.d0*pi)/Norm - PsiA
  else
    DeltaPsi = PsiB - PsiA
  end if

  rotvect_deltapsi=DeltaPsi
  return
 end function rotvect_deltapsi



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!->Function ROTVECT_LINTPSI.
!
!->Description.-
!
!   Interpolation of Cartesian Rotation vectors using linear interpolation of
!   the quaternions.
!
!->Remarks.-
!
!   1) The approach needs to be reviewed for cases where the two quaternions
!      go over the discontinuity at +/- pi. (RPN 25.05.2011)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 function rotvect_lintpsi (PsiA, PsiB, eta)

! I/O Variables.
  real(8),intent(in)   :: PsiA (3)
  real(8),intent(in)   :: PsiB (3)
  real(8),intent(in)   :: eta                ! Eta between 0 and 1.
  real(8),dimension(3) :: rotvect_lintpsi

! Local variables.
  real(8),dimension(4) :: QuatA,QuatB,QuatC     ! Quaternions.
  real(8),dimension(3) :: PsiC                  ! CRV at the interpolation point.
  real(8)              :: Norm
  real(8)              :: pi

! Initialize
  pi=4.d0*datan(1.d0)
  Norm=dsqrt(dot_product(PsiB,PsiB))

! Transform to quaternions. Check if the closest "distance" from the first to the second
! vector is obtained by adding or substracting the 2*pi. In FEM applications this should
! remove the discontinuity at +/- pi, but needs some extensive checking.
  QuatA=rotvect_psi2quat(PsiA)
  if (dot_product(PsiB*(Norm+2.d0*pi)/Norm-PsiA,PsiB*(Norm+2.d0*pi)/Norm-PsiA).lt.Pi*Pi) then
    QuatB=rotvect_psi2quat(PsiB*(Norm+2.d0*pi)/Norm)
  else if (dot_product(PsiB*(Norm-2.d0*pi)/Norm-PsiA,PsiB*(Norm-2.d0*pi)/Norm-PsiA).lt.Pi*pi) then
    QuatB=rotvect_psi2quat(PsiB*(Norm-2.d0*pi)/Norm)
  else
    QuatB=rotvect_psi2quat(PsiB)
  end if

  QuatC=rotvect_lintquat(QuatA,QuatB,eta)

! Transform back to CRV and force results within [-pi,+pi]
  PsiC=rotvect_quat2psi(QuatC)
  call rotvect_bounds (PsiC)

  rotvect_lintpsi=PsiC

  return
 end function rotvect_lintpsi


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!->Function ROTVECT_LINTQUAT.
!
!->Description.-
!
!   Linear interpolation of quaternions using the SLERP algorithm.
!
!->Remarks.-
!
!  1) See A.J. Hanson, "Visualizing quaternions", page 307.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 function rotvect_lintquat (QuatA, QuatB, eta)

! I/O Variables.
  real(8),intent(in)   :: QuatA (4)
  real(8),intent(in)   :: QuatB (4)
  real(8),intent(in)   :: eta              ! Interpolation param, between 0 and 1.
  real(8),dimension(4) :: rotvect_lintquat ! Interpolated rotation.

! Local variables.
  real(8):: QuatC(4)
  real(8):: NormQuat, Phi

! Minimum angle between A and B.
  Phi=dacos(dot_product(QuatA,QuatB))

! Compute interpolated rotation.
  if (abs(Phi).lt.rot_Epsilon) then
    QuatC=(1.d0-eta)*QuatA + eta*QuatB
  else
    QuatC=(dsin((1.d0-eta)*Phi)/dsin(Phi))*QuatA + &
&         (dsin(      eta *Phi)/dsin(Phi))*QuatB
  end if

! Enforce unit norm.
  NormQuat=dsqrt(dot_product(QuatC,QuatC))
  rotvect_lintquat =QuatC/NormQuat

  return
 end function rotvect_lintquat



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!->Function ROTVECT_PSI2THETA.
!
!->Description.-
!
!   Transform a rotation vector to the corresponding Rodrigues parameters.
!
!-> Remarks.-
!
!   1) It returns an error when Theta~=pi.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 function rotvect_psi2theta (Psi)

! I/O Variables.
  real(8),intent(in)  :: Psi(3)             ! Rotation Vector.
  real(8),dimension(3):: rotvect_psi2theta  ! Rodrigues Parameters.

! Local variables
  real(8):: NormPsi                         ! Norm of the Rotation Vector.
  real(8):: pi

! Initialize
  pi=4.d0*datan(1.d0)

! Compute normal vector.
  NormPsi=dsqrt(dot_product(Psi,Psi))

! The evaluation of the parameters depends on the magnitude of the rotation 
! vector.
  if (abs(NormPsi).le.rot_Epsilon) then
    rotvect_psi2theta=Psi
  else if (abs(NormPsi-pi).le.rot_EpsilonMin) then
    STOP 'ERROR: Rodrigues Parameters singularity was reached (11511)'
  else
    rotvect_psi2theta= (dtan(NormPsi/2.d0)/(NormPsi/2.d0)) * Psi
  end if

  return
 end function rotvect_psi2theta


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!->Function ROTVECT_THETA2PSI.
!
!->Description.-
!
!   Converts a vector of Rodrigues parameters into the corresponding Cartesian
!   rotation vector.
!
!-> Remarks.-
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 function rotvect_theta2psi (Theta)
 
! I/O Variables.
  real(8),intent(in)  :: Theta(3)       ! Rodrigues Parameters.
  real(8),dimension(3):: rotvect_theta2psi  ! Rotation Vector.

! Local variables
  real(8):: NormTheta      ! Norm of Theta.
  real(8):: Normal(3)      ! Unit normal vector to the rotation plane.
  real(8):: Psi   (3)      ! CRV
  real(8):: NormPsi        ! Norm of the CRV

  NormPsi = 0.d0    ! Rob: stops normpsi uninitialized warning flag.
                    ! Rob: the if statement below makes no sense.

! Compute normal vector.
  NormTheta=dsqrt(dot_product(Theta,Theta))
  if (abs(NormPsi).le.rot_Epsilon) then
    Normal=0.d0
  else
    Normal=Theta/NormPsi
  end if
  
! Compute rotation vector.
  Psi= 2.d0*datan(NormTheta/2.d0)*Normal

! Force bounds on the CRV
  call rotvect_bounds(Psi)
  rotvect_theta2psi=Psi

  return
 end function rotvect_theta2psi


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!->Function ROTVECT_PSI2QUAT.
!
!->Description.-
!
!   Transform a rotation vector to the corresponding quaternions.
!
!-> Remarks.-
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 function rotvect_psi2quat (Psi)

! I/O Variables.
  real(8),intent(in)  :: Psi(3)            ! Rotation Vector.
  real(8),dimension(4):: rotvect_psi2quat  ! Quaternions.

! Local variables
  real(8):: NormPsi        ! Norm of the Rotation Vector.

! Compute the norm.
  NormPsi=dsqrt(dot_product(Psi,Psi))

! Transform rotation vector to quaternions.
  rotvect_psi2quat(1)  =dcos(NormPsi/2.d0)

! Small angles approximation.
  if (abs(NormPsi).lt.rot_Epsilon) then
    rotvect_psi2quat(2:4)=(1.d0/2.d0-NormPsi*NormPsi/12.d0)*Psi

! Regular expression.
  else
    rotvect_psi2quat(2:4)=(dsin(NormPsi/2.d0)/NormPsi)*Psi

  end if

  return
 end function rotvect_psi2quat


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!->Function ROTVECT_QUAT2PSI.
!
!->Description.-
!
!   Transform quaternions to the corresponding rotation vector.
!
!-> Remarks.-
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 function rotvect_quat2psi (Quat)

! I/O Variables.
  real(8),intent(in)  :: Quat(4)            ! Quaternions.
  real(8),dimension(3):: rotvect_quat2psi   ! Rotation Vector.

! Local variables
  real(8):: NormPsi      ! Norm of the Rotation Vector.
  real(8):: Psi(3)       ! CRV


! Compute rotation angle (norm of the rotation vector). It is bounded in [0,2*pi].
  NormPsi = 2.d0 * dacos(max(-1.d0,min(Quat(1),1.d0)))
  
! Compute normal vector.
  if (abs(NormPsi).eq.0.d0) then
    Psi=0.d0
  else
    Psi=NormPsi*Quat(2:4)/dsin(NormPsi/2.d0)
  end if

! Set bounds in [-pi,pi].
  call rotvect_bounds (Psi)

! Evaluate rotation vector.
  rotvect_quat2psi= Psi

  return
 end function rotvect_quat2psi
 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!->Function ROTVECT_MAT2QUAT.
!
!->Description.-
!
!   Transform coordinate transformation matrix into quaternions form.
!
!-> Remarks.-
!
!   1) It uses the algorithm described in (Geradin & Cardona, 2001).
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 function rotvect_mat2quat (CoordTransMat)

! I/O Variables.
  real(8),intent(in)  :: CoordTransMat(3,3)   ! Quaternions.
  real(8),dimension(4):: rotvect_mat2quat     ! Rotation Vector.

! Local variables
  integer:: i             ! Counter.
  integer:: iMax          ! Index of SMax.
  real(8):: Quat(4)       ! Quaternions to be computed.
  real(8):: SMax          ! Max in the diagonal of S.
  real(8):: S(4,4)        ! Intermediate matrix.
  real(8):: R(3,3)        ! Rotation matrix corresponding to CoordTransMat. 

! Transpose to use standard expressions.
  R=transpose(CoordTransMat)

! Intermediate matrix S (see Geradin & Cardona, 2001)
  S(1,1)=1.d0+R(1,1)+R(2,2)+R(3,3)
  S(1,2)=R(3,2)-R(2,3)
  S(1,3)=R(1,3)-R(3,1)
  S(1,4)=R(2,1)-R(1,2)
  
  S(2,1)=R(3,2)-R(2,3)
  S(2,2)=1.d0+R(1,1)-R(2,2)-R(3,3)
  S(2,3)=R(1,2)+R(2,1)
  S(2,4)=R(1,3)+R(3,1)
  
  S(3,1)=R(1,3)-R(3,1)
  S(3,2)=R(2,1)+R(1,2)
  S(3,3)=1.d0-R(1,1)+R(2,2)-R(3,3)
  S(3,4)=R(2,3)+R(3,2)
  
  S(4,1)=R(2,1)-R(1,2)
  S(4,2)=R(1,3)+R(3,1)
  S(4,3)=R(2,3)+R(3,2)
  S(4,4)=1.d0-R(1,1)-R(2,2)+R(3,3)

! Find maximum term in diagonal.
  SMax=0.d0
  do i=1,4
    if (S(i,i).ge.SMax) then
      iMax=i
      SMax=S(i,i)
    end if
  end do
  
! Compute quaternion angles.
  Quat(iMax)=0.5d0*dsqrt(SMax)
  do i=1,4
    if (i.ne.iMax) Quat(i)=0.25d0*S(iMax,i)/Quat(iMax)
  end do

  rotvect_mat2quat=Quat
  return
 end function rotvect_mat2quat
 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!->Function ROTVECT_PSI2MAT.
!
!->Description.-
!
!   Transform a rotation vector to the corresponding coordinate transformation
!   matrix (transpose of the rotation matrix).
!
!-> Remarks.-
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 function rotvect_psi2mat (Psi)

! I/O Variables.
  real(8),intent(in)    :: Psi(3)           ! Rotation Vector.
  real(8),dimension(3,3):: rotvect_psi2mat  ! Coord transf matrix.

! Local variables
  real(8):: NormPsi        ! Norm of the Rotation Vector.
  real(8):: Normal(3)      ! Unit normal vector to the rotation plane.
  real(8):: RotMatrix(3,3) ! Rotation Matrix.
  real(8):: SkewNormal(3,3)! Skew-symmetric form of the normal vector.
  real(8):: SkewPsi(3,3)   ! Skew-symmetric form of the CRV.

! Compute norm.
  NormPsi=dsqrt(dot_product(Psi,Psi))
  
! Small angles.
  if (abs(NormPsi).le.rot_Epsilon) then
    SkewPsi=rot_skew(Psi)
    RotMatrix= Unit + SkewPsi + 0.5d0*matmul(SkewPsi,SkewPsi)

! Larger angles.
  else
    Normal=Psi/NormPsi
    SkewNormal=rot_skew(Normal)
    RotMatrix= Unit + dsin(NormPsi)*SkewNormal               &
&            + (1.d0-dcos(NormPsi))*matmul(SkewNormal,SkewNormal)
  end if

! The rotation matrix is the transpose of the coordinate transformation matrix.
  rotvect_psi2mat=transpose(RotMatrix)
!
  return
 end function rotvect_psi2mat


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!->Function ROTVECT_MAT2PSI.
!
!->Description.-
!
!   Transform coordinate transformation matrix into Cartesian rotation vector.
!
!-> Remarks.-
!
!  1) Quaternions are used as intermediate step. This is very efficient for large
!     rotations, but it yields zero for rotations below 1d-8 due to round-off
!     errors. It would be better in those cases to use the linear form of the
!     matrix.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 function rotvect_mat2psi (CoordTransMat)

! I/O Variables.
  real(8),intent(in)  :: CoordTransMat(3,3) ! Coordinate transformation matrix.
  real(8),dimension(3):: rotvect_mat2psi    ! Rotation Vector.

! Local variables.
  real(8),dimension(4):: Quat               ! Quaternions.  
  real(8),dimension(3):: Psi                ! Cartesian rotation matrix.

! Rotation Matrix -> Quaternions -> Cartesian rotation vector.
  Quat=rotvect_mat2quat(CoordTransMat)
  Psi= rotvect_quat2psi(Quat)

! For very small rotations, the previous formula fails.
  if (dot_product(Psi,Psi).eq.0.d0) then
    Psi(1)= CoordTransMat(2,3)
    Psi(2)= CoordTransMat(3,1)
    Psi(3)= CoordTransMat(1,2)
  end if

! Set bounds between [-pi,pi] and return the CRV.
  call rotvect_bounds (Psi)
  rotvect_mat2psi=Psi

  return
 end function rotvect_mat2psi
  
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!->Function ROTVECT_PSI2ROT.
!
!->Description.-
!
!   Obtain rotational (tangential) operator from Cartesian rotation vector.
!
!-> Remarks.-
!
!   1) The cartesian rotation vector from 1 to 2 (alpha_1^2) corresponds to the
!      coordinate transformation matrix from 1 to 2: C^2^1.
!
!      A vector given in its components in the base 1 (u_1), is expressed in its
!      components in a base 2 (u_2) as
!
!                  u_2= C^2^1*u_1
!
!  2) Equation (4.11) of (Cardona and Geradin, 2001).
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 function rotvect_psi2rot (Psi)

! I/O Variables. 
  real(8),intent(in)    :: Psi(3)            ! Cartesian rotation vector.
  real(8),dimension(3,3):: rotvect_psi2rot   ! Rotational operator.

! Local variables
  real(8):: K1,K2           ! Auxiliary variables.
  real(8):: NormPsi         ! Norm of Psi.
  real(8):: PsiSkew(3,3)    ! Skew-symmetric form of psi.

! Compute rotation angle.
  NormPsi=dsqrt(dot_product(Psi,Psi))
  PsiSkew= rot_skew(Psi)

! For small angles, use the linear approximation
  if (abs(NormPsi).le.rot_Epsilon) then
    K1=1.d0
    K2=1.d0/6.d0
  else
    K1=dsin(NormPsi/2.d0)/(NormPsi/2.d0)
    K2=(1.d0-dsin(NormPsi)/NormPsi)/(NormPsi*NormPsi)
  end if  

  rotvect_psi2rot= Unit - (0.5d0*K1*K1)*PsiSkew + K2*matmul(PsiSkew,PsiSkew)

  return
 end function rotvect_psi2rot


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!->Function ROTVECT_PSI2INTMAT.
!
!->Description.-
!
!   Transform a rotation vector to the integral of the rotation matrix for
!   an element of lenght 1.
!
!-> Remarks.-
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 function rotvect_psi2intmat (Psi)

! I/O Variables.
  real(8),intent(in)    :: Psi(3)              ! Rotation Vector.
  real(8),dimension(3,3):: rotvect_psi2intmat  ! Integral of coord transf matrix.

! Local variables
  real(8):: NormPsi        ! Norm of the Rotation Vector.
  real(8):: RotInt(3,3)    ! Integral of Rotation Matrix.
  real(8):: SkewPsi(3,3)   ! Skew-symmetric form of the CRV.

! Compute norm.
  NormPsi=dsqrt(dot_product(Psi,Psi))
  SkewPsi=rot_skew(Psi)

! Small angles.
  if (abs(NormPsi).le.rot_Epsilon) then
    RotInt= Unit + (1.d0/2.d0)*SkewPsi + (1.d0/6.d0)*matmul(SkewPsi,SkewPsi)

! Larger angles.
  else
    RotInt= Unit + ((1.d0-dcos(NormPsi))/(NormPsi*NormPsi))*SkewPsi     &
&                + ((NormPsi-dsin(NormPsi))/(NormPsi*NormPsi*NormPsi))  &
&                  *matmul(SkewPsi,SkewPsi)

  end if

! The rotation matrix is the transpose of the coordinate transformation matrix.
  rotvect_psi2intmat=transpose(RotInt)

  return
 end function rotvect_psi2intmat


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!->Function ROTVECT_DROTDPSIB.
!
!->Description.-
!
!   Compute the derivative of the tangential operator with the rotation vector
!   and multiply the result by a given vector, i.e., (dT/dPsi_i)*b_i
!
!-> Remarks.-
!
!   1) dT/dPsi is a third-order tensor. The first two indexes are the rows
!      and columns of the T matrix, while the third index goes along the 
!      dimensions of the Rotation vector. Vector b multiplies that third
!      dimension.
!   2) If b=Psi', the output is T'(Psi).
!   3) If b=PsiDot, the output is Tdot(Psi).
!   4) Note that droptdpsib(Psi,psiDot)*b= A1(Psi,b)*PsiDot, but A1 is not
!      the same as this function.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 function rotvect_drotdpsib (Psi,b)
 
! I/O Variables. 
  real(8),intent(in)    :: Psi(3)              ! Cartesian rotation vector.
  real(8),intent(in)    :: b(3)                ! Input Vector.
  real(8),dimension(3,3):: rotvect_drotdpsib   ! (dT/dPsi)*b

! Local variables
  real(8):: K1,K2           ! Auxiliary variables.
  real(8):: bSkew(3,3)      ! Skew-symmetric form of b.
  real(8):: NormalSkew(3,3) ! Skew-symmetric form of normal vector.
  real(8):: Normal(3)       ! Unit normal vector to the rotation plane.
  real(8):: NormPsi         ! Norm of Psi.
  real(8):: PsiSkew(3,3)    ! Skew-symmetric form of psi.

! Compute rotation angle.
  NormPsi=dsqrt(dot_product(Psi,Psi))
  bSkew=rot_skew(b)

! For small angles, use the linear approximation
  if (abs(NormPsi).le.rot_Epsilon) then
    PsiSkew=rot_skew(Psi)
    rotvect_drotdpsib= &
&     -0.5d0*bSkew + (1.d0/6.d0)*(matmul(PsiSkew,bSkew)+matmul(bSkew,PsiSkew))
    
! For larger angles.
  else
    Normal=Psi/NormPsi
    NormalSkew=rot_skew(Normal)
    K1= (1.d0-dcos(NormPsi))/(NormPsi*NormPsi)
    K2= dsin(NormPsi)/NormPsi

    rotvect_drotdpsib= &
&        (2.d0*K1-K2)*dot_product(b,Normal) * NormalSkew - K1*bSkew              &
&      + ((3.d0*K2-2.d0-dcos(NormPsi))/NormPsi) * dot_product(b,Normal)           &
&                                              * matmul(NormalSkew,NormalSkew)   &
&      + ((1.d0-K2)/NormPsi) * (matmul(NormalSkew,bSkew)+matmul(bSkew,NormalSkew))
  end if

  return
 end function rotvect_drotdpsib


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!->Function ROTVECT_DDROTDPSIB.
!
!->Description.-
!
!  Compute the derivative of the time or spatial derivative of the tangential
!  operator with the rotation vector and multiply the result by a given vector,
!  i.e., (dT'/dPsi)*b
!
!-> Remarks.-
!
!   1) dT'/dPsi is a third-order tensor. The first two indexes are the rows
!      and columns of the T' (or Tdot) matrix, while the third index goes along
!      the dimensions of the Rotation vector.
!
!   2) If b=Psi', the output is part of T"(Psi). T" is obtained as
!      T"=ddrotdpsib(Psi,PsiDot,PsiDot)+drotdpsib(Psi,PsiDotDot)
!
!   3) If b=PsiDot, the output is part of Tdotdot(Psi).
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 function rotvect_ddrotdpsib (Psi,PsiDot,b)

! I/O Variables.
  real(8),intent(in)    :: Psi(3)              ! Cartesian rotation vector.
  real(8),intent(in)    :: PsiDot(3)           ! Derivative of CRV.
  real(8),intent(in)    :: b(3)                ! Input Vector.
  real(8),dimension(3,3):: rotvect_ddrotdpsib  ! (dT'/dPsi)*b

! Local variables
  real(8):: bS(3,3)         ! Skew-symmetric form of b.
  real(8):: nS(3,3)         ! Skew-symmetric form of normal vector.
  real(8):: Normal(3)       ! Unit normal vector to the rotation plane.
  real(8):: NormPsi         ! Norm of Psi.
  real(8):: p1S(3,3)        ! Skew-symmetric form of PsiDot.
  real(8):: pS(3,3)         ! Skew-symmetric form of Psi.

  real(8):: C1,C2           ! Auxiliary variables.
  real(8):: K0,K1,K2
  real(8):: dnb,dnp1,dbp1

! Compute rotation angle.
  NormPsi=dsqrt(dot_product(Psi,Psi))

  bS  = rot_skew(b)
  pS  = rot_Skew(Psi)
  p1S = rot_skew(PsiDot)

! For small angles, use the linear approximation
  if (abs(NormPsi).le.rot_Epsilon) then
    rotvect_ddrotdpsib=                                      &
&               (1.d0/6.d0)*(matmul(bs,p1S)+matmul(p1S,bs))  &
&             + (dot_product(PsiDot,b)  /12.d0)*pS           &
&             + (dot_product(PsiDot,Psi)/12.d0)*bs           &
&             + (dot_product(Psi   ,b)  /12.d0)*p1S

! For larger angles.
  else
    Normal=Psi/NormPsi
    nS    =rot_skew(Normal)

    C1= 1.d0/NormPsi
    C2= C1*C1

    K0= dcos(NormPsi)
    K1= (1.d0-K0)/(NormPsi*NormPsi)
    K2= dsin(NormPsi)/NormPsi

    dnp1= dot_product(Normal,PsiDot)
    dnb = dot_product(Normal,b)
    dbp1= dot_product(b,PsiDot)

    rotvect_ddrotdpsib=                                        &
              - C1*(K0+8.d0*K1-5.d0*K2)*dnb*dnp1*nS            &
&             + C1*(2.d0*K1-K2)*(dbp1*nS+dnp1*bS+dnb*p1S)      &
&             - C2*((15.d0-NormPsi*NormPsi)*K2-8.d0-7.d0*K0)   &
&                 *dnb*dnp1*matmul(nS,nS)                      &
&             + C2*(3.d0*K2-2.d0-K0)                           &
&                 *( dbp1*matmul(nS,nS)+dnp1*matmul(bS,nS)     &
&                  + dnp1*matmul(nS,bS)+ dnb*matmul(nS,p1S)    &
&                  + dnb *matmul(p1S,nS))                      &
&             + C2*(1.d0-K2)*(matmul(bS,p1S)+matmul(p1S,bS))
  end if

  return
 end function rotvect_ddrotdpsib


 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!->Function ROTVECT_A1.
!
!->Description.-
!
!   A1 operator on the tangential operator.
!
!-> Remarks.-
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 function rotvect_a1 (Psi,b)
 
! I/O Variables. 
  real(8),intent(in)    :: Psi(3)              ! Cartesian rotation vector.
  real(8),intent(in)    :: b(3)                ! Input Vector.
  real(8),dimension(3,3):: rotvect_a1          ! (dT/dPsi)*b

! Local variables
  real(8):: K1,K2           ! Auxiliary variables.
  real(8):: bSkew(3,3)      ! Skew-symmetric form of b.
  real(8):: NormalMat(3,3)  ! Outer product of N by itself.
  real(8):: NormalSkew(3,3) ! Skew-symmetric form of normal vector.
  real(8):: Normal(3)       ! Unit normal vector to the rotation plane.
  real(8):: NormPsi         ! Norm of Psi.
  real(8):: PsiSkew(3,3)    ! Skew-symmetric form of psi.

! Compute rotation angle.
  NormPsi=dsqrt(dot_product(Psi,Psi))
  bSkew=rot_skew(b)

! For small angles, use the linear approximation
  if (abs(NormPsi).le.rot_Epsilon) then
    PsiSkew=rot_skew(Psi)
    rotvect_a1= 0.5d0*bSkew + (1.d0/6.d0)*matmul(bSkew,PsiSkew) &
&                           - (1.d0/3.d0)*matmul(PsiSkew,bSkew)

! For larger angles.
  else
    Normal=Psi/NormPsi
    NormalSkew=rot_skew(Normal)
    K1= (1.d0-dcos(NormPsi))/(NormPsi*NormPsi)
    K2= dsin(NormPsi)/NormPsi
    NormalMat=rot_outprod(Normal,Normal)

    rotvect_a1=                                                           &
&       -(2.d0*K1-K2)*matmul(bSkew,NormalMat) + K1*bSkew                  &
&       -((3.d0*K2-2.d0-dcos(NormPsi))/NormPsi) * matmul(NormalSkew,      &
&                                                matmul(bSkew,NormalMat)) &
&       +((1.d0-K2)/NormPsi) * (      matmul(bSkew,NormalSkew)            &
&                               -2.d0*matmul(NormalSkew,bSkew))
  end if

  return
 end function rotvect_a1

 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!->Function ROTVECT_A2.
!
!->Description.-
!
!   A2 operator on the tangential operator.
!
!-> Remarks.-
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 function rotvect_a2 (Psi,b)
 
! I/O Variables. 
  real(8),intent(in)    :: Psi(3)              ! Cartesian rotation vector.
  real(8),intent(in)    :: b(3)                ! Input Vector.
  real(8),dimension(3,3):: rotvect_a2          ! (dT/dPsi)*b

  rotvect_a2= -rotvect_a1(-Psi,b)

 end function rotvect_a2



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!->Function ROTVECT_B1.
!
!->Description.-
!
!   B1 operator on the tangential operator.
!
!-> Remarks.-
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 function rotvect_b1 (Psi,PsiDot,b)
 
! I/O Variables. 
  real(8),intent(in)    :: Psi(3)              ! Cartesian rotation vector.
  real(8),intent(in)    :: PsiDot(3)           ! Derivative of Psi.
  real(8),intent(in)    :: b(3)                ! Input Vector.
  real(8),dimension(3,3):: rotvect_b1          ! (dTdot/dPsi)*b

! Local variables
  real(8):: bS(3,3)         ! Skew-symmetric form of b.
  real(8):: nS(3,3)         ! Skew-symmetric form of normal vector.
  real(8):: Normal(3)       ! Unit normal vector to the rotation plane.
  real(8):: NormPsi         ! Norm of Psi.
  real(8):: p1S(3,3)        ! Skew-symmetric form of PsiDot.

  real(8):: C1,C2           ! Auxiliary variables.
  real(8):: K0,K1,K2,K3
  real(8),dimension(3,3):: nnT,bSnnT,nSbSnnT,np1T,p1nT,nSp1S2bnT

! Compute rotation angle.
  NormPsi=dsqrt(dot_product(Psi,Psi))

  bS  = rot_skew(b)
  p1S = rot_skew(PsiDot)

! For small angles, use the linear approximation
  if (abs(NormPsi).le.rot_Epsilon) then
    rotvect_b1= (1.d0/6.d0)*matmul(bS,p1S)                      &
&             - (1.d0/3.d0)*matmul(p1S,bS)                      &
&             - (1.d0/12.d0)*matmul(bS,rot_outprod(Psi,PsiDot)  &
&                                     +rot_outprod(PsiDot,Psi)) &
&             - (1.d0/12.d0)*dot_product(Psi,PsiDot)*bs

! For larger angles.
  else
    Normal=Psi/NormPsi
    nS    =rot_skew(Normal)

    C1= 1.d0/NormPsi
    C2= C1*C1
    
    K0= dcos(NormPsi)
    K1= (1.d0-K0)/(NormPsi*NormPsi)
    K2= dsin(NormPsi)/NormPsi
    K3= dot_product(Normal,PsiDot)

    nnT    = rot_outprod(Normal,Normal)
    bSnnT  = matmul(bS,nnT)
    nSbSnnT= matmul(nS,bSnnT)
    
    np1T   = rot_outprod(Normal,PsiDot)
    p1nT   = rot_outprod(PsiDot,Normal)
    nSp1S2bnT= matmul((matmul(nS,p1S)+matmul(p1S,nS)), &
&                     rot_outprod(b,Normal))

    rotvect_b1= C1*(K0+8.d0*K1-5.d0*K2)*K3*bSnnT                &
&             - C1*(2.d0*K1-K2)                                 &
&                 *(matmul(bS,np1T+p1nT)                        &
&                  +bs*dot_product(Normal,PsiDot))              &
&             + C2*((15.d0-NormPsi*NormPsi)*K2-8.d0-7.d0*K0)    &
&                 *K3*nSbSnnT                                   &
&             + C2*(3.d0*K2-2.d0-K0)                            &
&                 *(nSp1S2bnT-matmul(nS,matmul(bS,np1T))        &
&                   +K3*(matmul(bS,nS)-2.d0*matmul(nS,bS)))     &
&             + C2*(1.d0-K2)                                    &
&                 *(matmul(bS,p1S)-2.d0*matmul(p1S,bS))
  end if

  return
 end function rotvect_b1


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!->Function ROTVECT_B2.
!
!->Description.-
!
!   B2 operator on the tangential operator.
!
!-> Remarks.-
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 function rotvect_b2 (Psi,PsiDot,b)
 
! I/O Variables. 
  real(8),intent(in)    :: Psi(3)              ! Cartesian rotation vector.
  real(8),intent(in)    :: PsiDot(3)           ! Derivative of Psi.
  real(8),intent(in)    :: b(3)                ! Input Vector.
  real(8),dimension(3,3):: rotvect_b2          ! (dT/dPsi)*b

  rotvect_b2= rotvect_b1(-Psi,PsiDot,b)

 end function rotvect_b2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end module lib_rotvect
