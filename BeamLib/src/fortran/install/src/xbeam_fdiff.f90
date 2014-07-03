!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!-> Module.- XBEAM_FDIFF Henrik Hesse. 16/10/2010 - Last Update 07/01/2011
!
!-> Language: FORTRAN90, Free format.
!
!-> Description.-
!
!   Check linearisation of nonlinear eom of coupled system.
!
!-> Subroutines.-
!
!
!-> Remarks.-
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module xbeam_fdiff
  use lib_fem
  use lib_rot
  use lib_rotvect
  use lib_cbeam3
  use lib_xbeam
  use xbeam_shared

  implicit none

! Shared variables within the module.
  integer,private,parameter:: MaxNodCB3=3               ! Max number of nodes per element is 3.

  contains


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-> Subroutine FDIFF_CHECK
!
!-> Description:
!
!    Compute tangent matrices using finite differences.
!
!-> Remarks.-
!
!   1) This routine is only used to test the implementation of the
!      analytical matrices.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine fdiff_check (NumNE,rElem0,rElem,rElemDot,rElemDDot,Vrel,VrelDot,ElemMass,ElemStiff,NumGauss)
! I/O Variables.
  integer,intent(in) :: NumNE             ! Number of nodes in the element.
  real(8),intent(in) :: rElem0    (:,:)   ! Initial position/orientation of grid points.
  real(8),intent(in) :: rElem     (:,:)   ! Current position/orientation of grid points.
  real(8),intent(in) :: rElemDot  (:,:)   ! Current velocities of grid points.
  real(8),intent(in) :: rElemDDot (:,:)   ! Current velocities of grid points.
  real(8),intent(in) :: Vrel      (:)     ! Current velocities of a frame
  real(8),intent(in) :: VrelDot   (:)     ! Current velocities of a frame
  real(8),intent(in) :: ElemMass  (:,:)   ! Stiffness properties in the element.
  real(8),intent(in) :: ElemStiff (:,:)   ! Stiffness properties in the element.
  integer,intent(in) :: NumGauss          ! Number of Gauss points in the element..

! Local variables.
!  real(8)            :: DERIV(6*MaxNodCB3,6*MaxNodCB3), FDIFF(6*MaxNodCB3,6*MaxNodCB3)
!  real(8)            :: DERIV(6*MaxNodCB3,6), FDIFF(6*MaxNodCB3,6)
  real(8)            :: DERIV(6,6*MaxNodCB3), FDIFF(6,6*MaxNodCB3)
!  real(8)            :: DERIV(6,6), FDIFF(6,6)

! Initialize.
  DERIV=0.E0
  FDIFF=0.E0

! check if linearised form correct

!   KRS gyro
    call xbeam_kgyr   (NumNE,rElem0,rElem,rElemDot,rElemDDot,Vrel,VrelDot,ElemMass,DERIV,NumGauss)
    call fdiff_KRSgyr (NumNE,rElem0,rElem,rElemDot,rElemDDot,Vrel,VrelDot,ElemMass,FDIFF,NumGauss)

!   KRS mass
!    call xbeam_kmass    (NumNE,rElem0,rElem,rElemDDot,VrelDot,ElemMass,DERIV,NumGauss)
!    call fdiff_KRSmass  (NumNE,rElem0,rElem,rElemDot,rElemDDot,Vrel,VrelDot,ElemMass,FDIFF,NumGauss)

!   KSS gyro
!    call cbeam3_kgyr  (NumNE,rElem0,rElem,rElemDot,rElemDDot,Vrel,VrelDot,ElemMass,DERIV,NumGauss)
!    call fdiff_KSSgyr (NumNE,rElem0,rElem,rElemDot,rElemDDot,Vrel,VrelDot,ElemMass,FDIFF,NumGauss)

!   KSS stiff
!    call cbeam3_kmat     (NumNE,rElem0,rElem,ElemStiff,DERIV,NumGauss)
!    call cbeam3_kgeom    (NumNE,rElem0,rElem,ElemStiff,DERIV,NumGauss)
!    call fdiff_KSSstif   (NumNE,rElem0,rElem,ElemStiff,FDIFF,NumGauss)

!   CSR
!    call cbeam3_cvel (NumNE,rElem0,rElem,rElemDot,Vrel,ElemMass,DERIV,NumGauss)
!    call fdiff_CSR   (NumNE,rElem0,rElem,rElemDot,Vrel,ElemMass,FDIFF,NumGauss)

!   CRS
!    call xbeam_cgyr (NumNE,rElem0,rElem,rElemDot,Vrel,ElemMass,DERIV,NumGauss)
!    call fdiff_CRS  (NumNE,rElem0,rElem,rElemDot,Vrel,ElemMass,FDIFF,NumGauss)

!   CSS
!    call cbeam3_cgyr  (NumNE,rElem0,rElem,rElemDot,Vrel,ElemMass,DERIV,NumGauss)
!    call fdiff_CSS    (NumNE,rElem0,rElem,rElemDot,Vrel,ElemMass,FDIFF,NumGauss)

!   CRR
!    call xbeam_CRR  (NumNE,rElem0,rElem,rElemDot,Vrel,ElemMass,DERIV,NumGauss)
!    call fdiff_CRR  (NumNE,rElem0,rElem,rElemDot,Vrel,ElemMass,FDIFF,NumGauss)

!    open (unit=12,file='finite.txt',ACCESS='append')
!     write (12,'(1P18E12.3)') transpose(DERIV-FDIFF)/(maxval(abs(DERIV))-1)
!    close (12)

  return
 end subroutine fdiff_check

 subroutine fdiff_KRSmass (NumNE,rElem0,rElem,rElemDot,rElemDDot,Vrel,VrelDot,ElemMass,KRS,NumGauss)
! I/O Variables.
  integer,intent(in) :: NumNE              ! Number of nodes in the element.
  real(8),intent(in) :: rElem0   (:,:)     ! Initial position/orientation of grid points.
  real(8),intent(in) :: rElem    (:,:)     ! Current position/orientation of grid points.
  real(8),intent(in) :: rElemDot (:,:)     ! Current velocities of grid points.
  real(8),intent(in) :: rElemDDot(:,:)     ! Current velocities of grid points.
  real(8),intent(in) :: Vrel     (:)       ! Current velocities of a frame
  real(8),intent(in) :: VrelDot  (:)       ! Current velocities of a frame
  real(8),intent(in) :: ElemMass (:,:)     ! Stiffness properties in the element.
  real(8),intent(out):: KRS      (:,:)     ! Geometric tangent stiffness matrix.
  integer,intent(in) :: NumGauss           ! Number of Gauss points in the element..

! Local variables.
  integer            :: i,j,k           ! Counters.
  real(8),parameter  :: Delta=1.d-6     ! Delta for finite-differences.
  real(8)            :: drElem(MaxNodCB3,6), RiDDot(6*MaxNodCB3)
  real(8)            :: MRS(6,6*MaxNodCB3), MRR(6,6)
  real(8),allocatable:: Q0(:),Q1(:)     ! Values of the element stiffness functional.

! Initialize.
  allocate (Q0(6)); Q0=0.d0
  allocate (Q1(6)); Q1=0.d0

  k=0
  do i=1,NumNE
    do j=1,6
        k=k+1
        RiDDot(k)=rElemDDot(i,j)
    end do
  end do

! Compute force vector at reference point.
  call xbeam_mrs (NumNE,rElem0,rElem,ElemMass,MRS,NumGauss)
  call xbeam_mrr (NumNE,rElem0,rElem,ElemMass,MRR,NumGauss)
  Q0=matmul(MRS,RiDDot)+matmul(MRR,VrelDot)

! Loop in the element degrees of freedom.
  k=0
  do i=1,NumNE
    do j=1,6
      k=k+1
      drElem=0.d0
      drElem(i,j)=Delta
      Q1=0.d0

      call xbeam_mrs (NumNE,rElem0,rElem+drElem,ElemMass,MRS,NumGauss)
      call xbeam_mrr (NumNE,rElem0,rElem+drElem,ElemMass,MRR,NumGauss)
      Q1=matmul(MRS,RiDDot)+matmul(MRR,VrelDot)

      KRS(:,k)=(Q1-Q0)/Delta
    end do
  end do

  return
end subroutine fdiff_KRSmass


subroutine fdiff_KRSgyr (NumNE,rElem0,rElem,rElemDot,rElemDDot,Vrel,VrelDot,ElemMass,KRS,NumGauss)
! I/O Variables.
  integer,intent(in) :: NumNE              ! Number of nodes in the element.
  real(8),intent(in) :: rElem0   (:,:)     ! Initial position/orientation of grid points.
  real(8),intent(in) :: rElem    (:,:)     ! Current position/orientation of grid points.
  real(8),intent(in) :: rElemDot (:,:)     ! Current velocities of grid points.
  real(8),intent(in) :: rElemDDot(:,:)     ! Current velocities of grid points.
  real(8),intent(in) :: Vrel     (:)       ! Current velocities of a frame
  real(8),intent(in) :: VrelDot  (:)       ! Current velocities of a frame
  real(8),intent(in) :: ElemMass (:,:)     ! Stiffness properties in the element.
  real(8),intent(out):: KRS      (:,:)     ! Geometric tangent stiffness matrix.
  integer,intent(in) :: NumGauss           ! Number of Gauss points in the element..

! Local variables.
  integer            :: i,j,k           ! Counters.
  real(8),parameter  :: Delta=1.d-6     ! Delta for finite-differences.
  real(8)            :: drElem(MaxNodCB3,6)
  real(8),allocatable:: Q0(:),Q1(:)     ! Values of the element stiffness functional.

! Initialize.
  allocate (Q0(6)); Q0=0.d0
  allocate (Q1(6)); Q1=0.d0

! Compute force vector at reference point.
  call xbeam_fgyr (NumNE,rElem0,rElem,rElemDot,Vrel,ElemMass,Q0,NumGauss)

! Loop in the element degrees of freedom.
  k=0
  do i=1,NumNE
    do j=1,6
      k=k+1
      drElem=0.d0
      drElem(i,j)=Delta
      Q1=0.d0
      call xbeam_fgyr (NumNE,rElem0,rElem+drElem,rElemDot,Vrel,ElemMass,Q1,NumGauss)
      KRS(:,k)=(Q1-Q0)/Delta
    end do
  end do

  return
 end subroutine fdiff_KRSgyr



subroutine fdiff_KSSgyr (NumNE,rElem0,rElem,rElemDot,rElemDDot,Vrel,VrelDot,ElemMass,KSS,NumGauss)
! I/O Variables.
  integer,intent(in) :: NumNE              ! Number of nodes in the element.
  real(8),intent(in) :: rElem0   (:,:)     ! Initial position/orientation of grid points.
  real(8),intent(in) :: rElem    (:,:)     ! Current position/orientation of grid points.
  real(8),intent(in) :: rElemDot (:,:)     ! Current velocities of grid points.
  real(8),intent(in) :: rElemDDot(:,:)     ! Current velocities of grid points.
  real(8),intent(in) :: Vrel     (:)       ! Current velocities of a frame
  real(8),intent(in) :: VrelDot  (:)       ! Current velocities of a frame
  real(8),intent(in) :: ElemMass (:,:)     ! Stiffness properties in the element.
  real(8),intent(out):: KSS      (:,:)     ! Geometric tangent stiffness matrix.
  integer,intent(in) :: NumGauss           ! Number of Gauss points in the element..

! Local variables.
  integer            :: i,j,k           ! Counters.
  real(8),parameter  :: Delta=1.d-6     ! Delta for finite-differences.
  real(8)            :: drElem(MaxNodCB3,6)
  real(8),allocatable:: Q0(:),Q1(:)     ! Values of the element stiffness functional.

! Initialize.
  allocate (Q0(6*MaxNodCB3)); Q0=0.d0
  allocate (Q1(6*MaxNodCB3)); Q1=0.d0

! Compute force vector at reference point.
  call cbeam3_fgyr (NumNE,rElem0,rElem,rElemDot,Vrel,ElemMass,Q0,NumGauss)

! Loop in the element degrees of freedom.
  k=0
  do i=1,NumNE
    do j=1,6
      k=k+1
      drElem=0.d0
      drElem(i,j)=Delta
      Q1=0.d0
      call cbeam3_fgyr (NumNE,rElem0,rElem+drElem,rElemDot,Vrel,ElemMass,Q1,NumGauss)
      KSS(:,k)=(Q1-Q0)/Delta
    end do
  end do

  return
end subroutine fdiff_KSSgyr


subroutine fdiff_KSSstif (NumNE,rElem0,rElem,ElemStiff,KSS,NumGauss)
! I/O Variables.
  integer,intent(in) :: NumNE              ! Number of nodes in the element.
  real(8),intent(in) :: rElem0   (:,:)     ! Initial position/orientation of grid points.
  real(8),intent(in) :: rElem    (:,:)     ! Current position/orientation of grid points.
  real(8),intent(in) :: ElemStiff (:,:)    ! Stiffness properties in the element.
  real(8),intent(out):: KSS      (:,:)     ! Geometric tangent stiffness matrix.
  integer,intent(in) :: NumGauss           ! Number of Gauss points in the element..

! Local variables.
  integer            :: i,j,k           ! Counters.
  real(8),parameter  :: Delta=1.d-6     ! Delta for finite-differences.
  real(8)            :: drElem(MaxNodCB3,6)
  real(8),allocatable:: Q0(:),Q1(:)     ! Values of the element stiffness functional.

! Initialize.
  allocate (Q0(6*MaxNodCB3)); Q0=0.d0
  allocate (Q1(6*MaxNodCB3)); Q1=0.d0

! Compute force vector at reference point.
  call cbeam3_fstif (NumNE,rElem0,rElem,ElemStiff,Q0,NumGauss)

! Loop in the element degrees of freedom.
  k=0
  do i=1,NumNE
    do j=1,6
      k=k+1
      drElem=0.d0
      drElem(i,j)=Delta
      Q1=0.d0
      call cbeam3_fstif (NumNE,rElem0,rElem+drElem,ElemStiff,Q1,NumGauss)
      KSS(:,k)=(Q1-Q0)/Delta
    end do
  end do

  return
end subroutine fdiff_KSSstif


subroutine fdiff_CSS (NumNE,rElem0,rElem,rElemDot,Vrel,ElemMass,CSS,NumGauss)
! I/O Variables.
  integer,intent(in) :: NumNE               ! Number of nodes in the element.
  real(8),intent(in) :: rElem0   (:,:)      ! Initial position/orientation of grid points.
  real(8),intent(in) :: rElem    (:,:)      ! Current position/orientation of grid points.
  real(8),intent(in) :: rElemDot (:,:)      ! Current velocities of grid points.
  real(8),intent(in) :: Vrel     (:)        ! Current velocities of a frame
  real(8),intent(in) :: ElemMass (:,:)      ! Stiffness properties in the element.
  real(8),intent(out):: CSS      (:,:)      ! Geometric tangent stiffness matrix.
  integer,intent(in) :: NumGauss            ! Number of Gauss points in the element..

! Local variables.
  integer            :: i,j,k               ! Counters.
  real(8),parameter  :: Delta=1.E-6         ! Delta for finite-differences.
  real(8)            :: drElem(MaxNodCB3,6)
  real(8),allocatable:: Q0(:),Q1(:)         ! Values of the element stiffness functional.

! Initialize.
  allocate (Q0(6*MaxNodCB3)); Q0=0.d0
  allocate (Q1(6*MaxNodCB3)); Q1=0.d0

! Compute force vector at reference point.
  call cbeam3_fgyr (NumNE,rElem0,rElem,rElemDot,Vrel,ElemMass,Q0,NumGauss)

! Loop in the element degrees of freedom.
  k=0
  do i=1,NumNE
    do j=1,6
      k=k+1
      drElem=0.d0
      drElem(i,j)=Delta
      Q1=0.d0
      call cbeam3_fgyr (NumNE,rElem0,rElem,rElemDot+drElem,Vrel,ElemMass,Q1,NumGauss)
      CSS(:,k)=(Q1-Q0)/Delta
    end do
  end do

  return
end subroutine fdiff_CSS


subroutine fdiff_CSR (NumNE,rElem0,rElem,rElemDot,Vrel,ElemMass,CSR,NumGauss)
! I/O Variables.
  integer,intent(in) :: NumNE               ! Number of nodes in the element.
  real(8),intent(in) :: rElem0   (:,:)      ! Initial position/orientation of grid points.
  real(8),intent(in) :: rElem    (:,:)      ! Current position/orientation of grid points.
  real(8),intent(in) :: rElemDot (:,:)      ! Current velocities of grid points.
  real(8),intent(in) :: Vrel     (:)        ! Current velocities of a frame
  real(8),intent(in) :: ElemMass (:,:)      ! Stiffness properties in the element.
  real(8),intent(out):: CSR   (:,:)         ! Geometric tangent stiffness matrix.
  integer,intent(in) :: NumGauss            ! Number of Gauss points in the element..

! Local variables.
  integer            :: j               ! Counters.
  real(8),parameter  :: Delta=1.E-6     ! Delta for finite-differences.
  real(8)            :: dVrel(6)
  real(8),allocatable:: Q0(:),Q1(:)     ! Values of the element stiffness functional.

! Initialize.
  allocate (Q0(6*MaxNodCB3)); Q0=0.d0
  allocate (Q1(6*MaxNodCB3)); Q1=0.d0

! Compute force vector at reference point.
  call cbeam3_fgyr (NumNE,rElem0,rElem,rElemDot,Vrel,ElemMass,Q0,NumGauss)

! Loop in the element degrees of freedom.
  do j=1,6
    Q1=0.d0
    dVrel = 0.d0
    dVrel(j)=Delta
    call cbeam3_fgyr (NumNE,rElem0,rElem,rElemDot,Vrel+dVrel,ElemMass,Q1,NumGauss)
    CSR(:,j)=(Q1-Q0)/Delta
  end do

  return
end subroutine fdiff_CSR



subroutine fdiff_CRS (NumNE,rElem0,rElem,rElemDot,Vrel,ElemMass,CRS,NumGauss)
! I/O Variables.
  integer,intent(in) :: NumNE              ! Number of nodes in the element.
  real(8),intent(in) :: rElem0   (:,:)     ! Initial position/orientation of grid points.
  real(8),intent(in) :: rElem    (:,:)     ! Current position/orientation of grid points.
  real(8),intent(in) :: rElemDot (:,:)     ! Current velocities of grid points.
  real(8),intent(in) :: Vrel     (:)     ! Current velocities of a frame
  real(8),intent(in) :: ElemMass (:,:)   ! Stiffness properties in the element.
  real(8),intent(out):: CRS   (:,:)    ! Geometric tangent stiffness matrix.
  integer,intent(in) :: NumGauss        ! Number of Gauss points in the element..

! Local variables.
  integer            :: i,j,k           ! Counters.
  real(8),parameter  :: Delta=1.E-6     ! Delta for finite-differences.
  real(8)            :: drElem(MaxNodCB3,6)
  real(8),allocatable:: Q0(:),Q1(:)     ! Values of the element stiffness functional.

! Initialize.
  allocate (Q0(6)); Q0=0.d0
  allocate (Q1(6)); Q1=0.d0

! Compute force vector at reference point.
  call xbeam_fgyr (NumNE,rElem0,rElem,rElemDot,Vrel,ElemMass,Q0,NumGauss)

! Loop in the element degrees of freedom.
  k=0
  do i=1,NumNE
    do j=1,6
      k=k+1
      drElem=0.d0
      drElem(i,j)=Delta
      Q1=0.d0
      call xbeam_fgyr (NumNE,rElem0,rElem,rElemDot+drElem,Vrel,ElemMass,Q1,NumGauss)
      CRS(:,k)=(Q1-Q0)/Delta
    end do
  end do

  return
end subroutine fdiff_CRS



subroutine fdiff_CRR (NumNE,rElem0,rElem,rElemDot,Vrel,ElemMass,CRR,NumGauss)
! I/O Variables.
  integer,intent(in) :: NumNE               ! Number of nodes in the element.
  real(8),intent(in) :: rElem0   (:,:)      ! Initial position/orientation of grid points.
  real(8),intent(in) :: rElem    (:,:)      ! Current position/orientation of grid points.
  real(8),intent(in) :: rElemDot (:,:)      ! Current velocities of grid points.
  real(8),intent(in) :: Vrel     (:)        ! Current velocities of a frame
  real(8),intent(in) :: ElemMass (:,:)      ! Stiffness properties in the element.
  real(8),intent(out):: CRR   (:,:)         ! Geometric tangent stiffness matrix.
  integer,intent(in) :: NumGauss            ! Number of Gauss points in the element..

! Local variables.
  integer            :: j               ! Counters.
  real(8),parameter  :: Delta=1.E-6     ! Delta for finite-differences.
  real(8)            :: dVrel(6)
  real(8),allocatable:: Q0(:),Q1(:)     ! Values of the element stiffness functional.

! Initialize.
  allocate (Q0(6)); Q0=0.d0
  allocate (Q1(6)); Q1=0.d0

! Compute force vector at reference point.
  call xbeam_fgyr  (NumNE,rElem0,rElem,rElemDot,Vrel,ElemMass,Q0,NumGauss)

! Loop in the element degrees of freedom.
  do j=1,6
    Q1=0.d0
    dVrel = 0.d0
    dVrel(j)=Delta
    call xbeam_fgyr (NumNE,rElem0,rElem,rElemDot,Vrel+dVrel,ElemMass,Q1,NumGauss)
    CRR(:,j)=(Q1-Q0)/Delta
  end do

  return
end subroutine fdiff_CRR



end module xbeam_fdiff

