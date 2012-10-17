!->Module INTERFACE_LAPACK. Rafa Palacios. 04Oct2006
!
!->Description.-
!
!  This module includes interface for real*8 LAPACK Version 3 routines.
!
!->Subroutines:
!
!     -lapack_inv:         Inverse of general matrix.
!     -lapack_lufact:      LU factorization of general matrix.
!     -lapack_luback:      LU backsubstution of general matrix.
!     -lapack_band_lufact: LU factorization of banded matrix.
!     -lapack_band_luback: LU backsubstituion of banded matrix.
!     -lapack_band_inv:    Inverse of banded matrix.
!     -lapack_sparse:      Solve linear system defined with a sparse matrix.
!     -lapack_sparse_inv:  Sparse inverse of a sparse matrix.
!     -lapack_nonsymeigv:  Complex eigenvalues of generalized eigv problem.
!
!->Remarks.-
!
!  1) External routines are part of the LAPACK 3 library, which needs
!     to be compiled with this module.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module interface_lapack
 implicit none

!   (There are no public variables in this module).

 contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!->Subroutine LAPACK_INV
!
!->Description.-
!
!   Computes the inverse of a general matrix, using the LU factorization
!
!->Remarks.-
!
!  1) LAPACK routines used:
!
!      -dgetrf:
!        Computes an LU factorization of a general matrix, using partial
!        pivoting with row interchanges.
!
!      -dgetrI:
!         Computes the inverse of a general matrix, using the LU factorization
!         computed by DGETRF.  
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine lapack_inv (N,A,invA,Info)

  integer,intent(in) :: N         ! Size of the matrix.
  real(8),intent(in) :: A   (:,:) ! Matrix to be inverted.
  real(8),intent(out):: invA(:,:) ! A^-1.
  integer,intent(out):: Info      ! Error codes.

  integer,allocatable:: LUpivot(:)
  real(8),allocatable:: LUwork (:)

! Allocate auxiliary memory.

  allocate (LUpivot(N));  LUpivot=0
  allocate (LUwork(N*N)); LUwork =0.d0

! Perform LU decomposition.

  invA=A
  call DGETRF (N,N,invA,N,LUpivot,Info)
  if (Info.ne.0) return

! Perform matrix inversion.

  call DGETRI (N,invA(1:N,1:N),N,LUpivot(1:N),LUwork(1:N*N),N*N,Info)

  deallocate (LUpivot,LUwork)
  return
 end subroutine lapack_inv


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!->Subroutine LAPACK_LUFACT
!
!->Description.-
!
!   Computes the LU factorization of a general matrix.
!
!->Remarks.-
!
!  1) LAPACK routines used:
!
!      -dgetrf:
!        Computes an LU factorization of a general matrix, using partial
!        pivoting with row interchanges.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine lapack_lufact (N,A,LUpivot,Info)

  integer,intent(in)   :: N             ! Size of the matrix.
  real(8),intent(inout):: A   (1:N,1:N) ! LU factorization of A.
  integer,intent(out)  :: LUpivot(1:N)  ! LU pivoting index.
  integer,intent(out)  :: Info          ! Error codes.

! Perform LU decomposition.

  call DGETRF (N,N,A,N,LUpivot,Info)

  return
 end subroutine lapack_lufact


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!->Subroutine LAPACK_LUBACK
!
!->Description.-
!
!   Solves a general system of linear equations AX=B.
!
!->Remarks.-
!
!  1) LAPACK routines used:
!
!      -dgetrs:
!        Solves a general system of linear equations AX=B, A**T X=B
!        or A**H X=B, using the LU factorization computed by DGETRF.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine lapack_luback (N,A,LUpivot,B,Info)

  integer,intent(in)   :: N                ! Size of the matrix.
  real(8),intent(in)   :: A      (1:N,1:N) ! Matrix and its LU factorization.
  integer,intent(in)   :: LUpivot(1:N)     ! LU pivoting index.
  real(8),intent(inout):: B      (1:N)     ! B in input, X in output.
  integer,intent(out)  :: Info             ! Error codes.

! Local variables.
  real(8),allocatable  :: X(:,:)

  allocate (X(N,1))
  X(1:N,1)=B(1:N)

! Solve AX=B.
  call DGETRS ('N',N,1,A,N,LUpivot,X(:,1:1),N,Info)

  B(1:N)=X(1:N,1)

  deallocate (X)
  return
 end subroutine lapack_luback


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!->Subroutine LAPACK_BAND_LUFACT
!
!->Description.-
!
!   Computes the LU factorization of a banded matrix.
!
!->Remarks.-
!
!  1) LAPACK routines used:
!
!      -dgbtrf:
!        Computes an LU factorization of a banded matrix, using partial
!        pivoting with row interchanges.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine lapack_band_lufact (N,KL,KU,A,LUpivot,Info)

  integer,intent(in)   :: N         ! Size of the matrix.
  integer,intent(in)   :: KL        ! Number of subdiagonals in the band of B.
  integer,intent(in)   :: KU        ! Number of superdiagonals in the band of B.
  real(8),intent(inout):: A   (:,:) ! A; LU factorization of A.
  integer,intent(out)  :: LUpivot(:)! LU pivoting index.
  integer,intent(out)  :: Info      ! Error codes.

! Perform LU decomposition.

  call DGBTRF (N,N,KL,KU,A(1:2*KL+KU+1,1:N),2*KL+KU+1,LUpivot(1:N),Info)

  return
 end subroutine lapack_band_lufact


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!->Subroutine LAPACK_BAND_LUBACK
!
!->Description.-
!
!   Solves a system of linear equations AX=B with matrix in banded form.
!
!->Remarks.-
!
!  1) LAPACK routines used:
!
!      -dgbtrs:
!        Solves a general system of linear equations AX=B, A**T X=B
!        or A**H X=B, using the LU factorization computed by DGETRF.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine lapack_band_luback (N,KL,KU,A,LUpivot,B,Info)

  integer,intent(in)   :: N                ! Size of the matrix.
  integer,intent(in)   :: KL               ! Number of subdiagonals in the band of B.
  integer,intent(in)   :: KU               ! Number of superdiagonals in the band of B.
  real(8),intent(inout):: A   (:,:)        ! A; LU factorization of A.
  integer,intent(in)   :: LUpivot(1:N)     ! LU pivoting index.
  real(8),intent(inout):: B      (1:N)     ! B in input, X in output.
  integer,intent(out)  :: Info             ! Error codes.

! Local variables.
  real(8),allocatable  :: X(:,:)

  allocate (X(N,1))
  X(1:N,1)=B(1:N)

! Solve AX=B.

  call DGBTRS ('N',N,KL,KU,1,A,2*KL+KU+1,LUpivot,X(:,1:1),N,Info)

  B(1:N)=X(1:N,1)

  deallocate (X)
  return
 end subroutine lapack_band_luback

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!->Subroutine LAPACK_BAND_INV
!
!->Description.-
!
!   Computes the inverse of a banded matrix, using the LU factorization
!
!->Remarks.-
!
!  1) LAPACK routines used:
!
!      -dgbtrf:
!        Computes an LU factorization of a banded matrix, using partial
!        pivoting with row interchanges.
!
!      -dgbtri:
!         Computes the inverse of a banded matrix, using the LU factorization
!         computed by DGBTRF.  
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine lapack_band_inv (N,KL,KU,A,invA,Info)

  integer,intent(in) :: N         ! Size of the matrix.
  integer,intent(in) :: KL        ! Number of subdiagonals in the band of B.
  integer,intent(in) :: KU        ! Number of superdiagonals in the band of B.
  real(8),intent(in) :: A   (:,:) ! Banded atrix to be inverted.
  real(8),intent(out):: invA(:,:) ! A^-1 (full matrix)
  integer,intent(out):: Info      ! Error codes.

! Local variables.

  integer:: i
  real(8),allocatable:: LUdecom(:,:)
  integer,allocatable:: LUpivot(:)
  real(8),allocatable:: LUwork (:)

! Allocate auxiliary memory.

  allocate (LUpivot(N));  LUpivot=0
  allocate (LUwork(N*N)); LUwork =0.d0
  allocate (LUdecom(2*KL+KU+1,N)); LUdecom=A

! Perform LU decomposition.

  call DGBTRF (N,N,KL,KU,LUdecom,2*KL+KU+1,LUpivot,Info)
  if (Info.ne.0) return

! Perform matrix inversion.

  invA=0.d0
  do i=1,N
    invA(i,i)=1.d0
    call DGBTRS ('N',N,KL,KU,1,LUdecom,2*KL+KU+1,LUpivot,invA(i,:),N,Info)
  end do

  deallocate (LUpivot,LUwork,LUdecom)
  return
 end subroutine lapack_band_inv


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!->Subroutine LAPACK_SPARSE.
!
!->Description.-
!
!   Solve linear system defined by a sparse matrix and a forcing vector.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine lapack_sparse (DimArray,SprMat,b,X)
  use lib_sparse

! I/O Variables.
  integer,     intent(in) :: DimArray  ! Current storage dimension.
  type(sparse),intent(in) :: SprMat(:) ! Sparse matrix.
  real(8),     intent(in) :: b(:)      ! Forcing vector.
  real(8),     intent(out):: X(:)       ! Solution vector.

! Local variables.
  integer:: Info
  integer:: KL,KU                    ! Upper and lower bandwidth of matrix.
  integer:: N                        ! Size of the problem.
  real(8),allocatable:: Aband(:,:)   ! Matrix in banded form.
  integer,allocatable:: LUpoint(:)   ! Pointer in the LU decomposition.

! Size of the matrices and vectors.
  N=size(X)
  X=b
  
! Obtain banded matrix.
  call sparse_getbandwidth (DimArray,SprMat,KL,KU)
  allocate (LUPoint(N));  LUPoint=0
  allocate(Aband(2*KL+KU+1,N)); Aband  = 0.d0
  call sparse_copy2band (DimArray,SprMat,KL,KU,Aband)

! LU factorization of banded matrix.
  call lapack_band_lufact (N,KL,KU,Aband,LUPoint,Info)
  call lapack_band_luback (N,KL,KU,Aband,LUPoint,X,Info)

  if (Info.ne.0) STOP 'Error in the LU decomposition'

  deallocate (LUPoint,Aband)
  return 
end subroutine lapack_sparse


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!->Subroutine LAPACK_SPARSE_INV.
!
!->Description.-
!
!   Invert a sparse matrix.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine lapack_sparse_inv (N,DimA,A,DimAinv,Ainv)
  use lib_sparse

! I/O Variables.
  integer,     intent(in) :: N         ! Size of the matrix.
  integer,     intent(in) :: DimA      ! Current storage dimension.
  type(sparse),intent(in) :: A(:)      ! Sparse matrix.
  integer,     intent(out):: DimAinv   ! Current storage dimension.
  type(sparse),intent(out):: Ainv(:)   ! Sparse matrix.

! Local variables.
  integer:: Info,k
  integer:: KL,KU                    ! Upper and lower bandwidth of matrix.
  real(8),allocatable:: Aband(:,:)   ! Matrix in banded form.
  real(8),allocatable:: X(:,:)       ! Forcing vector.
  integer,allocatable:: LUpoint(:)   ! Pointer in the LU decomposition.

! Obtain banded matrix.
  call sparse_getbandwidth (DimA,A,KL,KU)
  allocate (LUPoint(N));  LUPoint=0
  allocate(Aband(2*KL+KU+1,N)); Aband  = 0.d0
  call sparse_copy2band (DimA,A,KL,KU,Aband)

! LU factorization of banded matrix.
  call lapack_band_lufact (N,KL,KU,Aband,LUPoint,Info)
  if (Info.ne.0) STOP 'Error in the LU decomposition'

! Backsubstitution with Unit vectors.
  allocate (X(N,1))
  do k=1,N
    X=0; X(k,1)=1.d0
    call lapack_band_luback (N,KL,KU,Aband,LUPoint,X(:,1),Info)
    call sparse_addmat (0,k-1,X,DimAinv,Ainv)
    if (Info.ne.0) STOP 'Error in the LU backsubstitution'
  end do

  deallocate (LUPoint,Aband,X)
  return
 end subroutine lapack_sparse_inv




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!->Subroutine LAPACK_NONSYMEIGV
!
!->Description.-
!
!   Computes the complex eigenvalues of the generalized problem AX= lambda*B*X.
!
!->Remarks.-
!
!  1) LAPACK routines used:
!
!      -DGGEV:
!          Computes for a pair of N-by-N complex nonsymmetric matrices
!          (A,B), the generalized eigenvalues, and optionally, the left and/or
!          right generalized eigenvectors.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine lapack_nonsymeigv (N,LenA,A,LenB,B,NumLambda,Lambda,Vectors,Info)
  use lib_sparse

  integer,     intent(in) :: N           ! Size of the matrix.
  integer,     intent(in) :: LenA        ! Size of the storage array for A matrix.
  type(sparse),intent(in) :: A (:)       ! Matrix A.
  integer,     intent(in) :: LenB        ! Size of the storage array for B matrix.
  type(sparse),intent(in) :: B (:)       ! Matrix B, definite positive.
  integer,     intent(in) :: NumLambda   ! Number of output eigenvalues.
  complex(8),  intent(out):: Lambda(:)   ! Complex Eigenvalues.
  complex(8),  intent(out):: Vectors(:,:)! Eigenvectors.
  integer,     intent(out):: Info        ! Error codes.

! Local variables.

  real(8),allocatable:: AlphaI(:) ! Eigenvalues are Alpha/Beta
  real(8),allocatable:: AlphaR(:) ! Eigenvalues are Alpha/Beta
  real(8),allocatable:: Afull(:,:)! A matrix in NxN form.
  real(8),allocatable:: Bfull(:,:)! B matrix in NxN form.
  real(8),allocatable:: Beta (:)
  real(8)            :: DummyL
  logical            :: Flag
  integer            :: i,j,k     ! Counters.
  complex(8)         :: iu        ! Imaginary Unit.
  character(len=1)   :: JobVL     ! = 'N':  do not compute the left generalized eigenvectors;
                                  ! = 'V':  compute the left generalized eigenvectors.
  character(len=1)   :: JobVR     ! = 'N':  do not compute the right generalized eigenvectors;
                                  ! = 'V':  compute the right generalized eigenvectors.
  integer,allocatable:: List(:)   ! Ordered list of eigenvalues.
  integer            :: LWork     ! The dimension of the array WORK.
  real(8)            :: MaxIm     ! Maximum of imaginary values.
  real(8),allocatable:: VFull(:,:)! All eigenvectors.
  real(8),allocatable:: Work (:)  ! Complex workspace.

! Initialize.

  allocate(Afull(N,N)); Afull=0.d0
  allocate(Bfull(N,N)); Bfull=0.d0
  allocate(Vfull(N,N)); Afull=0.d0

  allocate(AlphaI(N)); AlphaI=0.d0
  allocate(AlphaR(N)); AlphaR=0.d0
  allocate(Beta(N));   Beta  =0.d0

  LWork=16*N
  allocate(List(NumLambda)); List  =0
  allocate(Work(LWork))

  iu= cdsqrt(dcmplx(-1.D0))

! Convert input matrix to full form.

  do k=1,LenA
    Afull(A(k)%i,A(k)%j)= A(k)%a
  end do
  do k=1,LenB
    Bfull(B(k)%i,B(k)%j)= B(k)%a
  end do

! Call LAPACK Routine.

  JobVL='N'
  JobVR='V'

  call DGGEV (JobVL,JobVR,N,Afull,N,Bfull,N,AlphaR,AlphaI,Beta, &
&             DummyL,1,Vfull,N,Work,LWork,Info)

  if (Info.ne.0) return

! Identify the larger M eigenvalues.

  Vectors=(0.d0,0.d0)
  do i=1,NumLambda
    MaxIm=0.d0
    do j=1,N

      Flag=.true.
      do k=1,i-1
        if (List(k).eq.j) Flag=.false.
      end do

      if (Flag) then
        if (abs(AlphaI(j)/ Beta(j)).gt.MaxIm) then
          MaxIm=abs(AlphaI(j)/Beta(j))
          List(i)=j
        end if
      end if
    end do

! Store eigenvalues.

    Lambda(i)=   dcmplx(AlphaR(List(i))/Beta(List(i))) + &
&             iu*dcmplx(AlphaI(List(i))/Beta(List(i)))
  end do

! Store eigenvectors.

  i=1
  do while (i.le.NumLambda)
    Vectors(:,i)  =dcmplx(Vfull(:,List(i)))+iu*dcmplx(Vfull(:,List(i+1)))
    Vectors(:,i+1)=dcmplx(Vfull(:,List(i)))-iu*dcmplx(Vfull(:,List(i+1)))
    i=i+2
  end do
  
  deallocate(Afull,Bfull,Vfull,AlphaI,AlphaR,Beta,List,Work)
  return
  end subroutine lapack_nonsymeigv

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end module interface_lapack
 
