!->Copyright by The University of Michigan, Aerospace Department. 2003
!
!->Module LIB_LU. Rafa Palacios. 22Jun2003
!
!->Description.-
!
!  Matrix operations using LU decomposition.
!
!->Subroutines.-
!
!   lu_bksubs:    LU back-substitution.
!   lu_decomp:    LU decomposition of a general square matrix.
!   lu_determ:    Determinant of a square matrix.
!   lu_invers:    Inversion of a square matrix.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module lib_lu
 implicit none
!
!   (There are no public variables in this module).
!
 contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!->Subroutine LU_BKSUBS.
!
!->Description.-
!
!   Solves the set of N linear equations A�X=B. A is given by its LU decomposition.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine lu_bksubs (A, INDX, B)
!
!-> I/O Variables.
!
  real(8), intent(in)   :: A(:,:)
  integer, intent(in)   :: INDX(:)
  real(8), intent(inout):: B(:)
!
!-> Local Variables.
!
  integer :: N               ! Size of the matrices.
  integer :: I, II, LL, J
  real(8) :: Sum
!
      n= size(A,1)
!
      II = 0
      DO I = 1, N
        LL = INDX(I)
        SUM = B(LL)
        B(LL) = B(I)
        IF(II.NE.0)THEN
          DO J = II, I-1
            SUM = SUM - A(I, J)*B(J)
          END DO
        ELSEIF(SUM.NE.0)THEN
          II = I
        ENDIF
        B(I) = SUM
      END DO
!
      DO I = N, 1, -1
        SUM = B(I)
        IF(I.LT.N)THEN
          DO J = I+1, N
            SUM = SUM - A(I, J)*B(J)
          END DO
        ENDIF
        B(I) = SUM/A(I, I)
      END DO
!
  return
 end subroutine lu_bksubs

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!->Subroutine LU_DECOMP.
!
!->Description.-
!
!   LU Decomposition of a NxN Matrix.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine lu_decomp (A, INDX, D)
!
!-> I/O Variables.
!
  real(8), intent(inout):: A(:,:)
  integer, intent(out)  :: INDX(:)
  real(8), intent(out)  :: D
!
!-> Local Variables.
!
  integer :: i,j,k                   ! Counters.
  integer :: imax                    ! Max Counter.
  integer :: N                       ! Size of the matrices.
  real(8) :: Aamax, Sum, Dum
  real(8), allocatable :: VV (:)     ! Stores the implicit scaling of each row.
!
  real(8), parameter :: Tiny=1.0d-20 ! Small number.
!
!
      n= size(A,1)
      allocate (VV(N))
!
      D = 1.0
      DO I = 1, N
        AAMAX = 0.0
        DO J = 1, N
          IF(ABS(A(I,J)).GT.AAMAX) AAMAX=ABS(A(I,J))
        END DO
        IF(AAMAX.EQ.0.0) STOP 'Singular matrix'
        VV(I:I) = 1.d0/AAMAX
      END DO
!
      DO J = 1, N
        DO I = 1, J-1
          SUM = A(I, J)
          DO K = 1, I-1
            SUM = SUM - A(I, K)*A(K, J)
          END DO
          A(I, J) = SUM
        END DO
!
        AAMAX = 0.0
        DO I = J, N
          SUM = A(I, J)
          DO K = 1, J-1
            SUM = SUM - A(I, K)*A(K, J)
          END DO
          A(I, J) = SUM
          DUM = VV(I)*ABS(SUM)
          IF(DUM.GE.AAMAX)THEN
            IMAX = I
            AAMAX = DUM
          ENDIF
        END DO
        IF(J.NE.IMAX)THEN
          DO K = 1, N
            DUM = A(IMAX, K)
            A(IMAX, K) = A(J, K)
            A(J, K) = DUM
          END DO
          D = -D
          VV(IMAX) = VV(J)
        ENDIF
        INDX(J) = IMAX
        IF(A(J, J).EQ.0.0) A(J, J) = TINY
        IF(J.NE.N)THEN
          DUM = 1.0/A(J, J)
          DO I = J+1, N
            A(I, J) = A(I, J)*DUM
          END DO
        ENDIF
      END DO
!
  deallocate (VV)
  return
 end subroutine lu_decomp
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!->Subroutine LU_INVERS.
!
!->Description.-
!
!   Invert a matrix using LU Decomposition. It can also return the determinant.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine lu_invers (Matrix, InvMatrix, Det)
!
!-> I/O Variable.
! 
  real(8),intent(in) :: Matrix(:,:)
  real(8),intent(out):: InvMatrix(:,:)
  real(8),optional,intent(out):: Det
!
!-> Local variables.
!
  integer::i,N
  real(8) :: D
  real(8), allocatable :: MatrixLU(:,:)
  integer, allocatable :: Indexes(:)
!
! Set dimension and the auxiliary vector.
!
  N= size(Matrix,1)
  allocate (Indexes(N))
  allocate (MatrixLU(N,N))
!
! Initialize InvMatrix as the identity matrix and MatrixLU as Matrix.
!
  MatrixLU=Matrix
!
  InvMatrix=0.d0
  do i=1,N
    InvMatrix(i,i)=1.d0
  end do
!
! LU Decomposition. 
!
  call lu_decomp (MatrixLU,Indexes,D)
!
! The inverse is the solution to units vectors.
!
  do i=1,N
    call lu_bksubs (MatrixLU,Indexes,InvMatrix(:,i))
  end do
!
! Return the determinant of the matrix.
!
  if (present(Det)) then
    Det=D
    do i=1,N
      Det=Det*MatrixLU(i,i)
    end do
  end if
!
  return
 end subroutine lu_invers
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!->Subroutine LU_DETERM.
!
!->Description.-
!
!   Invert a matrix using LU Decomposition.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine lu_determ (Matrix, Det)
!
!-> I/O Variable.
! 
  real(8), intent(in) :: Matrix(:,:)
  real(8), intent(out):: Det
!
!-> Local variables.
!
  integer::i,N
  real(8), allocatable :: MatrixLU(:,:)
  integer, allocatable :: Indexes(:)
!
! Set dimension and the auxiliary vector.
!
  N= size(Matrix,1)
  allocate (Indexes(N))
  allocate (MatrixLU(N,N)); MatrixLU=Matrix
!
! LU Decomposition. 
!
  call lu_decomp (MatrixLU,Indexes,Det)
!
! The determinant of a LU decomposed matrix is just the product of the
! diagonal elements.
!
  do i=1,N
    Det=Det*MatrixLU(i,i)
  end do
!
  return
 end subroutine lu_determ
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end module lib_lu
