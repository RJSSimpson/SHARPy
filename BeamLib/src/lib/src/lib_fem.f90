!->Module LIB_FEM. Rafa Palacios. Last Update 15Jul2008.
!
!->Description.-
!
!    This module includes specific subroutines for Finite Element solution
!    with banded matrices.
!
!->Subroutines.-
!
!  fem_bandmat:         Multiply a banded matrix by a full matrix.
!  fem_bandvec:         Multiply a banded matrix by a vector.
!  fem_choleski_dcmp:   Choleski decomposition of a symmetric banded matrix.
!  fem_choleski_solv:   Solution of linear system with Choleski decomposition.
!  fem_1d_gauss_val:    Database of 1-D Gauss point coordinates and weights.
!  fem_2d_gauss_num:    Determines the number of gaussian Points in a 2-D element.
!  fem_2d_gauss_val:    Database of 2-D Gauss point coordinates and weights.
!  fem_1d_shapefun:     1-D isoparametric shape functions.
!  fem_2d_shapefun:     2-D isoparametric Shape Functions.
!  fem_elmdof:          Get index of element d.o.f. in the model.
!  fem_glob2loc_map:    Get dof data from global to local matrix.
!  fem_glob2loc_extract: Get data from global vector into element grids.
!  fem_v2m:             Convert vector to matrix.
!  fem_m2v:             Convert matrix to vector array
!
!->Remarks.-
!
!  1) Compact storage of (anti)symmetric banded matrices:
!
!       Actual matrix:
!
!        -------------------------
!        |a11  a12  a13  0   0   |
!        |a12  a22  a23  a24 0   |
!        |a13  a23  a33  a34 a35 |
!        |0    a24  a34  a44 a45 |
!        |0    0    a35  a45 a55 |
!        -------------------------
!
!      Stored matrix:
!
!        ----------------
!        |a11  a12  a13 |
!        |a22  a23  a24 |
!        |a33  a34  a35 |
!        |a44  a45  0   |
!        |a55  0    0   |
!        ----------------
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module lib_fem
 implicit none
!
! (This module does not contain Public Variables).
!

!-> Private variables.
!
! This parameter to simulate the machine epsilon for 
! a very small number which will cause singularity
 real(8), private, parameter:: Epsilon = 1.0d-13


 contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!-> Function FEM_BANDMAT
!
!-> Description.-
!
!   Multiply a banded matrix by a full matrix.
!
!->Remarks:
!
!   1) The banded matrix is defined by its upper triangular submatrix.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 function fem_bandmat (Symmetry,A,B)
!
!-> I/O Variables.
!
  character(len=1),intent(in) :: Symmetry               ! Symmetry properties.
  real(8),intent(in) :: A(:,:)                          ! Banded Matrix, A.
  real(8),intent(in) :: B(:,:)                          ! Regular Matrix,B.
  real(8) :: fem_bandmat(size(A,DIM=1), size(B,DIM=2))  ! A�B
!
!-> Local Variables.
!
  integer :: iCol       ! Counters on columns of matrix B.
!
! Multiply each column of B using the vector multiplication.
!
  do iCol=1, size(B,DIM=2)
    fem_bandmat(:,iCol)= fem_bandvec (Symmetry,A,B(:,iCol))
  end do
!
  return
 end function fem_bandmat


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!-> Function FEM_BANDVEC
!
!-> Description.-
!
!   Multiply a banded matrix by a vector.
!
!->Remarks:
!
!   1) The banded matrix is defined by its upper triangular submatrix.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 function fem_bandvec (Symmetry,A,b)
!
!-> I/O Variables.
!
  character(len=1), intent(in) :: Symmetry    ! =S or A: sym or ant matrix.
  real(8),intent(in) :: A(:,:)                ! Banded Matrix, A.
  real(8),intent(in) :: b(:)                  ! Regular vector,b.
  real(8)            :: fem_bandvec(size(b))  ! A�b
!
!-> Local Variables.
!
  real(8) :: Delta      ! +1 or -1 depending of symmetry properties.
  integer :: DimBand    ! Length of the band for the current row.
  integer :: iRow       ! Counter on rows.
  integer :: iCol       ! Counters on columns.
  integer :: NumRows    ! Number of rows.
  integer :: NumCols    ! Number of columns in the matrix.
!
! Initialization of variables.
!
  fem_bandvec=0.d0  
  if (Symmetry.eq.'S') Delta= 1.d0
  if (Symmetry.eq.'A') Delta=-1.d0
  NumRows= size(b)
  NumCols= size(A,DIM=2)
!
! Multiply the upper triangular matrix by the vector.
!
  DimBand=NumCols
  do iRow=1,NumRows
    if (iRow.gt.NumRows-NumCols+1) DimBand=DimBand-1
    do iCol=1,DimBand
      fem_bandvec(iRow)= fem_bandvec(iRow) + A(iRow,iCol)*b(iRow+iCol-1)
    end do
  end do
!
! Multiply the lower triangular matrix by the vector.
!
  DimBand=NumCols
  do iRow=1,NumRows-1
    if (iRow.gt.NumRows-NumCols+1) DimBand=DimBand-1
    do iCol=2,DimBand
      fem_bandvec(iRow+iCol-1)= fem_bandvec(iRow+iCol-1) + Delta*A(iRow,iCol)*b(iRow)
    end do
  end do
!
  return
 end function fem_bandvec

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!->Subroutine FEM_CHOLESKI_DCMP.
!
!->Description.-
!
!    Cholesky decomposition for a banded symmetric positive definite 
!    coefficient matrix.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine fem_choleski_dcmp (A,Error)
!
!-> I/O Variables.
!
  real(8), intent(inout)::A (:,:)  ! Input and output matrices.
  integer, intent(out) :: Error    ! Error code.
! 
!-> Local variables.
!
  integer :: i, j, k    ! Counters.
  integer :: NumBand    ! Half-Bandwidth of the problem.
  real(8) :: Sum        ! Partial sum of coefficient rates.
!
! Initialize variables.
!
  Error=0
  NumBand= size(A,DIM=2)  
!
! Loop in the coefficients of the matrix.
!
  do i=1,size(A,DIM=1)
    do j=1,NumBand-1
!
! For each element, get the corresponding numerator in the decomposition
! and store in the original matrix..
!
      Sum=0.0d0
      do k= max(i+j - NumBand,1) , i-1
        Sum= Sum + A(k,1) * A(k,i+1-k) * A(k,i+j-k)
      enddo  
      A(i,j)= A(i,j) - Sum
    enddo
!
! Update the element value in the decomposition and give an error if
! the matrix is ill-conditioned.
! 
    if(abs(A(i,1)).le.Epsilon) then
      Error=-1
      return
    else
      do j = 2, NumBand
        A(i,j)= A(i,j) / A(i,1)
      enddo
    end if
  enddo
!
  return
 end subroutine fem_choleski_dcmp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!->Subroutine FEM_CHOLESKI_SOLV
!
!->Description.-
!
!    Solve a linear system Ax=b, where A is given in its Choleski
!    factorization. The solution is performed through forward and
!    back substitution.
!
!->Remarks.-
!
!   1) Algorithms can be found in Belengundu & Chandrupatla's book
!      "Finite Elements in Engineering".
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine fem_choleski_solv (A,b)
!
!-> I/O Variables.
!
  real(8), intent(in)::     A(:,:)  ! Matrix in Choleski form.
  real(8), intent(inout) :: b(:)    ! Right Hand Side.
!
!-> Local variables.
!  
  integer :: i, j      ! Counters on rows and columns.
  integer :: NumCols   ! Half-Bandwidth of the matrix.
  integer :: NumRows   ! Number of equations in the system.
  real(8) :: Sum       ! Value to store sums.

  NumRows= size(A,DIM=1)
  NumCols= size(A,DIM=2)

! Forward Substitution.
  do i=2,NumRows
    Sum= 0.0d0
    do j= max(1,i+1-NumCols), i-1 
      Sum= Sum + A(j,i-j+1) * b(j)
    enddo
    b(i)= b(i)-Sum
  enddo

! Backwards Substitution.
  b(NumRows) = b(NumRows)/A(NumRows,1)

  do i = NumRows-1, 1, -1
    Sum = 0.0d0
    do j= 2, min(NumCols, NumRows - i + 1)
      Sum = Sum + A(i,j) * b(j + i - 1)
    enddo
    b(i) = b(i)/A(i,1)- Sum
  enddo

  return
 end subroutine fem_choleski_solv


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!-> Subroutine FEM_1D_GAUSS_VAL
!
!-> Description.-
!
!   Database of Weights and Coordinates for Gaussian Quadrature for 1-D
!   integration.
!
!->Remarks:
!
!   1) Coefficients from Abramowitz and Stegun.
! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine fem_1d_gauss_val (NumGauss,Coords,Weight)

!-> I/O Variables.

  integer,intent(in) :: NumGauss        ! Number of Gaussian points.
  real(8),intent(out):: Coords(:)       ! Local Coordinates of the Points.
  real(8),intent(out):: Weight(:)       ! Weights of the Gaussian Points.

! Default output values.
  Coords=0.d0
  Weight=0.d0

! Database of weights and location of gaussian points.
  select case(NumGauss)
  case(1)
    Coords(1)= 0.d0
    Weight(1)= 2.d0
  case(2)
    Coords(1)=-0.577350269189626d0
    Coords(2)= 0.577350269189626d0
    Weight(1)= 1.d0
    Weight(2)= 1.d0
  case(3)
    Coords(1)=-0.774596669241483d0
    Coords(2)= 0.d0
    Coords(3)= 0.774596669241483d0
    Weight(1)= 0.888888888888889d0
    Weight(2)= 0.555555555555556d0
    Weight(3)= 0.888888888888889d0
  case(4)
    Coords(1)=-0.861136311594053d0
    Coords(2)=-0.339981043584856d0
    Coords(3)= 0.339981043584856d0
    Coords(4)= 0.861136311594053d0
    Weight(1)= 0.347854845137454d0
    Weight(2)= 0.652145154862546d0
    Weight(3)= 0.652145154862546d0
    Weight(4)= 0.347854845137454d0
  case(5)
    Coords(1)=-0.906179845938664d0
    Coords(2)=-0.538469310105683d0
    Coords(3)= 0.d0
    Coords(4)= 0.538469310105683d0
    Coords(5)= 0.906179845938664d0
    Weight(1)= 0.236926885056189d0
    Weight(2)= 0.478628670499366d0
    Weight(3)= 0.568888888888889d0
    Weight(4)= 0.478628670499366d0
    Weight(5)= 0.236926885056189d0
  case(6)
    Coords(1)=-0.932469514203152d0
    Coords(2)=-0.661209386466265d0
    Coords(3)=-0.238619186083197d0
    Coords(4)= 0.238619186083197d0
    Coords(5)= 0.661209386466265d0
    Coords(6)= 0.932469514203152d0
    Weight(1)= 0.171324492379170d0
    Weight(2)= 0.360761573048139d0
    Weight(3)= 0.467913934572691d0
    Weight(4)= 0.467913934572691d0
    Weight(5)= 0.360761573048139d0
    Weight(6)= 0.171324492379170d0

  case default
    STOP 'Error: Number of Gauss points in 1-D elements larger than 6.'
  end select

  return
 end subroutine fem_1d_gauss_val



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!-> Subroutine FEM_1D_SHAPEFUN
!
!-> Description.-
!
!   Calculate the value of 1-D isoparametric shape function and its
!   derivative at a inner point in the element defined by local coordinates.
!
!->Remarks:
!
!   1) The isoparametric shape functions were extracted from Chandrupatla and
!      Belengundu's book "Finite Elements in Engineering". All conventions in
!      the book are used here.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine fem_1d_shapefun (NumNodeElem, z, Fun, Der)
!
!-> I/O Variables.
  integer,intent(in) :: NumNodeElem   ! Number of nodes in the element.
  real(8),intent(in) :: z             ! Local Coordinate of the point.
  real(8),intent(out):: Fun(:)        ! Shape Function.
  real(8),intent(out):: Der(:)        ! Derivative of the Shape Function.

! Initialization.
  Fun=0.d0
  Der=0.d0

  select case (NumNodeElem)

! Three-noded beam:
  case (3)
    Fun(1)=-0.5d0*z*(1.d0-z)
    Fun(2)= 0.5d0*z*(1.d0+z)
    Fun(3)= 1.0d0-z*z
    Der(1)=-0.5d0*(1.d0-2.d0*z)
    Der(2)= 0.5d0*(1.d0+2.d0*z)
    Der(3)=-2.0d0*z

! Two-noded beam:
  case (2)
    Fun(1)= 0.5d0*(1.d0-z)
    Fun(2)= 0.5d0*(1.d0+z)
    Der(1)=-0.5d0
    Der(2)= 0.5d0

  case default
    STOP 'Error: 1-D elements must have either 2 or 3 nodes'
  end select

  return
 end subroutine fem_1d_shapefun



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!-> Subroutine FEM_2D_GAUSS_NUM
!
!-> Description.-
!
!   Determine the number of gauss points for the current type of element.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine fem_2d_gauss_num (Flag, NumNode_Elem, NumGauss)
!
!-> I/O Variables
!
  integer,intent(in) :: NumNode_Elem  ! # nodes in the element.
  integer,intent(out):: NumGauss      ! # Gaussian points.
  integer,intent(in) :: Flag          ! =0 for CTRIA, /=0 for CQUAD.
!
! Set the number of Gauss points depending of the type of element (CTRIA or
! CQUAD) and the number of nodes of the element.
!
  if (Flag.eq.0) then
    if (NumNode_Elem.gt.3) then 
      NumGauss=7
    else
      NumGauss=3
    end if
  else
    if (NumNode_Elem.gt.4) then 
      NumGauss=9
    else
      NumGauss=4
    end if
  endif
!
  return
 end subroutine fem_2d_gauss_num


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!-> Subroutine FEM_2D_GAUSS_VAL
!
!-> Description.-
!
!   Database of Weights and Coordinates for Gaussian Quadrature.
!
!->Remarks:
!
!   1) Coefficients are extracted from Zienkiewitz' book "The Finite
!      Element Method".
! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine fem_2d_gauss_val (Flag, NumGauss, Coords, Weight, iErr)
!
!-> I/O Variables.
!
  integer,intent(in) :: Flag              ! =0 for CTRIA, /=0 for CQUAD.
  integer,intent(in) :: NumGauss          ! Number of Gaussian points.
  real(8),intent(out):: Coords(:,:)       ! Local Coordinates of the Points.
  real(8),intent(out):: Weight(:)         ! Weights of the Gaussian Points.
  integer,intent(out):: iErr              ! Error Code.
!
!-> Local variables.
!
  integer, parameter:: MaxPoints=10 ! Maximum # of linear points.
  integer :: iGauss                 ! Counter in the Gaussian points.
  integer :: i1,i2                  ! Counters in the linear points.
  integer :: NumPoints              ! Number of linear points.
  real(8) :: CoordsLin(MaxPoints)   ! Local Coordinates for [-1,1] integration.
  real(8) :: WeightLin(MaxPoints)   ! Local Coordinates for [-1,1] integration.

! Default output values.
  iErr=0
  Coords=0.d0
  Weight=0.d0

! Quadrilateral Elements.
  if (Flag.ne.0) then
    NumPoints=nint(sqrt(float(NumGauss)))
    if ((NumGauss-NumPoints*NumPoints).ne.0) NumPoints=0

    select case(NumPoints)
    case(1)
      CoordsLin(1)=0.d0
      WeightLin(1)=0.d0
    case(2)
      CoordsLin(1)=-1.d0/dsqrt(3.d0)
      CoordsLin(2)= 1.d0/dsqrt(3.d0)
      WeightLin(1)= 1.d0
      WeightLin(2)= 1.d0
    case(3)
      CoordsLin(1)=-dsqrt(.6d0)
      CoordsLin(2)= 0.d0
      CoordsLin(3)= dsqrt(.6d0)
      WeightLin(1)= 5.d0/9.d0
      WeightLin(2)= 8.d0/9.d0
      WeightLin(3)= 5.d0/9.d0
    case(4)
      CoordsLin(1)=-0.861136311594053d0
      CoordsLin(2)=-0.339981043584856d0
      CoordsLin(3)= 0.339981043584856d0
      CoordsLin(4)= 0.861136311594053d0
      WeightLin(1)= 0.347854845137454d0
      WeightLin(2)= 0.652145154862546d0
      WeightLin(3)= 0.652145154862546d0
      WeightLin(4)= 0.347854845137454d0
    case(5)
      CoordsLin(1)=-0.906179845938664d0
      CoordsLin(2)=-0.538469310105683d0
      CoordsLin(3)= 0.d0
      CoordsLin(4)= 0.538469310105683d0
      CoordsLin(5)= 0.906179845938664d0
      WeightLin(1)= 0.236926885056189d0
      WeightLin(2)= 0.478628670499366d0
      WeightLin(3)= 0.568888888888889d0
      WeightLin(4)= 0.478628670499366d0
      WeightLin(5)= 0.236926885056189d0
    case default
      iErr=-1
    end select

    iGauss=0
    do i1=1,NumPoints
      do i2=1,NumPoints
        iGauss= iGauss+1
        Coords(iGauss,1)= CoordsLin(i1)
        Coords(iGauss,2)= CoordsLin(i2)
        Weight(iGauss)  = WeightLin(i1)*WeightLin(i2)
      end do
    end do

! Triangular elements.
  else if (Flag.eq.0) then
    NumPoints=NumGauss
    select case(NumPoints)
    case(3)
      Coords(1,1)=1.d0/6.d0
      Coords(1,2)=1.d0/6.d0
      Coords(2,1)=2.d0/3.d0
      Coords(2,2)=1.d0/6.d0
      Coords(3,1)=1.d0/6.d0
      Coords(3,2)=2.d0/3.d0
      Weight(1:3)=1.d0/6.d0
    case(7)
      Coords(1,1)= 0.101286507323456d0
      Coords(1,2)= 0.101286507323456d0
      Weight(1)  = 0.0629695902724d0
      Coords(2,1)= 0.797426985353087d0
      Coords(2,2)= 0.101286507323456d0
      Weight(2)  = 0.0629695902724d0
      Coords(3,1)= 0.101286507323456d0
      Coords(3,2)= 0.797426985353087d0
      Weight(3)  = 0.0629695902724d0
      Coords(4,1)= 0.470142064105115d0
      Coords(4,2)= 0.059715871789770d0
      Weight(4)  = 0.0661970763943d0
      Coords(5,1)= 0.470142064105115d0
      Coords(5,2)= 0.470142064105115d0
      Weight(5)  = 0.0661970763943d0
      Coords(6,1)= 0.059715871789770d0
      Coords(6,2)= 0.470142064105115d0
      Weight(6)  = 0.0661970763943d0
      Coords(7,1)= 0.333333333333333d0
      Coords(7,2)= 0.333333333333333d0
      Weight(7)  = 0.1125d0
    case default
      iErr=-1
    end select
  end if

  return
 end subroutine fem_2d_gauss_val


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!-> Subroutine FEM_2D_SHAPEFUN
!
!-> Description.-
!
!   Calculate the value of the 2-D isoparametric shape function and its
!   derivatives at a inner point in the element defined by local coordinates.
!
!->Remarks:
!
!   1) The isoparametric shape functions were extracted from Chandrupatla and
!      Belengundu's book "Finite Elements in Engineering". All conventions in
!      the book are used here.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine fem_2d_shapefun (Flag, Node, z1, z2, Fun, Der)
!
!-> I/O Variables.
!
  integer,intent(in) :: Flag                   ! =0: for CTRIA, /=0: for CQUAD.
  integer,intent(in) :: Node(:)                ! Array on nodes (=0: node not present).
  real(8),intent(in) :: z1,z2         ! Local Coordinates of the Points.
  real(8),intent(out):: Fun(:)        ! Shape Function.
  real(8),intent(out):: Der(:,:)      ! Derivative of the Shape Function.
!
!-> Local variables.
!
  real(8) :: z3                       ! Third variable for triangular elements.
!
! Initialization.
!
  Fun=0.d0
  Der=0.d0
!
!
!-> Quadrilateral Elements:
!   ----------------------
!
  if (Flag.ne.0) then
!
! Shape Functions.
!
    if (Node(9).ne.0) Fun(9)=        (1.d0-z1*z1)*(1.d0-z2*z2)
    if (Node(5).ne.0) Fun(5)=.5d0 * ((1.d0-z1*z1)*(1.d0-z2)-Fun(9))
    if (Node(7).ne.0) Fun(7)=.5d0 * ((1.d0-z1*z1)*(1.d0+z2)-Fun(9))
    if (Node(6).ne.0) Fun(6)=.5d0 * ((1.d0-z2*z2)*(1.d0+z1)-Fun(9))
    if (Node(8).ne.0) Fun(8)=.5d0 * ((1.d0-z2*z2)*(1.d0-z1)-Fun(9))
!
    Fun(1)=.25d0 * ((1.d0-z1)*(1.d0-z2)-Fun(9)) - .5d0*(Fun(8)+Fun(5))
    Fun(2)=.25d0 * ((1.d0+z1)*(1.d0-z2)-Fun(9)) - .5d0*(Fun(5)+Fun(6))
    Fun(3)=.25d0 * ((1.d0+z1)*(1.d0+z2)-Fun(9)) - .5d0*(Fun(6)+Fun(7))
    Fun(4)=.25d0 * ((1.d0-z1)*(1.d0+z2)-Fun(9)) - .5d0*(Fun(7)+Fun(8))
!
! Derivatives with z1.
!
    if (Node(9).ne.0) Der(9,1)=-2.d0*z1*(1.d0-z2*z2)
    if (Node(5).ne.0) Der(5,1)=-z1*(1.d0-z2) -.5d0*Der(9,1)
    if (Node(7).ne.0) Der(7,1)=-z1*(1.d0+z2) -.5d0*Der(9,1)
    if (Node(6).ne.0) Der(6,1)=.5d0*( (1.d0-z2*z2)-Der(9,1))
    if (Node(8).ne.0) Der(8,1)=.5d0*(-(1.d0-z2*z2)-Der(9,1))
!
    Der(1,1)=.25d0*(-(1.d0-z2)-Der(9,1)) - .5d0*(Der(8,1)+Der(5,1))
    Der(2,1)=.25d0*( (1.d0-z2)-Der(9,1)) - .5d0*(Der(5,1)+Der(6,1))
    Der(3,1)=.25d0*( (1.d0+z2)-Der(9,1)) - .5d0*(Der(6,1)+Der(7,1))
    Der(4,1)=.25d0*(-(1.d0+z2)-Der(9,1)) - .5d0*(Der(7,1)+Der(8,1))
!
! Derivatives with z2.
!      
    if (Node(9).ne.0) Der(9,2)=-2.d0*z2*(1.d0-z1*z1)
    if (Node(5).ne.0) Der(5,2)=.5d0*(-(1.d0-z1*z1)-Der(9,2))
    if (Node(7).ne.0) Der(7,2)=.5d0*( (1.d0-z1*z1)-Der(9,2))
    if (Node(6).ne.0) Der(6,2)=-z2*(1.d0+z1) -.5d0*Der(9,2)
    if (Node(8).ne.0) Der(8,2)=-z2*(1.d0-z1) -.5d0*Der(9,2)
!
    Der(1,2)=.25d0*(-(1.d0-z1)-Der(9,2)) - .5d0*(Der(8,2)+Der(5,2))
    Der(2,2)=.25d0*(-(1.d0+z1)-Der(9,2)) - .5d0*(Der(5,2)+Der(6,2))
    Der(3,2)=.25d0*( (1.d0+z1)-Der(9,2)) - .5d0*(Der(6,2)+Der(7,2))
    Der(4,2)=.25d0*( (1.d0-z1)-Der(9,2)) - .5d0*(Der(7,2)+Der(8,2))
!
!
!-> Triangular Elements:
!   -------------------
!
  else
    z3=1.d0-z1-z2
!
! Shape Functions.
!
    if (Node(5).ne.0) Fun(5)= 4.d0*z3*z1
    if (Node(6).ne.0) Fun(6)= 4.d0*z1*z2
    if (Node(7).ne.0) Fun(7)= 4.d0*z2*z3
!
    Fun(1)= z3 - .5d0*(Fun(5)+Fun(7))
    Fun(2)= z1 - .5d0*(Fun(6)+Fun(5))
    Fun(3)= z2 - .5d0*(Fun(7)+Fun(6))
!
! Derivatives with z1.
!
    if (Node(5).ne.0) Der(5,1)= 4.d0*(z3-z1)
    if (Node(6).ne.0) Der(6,1)= 4.d0*z2
    if (Node(7).ne.0) Der(7,1)=-4.d0*z2
!
    Der(1,1)=-1.d0 - .5d0*(Der(5,1)+Der(7,1))
    Der(2,1)= 1.d0 - .5d0*(Der(6,1)+Der(5,1))
    Der(3,1)=      - .5d0*(Der(7,1)+Der(6,1))
!
! Derivatives with z2.
!
    if (Node(5).ne.0) Der(5,2)=-4.d0*z1
    if (Node(6).ne.0) Der(6,2)= 4.d0*z1
    if (Node(7).ne.0) Der(7,2)= 4.d0*(z3-z2)
!
    Der(1,2)=-1.d0 - .5d0*(Der(5,2)+Der(7,2))
    Der(2,2)=      - .5d0*(Der(6,2)+Der(5,2))
    Der(3,2)= 1.d0 - .5d0*(Der(7,2)+Der(6,2))
!
  end if
!
  return
 end subroutine fem_2d_shapefun

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!->Subroutine FEM_ELMDOF.
!
!->Description.-
!
!   Find the corresponding d.o.f's for an element and store them in
!   a steering vector.  
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine fem_elmdof(NumDofNode,NodeList,g)
!
!-> I/O Variables.
!
  integer,intent(in) :: NumDofNode    ! Number of d.o.f.'s in each node.
  integer,intent(in) :: NodeList(:)   ! List of nodes in the element.
  integer,intent(out):: g(:)          ! Global d.o.f's in the element.
!
!-> Local variables.
!
  integer::i,j,m
!
! Initialize.
  g= 0;
  m= 0;  
!
! Loop in the nodes.
  do i=1,size(NodeList)
    if(NodeList(i).ne.0) then
      do j=1,NumDofNode
        m= m+1
        g(m)= NumDofNode*(NodeList(i)-1)+j
      end do
   end if
  end do
!
  return
 end subroutine fem_elmdof

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!->Subroutine FEM_GLOB2LOC_MAP.
!
!->Description.-
!
!    Extract the elements corresponding to a given dof from a global
!    matrix and add them to a given local matrix.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine fem_glob2loc_map (global,local,g)

  ! I/O Variables.
  real(8),intent(in) :: global(:,:)
  real(8),intent(out) :: local(:,:)
  integer, intent(in) :: g(:)
!
!-> Local Variables.
!
  integer:: i
!
  do i=1,size(local,1) 
    local(i,:)=local(i,:)+global(g(i),:)
  end do
!
  return
 end subroutine fem_glob2loc_map


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!->Subroutine FEM_GLOB2LOC_EXTRACT.
!
!->Description.-
!
!    Extract the values corresponding to given nodes from a vector of global
!    data.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine fem_glob2loc_extract (ElemNodes,GlobVector,LocVector,NumNodeElem)

! I/O variables.
  integer,intent(in) :: ElemNodes (:)      ! Nodes in the current element.
  real(8),intent(in) :: GlobVector(:,:)    ! Vector of global values.
  real(8),intent(out):: LocVector (:,:)    ! Vector of local values.
  integer,intent(out):: NumNodeElem        ! Actual number of nodes in the element.

! Local variables.
  integer:: i

! Compute number of non-zero connectivies in the element.
  NumNodeElem=0
  do i=1,size(ElemNodes)
    if (ElemNodes(i).gt.0) NumNodeElem=NumNodeElem+1
  end do

! Extract data from global vector.
  LocVector=0.d0
  LocVector(1:NumNodeElem,:)=GlobVector(ElemNodes(1:NumNodeElem),:)

  return
 end subroutine fem_glob2loc_extract


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!->Subroutine FEM_V2M.
!
!->Description.-
!
!    Convert array of degrees of freedom to matrix form.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 function fem_v2m (Vector,Nrow,Ncol)

! I/O variables.
  real(8),intent(in) :: Vector(:)          ! Vector.
  integer,intent(in) :: Nrow               ! Number of rows in the matrix.
  integer,intent(in) :: Ncol               ! Number of columns in the matrix.
  real(8)            :: fem_v2m(Nrow,Ncol) ! Matrix.

! Local variables.
  integer:: i,iCol,iRow    ! Counters.

  iRow=1
  iCol=1
  do i=1,size(Vector)

    fem_v2m(iRow,iCol)=Vector(i)
    
    if (iCol.eq.NCol) then
      iRow=iRow+1
      iCol=1
    else
      iCol=iCol+1
    end if
  end do

  return
 end function fem_v2m


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!->Subroutine FEM_M2V.
!
!->Description.-
!
!    Convert matrix to vector array.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 function fem_m2v (Matrix,N,Filter)

! I/O variables.
  real(8),intent(in)         :: Matrix(:,:)      ! Matrix.
  integer,intent(in)         :: N                ! Size of the vector.
  integer,optional,intent(in):: Filter(:)        ! Filter on rows of the matrix.
  real(8)                    :: fem_m2v(N)       ! Vector.

! Local variables.
  integer:: i,j,k      ! Counters.

! Initialize
  k=0
  fem_m2v=0.d0
  
  do i=1,size(Matrix,DIM=1)
    do j=1,size(Matrix,DIM=2)

      if (present(Filter)) then
        if (Filter(i).ne.0) then
          k=k+1
          fem_m2v(k)=Matrix(i,j)
        end if

      else
        k=k+1
        fem_m2v(k)=Matrix(i,j)
      end if
    end do
  end do

  return
 end function fem_m2v


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!->Subroutine FEM_N2D.
!
!->Description.-
!
!    Convert node array to degrees of freedom.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 function fem_n2d (NumNodes,NumDofNode,NodeList)

! I/O variables.
   integer,intent(in):: NumNodes        ! Number of degrees of freedom in the model.
   integer,intent(in):: NumDofNode      ! Number of degrees of freedom on each node.
   integer,intent(in):: NodeList(:)     ! Nodes in the global array.
   integer           :: fem_n2d(NumNodes*NumDofNode) ! Vector.

! Local variables.
  integer:: i,j,k      ! Counters.

! Initialize
  k=0
  fem_n2d=0

  do i=1,size(NodeList)
    do j=1,NumDofNode
      k=k+1
      if (NodeList(i).ne.0) fem_n2d(k)= NumDofNode*(NodeList(i)-1)+j
    end do
  end do

  return
 end function fem_n2d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end module lib_fem