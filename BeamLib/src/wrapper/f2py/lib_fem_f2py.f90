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
  real(8),intent(out):: Fun(NumNodeElem)        ! Shape Function.
  real(8),intent(out):: Der(NumNodeElem)        ! Derivative of the Shape Function.

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
