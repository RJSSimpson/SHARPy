!-> Module.- XBEAM_UNDEF Rafa Palacios. 15Jul2008
!
!-> Language: FORTRAN90, Free format.
!
!-> Description.-
!
!  Compute properties of undeformed structure.
!
!-> Subroutines.-
!
!    -xbeam_undef_geom      : Compute geometric properties of elements.
!    |-xbeam_undef_relatree : Compute tree of nodes shared between elements.
!                              Identify master and slave nodes.
!    -xbeam_undef_dofs      : Identify degrees of freedom in the problem.
!    |-xbeam_undef_nodeindep: Identify independent nodes.
!
!-> Remarks.-
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module xbeam_undef
  use xbeam_shared
  implicit none

 contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-> Subroutine XBEAM_UNDEF_GEOM
!
!-> Description:
!
!    Compute basic geometric parameters on the undeformed configuration.
!
!-> Remarks.-
!
!  1) The undeformed frame (b) is defined within each element. For nodes
!     belonging to two or more elements, there will be as many undeformed
!     frames.
!
!  2) The displaced coordinate system (a) is given as the cartesian rotation
!     vector that defines its orientation with respect to the global frame.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine xbeam_undef_geom (Elem,Coords,PhiNodes,Psi0,Options)
  use lib_fem
  use lib_bgeom

! I/O Variables.
  type(xbelem),intent(inout):: Elem     (:)    ! Element information.
  real(8),      intent(in)   :: Coords   (:,:)  ! Initial coordinates of the grid points.
  real(8),      intent(in)   :: PhiNodes (:)    ! Pretwist at the nodes.
  real(8),      intent(out)  :: Psi0     (:,:,:)! Initial CRV at element nodes.
  type(xbopts),intent(in)   :: Options         ! Solver parameters.

! Local variables.
  real(8)            :: ElemPhi(MaxElNod)      ! Pretwist at the nodes of the element.
  integer            :: j                      ! Counter on the elements in the model.
  real(8)            :: LocCoords(MaxElNod,3)  ! Global coordinates of nodes in the element.
 
! Loop in all elements in the model.
  do j=1,size(Elem)

! Extract element information.
    call fem_glob2loc_extract (Elem(j)%Conn,Coords,LocCoords,Elem(j)%NumNodes)

    ElemPhi(1:Elem(j)%NumNodes)= PhiNodes(Elem(j)%Conn(1:Elem(j)%NumNodes))

! Compute initial displaced coordinate system (a).
    call bgeom_elemframe (LocCoords(1:2,:),Elem(j)%Vector,Elem(j)%Psi, &
&                         Options%DeltaCurved)

! Compute undeformed frame (b) at all nodes of each element. Given by CRV in Psi0.
    call bgeom_nodeframe (Elem(j)%NumNodes,LocCoords,ElemPhi,Elem(j)%Vector,&
&                         Psi0(j,:,:),Options%DeltaCurved)

! Compute element length.
    call bgeom_elemlength (Options%NumGauss,Elem(j)%NumNodes,LocCoords,Elem(j)%Length)

! Compute curvature in the initial configuration.
    if (Elem(j)%NumNodes.eq.2) then
      call bgeom_elemcurv2 (LocCoords(1:2,:),Psi0(j,:,:),Elem(j)%PreCurv)
    end if

  end do

! Determine tree of connectivities. For each element determine if their nodes are master or
! slaves. In the last case, determine its master node.
  call xbeam_undef_relatree (Elem)

  return
 end subroutine xbeam_undef_geom


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-> Subroutine XBEAM_UNDEF_DOFS
!
!-> Description:
!
!    Initialize the solution process.
!
!-> Remarks.-
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine xbeam_undef_dofs (Elem,BoundConds,Node,NumDof)
  use lib_fem
  use lib_bgeom

! I/O Variables.
  type(xbelem),intent(in)   :: Elem      (:)      ! Element information.
  integer,      intent(in)   :: BoundConds(:)      ! Boundary conditions.
  type(xbnode),intent(inout):: Node      (:)      ! Nodal information.
  integer,      intent(out)  :: NumDof             ! Number of independent degrees of freedom.

! Local variables.
  integer:: NumE                           ! Number of elements in the model.
  integer:: NumN                           ! Number of nodes in the model.
  integer :: j,k                           ! Counter on the elements.
  integer,allocatable:: ListFr    (:)      ! List of non-free nodes.
  integer,allocatable:: ListIN    (:)      ! List of independent nodes.

! Initialize
  NumE=size(Elem)
  NumN=size(Node)

! Determine ID of master for each node in the model.
  do j=1,NumE
    do k=1,Elem(j)%NumNodes
        if (Elem(j)%Master(k,1).eq.0) then
          Node(Elem(j)%Conn(k))%Master(1)= j
          Node(Elem(j)%Conn(k))%Master(2)= k
        end if
    end do
  end do

! Get list of independent nodes.
  allocate (ListIN(NumN)); ListIN= 0
  allocate (ListFr(NumN)); ListFr= 0
  call xbeam_undef_nodeindep (NumN,BoundConds,NumDof,ListIN,ListFr)
  NumDof=NumDof*6

  do j=1,NumN
    Node(j)%Vdof=ListIN(j)
    Node(j)%Fdof=ListFr(j)
  end do

  deallocate (ListIn,ListFr)
  return
 end subroutine xbeam_undef_dofs


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-> Subroutine XBEAM_UNDEF_RELATREE
!
!-> Description:
!
!    Compute tree of nodes shared between elements.
!
!-> Remarks.-
!
!  1) Elem(i)%Master(j,1)=0 -> Node j of element i is a master one.
!
!     Elem(i)%Master(j,1)=i0 & Elem(i)%Master(j,2)=j0 -> Node j of element i
!     is slave node of node j0 of element i0.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine xbeam_undef_relatree (Elem)

! I/O Variables.
  type(xbelem),intent(inout):: Elem (:)     ! Element information

! Local variables.
  integer :: iElem,jElem        ! Counters on the elements in the model.
  integer :: iNode,jNode        ! Counters on the nodes in the elements.

! Loop in the nodes of all elements.
  do iElem=1,size(Elem)

! Initalize to zero (i.e., by default all nodes are master nodes).
    Elem(iElem)%Master=0

! Loop in the nodes within the element
    do iNode=1,Elem(iElem)%NumNodes

! For each node, see if it has already appeared. For that purpose, scan
! all elements from 1 to iElem, and compare their connectivites to the current
! node.
      jElem=1
      do while ((Elem(iElem)%Master(iNode,1).eq.0).and.(jElem.lt.iElem))
        do jNode=1,Elem(iElem)%NumNodes
          if (Elem(jElem)%Conn(jNode) .eq. Elem(iElem)%Conn(iNode)) then
            Elem(iElem)%Master(iNode,1)=jElem
            Elem(iElem)%Master(iNode,2)=jNode
          end if
        end do
        jElem=jElem+1
      end do
    end do
  end do

  return
 end subroutine xbeam_undef_relatree


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-> Subroutine XBEAM_UNDEF_NODEINDEP
!
!-> Description:
!
!    Identify independent nodes in the problem.
!
!-> Remarks.-
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine xbeam_undef_nodeindep (NumN,BoundCond,NumIN,ListIN,ListFr)

! I/O Variables.
  integer,intent(in)::  NumN               ! Number of nodes in the model.
  integer,intent(in)::  BoundCond(:)       ! Boundary conditions: 1: Clamped; -1; free.
  integer,intent(out):: NumIN              ! Number of independent nodes.
  integer,intent(out):: ListIN (:)         ! List of independent nodes.
  integer,intent(out):: ListFr (:)         ! List of nodes with independent force vector (not free nodes).

! Local variables.
  integer:: iNode                          ! Counter on the nodes.
  integer:: NumFr                          ! Counter on the force vector.

! Loop on the nodes and remove then from the final list if they are constrained.
  NumIN=0
  NumFr=0
  ListIN=0
  ListFr=0

  do iNode=1,NumN
    select case (BoundCond(iNode))
      case (0)
        NumIN=NumIN+1
        NumFr=NumFr+1
        ListIN(iNode)=NumIN
        ListFr(iNode)=NumFr
      case (-1)
        NumIN=NumIN+1
        ListIN(iNode)=NumIN
      case (1)
        NumFr=NumFr+1
        ListFr(iNode)=NumFr
    end select
  end do

  return
 end subroutine xbeam_undef_nodeindep

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end module xbeam_undef
