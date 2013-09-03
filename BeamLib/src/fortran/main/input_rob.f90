!-> Module.- INPUT Rob Simpson 22/11/2012
!
!-> Author:Rob Simpson copied from input_andrea.f90 (R Palacios and H Hesse)
!
!-> Language: FORTRAN90, Free format.
!
!-> Description.-
!
!  Input data for XBeam test cases (test case x000).
!
!-> Subroutines.-
!
!    -input_setup:    Setup test case.
!    -input_elem :    Compute element information.
!    -input_node :    Compute nodal information.
!    -input_modal:    Compute natural vibration modes.
!    -input_dynsetup: Setup parameters for dynamic solution
!    -input_dynforce: Define time history of applied forces.
!    -input_foredvel: Define time-varying forced velocities.
!    -output_elem:    Write element information in output file.
!
!-> Remarks.-
!
!  1) Test cases from Geradin & Cardona's book.
!  2) 1 point for 2-node beam, to prevent shear locking.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module input
 use xbeam_shared
 implicit none

 ! Shared variables.
 real(8),private,save:: BeamLength1, BeamLength2     ! Beam defining lengths.
 real(8),private,save:: BeamStiffness(6,6)           ! Beam element stiffness matrix (assumed constant).
 real(8),private,save:: BeamMass(6,6)                ! Beam element mass matrix (assumed constant).
 real(8),private,save:: ExtForce(3)                  ! Applied forces at the tip.
 real(8),private,save:: ExtMomnt(3)                  ! Applied moments at the tip.
 integer,private,save:: NumNodesElem                 ! Number of nodes on each element.
 real(8),private,save:: SectWidth,SectHeight         ! Height and width of the cross section.
 real(8),private,save:: ThetaRoot=0.d0               ! Pretwist angle at root.
 real(8),private,save:: ThetaTip =0.d0               ! Pretwist angle at tip.
 real(8),private,save:: TipMass  =0.d0               ! Mass at the beam tip.
 real(8),private,save:: TipMassY =0.d0               ! Y offset of the tip mass.
 real(8),private,save:: TipMassZ =0.d0               ! Z offset of the tip mass.

 real(8),private,save:: Omega   =0.d0                ! Frequency of oscillatory motions.

 character(len=4),private,save:: ElemType            ! ='STRN','DISP'
 character(len=4),private,save:: TestCase            ! Define the test case (CANT,ARCH).
 character(len=2),private,save:: BConds              ! ='CC': Clamped-Clamped
                                                     ! ='CF': Clamped-Free'
 contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-> Subroutine INPUT_SETUP
!
!-> Description:
!
!    Setup test case.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine input_setup (NumElems,OutFile,Options)

! I/O Variables.
  integer,          intent(out):: NumElems       ! Number of elements in the model.
  character(len=25),intent(out):: OutFile        ! Output file.
  type(xbopts),     intent(out):: Options        ! Solution options.

! Local variables.
  real(8) :: E=0.d0,G=0.d0,rho=0.d0
  real(8) :: sigma

  TestCase='TPY0'
  Options%Solution=312     ! See xbeam_shared for options.

! Default values.
  ExtForce(1:3)=0.d0
  ExtMomnt(1:3)=0.d0

! Select test case.
  select case (trim(TestCase))

  case ('HALE')
  ! half wing of HALE model aircraft presented in Table 4 in Murua et al (2011), AIAA-2011-1915
  ! Rob: modified for SharPy PyBeam NonLinearStatic test case 0
    NumElems    = 8
    NumNodesElem= 2
    ThetaRoot   = 0.d0
    ThetaTip    = 0.d0
    BeamLength1 = 16.0d0

    sigma=1.0d0     ! Parameter to change torsion and bending stiffness

    TipMass =0.0d0;
    TipMassY=0.0d0;
    TipMassZ=0.0d0;

    BConds  ='CF'
    ExtForce=(/0.d0,0.d0,800.d0/)*1.d0
    ExtMomnt=(/0.d0,0.d0,800.d0/)*1.d0

    Options%FollowerForce    = .true.
    Options%FollowerForceRig = .true.
    Options%OutInaframe      = .false.
    Options%NumLoadSteps  = 10
    Options%MinDelta      = 1.d-4
    Options%MaxIterations = 99

  case ('GOLD')
  ! half Goland wing taken from Table 3 in Murua et al (2011), AIAA-2011-1915
    NumElems    = 6
    NumNodesElem= 3
    ThetaRoot   = 0.d0
    ThetaTip    = 0.d0
    BeamLength1 = 6.096d0

    TipMass =0.0d0;
    TipMassY=0.0d0;
    TipMassZ=0.0d0;

    BConds  ='CF'
    ExtForce=(/0.d0,0.d0,0.d0/)*1.d5
    ExtMomnt=(/0.d0,0.d0,0.d0/)*1.d0

    Options%FollowerForce    = .false.
    Options%FollowerForceRig = .true.
    Options%OutInaframe      = .false.
    Options%NumLoadSteps  = 10
    Options%MinDelta      = 1.d-4
    Options%MaxIterations = 99

  case ('TPY0')
  ! half wing of HALE model aircraft presented in Table 4 in Murua et al (2011), AIAA-2011-1915
  ! Rob: modified for SharPy PyBeam NonLinearStatic test case 0
    NumElems    = 8
    NumNodesElem= 3
    ThetaRoot   = 0.d0
    ThetaTip    = 0.d0
    BeamLength1 = 16.0d0

    sigma=1.0d0     ! Parameter to change torsion and bending stiffness

    TipMass =0.0d0;
    TipMassY=0.0d0;
    TipMassZ=0.0d0;

    BConds  ='CF'
    ExtForce=(/0.d0,0.d0,80.d0/)*1.d0
    ExtMomnt=(/0.d0,0.d0,0.d0/)*1.d0

    Options%FollowerForce    = .false.
    Options%FollowerForceRig = .true.
    Options%OutInaframe      = .false.
    Options%NumLoadSteps  = 10
    Options%MinDelta      = 1.d-4
    Options%MaxIterations = 99

  end select

! Stiffness and mass properties.
  BeamStiffness=0.d0
  BeamMass     =0.d0

  select case (trim(TestCase))

  case ('GOLD')
    BeamMass(1,1)=35.709121d0
    BeamMass(2,2)=BeamMass(1,1)
    BeamMass(3,3)=BeamMass(1,1)
    BeamMass(4,4)=8.6405832d0
    BeamMass(5,5)=1.d-3
    BeamMass(6,6)=1.d-3

!    BeamMass(1,6)= 6.53048405d0
    BeamMass(3,4)=-BeamMass(1,6)
    BeamMass(4:6,1:3)=transpose(BeamMass(1:3,4:6))

    BeamStiffness(1,1)=1.0d9
    BeamStiffness(2,2)=BeamStiffness(1,1)
    BeamStiffness(3,3)=BeamStiffness(1,1)
    BeamStiffness(4,4)=0.987581d6
    BeamStiffness(5,5)=9.77221d6
    BeamStiffness(6,6)=9.77221d8

  case ('HALE','TPY0')
    BeamMass(1,1)=0.75d0
    BeamMass(2,2)=BeamMass(1,1)
    BeamMass(3,3)=BeamMass(1,1)
    BeamMass(4,4)=1.d-1
    BeamMass(5,5)=1.d-3
    BeamMass(6,6)=1.d-3

    BeamStiffness(1,1)=1.0d9
    BeamStiffness(2,2)=BeamStiffness(1,1)
    BeamStiffness(3,3)=BeamStiffness(1,1)
    BeamStiffness(4,4)=1.0d4
    BeamStiffness(5,5)=2.0d4
    BeamStiffness(6,6)=4.0d6

    BeamStiffness=BeamStiffness*sigma
  end select

! Set name for output file.
  OutFile(1:8)=trim(TestCase)//'_SOL'
  write (OutFile(9:11),'(I3)') Options%Solution

! Solver options.
  select case (Options%Solution)
  case (102,112,142,202,212,302,312,322,900,902,910,912,922,952)
    ElemType= 'DISP'
  case default
    STOP 'Error: Wrong solution code (51707)'
  end select

! Define number of Gauss points (2-noded displacement-based element needs reduced integration).
  select case (NumNodesElem)
    case (2)
      Options%NumGauss=1
    case (3)
      Options%NumGauss=2
  end select

! Minimum angle for two unit vectors to be parallel.
  Options%DeltaCurved=1.d-5

  return
 end subroutine input_setup

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-> Subroutine INPUT_ELEM
!
!-> Description:
!
!    Define element properties.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine input_elem (NumElems,NumNodes,Elem)
 use lib_lu
 use lib_rot

! I/O Variables.
  integer,intent(in) :: NumElems                ! Number of elements in the model.
  integer,intent(out):: NumNodes                ! Number of nodes in the model.
  type(xbelem),intent(out):: Elem(:)            ! Element information

! Local variables.
  integer:: i                            ! Counter.
  integer,save :: fl=2,fw=2   ! Multiplier.
  real(8):: BeamInvStiffness(6,6)        ! Inverse of the stiffness matrix.
  real(8):: LocPos(3)                    ! Local position vector of the lumped mass.

! Connectivies.
  select case (trim(TestCase))

    case default
      select case (NumNodesElem)
      case (2)
        do i=1,NumElems
            Elem(i)%Conn=0
            Elem(i)%Conn(1)=i
            Elem(i)%Conn(2)=i+1
            Elem(i)%NumNodes=2
        end do
        NumNodes=NumElems+1

      case (3)
        do i=1,NumElems
            Elem(i)%Conn=0
            Elem(i)%Conn(1)=2*(i-1)+1
            Elem(i)%Conn(2)=2*(i-1)+3
            Elem(i)%Conn(3)=2*(i-1)+2
            Elem(i)%NumNodes=3
        end do
        NumNodes=2*NumElems+1
    end select
  end select

! Store element stiffness/mass (constant)
  BeamInvStiffness=0.d0
  call lu_invers (BeamStiffness, BeamInvStiffness)
  do i=1,NumElems
    Elem(i)%Stiff   = BeamStiffness
    Elem(i)%Mass    = BeamMass
    Elem(i)%InvStiff= BeamInvStiffness
  end do

! Define lumped masses at element nodes.
  do i=1,NumElems
    Elem(i)%RBMass = 0.d0
  end do

  select case (trim(TestCase))

  case default
        LocPos(1)=0.d0
        LocPos(2)=TipMassY
        LocPos(3)=TipMassZ

!       Elem(NumElems/2)%RBMass(2,1:3,1:3)= TipMass*Unit
!       Elem(NumElems/2)%RBMass(2,1:3,4:6)=-TipMass*rot_skew(LocPos)
!       Elem(NumElems/2)%RBMass(2,4:6,1:3)= TipMass*rot_skew(LocPos)
!       Elem(NumElems/2)%RBMass(2,4:6,4:6)=-TipMass*matmul(rot_skew(LocPos),rot_skew(LocPos))
!
!       Elem(NumElems)%RBMass(1,:,:)= Elem(NumElems/2)%RBMass(2,:,:)
  end select

! Element orientation.
  do i=1,NumElems
    select case (trim(TestCase))
    case ('HALE','GOLD','TPY0')
      Elem(i)%Vector(2)= 1.d0
    end select
  end do

! Define element types.
  select case (ElemType)
  case ('DISP')
    Elem(1:NumElems)%MemNo=0
  end select

  return
 end subroutine input_elem


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-> Subroutine INPUT_NODE
!
!-> Description:
!
!    Define nodal properties.
!
!-> Remarks.-
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine input_node (NumNodes,Elem,BoundConds,Coords,Forces,PhiNodes)

! I/O Variables.
  integer,      intent(in) :: NumNodes            ! Number of nodes in the model.
  type(xbelem), intent(in) :: Elem      (:)       ! Element information
  integer,      intent(out):: BoundConds(:)       ! =0 on free nodes; =1 on clamped nodes; =-1 on free nodes.
  real(8),      intent(out):: Coords    (:,:)     ! Initial nodal Coordinates.
  real(8),      intent(out):: Forces    (:,:)     ! Applied nodal forces.
  real(8),      intent(out):: PhiNodes  (:)       ! Initial twist at grid points.

! Local variables.
  integer      :: i           ! Counters.
  integer,save :: fl=2,fw=2   ! Multiplier.
  real(8)      :: Theta       ! Parameter on curved beams.

! Initial position vector of grid points.
  Coords= 0.d0
  do i=1,NumNodes
    select case (trim(TestCase))
    case ('HALE','GOLD','TPY0')
      Coords(i,1)=BeamLength1*(dble(i-1)/dble(NumNodes-1))
    end select
  end do

! Initial pretwist angle.
  do i=1,NumNodes
    PhiNodes(i)=ThetaRoot+(ThetaTip-ThetaRoot)*(dble(i-1)/dble(NumNodes-1))
  end do

! Static point forces.
  Forces=0.d0
  select case (trim(TestCase))

  case default
    Forces(1,1:3)=-ExtForce
    Forces(1,4:6)=-ExtMomnt
    Forces(NumNodes,1:3)=ExtForce
    Forces(NumNodes,4:6)=ExtMomnt
  end select

! Boundary conditions
  BoundConds=0
  select case (trim(TestCase))

  case default
    if (BConds(1:1).eq.'C') BoundConds(1)       = 1
    if (BConds(2:2).eq.'C') BoundConds(NumNodes)= 1
    if (BConds(2:2).eq.'F') BoundConds(NumNodes)=-1
  end select

  return
 end subroutine input_node


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-> Subroutine INPUT_DYNSETUP
!
!-> Description:
!
!    Setup test case.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine input_dynsetup (NumSteps,t0,dt,Options)

! I/O Variables.
  integer,intent(out):: NumSteps             ! Number of time steps.
  real(8),intent(out):: t0                   ! Initial time.
  real(8),intent(out):: dt                   ! Initial, final and delta time.
  type(xbopts),intent(inout):: Options      ! Solution options.

! Local variables.
  real(8):: tfin                             ! Final time.

  select case (trim(TestCase))

  case ('GOLD')
    t0  = 0.0d0   ! 0.d0
    dt  = 1.0d-3
    tfin= 1.0d0
    NumSteps= ceiling((tfin-t0)/dt)
    Options%NewmarkDamp=0.01d0
    Omega=20.d0    ! rad/s

  case ('HALE')
    t0  = 0.0d0   ! 0.d0
    dt  = 1.0d-2
    tfin= 10.0d0
    NumSteps= ceiling((tfin-t0)/dt)
    Options%NewmarkDamp=0.05d0
    Omega=20.d0    ! rad/s

  case ('TPY0')
    t0  = 0.0d0   ! 0.d0
    dt  = 1.0d-2
    tfin= 1.0d0
    NumSteps= ceiling((tfin-t0)/dt) + 2
    Options%NewmarkDamp=0.01d0
    Omega=20.d0    ! rad/s

  case default
    STOP 'Error: Input data for dynamic analysis not defined.'
  end select

  return
 end subroutine input_dynsetup


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-> Subroutine INPUT_DYNFORCE
!
!-> Description:
!
!    Define time-varying forcing terms.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine input_dynforce (NumNodes,Time,ForceStatic,ForceDynAmp,ForceTime)

! I/O Variables.
  integer,intent(in) :: NumNodes            ! Number of nodes in the model.
  real(8),intent(in) :: Time(:)             ! Time history.
  real(8),intent(in) :: ForceStatic(:,:)    ! Static force.
  real(8),intent(out):: ForceDynAmp(:,:)    ! Nodal force amplitudes.
  real(8),intent(out):: ForceTime  (:)      ! Time history of the applied forces.

! Local variables.
  integer:: i             ! Counter.
  integer:: NumSteps
  real(8):: Time0         ! Time for discontinuity in function.

! Initialize
  NumSteps=size(Time)-1

  select case (trim(TestCase))

  case default
    !ForceDynAmp = 2.d0*ForceStatic
    ForceTime(:)   = 1.d0
    ForceDynAmp(NumNodes,2) = 160.d0

! Ramped harmonic load.
    if (.false.) then
      Time0= Time(NumSteps)/2.d0
      ForceTime=0.d0
      do i=1,NumSteps
        ForceTime(i+1)=sin(Omega*Time(i+1))
        if (Time(i+1) < Time0) ForceTime(i+1)=ForceTime(i+1)*Time(i+1)/Time0
      end do
    end if

! Initial Load
    if (.false.) then
      ForceTime(:) = 0.d0
      ForceTime(1)= 1.d0
      ForceDynAmp(NumNodes,3) = 100.d3
    end if

! 1-Cos load.
    if (.false.) then
     ForceDynAmp(NumNodes,3) = 1.d03
     do i=1,NumSteps+1
         if ((Time(i).ge.0.d0).and.(Time(i).le.1.d-2)) then
             ForceTime(i)=(1.d0-cos(2*Pi*dble(Time(i)-0.d0)/dble(1.d-2)))/2.d0
         end if
     end do
    end if

! Ramped load.
    if (.false.) then
      ForceTime=1.d0
      do i=1,NumSteps+1
        if ((Time(i).ge.0.d0).and.(Time(i).le.2.5d0)) then
            ForceTime(i)=Time(i)/2.5d0
        end if
      end do
    end if

! sinusiodal load.
    if (.false.) then
     do i=1,NumSteps+1
        ! ForceTime(i) = sin(2.*1.*Time(i))
        ForceTime(i)=sin(Pi*dble(Time(i)-0.d0)/dble(0.5))
     end do
    end if
  end select
print *, ForceDynAmp
  return
 end subroutine input_dynforce


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-> Subroutine INPUT_FORCEDVEL
!
!-> Description:
!
!    Define time-varying forcing velocity.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine input_forcedvel (NumNodes,Time,ForcedVel,ForcedVelDot)

! I/O Variables.
  integer,intent(in) :: NumNodes           ! Number of nodes in the model.
  real(8),intent(in) :: Time(:)            ! Time history.
  real(8),intent(out):: ForcedVel(:,:)     ! Forced root velocities.
  real(8),intent(out):: ForcedVelDot(:,:)  ! Time derivatives of the forced root velocities.

! Local variables.
  integer:: i             ! Counter.
  integer:: NumSteps
  real(8):: VelAmp(6)     ! Amplitude of the velocities.

! Initialize.
  NumSteps=size(Time)-1
  ForcedVel   =0.d0
  ForcedVelDot=0.d0
  VelAmp=0.d0

  return
 end subroutine input_forcedvel

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-> Subroutine OUTPUT_ELEMS
!
!-> Description:
!
!    Write information in elements.
!
!-> Remarks.-
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine output_elems (iuOut,Elem,Coords,Psi)
  use lib_fem

! I/O Variables.
  integer,      intent(in)   :: iuOut             ! Output file.
  type(xbelem),intent(in)   :: Elem   (:)        ! Element information.
  real(8),      intent(in)   :: Coords (:,:)      ! Coordinates of the grid points.
  real(8),      intent(in)   :: Psi    (:,:,:)    ! CRV of the nodes in the elements.

! Local variables.
  integer:: i                    ! Counter.
  integer:: iElem                ! Counter on the finite elements.
  integer:: NumE                 ! Number of elements in the model.
  integer:: NumNE                ! Number of nodes in an element.
  real(8):: PosElem (MaxElNod,3) ! Coordinates/CRV of nodes in the element.

  NumE=size(Elem)

! Loop in the elements in the model.
  do iElem=1,NumE
    call fem_glob2loc_extract (Elem(iElem)%Conn,Coords,PosElem,NumNE)

    do i=1,NumNE
      write (iuOut,'(2I4,1P6E15.6)') iElem,i,PosElem(i,:),Psi(iElem,i,:)
    end do
  end do

end subroutine output_elems

end module input
