!-> Module.- INPUT Henrik Hesse & Rafa Palacios. 07/01/2011
!
!-> Language: FORTRAN90, Free format.
!
!-> Description.-
!
!  Input data for XBeam test cases (test case x000)
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
 real(8),private,save:: BeamLength1=0.d0             ! Beam defining lengths.
 real(8),private,save:: BeamLength2=0.d0             ! Beam defining lengths.
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

  TestCase='NCB1'
  Options%Solution=112     ! See xbeam_shared for options.

! Default values.
  ExtForce(1:3)=0.d0
  ExtMomnt(1:3)=0.d0
  BeamStiffness=0.d0
  BeamMass     =0.d0

! Select test case.
  select case (trim(TestCase))
   case ('NCB1')
  ! Cantilever beam in Geradin & Cardona, p 133.
  ! Analytical NatFreqs (rad/s): 43.0(B1), 99.3(T1), 269.4(B2), 298.0(T2), 496.7(T3), 699.9(A1), 759.5(B3).
    NumElems    = 20
    NumNodesElem= 2
    ThetaRoot   = 0.d0
    ThetaTip    = 0.d0
    BeamLength1 = 5.0d0

    BConds  ='CF'
    ExtForce=(/0.d0,0.d0,600.d3/)
    ExtMomnt=(/0.d0,0.d0,0.d0/)
    Options%FollowerForce=.false.

    BeamStiffness(1,1)= 4.8d8
    BeamStiffness(2,2)= 3.231d8
    BeamStiffness(3,3)= 3.231d8
    BeamStiffness(4,4)= 1.d6
    BeamStiffness(5,5)= 9.346d6
    BeamStiffness(6,6)= 9.346d6

    BeamMass(1,1)=100.d0
    BeamMass(2,2)=BeamMass(1,1)
    BeamMass(3,3)=BeamMass(1,1)
    BeamMass(4,4)=10.d0
    BeamMass(5,5)=10.d0
    BeamMass(6,6)=10.d0

    Options%NumLoadSteps  = 5
    Options%MinDelta      = 1.d-5
    Options%MaxIterations = 99

  case ('WINF')
    NumElems   =  6

    ThetaRoot  =  0.d0
    ThetaTip   =  0.d0   ! Pi/6.d0
    BeamLength1= 10.d0
    BeamLength2=  5.d0

    SectWidth =1.d-1
    SectHeight=5.d-2
    E  = 70.d11
    G  = E/(2.d0*(1.d0+0.30d0))
    rho= 2700.d0

    BConds  ='CF'
    ExtForce=(/0.d0,-1.d0,0.d0/)*1.d3
    ExtMomnt=(/-5.d0,0.d0,0.d0/)*1.d2

    Options%FollowerForce    = .false.
    Options%FollowerForceRig = .false.
    Options%OutInaframe      = .false.
    Options%NumLoadSteps     = 10
    Options%MinDelta         = 1.d-2
    Options%MaxIterations    = 99

  case ('CANT')
    NumElems   =  2

    ThetaRoot  =  0.d0
    ThetaTip   =  0.d0   ! Pi/6.d0
    BeamLength1=  1.D0   ! 16.d0
    BeamLength2=  5.d0

    SectWidth =1.d-1
    SectHeight=5.d-2
    E  = 70.d9
    G  = E/(2.d0*(1.d0+0.30d0))
    rho= 2700.d0

    BConds  ='CF'
    ExtForce=(/0.d0,0.d0,1.d0/)*1.d8
    ExtMomnt=(/0.d0,0.d0,0.d0/)*1.d2

    Options%FollowerForce    = .true.
    Options%FollowerForceRig = .true.
    Options%OutInaframe      = .false.
    Options%NumLoadSteps  = 1
    Options%MinDelta      = 1.d-3
    Options%MaxIterations = 99

  case ('PATE')
    ! Follo Patil and Hodges Flying Wing (2006)
    NumElems   = 24

    ThetaRoot  = 0.d0
    ThetaTip   = 0.d0       ! Pi/6.d0
    BeamLength1= 80.d0
    BeamLength2= 40.d0

    BConds  ='CF'
    ExtForce=(/0.d0,0.d0,5.2d0/)*1.d3
    ExtMomnt=(/0.d0,0.d0,0.d0/)*1.d0
    Options%FollowerForce    = .false.
    Options%FollowerForceRig = .false.
    Options%OutInaframe      = .false.
    Options%NumLoadSteps  = 20
    Options%MinDelta      = 1.d-4
    Options%MaxIterations = 999

  case ('WING')
! This is a full Aluminium wing with tip vertical forces of opposite sign at both ends.
! This test case was used to test the routines for lumped masses for SOL 202/212.
    NumElems   = 50
    BeamLength1= 10.d0

    SectWidth =1.d-1
    SectHeight=5.d-2
    E  = 70.d9
    G  = E/(2.d0*(1.d0+0.30d0))
    rho= 2700.d0
		TipMass =0.0d0;
		TipMassY=0.0d0;
		TipMassZ=0.0d0;

    BConds  ='CF'
    ExtForce=(/0.d0,0.d0,100.d0/)*1.d0
    ExtMomnt=(/100.d0,0.d0,0.d0/)*0.d0

    Options%FollowerForce    = .true.
    Options%FollowerForceRig = .true.
    Options%OutInaframe      = .false.
    Options%NumLoadSteps  = 1
    Options%MinDelta      = 1.d-6
    Options%MaxIterations = 99

  ! L-frame from previous case.
  case ('FRAM')
    NumElems=20
    BeamLength1=5.0d0
    BeamLength2=5.0d0

    BConds='CF'
    ExtForce=(/ 0.d0, 0.d0, -1.0d0/)*1.0d5
    ExtMomnt=(/ 0.d0, 0.d0, 0.d0/)
    Options%FollowerForce=.true.
    Options%NumLoadSteps = 1
    Options%MinDelta     = 1.d-7

  ! Cantilever beam as in Eric Brown's thesis (p. 84).
  case ('CB2F')
    NumElems   = 40
    ThetaRoot  = 0.d0
    ThetaTip   = 0.d0   ! Pi/6.d0
    BeamLength1= 1.0d0

    ExtForce=(/0.d0,0.d0,10.d0/)
    ExtMomnt=(/0.d0,0.d0,0.d0/)
    Options%FollowerForce=.false.
    BConds  ='CF'
    Options%MinDelta     = 1.d-5

  ! 45-degree bend in Geradin & Cardona, p 133.
  case ('BEND')
    NumElems= 8
    BeamLength2=100.d0

    SectWidth= 1.d0
    SectHeight=1.d0
    E  = 1.0d7
    G  = 0.5d7
    rho= 1.

    ExtForce=(/0.d0,0.d0,200.d0/)
    ExtMomnt=(/0.d0,0.d0,0.d0/)
    Options%FollowerForce=.true.
    BConds  ='CF'

    Options%MinDelta     = 1.d-4     ! 1.d-4
    Options%NumLoadSteps =  2        ! 1
    Options%MaxIterations = 999

  ! Results in Pai's book (p. 392). X axis goes along the vertical beam.
  case ('FPAI')
    NumElems   = 40
    ThetaRoot  = 0.d0
    ThetaTip   = 0.d0
    BeamLength1= 479.d-3

    SectWidth= 50.8d-3
    SectHeight=0.45d-3
    E  = 1.27d11
    G  = E/(2.d0*(1.d0+0.36d0))
    rho= 4.43d3

    ExtForce=(/-9.81d0,0.d0,0.d0/)*rho*SectWidth*SectHeight*BeamLength1/dble(NumElems)
    Options%FollowerForce=.false.
    Options%NumLoadSteps = 10
    BConds ='CF'
    Options%MinDelta = 1.d-5

  case ('ARCH')
    NumElems   = 100
    BeamLength2= 100.d0
    ThetaRoot  = 0.d0
    ThetaTip   = 0.D0   ! Pi/6.d0
    SectWidth  = 2.d0
    SectHeight = 1.d0

    ExtForce= (/0.d0,0.d0,5.d0/)
    ExtMomnt=   0.d0
    Options%FollowerForce=.false.
    Options%NumLoadSteps = 20
    BConds  ='CF'

  end select



  select case (trim(TestCase))

  case ('PATE','CANT')
    BeamMass(1,1)=10.d0
    BeamMass(2,2)=BeamMass(1,1)
    BeamMass(3,3)=BeamMass(1,1)
    BeamMass(4,4)=1.d-3
    BeamMass(5,5)=1.d-3
    BeamMass(6,6)=1.d-3
    BeamStiffness(1,1)=1.0d6
    BeamStiffness(2,2)=1.0d6
    BeamStiffness(3,3)=1.0d6
    BeamStiffness(4,4)=1.0d6
    BeamStiffness(5,5)=1.0d6
    BeamStiffness(6,6)=1.0d6


  case ('WING')
    BeamMass(1,1)=rho*SectWidth*SectHeight
    BeamMass(2,2)=rho*SectWidth*SectHeight
    BeamMass(3,3)=rho*SectWidth*SectHeight
    BeamMass(5,5)=rho*SectWidth*(SectHeight**3.d0)/12.d0
    BeamMass(6,6)=rho*SectHeight*(SectWidth**3.d0)/12.d0
    BeamMass(4,4)=0.5d0*(BeamMass(5,5)+BeamMass(6,6))

! Following Geradin and Cardona (p.119)
    BeamStiffness(1,1)=SectWidth*SectHeight*E
    BeamStiffness(2,2)=(5.d0/6.d0)*SectWidth*SectHeight*G
    BeamStiffness(3,3)=(5.d0/6.d0)*SectWidth*SectHeight*G
    BeamStiffness(5,5)=E*(SectWidth*(SectHeight**3.d0)/12.d0)
    BeamStiffness(6,6)=E*(SectHeight*(SectWidth**3.d0)/12.d0)
    BeamStiffness(4,4)=G*0.228d0*SectWidth*(SectHeight**3.d0)

  case ('WINF')
    BeamMass(1,1)=rho*SectWidth*SectHeight
    BeamMass(2,2)=BeamMass(1,1)
    BeamMass(3,3)=BeamMass(1,1)
    BeamMass(4,4)=rho*SectWidth*(SectHeight**3.d0)/12.d0 + rho*SectHeight*(SectWidth**3.d0)/12.d0
    BeamMass(5,5)=rho*(SectHeight**3.d0)/12.d0
    BeamMass(6,6)=(BeamMass(4,4)+BeamMass(5,5))

! Following Geradin and Cardona (p.119)
    BeamStiffness(1,1)=SectWidth*SectHeight*E
    BeamStiffness(2,2)=(5.d0/6.d0)*SectWidth*SectHeight*G
    BeamStiffness(3,3)=(5.d0/6.d0)*SectWidth*SectHeight*G
    BeamStiffness(5,5)=E*(SectWidth*(SectHeight**3.d0)/12.d0)
    BeamStiffness(6,6)=E*(SectHeight*(SectWidth**3.d0)/12.d0)
    BeamStiffness(4,4)=G*0.228d0*SectWidth*(SectHeight**3.d0)

  case ('FRAM')
    BeamStiffness(1,1)= 4.8d8
    BeamStiffness(2,2)= 3.231e8
    BeamStiffness(3,3)= 3.231e8
    BeamStiffness(4,4)= 1e6
!    BeamStiffness(5,5)= 9.346e6
!    BeamStiffness(6,6)= 9.346e6
    BeamStiffness(5,5)= 23.e6
    BeamStiffness(6,6)= 37.e10

!    BeamStiffness=1e-4*BeamStiffness

    BeamMass(1,1)=10.d0
    BeamMass(2,2)=BeamMass(1,1)
    BeamMass(3,3)=BeamMass(1,1)
    BeamMass(4,4)=1.d0
    BeamMass(5,5)=1.d0
    BeamMass(6,6)=1.d0

  case ('BEND','FPAI')
    BeamStiffness(1,1)=SectWidth*SectHeight*E
    BeamStiffness(2,2)=(5.d0/6.d0)*SectWidth*SectHeight*G
    BeamStiffness(3,3)=(5.d0/6.d0)*SectWidth*SectHeight*G
    BeamStiffness(5,5)=SectWidth*(SectHeight**3.d0)*E/12.d0
    BeamStiffness(6,6)=SectHeight*(SectWidth**3.d0)*E/12.d0
    BeamStiffness(4,4)=0.5d0*(BeamStiffness(5,5)+BeamStiffness(6,6))

    BeamMass(1,1)= rho*SectWidth*SectHeight
    BeamMass(2,2)=BeamMass(1,1)
    BeamMass(3,3)=BeamMass(1,1)
    BeamMass(5,5)=rho*SectWidth*(SectHeight**3.d0)/12.d0
    BeamMass(6,6)=rho*SectHeight*(SectWidth**3.d0)/12.d0
    BeamMass(4,4)=0.5d0*(BeamMass(5,5)+BeamMass(6,6))

  case ('CB2F')
    BeamStiffness(1,1)=1.0d6
    BeamStiffness(2,2)=1.0d10
    BeamStiffness(3,3)=1.0d10
    BeamStiffness(4,4)=5.0d1
    BeamStiffness(5,5)=5.0d1
    BeamStiffness(6,6)=1.0e3

    BeamMass(1,1)=0.2d0
    BeamMass(2,2)=BeamMass(1,1)
    BeamMass(3,3)=BeamMass(1,1)
    BeamMass(4,4)=1.0d-4
    BeamMass(5,5)=1.0d-6
    BeamMass(6,6)=1.0d-4

  case ('ARCH')
    BeamStiffness(1,1)=SectWidth*SectHeight*1.0d7
    BeamStiffness(2,2)=(5.d0/6.d0)*SectWidth*SectHeight*0.5d7
    BeamStiffness(3,3)=(5.d0/6.d0)*SectWidth*SectHeight*0.5d7
    BeamStiffness(5,5)=SectWidth*(SectHeight**3.d0)*1.0d7/12.d0
    BeamStiffness(6,6)=SectHeight*(SectWidth**3.d0)*1.0d7/12.d0
    BeamStiffness(4,4)=0.5d0*(BeamStiffness(5,5)+BeamStiffness(6,6))
  end select

! Set name for output file.
  OutFile(1:8)=trim(TestCase)//'_SOL'
  write (OutFile(9:11),'(I3)') Options%Solution

! Solver options.
  select case (Options%Solution)
  case (101,111,201,211,311)
    ElemType='STRN'
    NumNodesElem=2
  case (102,112,142,202,212,302,312,322,900,902,910,912,922)
    ElemType= 'DISP'
    NumNodesElem=2
  case default
    STOP 'Error: Wrong solution code (51707)'
  end select

! Define number of Gauss points (2-noded displacement-based element needs reduced integration).
  select case (NumNodesElem)
    case (2)
      Options%NumGauss=1
    case (3)
      Options%NumGauss=5
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
  real(8):: BeamInvStiffness(6,6)        ! Inverse of the stiffness matrix.
	real(8):: LocPos(3)                    ! Local position vector of the lumped mass.

! Connectivies.
  select case (trim(TestCase)) 
    case ('WINF')
      select case (NumNodesElem)
      case (2)
        NumNodes=NumElems+1
        do i=1,NumElems
            Elem(i)%Conn=0
            Elem(i)%Conn(1)=i
            Elem(i)%Conn(2)=i+1
            Elem(i)%NumNodes=2
        end do

        Elem(NumElems/3+1)%Conn(1)=1
        Elem(2*NumElems/3+1)%Conn(1)=(NumNodes-1)/3+1
        Elem(5*NumElems/6+1)%Conn(1)=2*(NumNodes-1)/3+1

      case (3)
        NumNodes=2*NumElems+1
        do i=1,NumElems/2
            Elem(i)%Conn=0
            Elem(i)%Conn(1)=2*(i-1)+1
            Elem(i)%Conn(2)=2*(i-1)+3
            Elem(i)%Conn(3)=2*(i-1)+2
            Elem(i)%NumNodes=3
        end do
        do i=NumElems/2+1,NumElems
            Elem(i)%Conn=0
            Elem(i)%Conn(1)=2*(i-1)+3
            Elem(i)%Conn(2)=2*(i-1)+1
            Elem(i)%Conn(3)=2*(i-1)+2
            Elem(i)%NumNodes=3
        end do
        Elem(NumElems/2+1)%Conn(2)=1
      end select

    case ('WING','PATE')
      select case (NumNodesElem)
      case (2)
        NumNodes=NumElems+1
        do i=1,NumElems/2
            Elem(i)%Conn=0
            Elem(i)%Conn(1)=i
            Elem(i)%Conn(2)=i+1
            Elem(i)%NumNodes=2
        end do
        do i=NumElems/2+1,NumElems
            Elem(i)%Conn=0
            Elem(i)%Conn(1)=i+1
            Elem(i)%Conn(2)=i
            Elem(i)%NumNodes=2
        end do
        Elem(NumElems/2+1)%Conn(2)=1

      case (3)
        NumNodes=2*NumElems+1
        do i=1,NumElems/2
            Elem(i)%Conn=0
            Elem(i)%Conn(1)=2*(i-1)+1
            Elem(i)%Conn(2)=2*(i-1)+3
            Elem(i)%Conn(3)=2*(i-1)+2
            Elem(i)%NumNodes=3
        end do
        do i=NumElems/2+1,NumElems
            Elem(i)%Conn=0
            Elem(i)%Conn(1)=2*(i-1)+3
            Elem(i)%Conn(2)=2*(i-1)+1
            Elem(i)%Conn(3)=2*(i-1)+2
            Elem(i)%NumNodes=3
        end do
        Elem(NumElems/2+1)%Conn(2)=1
      end select

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
  case ('WING')
		LocPos(1)=0.d0
		LocPos(2)=TipMassY
		LocPos(3)=TipMassZ
		
		Elem(NumElems/2)%RBMass(2,1:3,1:3)= TipMass*Unit
		Elem(NumElems/2)%RBMass(2,1:3,4:6)=-TipMass*rot_skew(LocPos)
		Elem(NumElems/2)%RBMass(2,4:6,1:3)= TipMass*rot_skew(LocPos)
		Elem(NumElems/2)%RBMass(2,4:6,4:6)=-TipMass*matmul(rot_skew(LocPos),rot_skew(LocPos))

		Elem(NumElems)%RBMass(1,:,:)= Elem(NumElems/2)%RBMass(2,:,:)
	end select

! Element orientation.
  do i=1,NumElems
		Elem(i)%Vector=0.d0
    select case (trim(TestCase))
    case ('WING')
      if (i.le.NumElems/2) then
        Elem(i)%Vector(2)= 1.d0
      else
        Elem(i)%Vector(2)= 1.d0
      end if

    case ('BEND')
      Elem(i)%Vector= 0.d0
      Elem(i)%Vector(1)=-dcos((Pi/4.d0)*(dble(i-1)/dble(NumElems)))
      Elem(i)%Vector(2)= dsin((Pi/4.d0)*(dble(i-1)/dble(NumElems)))

    case ('WINF')
      if      (i.le.NumElems/3) then
        Elem(i)%Vector(2)= 1.d0
      elseif  (i.le.2*NumElems/3) then
        Elem(i)%Vector(2)= 1.d0
      elseif  (i.le.5*NumElems/6) then
        Elem(i)%Vector(2)= 1.d0
      else
        Elem(i)%Vector(2)= 1.d0
      end if
      
    case default
      Elem(i)%Vector(2)= 1.d0
    end select
  end do


! Define element types.
  select case (ElemType)
  case ('DISP')
    Elem(1:NumElems)%MemNo=0
  case ('STRN')
    Elem(1:NumElems)%MemNo=1
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
    case ('ZEPH')
      ! ZEPHYR REPRESENTATION
      ! =====================
      if      (i.le.(fl*BeamLength2+1)) then
      ! FUSELAGE
          Coords(i,2)=BeamLength2/4 - BeamLength2*(dble(i-1)/dble(fl*BeamLength2))
      elseif  (i.le.((fl*BeamLength2+1)+(fw*BeamLength1))) then
      ! Wing 1
          Coords(i,1)=-BeamLength1*(1-(dble(i-1-(fl*BeamLength2+1))/dble(fw*BeamLength1)))
      elseif  (i.le.((fl*BeamLength2+1)+(fw*2*BeamLength1))) then
      ! Wing 2
          Coords(i,1)=BeamLength1*(dble(i-(fl*BeamLength2+1)-(fw*BeamLength1))/dble(fw*BeamLength1))
      elseif  (i.le.((fl*BeamLength2+1)+(fw*2*BeamLength1)+(fw*(BeamLength1/5)))) then
      ! TAIL Wing 1
          Coords(i,2)=-3*BeamLength2/4
          Coords(i,1)=-(BeamLength1/5)*(1-(dble(i-1-(fl*BeamLength2+1)-(fw*2*BeamLength1))/dble(fw*BeamLength1/5)))
      elseif  (i.le.((fl*BeamLength2+1)+(fw*2*BeamLength1)+(fw*2*(BeamLength1/5)))) then
      ! TAIL Wing 2
          Coords(i,2)=-3*BeamLength2/4
          Coords(i,1)=(BeamLength1/5)*(dble(i-(fl*BeamLength2+1)-(fw*2*BeamLength1)-(fw*(BeamLength1/5)))/dble(fw*BeamLength1/5))
      else
      ! TAIL BOOM
          Coords(i,2)=-3*BeamLength2/4
          Coords(i,3)=(BeamLength1/5)*  &
&                      (dble(i-(fl*BeamLength2+1)-(fw*2*BeamLength1)-(fw*2*(BeamLength1/5)))/dble(fl*BeamLength1/5))
      end if
      
    case ('WINF')
      select case (NumNodesElem)
      case (2)
        if      (i.le.(NumNodes-1)/3+1) then
           Coords(i,1)=BeamLength1*(dble(i-1)/dble((NumNodes-1)/3))
        elseif  (i.le.2*(NumNodes-1)/3+1) then
           Coords(i,1)=-BeamLength1*(dble(i-1-(NumNodes-1)/3)/dble((NumNodes-1)/3))
        elseif  (i.le.5*(NumNodes-1)/6+1) then
           Coords(i,1)=BeamLength1
           Coords(i,3)=BeamLength2*(dble(i-1-2*(NumNodes-1)/3)/dble((NumNodes-1)/6))
        else
           Coords(i,1)=-BeamLength1
           Coords(i,3)=BeamLength2*(dble(i-1-5*(NumNodes-1)/6)/dble((NumNodes-1)/6))
        end if
      case (3)
        if      (i.le.(NumNodes-1)/3+1) then
           Coords(i,1)=BeamLength1*(dble(i-1)/dble((NumNodes-1)/3))
        elseif  (i.le.(NumNodes+1)/2) then
           Coords(i,1)=BeamLength1
           Coords(i,3)=BeamLength2*(dble(i-1-(NumNodes-1)/3)/dble((NumNodes-1)/6))
        elseif  (i.le.5*(NumNodes-1)/6+1) then
           Coords(i,1)=-BeamLength1*(dble(i-(NumNodes+1)/2)/dble((NumNodes-1)/3))
        else
           Coords(i,1)=-BeamLength1
           Coords(i,3)=BeamLength2*(dble(i-1-5*(NumNodes-1)/6)/dble((NumNodes-1)/6))
        end if
      end select

    case ('PATE')
      select case (NumNodesElem)
      case (2)
        if      (i.le.(NumNodes-1)/3+1) then
           Coords(i,1)=BeamLength1*(dble(i-1)/dble((NumNodes-1)/3))
        elseif  (i.le.(NumNodes+1)/2) then
           Coords(i,1)=BeamLength1 + BeamLength2*(dble(i-1-(NumNodes-1)/3)/dble((NumNodes-1)/6))*cos(0.174532925)
           Coords(i,3)=BeamLength2*(dble(i-1-(NumNodes-1)/3)/dble((NumNodes-1)/6))*sin(0.174532925)
        elseif  (i.le.5*(NumNodes-1)/6+1) then
           Coords(i,1)=-BeamLength1*(dble(i-(NumNodes+1)/2)/dble((NumNodes-1)/3))
        else
           Coords(i,1)=-BeamLength1 - BeamLength2*(dble(i-1-5*(NumNodes-1)/6)/dble((NumNodes-1)/6))*cos(0.174532925)
           Coords(i,3)= BeamLength2*(dble(i-1-5*(NumNodes-1)/6)/dble((NumNodes-1)/6))*sin(0.174532925)
        end if
      case (3)
        if      (i.le.(NumNodes-1)/3+1) then
           Coords(i,1)=BeamLength1*(dble(i-1)/dble((NumNodes-1)/3))
        elseif  (i.le.(NumNodes+1)/2) then
           Coords(i,1)=BeamLength1
           Coords(i,3)=BeamLength2*(dble(i-1-(NumNodes-1)/3)/dble((NumNodes-1)/6))
        elseif  (i.le.5*(NumNodes-1)/6+1) then
           Coords(i,1)=-BeamLength1*(dble(i-(NumNodes+1)/2)/dble((NumNodes-1)/3))
        else
           Coords(i,1)=-BeamLength1
           Coords(i,3)=BeamLength2*(dble(i-1-5*(NumNodes-1)/6)/dble((NumNodes-1)/6))
        end if
      end select

    case ('WING')
      select case (NumNodesElem)
      case (2)
        if (i.le.NumNodes/2+1) then
           Coords(i,1)= BeamLength1*(dble(i-1)/dble(NumNodes/2))
        else
           Coords(i,1)=-BeamLength1*(dble(i-1-NumNodes/2)/dble(NumNodes/2))
        end if
      case (3)
        if (i.le.(NumNodes+1)/2) then
           Coords(i,1)=BeamLength1*(dble(i-1)/dble((NumNodes+1)/2-1))
        else
           Coords(i,1)=-BeamLength1*(dble(i-(NumNodes+1)/2)/dble((NumNodes+1)/2-1))
        end if
      end select


    case ('FRAM')
      if (i.le.NumNodes/2+1) then
        Coords(i,3)=BeamLength1*(dble(i-1)/dble(NumNodes/2))
      else
        Coords(i,1)=BeamLength2*(dble(i-1-NumNodes/2)/dble(NumNodes/2))
        Coords(i,3)=BeamLength1
      end if

    case ('BEND')
      Theta=(Pi/4.d0)*(dble(i-1)/dble(NumNodes-1))
      Coords(i,1)=BeamLength2*(1.d0-dcos(Theta))
      Coords(i,2)=BeamLength2*(     dsin(Theta))

	  case ('ARCH')
      Theta=dble(i-1)*(Pi/dble(NumNodes-1))
      Coords(i,1)=BeamLength2*(1.d0-dcos(Theta))
      Coords(i,3)=BeamLength2*(     dsin(Theta))

    case default
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
  case ('WINF')
    Forces(1,1:3)       = ExtForce
    Forces(1,4:6)       = ExtMomnt
!    Forces(5*(NumNodes-1)/6+1,1:3)      = ExtForce
!    Forces(NumNodes,1:3)                = ExtForce
    do i=2,NumNodes
!        Forces(i,1:3)=ExtForce/(NumNodes-1)
    end do

  case ('WING')
    Forces( NumNodes,       1:3)=-ExtForce
    Forces((NumNodes-1)/2+1,1:3)= ExtForce
    Forces( NumNodes,       4:6)=-ExtMomnt
    Forces((NumNodes-1)/2+1,4:6)= ExtMomnt

  case ('FPAI')
    do i=2,NumNodes
      Forces(i,1:3)=ExtForce
      Forces(i,4:6)=0.0d0
    end do

  case default
    Forces(NumNodes,1:3)=ExtForce
    Forces(NumNodes,4:6)=ExtMomnt
  end select

! Boundary conditions
  BoundConds=0
  select case (trim(TestCase))
  case ('WING','PATE')
    BoundConds(1)               = 1
    BoundConds((NumNodes-1)/2+1)=-1
    BoundConds(NumNodes)        =-1
  case ('WINF')
    BoundConds(1)                 = 1
    BoundConds(5*(NumNodes-1)/6+1)=-1
    BoundConds(NumNodes)          =-1
  case default
    if (BConds(1:1).eq.'C') BoundConds(1)=1
    if (BConds(2:2).eq.'C') BoundConds(NumNodes)=1
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

  case ('PATE')
    t0  =  0.d0    ! 0.d0
    dt  =  2.d-2
    tfin=  4.d+0
    NumSteps= ceiling((tfin-t0)/dt)
    Options%NewmarkDamp=0.1d0
    Omega=20.d0    ! rad/s

  case ('WINF')
    t0  = 0.0d0   ! 0.d0
    dt  = 2.0d-1
    tfin= 10.d0
    NumSteps= ceiling((tfin-t0)/dt)
    Options%NewmarkDamp=0.05d0
    Omega=20.d0    ! rad/s

  case ('WING')
    t0  = 0.0
    dt  = 1.0d-2
    tfin= 10.d0
    NumSteps= ceiling((tfin-t0)/dt)
    Options%NewmarkDamp=0.0d0
    Omega=10.d0    ! rad/s

  case ('CANT')
    t0  = 0.0d0   ! 0.d0
    dt  = 1.0d-1
    tfin= 10.d0
    NumSteps= ceiling((tfin-t0)/dt)
    Options%NewmarkDamp=0.0d0
    Omega=20.d0    ! rad/s

  case ('CB2F')
    t0      = 0.d0
    dt      = 1.d-3
    tfin    = 2.0d0
    NumSteps= ceiling((tfin-t0)/dt)
    Options%NewmarkDamp=0.005d0
    Omega= 50.0d0   ! rad/s

  case ('FPAI')
    Omega=  2.d0*Pi*9.d0           !*32.d0
    t0   =-(2.d0*Pi/Omega)*250.d0  !*500.d0    ! Ramp up to t=0.
    dt   = (2.d0*Pi/Omega)/100.d0
    tfin =  2.d0
    NumSteps= ceiling((tfin-t0)/dt)
    Options%NewmarkDamp=5.d-3

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

! Harmonic load with no ramping.
  case ('CB2F')
    ForceTime=0.d0
    ForceDynAmp=ForceStatic
    do i=1,NumSteps
      ForceTime(i+1)=sin(Omega*Time(i+1))
    end do

  case ('FPAI')
     ForceTime=0.d0
     ForceDynAmp=0.d0

  case default
    ForceDynAmp = 1.d0*ForceStatic
    ForceTime   = 1.d0

! Ramped harmonic load.
    if (.false.) then
      Time0= Time(NumSteps)/2.d0
      ForceTime=0.d0
      do i=1,NumSteps
        ForceTime(i+1)=sin(Omega*Time(i+1))
        if (Time(i+1) < Time0) ForceTime(i+1)=ForceTime(i+1)*Time(i+1)/Time0
      end do
    end if

! 1-Cos load.
    if (.false.) then
     do i=1,NumSteps+1
         if ((Time(i).ge.0.d0).and.(Time(i).le.1.d1)) then
             ForceTime(i)=(1.d0-cos(Pi*dble(Time(i)-0.d0)/dble(10.d0)))/2.d0
         end if
     end do
    end if

! Ramped load.
    if (.false.) then
      ForceTime=1.d0
      do i=1,NumSteps+1
        if ((Time(i).ge.0.d0).and.(Time(i).le.10.d0)) then
            ForceTime(i)=Time(i)/10.d0
        end if
      end do
    end if

! sinusoidal load.
    if (.true.) then
     do i=1,NumSteps+1
         ForceTime(i)=sin(Pi*dble(Time(i)))
     end do
    end if
  end select

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

! Ramped harmonic motions of the base.
  select case (trim(TestCase))
  case ('WINF','WING')
    VelAmp(3)= 0*10.d0

    do i=1,NumSteps+1
      ForcedVelDot(i,:)   = VelAmp(:)*Time(i)/1.d0
    end do

  case ('CANT')
    VelAmp(3)= 0.d0    ! 10.d0

    do i=1,NumSteps
      ForcedVel(i+1,:)   = VelAmp(:)*sin(Omega*Time(i+1)-Time(1))
      ForcedVelDot(i+1,:)= VelAmp(:)*cos(Omega*Time(i+1)-Time(1))*Omega

      if (Time(i+1) < 0.d0) then
        ForcedVelDot(i+1,:)= ForcedVelDot(i+1,:)*(Time(i+1)-Time(1))/(-Time(1)) + ForcedVel(i+1,:)
        ForcedVel   (i+1,:)= ForcedVel   (i+1,:)*(Time(i+1)-Time(1))/(-Time(1))
      end if
    end do


    case ('FPAI')
    VelAmp(3)= 0.1399d0 !0.3414d0

    do i=1,NumSteps
      ForcedVel(i+1,:)   = VelAmp(:)*sin(Omega*Time(i+1)-Time(1))
      ForcedVelDot(i+1,:)= VelAmp(:)*cos(Omega*Time(i+1)-Time(1))*Omega

      if (Time(i+1) < 0.d0) then
        ForcedVelDot(i+1,:)= ForcedVelDot(i+1,:)*(Time(i+1)-Time(1))/(-Time(1)) + ForcedVel(i+1,:)
        ForcedVel   (i+1,:)= ForcedVel   (i+1,:)*(Time(i+1)-Time(1))/(-Time(1))
      end if
    end do
  end  select

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