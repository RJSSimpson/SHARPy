!-> Program.- MAIN - 06Jan2011 Updated: 22/21/2012
!
!-> Author.- Henrik Hesse  (h.hesse09@imperial.ac.uk)
!            Rafa Palacios (rpalacio@imperial.ac.uk)
!            Rob Simpson (rjs10@imperial.ac.uk) copied from main_andrea.f90
!
!-> Language.- FORTRAN90, Free format.
!
!-> Description.-
!
!   Main programme for the core routines of the multibeam assembly.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program main
 use xbeam_shared
 use xbeam_undef
 use cbeam3_solv
 use xbeam_solv
 use xbeam_perturb
 use input
 use lib_out


 implicit none

 real(8):: t0,dt                               ! Initial time and time step.
 integer:: i,j                                 ! Counter.
 integer:: NumElems,NumNodes                   ! Number of elements/nodes in the model.
 integer:: NumSteps                            ! Number of time steps.
 integer:: NumDof                              ! Number of independent degrees of freedom (2nd-order formulation).
 type(xbopts)            :: Options           ! Solution options (structure defined in xbeam_shared).
 type(xbelem),allocatable:: Elem(:)           ! Element information.
 type(xbnode),allocatable:: Node(:)           ! Nodal information.
 integer,      allocatable:: BoundConds(:)     ! =0 on free nodes; =1 on clamped nodes.
 real(8),      allocatable:: ForceStatic (:,:) ! Applied static nodal forces.
 real(8),      allocatable:: ForceDynAmp (:,:) ! Amplitude of the applied dynamic nodal forces.
 real(8),      allocatable:: ForceTime   (:)   ! Time history of the dynamic nodal forces.
 real(8),      allocatable:: ForcedVel   (:,:) ! Forced velocities at the support.
 real(8),      allocatable:: ForcedVelDot(:,:) ! Derivatives of the forced velocities at the support.
 real(8),      allocatable:: PhiNodes (:)      ! Initial twist at grid points.
 real(8),      allocatable:: InternalForces(:,:)  ! Internal force/moments at nodes.
 logical,      allocatable:: OutGrids(:)       ! Grid nodes where output is written.
 character(len=25)        :: OutFile           ! Output file.

 real(8),      allocatable:: PosIni   (:,:)    ! Initial nodal Coordinates.
 real(8),      allocatable:: PsiIni (:,:,:)    ! Initial element orientation vectors.
 real(8),      allocatable:: PosDef (:,:)      ! Current nodal position vector.
 real(8),      allocatable:: PsiDef (:,:,:)    ! Current element orientation vectors.
 real(8),      allocatable:: PosDotDef (:,:)   ! Current nodal position vector.
 real(8),      allocatable:: PsiDotDef (:,:,:) ! Current element orientation vectors.

 real(8),      allocatable:: PosPsiTime(:,:)   ! Position vector/rotation history at beam tip.
 real(8),      allocatable:: ForcesTime(:,:)   ! History of the force/moment vector at the beam root element.
 real(8),      allocatable:: VelocTime(:,:)    ! History of velocities.
 real(8),      allocatable:: Time(:)           ! Discrete time vector in the dynamic simulation.

 ! Rigid-body variables
 real(8),      allocatable:: RefVel   (:,:)    ! Velocities of reference frame at the support (rigid body).
 real(8),      allocatable:: RefVelDot(:,:)    ! Derivatives of the velocities of reference frame a.
 real(8),      allocatable:: Quat     (:)      ! Quaternions to describe propagation of reference frame a.
 real(8),      allocatable:: DynOut   (:,:)    ! Position of all nodes wrt to global frame a for each time step

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Read input data.
 call input_setup (NumElems,OutFile,Options)

 allocate (Elem(NumElems))
 call input_elem (NumElems,NumNodes,Elem)

 allocate(PosIni     (NumNodes,3)); PosIni     = 0.d0
 allocate(ForceStatic(NumNodes,6)); ForceStatic= 0.d0
 allocate(PhiNodes   (NumNodes));   PhiNodes   = 0.d0
 allocate(BoundConds (NumNodes));   BoundConds = 0
 call input_node (NumNodes,Elem,BoundConds,PosIni,ForceStatic,PhiNodes)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Open main output file and select grid points where output will be written.
 open (unit=12,file=OutFile(1:11)//'.mrb',status='replace')
 allocate(OutGrids(NumNodes))
 OutGrids          = .false.
 OutGrids(NumNodes)= .true.

 call out_title (12,'GLOBAL CONSTANTS IN THE MODEL:')
 write (12,'(14X,A,I12)')    'Number of Beam DOFs:    ', 6
 call out_title (12,'OUTPUT OPTIONS:')
 write (12,'(14X,A,I12)')    'Number of Output Nodes: ', 1
 write (12,'(14X,A,I12)')    'Print Displacements:    ', 1
 write (12,'(14X,A,I12)')    'Print Velocities:       ', 1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Compute initial (undeformed) geometry.
 allocate(PsiIni(NumElems,MaxElNod,3)); PsiIni=0.d0
 call xbeam_undef_geom (Elem,PosIni,PhiNodes,PsiIni,Options)

! Store undeformed geometry in external text file.
 open (unit=11,file=OutFile(1:11)//'_und.txt',status='replace')
   call output_elems (11,Elem,PosIni,PsiIni)
 close (11)

! Identify nodal degrees of freedom.
 allocate (Node(NumNodes))
 call xbeam_undef_dofs (Elem,BoundConds,Node,NumDof)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Static solution.

! Allocate memory for problem variables.
 allocate (PosDef(NumNodes,3));          PosDef= PosIni
 allocate (PsiDef(NumElems,MaxElNod,3)); PsiDef= PsiIni
 allocate (InternalForces(NumNodes,6));  InternalForces= 0.d0

 select case (Options%Solution)
 case (102,302)
   call cbeam3_solv_linstatic (NumDof,Elem,Node,ForceStatic,PosIni,PsiIni, &
&                              PosDef,PsiDef,Options)

 case (112,142,312,322)
   call cbeam3_solv_nlnstatic (NumDof,Elem,Node,ForceStatic,PosIni,PsiIni, &
&                              PosDef,PsiDef,Options)

 case default
   print *, 'No static solution'
 end select

! Store results in text file.
 open (unit=11,file=OutFile(1:11)//'_def.txt',status='replace')
 call output_elems (11,Elem,PosDef,PsiDef)
 close (11)
! write (*,'(1P6E12.4)') PosDef(NumNodes,:),PsiDef(NumElems,2,:)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Vibration analysis.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  allocate (ForcedVel(1,6)); ForcedVel   = 0.d0

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! CBEAM3: Tangent linear vibration analysis (around the current deformed beam).
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  select case (Options%Solution)
  case (142)
     call cbeam3_solv_modal (12,NumDof,Elem,Node,ForcedVel,PosIni,PsiIni,     &
&                            PosDef,PsiDef,Options)
  end select
  deallocate(ForcedVel)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Dynamic solution.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if ((Options%Solution.ge.200).and.(Options%Solution.le.952)) then

    ! Input data for transient dynamic solution.
    call input_dynsetup (NumSteps,t0,dt,Options)
    allocate (Time(NumSteps+1))
    do i=1,NumSteps+1
      Time(i)=t0+dt*dble(i-1)
    end do

    ! Force or velocity input.
    allocate (ForceTime   (NumSteps+1));   ForceTime   = 0.d0
    allocate (ForceDynAmp (NumNodes,6));   ForceDynAmp = 0.d0
    allocate (ForcedVel   (NumSteps+1,6)); ForcedVel   = 0.d0
    allocate (ForcedVelDot(NumSteps+1,6)); ForcedVelDot= 0.d0

    call input_dynforce  (NumNodes,Time,ForceStatic,ForceDynAmp,ForceTime)
    call input_forcedvel (NumNodes,Time,ForcedVel,ForcedVelDot)

    open (unit=11,file=OutFile(1:11)//'_force.txt',status='replace')
      do i=1,NumSteps
        write (11,'(1X,1P14E13.5)') Time(i), ForceTime(i), ForcedVel(i,:), ForcedVelDot(i,:)
      end do
    close (11)

    allocate (PosDotDef(NumNodes,3));           PosDotDef= 0.d0
    allocate (PsiDotDef(NumElems,MaxElNod,3));  PsiDotDef= 0.d0
    allocate (PosPsiTime(NumSteps+1,6));        PosPsiTime=0.d0
    allocate (VelocTime(NumSteps+1,NumNodes));  VelocTime= 0.d0
    allocate (DynOut((NumSteps+1)*NumNodes,3)); DynOut=0.d0

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Structural dynamic analysis only
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    select case (Options%Solution)
    ! CBEAM3: Tangent linear dynamic (around the current deformed beam).
    case (202,322)
      call cbeam3_solv_lindyn (12,NumDof,Time,Elem,Node,ForceStatic*0.d0,ForceDynAmp,         &
&                              ForceTime,ForcedVel,ForcedVelDot,PosIni,PsiIni,                &
&                              PosDef,PsiDef,PosDotDef,PsiDotDef,PosPsiTime,VelocTime,DynOut, &
&                              OutGrids,Options)

    ! CBEAM3: Linear static + linear dynamic.
    case (302)
!      PosDef=PosIni
!      PsiDef=PsiIni
      call cbeam3_solv_lindyn (12,NumDof,Time,Elem,Node,ForceStatic,ForceDynAmp*0.d0,             &
&                              ForceTime*0.d0,ForcedVel*0.d0,ForcedVelDot*0.d0,PosIni,PsiIni,     &
&                              PosDef,PsiDef,PosDotDef,PsiDotDef,PosPsiTime,VelocTime,DynOut,     &
&                              OutGrids,Options)

    ! CBEAM3: Nonlinear dynamic (around the current deformed beam).
    case (212,312)
      call cbeam3_solv_nlndyn (12,NumDof,Time,Elem,Node,ForceStatic,ForceDynAmp,             &
&                              ForceTime,ForcedVel,ForcedVelDot,PosIni,PsiIni,                    &
&                              PosDef,PsiDef,PosDotDef,PsiDotDef,PosPsiTime,VelocTime,DynOut,     &
&                              OutGrids,Options)

    end select

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Coupled analysis of structural and rigid-body dynamics.
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if ((Options%Solution.ge.900).and.(Options%Solution.le.952)) then

      ! Initialize
      allocate (RefVel   (NumSteps+1,6));      RefVel   =ForcedVel;       ! RefVel(1,5)=0.5d0
      allocate (RefVelDot(NumSteps+1,6));      RefVelDot=ForcedVelDot
      allocate (Quat     (4));                 Quat     =(/1.d0,0.d0,0.d0,0.d0/)

      select case (Options%Solution)
      ! linear rigid body dynamics only
      case (900)
        call xbeam_solv_rigidlndyn (12,NumDof,Time,Elem,Node,ForceStatic*0.d0,ForceDynAmp,      &
&                                   ForceTime,RefVel,RefVelDot,Quat,PosIni,PsiIni,              &
&                                   PosDef,PsiDef,PosDotDef,PsiDotDef,Options)

      ! nonlinear rigid body dynamics only
      case (910)
        call xbeam_solv_rigidnlndyn (12,NumDof,Time,Elem,Node,ForceStatic*0.d0,ForceDynAmp,     &
&                                    ForceTime,RefVel,RefVelDot,Quat,PosIni,PsiIni,             &
&                                    PosDef,PsiDef,PosDotDef,PsiDotDef,DynOut,Options)

      ! coupled linear rigid body dynamics and linear structural dynamics
      case (902)
        call xbeam_solv_coupledlindyn (12,NumDof,Time,Elem,Node,ForceStatic*0.d0,ForceDynAmp,   &
&                                      ForceTime,RefVel,RefVelDot,Quat,PosIni,PsiIni,           &
&                                      PosDef,PsiDef,PosDotDef,PsiDotDef,DynOut,Options)

      ! coupled nonlinear rigid body dynamics and nonlinear structural dynamics
      case (912)
        call xbeam_solv_couplednlndyn (12,NumDof,Time,Elem,Node,ForceStatic*0.d0,ForceDynAmp,   &
&                                      ForceTime,RefVel,RefVelDot,Quat,PosIni,PsiIni,           &
&                                      PosDef,PsiDef,PosDotDef,PsiDotDef,DynOut,Options)

      ! static and then nonlinear rigid body dynamics and nonlinear structural dynamics
      case (922)
        call cbeam3_solv_nlnstatic (NumDof,Elem,Node,ForceStatic,PosIni,PsiIni,PosDef,PsiDef,Options)

        PosIni = PosDef
        PsiIni = PsiDef

        call xbeam_solv_couplednlndyn (12,NumDof,Time,Elem,Node,ForceStatic,ForceDynAmp,   &
&                                      ForceTime,RefVel,RefVelDot,Quat,PosIni,PsiIni,      &
&                                      PosDef,PsiDef,PosDotDef,PsiDotDef,DynOut,Options)

      ! coupled linear-elastic dynamics
      case (952)
        call xbeam_perturb_solv (12,NumDof,Time,Elem,Node,ForceStatic*0.d0,ForceDynAmp,    &
&                                ForceTime,RefVel,RefVelDot,Quat,PosIni,PsiIni,            &
&                                PosDef,PsiDef,PosDotDef,PsiDotDef,DynOut,Options)

      end select

      ! Store rigid body velocities and accelerations of global reference frame
      open (unit=11,file=OutFile(1:11)//'_rigid.txt',status='replace')
        do i=1,NumSteps+1;  write (11,'(1X,1P7E15.6)') Time(i),RefVel(i,:);  end do
      close (11)

      open (unit=11,file=OutFile(1:11)//'_vreldot.txt',status='replace')
        do i=1,NumSteps+1;  write (11,'(1X,1P7E15.6)') Time(i),RefVelDot(i,:);  end do
      close (11)

    end if      ! Coupled analysis


! Store results for general dynamic analysis.
    ! Dynamic response of specified node with respect to global frame a.
    open (unit=11,file=OutFile(1:11)//'_dyn.txt',status='replace')
      do i=1,NumSteps-1;  write (11,'(1X,1P7E15.6)') Time(i),PosPsiTime(i,:);  end do
    close (11)

    open (unit=11,file=OutFile(1:11)//'_vel.txt',status='replace')
      do i=1,NumSteps-1;
        if (Time(i).ge.0.d0) then
          do j=1,NumNodes
            write (11,'(1X,2I8,1PE15.6)') i,j,VelocTime(i,j)
          end do
        end if
      end do
    close (11)

    ! Position vector of every node wrt global frame a at each time step.
    open (unit=11,file=OutFile(1:11)//'_shape.txt',status='replace')
      do i=1,NumSteps-1;
        do j=1,NumNodes
           write (11,'(1X,1P7E15.6)') Time(i), DynOut((i-1)*NumNodes+j,:);
        end do
      end do
    close (11)

    ! Screen output position and CRV of specified node at last time step.
    write (*,'(1P6E12.4)') PosDef(NumNodes,:),PsiDef(NumElems,2,:)

  end if     ! Dynamic analysis

end program main
