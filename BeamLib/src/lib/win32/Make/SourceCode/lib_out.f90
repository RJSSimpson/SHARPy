!->Copyright by Imperial College London, Aeronautics, 2009.
!
!->Module LIB_OUT. Rafa Palacios. 2Oct2009
!
!->Description.-
!
!  This module includes different elementary tools.
!
!->Subroutines:
!
!   out_outgrid:   Write informaiton in output grid points.
!   out_newline:   Write a blank line in a file.
!   out_underline: Echo of computations in output file.
!   out_comment:   Write commented line in output file.
!   out_header:    Write header in columns.
!   out_title:     Write main title.
!   out_time:      Create line with current time step.
!   out_final:     Final output comments.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module lib_out
 implicit none

! Output options (with default values).
 type outopts
  logical:: PrintDispl    = .true.        ! Print nodal displacements.
  logical:: PrintInflow   = .false.       ! Print inflow at elements.
  logical:: PrintVeloc    = .false.       ! Print nodal velocities.
  logical:: PrintExtForces= .false.       ! Print external forces.
  logical:: PrintIntForces= .false.       ! Print internal forces.
  logical:: PrintResForces= .false.       ! Print internal forces.
  real(8):: ScaleForce    = 1.d0          ! Scale factor for forces.
  real(8):: ScaleLength   = 1.d0          ! Scale factor for lengths.
  real(8):: ScaleMass     = 1.d0          ! Scale factor for mass.
  real(8):: ScaleMoment   = 1.d0          ! Scale factor for moments.
  real(8):: ScaleTime     = 1.d0          ! Scale factor for time.
 end type


 contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!->Subroutine OUT_OUTGRID.
!
!->Description.-
!
!     Write time evolution at selected nodes.
!
!-> Remarks.-
!
!   1) Time step information is given in the global coordinate system.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine out_outgrid (iOut,Node_or_Elem,Options,NumMembers,NumElems,NumDofs, &
&                        OutGrids,DISPL,VELOC,FORCEI,FORCEE)

!--> I/O variables.
  integer,         intent(in):: iOut         ! Output file.
  character(len=4),intent(in):: Node_or_Elem !='NODE': Output in Nodes
                                             !='ELEM': Output in Elements.
  type(outopts),   intent(in):: Options      ! Output parameters.
  integer,         intent(in):: NumMembers   ! Number of members in the model.
  integer,         intent(in):: NumElems(:)  ! Number of elements in each member.
  integer,         intent(in):: NumDofs      ! Number of degrees of freedom on each node.
  logical,         intent(in):: OutGrids(:)  ! =T for grids where output is written.
  real(8),optional,intent(in):: DISPL (:,:)  ! Displacements.
  real(8),optional,intent(in):: VELOC (:,:)  ! Velocities.
  real(8),optional,intent(in):: FORCEI(:,:)  ! Internal forces.
  real(8),optional,intent(in):: FORCEE(:,:)  ! External (applied) forces.

!-> Local variables.
  integer:: AddVal      ! Either 1 ('NODE') or 0 ('ELEM').
  integer:: j           ! Counters.
  integer:: iMember     ! Counter on the members.
  integer:: iNode       ! Counter on the nodes.
  integer:: jNode       ! Index of the node in the overall model.
  real(8):: OutForce(3) ! Output forces.
  real(8),         allocatable:: ResForce(:) ! Resultant forces on members.
  character(len=8),allocatable:: Labels  (:) ! Labels in the headers.
!
  allocate(Labels(2+NumDofs))

  select case (Node_or_Elem)
  case ('NODE')
    AddVal=1
  case ('ELEM')
    AddVal=0
  end select

!-> Write time evolution of displacements at Output nodes.
  if (Options%PrintDispl) then

! Write labels.
    Labels(1:8)=(/' MEMBER ','  '//Node_or_Elem//'  ', &
&                 '   D1   ','   D2   ','   D3   ',     &
&                 '   R1   ','   R2   ','   R3   '/)
    do j=1,NumDofs-6
      write (Labels(8+j),'(A,I2.2,A)') '   q',j,'  '
    end do

    call out_header (iOut,NumDofs+2,Labels)
!
! Loop in in the set of nodes defined in OUTGRID.
!
    jNode=0
    do iMember=1,NumMembers
      do iNode=1,NumElems(iMember)+AddVal
        jNode=jNode+1
        if ((Node_or_Elem.eq.'ELEM').or.OutGrids(jNode)) then
!
! Write data in output file.
!
          write (iOut,'(3X,2(4X,I8,4X),$)') iMember, iNode

          do j=1,3
             write (iOut,'(1PE16.8,$)') Displ(jNode,j)*Options%ScaleLength
          end do

          do j=4,NumDofs
            write (iOut,'(1PE16.8,$)') Displ(jNode,j)
          end do

          write (iOut,'(A1)') ' '

        end if
      end do
    end do
  end if

!
!-> Write time evolution of velocities at Output nodes.
!
  if (present(Veloc).and.Options%PrintVeloc) then
!
! Write labels.
!
    Labels(1:8)=(/' MEMBER ','  '//Node_or_Elem//'  ', &
&                 '   V1   ','   V2   ','   V3   ',     &
&                 '   W1   ','   W2   ','   W3   '/)
    do j=1,NumDofs-6
      write (Labels(8+j),'(A,I2.2,A)') ' qdot',j,' '
    end do
!
    call out_header (iOut,NumDofs+2,Labels)
!
! Loop in in the set of nodes defined in OUTGRID.
!
    jNode=0
    do iMember=1,NumMembers
      do iNode=1,NumElems(iMember)+AddVal
        jNode=jNode+1
        if ((Node_or_Elem.eq.'ELEM').or.OutGrids(jNode)) then
!
! Write data in output file.
!
          write (iOut,'(3X,2(4X,I8,4X),$)') iMember, iNode

          do j=1,3
             write (iOut,'(1PE16.8,$)') Veloc(jNode,j)*Options%ScaleLength
          end do

          do j=4,NumDofs
            write (iOut,'(1PE16.8,$)') Veloc(jNode,j)
          end do

          write (iOut,'(A1)') ' '

        end if
      end do
    end do
  end if

!
!-> Write internal Forces on Output nodes.
  if (present(ForceI).and.Options%PrintIntForces) then

! Write labels.

    Labels(1:8)=(/' MEMBER ','  '//Node_or_Elem//'  ', &
&                 '  FI1   ','  FI2   ','  FI3   ',     &
&                 '  MI1   ','  MI2   ','  MI3   '/)

    do j=1,NumDofs-6
      write (Labels(j+8),'(A,I2.2,A)') '  QI',j,'  '
    end do

    call out_header (iOut,NumDofs+2,Labels)
!
! Loop in in the set of nodes defined in OUTGRID.
    jNode=0
    do iMember=1,NumMembers
      do iNode=1,NumElems(iMember)+AddVal
        jNode=jNode+1
        if ((Node_or_Elem.eq.'ELEM').or.OutGrids(jNode)) then
!
! Write data in output file.
          write (iOut,'(3X,2(4X,I8,4X),$)') iMember, iNode

          do j=1,3
            write (iOut,'(1PE16.8,$)') ForceI(jNode,j)*Options%ScaleForce
          end do

          do j=4,6
            write (iOut,'(1PE16.8,$)') ForceI(jNode,j)*Options%ScaleMoment
          end do

          do j=7,NumDofs
            write (iOut,'(1PE16.8,$)') ForceI(jNode,j)
          end do

          write (iOut,'(A1)') ' '
!
        end if
      end do
    end do
  end if
!

!
!-> Write applied Forces on Output nodes at current time step.
!
  if (present(ForceE).and.Options%PrintExtForces) then

! Write labels.

    Labels(1:8)=(/' MEMBER ','  '//Node_or_Elem//'  ', &
&                 '  FE1   ','  FE2   ','  FE3   ',     &
&                 '  ME1   ','  ME2   ','  ME3   '/)
    do j=1,NumDofs-6
      write (Labels(j+8),'(A,I2.2,A)') '  QE',j,'  '
    end do

    call out_header (iOut,NumDofs+2,Labels)
!
! Loop in in the set of nodes defined in OUTGRID.
!
    jNode=0
    do iMember=1,NumMembers
      do iNode=1,NumElems(iMember)+AddVal
        jNode=jNode+1
        if ((Node_or_Elem.eq.'ELEM').or.OutGrids(jNode)) then
!
! Write data in output file.
          write (iOut,'(3X,2(4X,I8,4X),$)') iMember, iNode

          OutForce=ForceE(jNode,1:3)
          !OutForce=matmul(transpose(MemRotMat(iMember,:,:)),ForceE(jNode,1:3))
          do j=1,3
            write (iOut,'(1PE16.8,$)') OutForce(j)*Options%ScaleForce
          end do

          OutForce=ForceE(jNode,4:6)
          !OutForce=matmul(transpose(MemRotMat(iMember,:,:)),ForceE(jNode,4:6))
          do j=1,3
            write (iOut,'(1PE16.8,$)') OutForce(j)*Options%ScaleMoment
          end do

          do j=7,NumDofs
            write (iOut,'(1PE16.8,$)') ForceE(jNode,j)
          end do

          write (iOut,'(A1)') ' '
        end if
      end do
    end do
  end if

!
!-> Write resultant forces on Members at current time step.
!
  if (present(ForceE).and.Options%PrintResForces) then
    allocate(ResForce(NumDofs))

! Write labels.
    Labels(1:8)=(/' MEMBER ','        ','  FR1   ','  FR2   ',&
&                 '  FR3   ','  MR1   ','  MR2   ','  MR3   '/)
    do j=1,NumDofs-6
      write (Labels(j+8),'(A,I2.2,A)') '  QR',j,'  '
    end do

    call out_header (iOut,NumDofs+2,Labels)
!
! Loop in the members and in all nodes for each member.
!
    jNode=0
    do iMember=1,NumMembers
      ResForce=0.d0
      do iNode=1,NumElems(iMember)+AddVal
        jNode=jNode+1
        ResForce=ResForce+ForceE(jNode,1:NumDofs)
      end do
!
! Write data in output file.
!
      write (iOut,'(3X,4X,I8,4X,A12,$)') iMember, ' '

      OutForce=ResForce(1:3)
      !OutForce=matmul(transpose(MemRotMat(iMember,:,:)),ResForce(1:3))
      do j=1,3
        write (iOut,'(1PE16.8,$)') OutForce(j)*Options%ScaleForce
      end do

      OutForce=ResForce(4:6)
      !OutForce=matmul(transpose(MemRotMat(iMember,:,:)),ResForce(4:6))
      do j=1,3
        write (iOut,'(1PE16.8,$)') OutForce(j)*Options%ScaleMoment
      end do

      do j=7,NumDofs
        write (iOut,'(1PE16.8,$)') ResForce(j)
      end do

      write (iOut,'(A1)') ' '
    end do

    deallocate(ResForce)
  end if

  deallocate (Labels)
  return
 end subroutine out_outgrid




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!->Subroutine OUT_NEWLINE.
!
!->Description.-
!
!   Output a new line to a file connected to unit iOut.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine out_newline(iOut)
!
  integer,intent(in)::iOut
  write(iOut,'(1x,/)')
  return
 end subroutine out_newline



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!->Subroutine OUT_UNDERLINE.
!
!->Description.-
!
!   Output a line with a bunch of '=' to iOut
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine out_underline(iOut)
!  
  integer,intent(in)::iOut
  write(iOut,'(1x,50("="))')
!  
  return
 end subroutine out_underline






!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!->Subroutine OUT_COMMENT.
!
!->Description.-
!
!   Write a commented line in the given output file.
!   Commented lines start with the $ symbol.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine out_comment (iOut,Text,Line,Date)
!
!-> I/O Variables.
!
  integer, intent(in) ::           iOut     ! Logical unit of the file.
  character(len=*), intent(in) ::  Text     ! Text to be written.
  logical, intent(in), optional :: Line     ! Option to write a line of '$'.
  logical, intent(in), optional :: Date     ! Option to write the date.
!
  write (iOut,'(A)') '$'
  write (iOut,'(A)') '$ '//Text
  write (iOut,'(A)') '$'
  if (present(Line)) then
    if (Line) write (iOut,'(A80)') '$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'// &
&                                   '$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'
  end if
!
  return
 end subroutine out_comment
 
 

 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!-> Subroutine OUT_HEADER.
!
!-> Description.-
!
!     Writes headers for the columns in the output file iOut.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine out_header (iOut, NumCols, Headers)
!
!-> Input Variables.
!
  integer      iOut        ! Output file.
  integer      NumCols     ! Number of columns.
  character*8  Headers(*)  ! Headers for the columns.
!
!-> Local variables.
!
  integer      i
!
  write (iOut, '(/16(4X,A8,4X))') (Headers(i), i=1, NumCols)
!
  return
 end subroutine out_header


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!-> Subroutine OUT_TITLE.
!
!-> Description.-
!
!     Writes a main title in a given file iOut.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine out_title (iOut, Title)
!
!-> Input Variables.
  integer      iOut        ! Output file.
  character(*) Title       ! Title to be written.
!
!-> Local variables.
!
  integer      i
!
!
  write (iOut,'(//10X,A)') Title
  write (iOut,'(10X,100A1)') ('=',i=1,len_trim(Title)+1)
!
  return
 end subroutine out_title

 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!-> Subroutine OUT_TIME.
!
!-> Description.-
!
!     Writes a main title in a given file iOut.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine out_time (iStep,Time,Text)

! I/O Variables.
  integer,          intent(in) :: iStep  ! Current time step.
  real(8),          intent(in) :: Time   ! Current time.
  character(len=80),intent(out):: Text   ! Text to be written.

  Text=' '
  write (Text,'(A,I7,A,1PE14.6)') 'GLOBAL TIME STEP:',iStep,'. TIME:',Time

  return
 end subroutine out_time


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!-> Subroutine OUT_FINAL.
!
!-> Description.-
!
!     Manages output at final time step.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine out_final (iOut)

! I/O Variables.
  integer,intent(in)::      iOut        ! Output file.
!
! End of program.
  write (iOut,'(//10X,A)') '*** End of NLABS ***'
  close (iOut)

  write (*,'(/A/)') ' End of NLABS OK.'
  return
 end subroutine out_final
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end module lib_out
 
