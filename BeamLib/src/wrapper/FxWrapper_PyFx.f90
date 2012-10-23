
!-------------------------------------------------------------------------------
! FlexFlight project
!   A.Da Ronch, email: andreadr@liverpool.ac.uk
!
!-------------------------------------------------------------------------------
!
! DESCRIPTION
! Wrapper file to interface Python scripting and Fortran 90 code.
! The file contains various F90 routines for Python-driven analyses and 
! circumvents the issue of derived types by passing each member of the structure 
! separately as function arguments.
! 
!-------------------------------------------------------------------------------
!
! DATE     VERSION PROGRAMMER  DESCRIPTION
! 20120205 1.0     A.Da Ronch  Created/allocated memory set to zero
! 20120210 1.0     A.Da Ronch  Function wrap_cbeam3_asbly_dynamic modified
! 20120222 1.0     A.Da Ronch  Function wrap_cbeam3_solv_state2disp modified
! 20120929 1.0     A.Da Ronch  Added function wrap_fem_glob2loc_extract
! 20121011 1.1     R. Simpson  Added wrapper for F90 solvers...
!                               wrap_cbeam3_solv_nlndyn and...
!                                wrap_cbeam3_solv_nlnstatic
! 
!-------------------------------------------------------------------------------
! 
! To be done
! - 
!
!-------------------------------------------------------------------------------
module test
    
    use xbeam_shared
    use input
    use xbeam_undef
    use lib_sparse
    use cbeam3_solv
    use cbeam3_asbly
    use lib_fem
    use lib_xbeam
    use lib_out
    use xbeam_asbly
    
    implicit none
    
    contains
    
    !---------------------------------------------------------------------------
    ! Collect separate variables into xbelem fields
    ! Memory allocation is done before the call
    !---------------------------------------------------------------------------
    subroutine pack_xbelem(NumElems,Elem,NumNodes,MemNo,Conn,Master,Length,&
    &           PreCurv,Psi,Vector,Mass,Stiff,InvStiff,RBMass)
        
        integer     ,intent(in) :: NumElems
        type(xbelem),intent(out):: Elem(NumElems)
        integer     ,intent(in) :: NumNodes(NumElems)
        integer     ,intent(in) :: MemNo(NumElems)
        integer     ,intent(in) :: Conn(MaxElNod*NumElems)
        integer     ,intent(in) :: Master(MaxElNod*NumElems,2)
        real(8)     ,intent(in) :: Length(NumElems)
        real(8)     ,intent(in) :: PreCurv(3*NumElems)
        real(8)     ,intent(in) :: Psi(3*NumElems)
        real(8)     ,intent(in) :: Vector(3*NumElems)
        real(8)     ,intent(in) :: Mass(6*NumElems,6)
        real(8)     ,intent(in) :: Stiff(6*NumElems,6)
        real(8)     ,intent(in) :: InvStiff(6*NumElems,6)
        real(8)     ,intent(in) :: RBMass(MaxElNod*NumElems,6,6)
        
        integer:: i,i3,i6,iMaxElNod
        
        do i=1,NumElems
            i3 = (i-1)*3; i6 = i3*2
            iMaxElNod = (i-1)*MaxElNod
            
            Elem(i)%NumNodes = NumNodes(i)
            Elem(i)%MemNo = MemNo(i)
            Elem(i)%Conn = Conn(1+iMaxElNod:MaxElNod+iMaxElNod)
            Elem(i)%Master = Master(1+iMaxElNod:MaxElNod+iMaxElNod,1:2)
            Elem(i)%Length = Length(i)
            Elem(i)%PreCurv = PreCurv(1+i3:3+i3)
            Elem(i)%Psi = Psi(1+i3:3+i3)
            Elem(i)%Vector = Vector(1+i3:3+i3)
            Elem(i)%Mass = Mass(1+i6:6+i6,:)
            Elem(i)%Stiff = Stiff(1+i6:6+i6,:)
            Elem(i)%InvStiff = InvStiff(1+i6:6+i6,:)
            Elem(i)%RBMass = RBMass(1+(i-1)*MaxElNod:3+i*MaxElNod,:,:)
        end do
        
        return
        
    end subroutine pack_xbelem
    
    !---------------------------------------------------------------------------
    ! Collect separate variables into xbnode fields
    ! Memory allocation is done before the call
    !---------------------------------------------------------------------------
    subroutine pack_xbnode(NumNodes,Node,Master,Vdof,Fdof)
        
        integer     ,intent(in) :: NumNodes
        type(xbnode),intent(out):: Node(NumNodes)
        integer     ,intent(in) :: Master(2*NumNodes)
        integer     ,intent(in) :: Vdof(NumNodes)
        integer     ,intent(in) :: Fdof(NumNodes)
        
        integer:: i,i2
        
        do i=1,NumNodes
            i2 = (i-1)*2
            Node(i)%Master(1) = Master(1+i2)
            Node(i)%Master(2) = Master(2+i2)
            Node(i)%Vdof = Vdof(i)
            Node(i)%Fdof = Fdof(i)
        end do
        
        return
        
    end subroutine pack_xbnode
    
    !---------------------------------------------------------------------------
    ! Collect separate variables into xbopts fields
    ! Memory allocation is done before the call
    !---------------------------------------------------------------------------
    subroutine pack_xbopts(Options,FollowerForce,FollowerForceRig,PrintInfo,&
    &           OutInBframe,OutInaframe,ElemProj,MaxIterations,NumLoadSteps,&
    &           NumGauss,Solution,DeltaCurved,MinDelta,NewmarkDamp)
        
        type(xbopts),intent(out):: Options
        logical     ,intent(in) :: FollowerForce
        logical     ,intent(in) :: FollowerForceRig
        logical     ,intent(in) :: PrintInfo
        logical     ,intent(in) :: OutInBframe
        logical     ,intent(in) :: OutInaframe
        integer     ,intent(in) :: ElemProj
        integer     ,intent(in) :: MaxIterations
        integer     ,intent(in) :: NumLoadSteps
        integer     ,intent(in) :: NumGauss
        integer     ,intent(in) :: Solution
        real(8)     ,intent(in) :: DeltaCurved
        real(8)     ,intent(in) :: MinDelta
        real(8)     ,intent(in) :: NewmarkDamp
        
        Options%FollowerForce = FollowerForce
        Options%FollowerForceRig = FollowerForceRig
        Options%PrintInfo = PrintInfo
        Options%OutInBframe = OutInBframe
        Options%OutInaframe = OutInaframe
        Options%ElemProj = ElemProj
        Options%MaxIterations = MaxIterations
        Options%NumLoadSteps = NumLoadSteps
        Options%NumGauss = NumGauss
        Options%Solution = Solution
        Options%DeltaCurved = DeltaCurved
        Options%MinDelta = MinDelta
        Options%NewmarkDamp = NewmarkDamp
        
        return
        
    end subroutine pack_xbopts
    
    !---------------------------------------------------------------------------
    ! Collect separate variables into sparse fields
    ! Memory allocation is done before the call
    !---------------------------------------------------------------------------
    subroutine pack_sparse(DimSparse,IArray,JArray,Vec,Vec1)
        
        integer     ,intent(in) :: DimSparse
        integer     ,intent(in) :: IArray(DimSparse)
        integer     ,intent(in) :: JArray(DimSparse)
        real(8)     ,intent(in) :: Vec(DimSparse)
        type(sparse),intent(out):: Vec1(DimSparse)
        
        integer:: k
        
        do k=1,DimSparse
            Vec1(k)%i = IArray(k)
            Vec1(k)%j = JArray(k)
            Vec1(k)%a = Vec(k)
        end do
        
        return
        
    end subroutine pack_sparse
    
    !---------------------------------------------------------------------------
    ! Extract xbelem fields to separate variables
    ! Memory allocation is done before the call
    !---------------------------------------------------------------------------
    subroutine unpack_xbelem(NumElems,Elem,NumNodes,MemNo,Conn,Master,Length,&
    &           PreCurv,Psi,Vector,Mass,Stiff,InvStiff,RBMass)
        
        integer     ,intent(in) :: NumElems
        type(xbelem),intent(in) :: Elem(NumElems)
        integer     ,intent(out):: NumNodes(NumElems)
        integer     ,intent(out):: MemNo(NumElems)
        integer     ,intent(out):: Conn(MaxElNod*NumElems)
        integer     ,intent(out):: Master(MaxElNod*NumElems,2)
        real(8)     ,intent(out):: Length(NumElems)
        real(8)     ,intent(out):: PreCurv(3*NumElems)
        real(8)     ,intent(out):: Psi(3*NumElems)
        real(8)     ,intent(out):: Vector(3*NumElems)
        real(8)     ,intent(out):: Mass(6*NumElems,6)
        real(8)     ,intent(out):: Stiff(6*NumElems,6)
        real(8)     ,intent(out):: InvStiff(6*NumElems,6)
        real(8)     ,intent(out):: RBMass(MaxElNod*NumElems,6,6)
        
        integer:: i,i3,i6,iMaxElNod
        
        do i=1,NumElems
            i3 = (i-1)*3
            i6 = (i-1)*6
            iMaxElNod = (i-1)*MaxElNod
            
            NumNodes(i) = Elem(i)%NumNodes
            MemNo(i) = Elem(i)%MemNo
            Conn(1+iMaxElNod:MaxElNod+iMaxElNod) = Elem(i)%Conn
            Master(1+iMaxElNod:MaxElNod+iMaxElNod,:) = Elem(i)%Master
            Length(i) = Elem(i)%Length
            PreCurv(1+i3:3+i3) = Elem(i)%PreCurv
            Psi(1+i3:3+i3) = Elem(i)%Psi
            Vector(1+i3:3+i3) = Elem(i)%Vector
            Mass(1+i6:6+i6,:) = Elem(i)%Mass
            Stiff(1+i6:6+i6,:) = Elem(i)%Stiff
            InvStiff(1+i6:6+i6,:) = Elem(i)%InvStiff
            RBMass(1+(i-1)*MaxElNod:3+i*MaxElNod,:,:) = Elem(i)%RBMass
        end do
        
        return
        
    end subroutine unpack_xbelem
    
    !---------------------------------------------------------------------------
    ! Extract xbnode fields to separate variables
    ! Memory allocation is done before the call
    !---------------------------------------------------------------------------
    subroutine unpack_xbnode(NumNodes,Node,Master,Vdof,Fdof)
        
        integer     ,intent(in) :: NumNodes
        type(xbnode),intent(in) :: Node(NumNodes)
        integer     ,intent(out):: Master(2*NumNodes)
        integer     ,intent(out):: Vdof(NumNodes)
        integer     ,intent(out):: Fdof(NumNodes)
        
        integer:: i,i2
        
        do i=1,NumNodes
            i2 = (i-1)*2
            Master(1+i2:2+i2) = Node(i)%Master
            Vdof(i) = Node(i)%Vdof
            Fdof(i) = Node(i)%Fdof
        end do
        
        return
        
    end subroutine unpack_xbnode
    
    !---------------------------------------------------------------------------
    ! Extract xbopts fields to separate variables
    ! Memory allocation is done before the call
    !---------------------------------------------------------------------------
    subroutine unpack_xbopts(Options,FollowerForce,FollowerForceRig,PrintInfo,&
    &           OutInBframe,OutInaframe,ElemProj,MaxIterations,NumLoadSteps,&
    &           NumGauss,Solution,DeltaCurved,MinDelta,NewmarkDamp)
        
        type(xbopts),intent(in) :: Options
        logical     ,intent(out):: FollowerForce
        logical     ,intent(out):: FollowerForceRig
        logical     ,intent(out):: PrintInfo
        logical     ,intent(out):: OutInBframe
        logical     ,intent(out):: OutInaframe
        integer     ,intent(out):: ElemProj
        integer     ,intent(out):: MaxIterations
        integer     ,intent(out):: NumLoadSteps
        integer     ,intent(out):: NumGauss
        integer     ,intent(out):: Solution
        real(8)     ,intent(out):: DeltaCurved
        real(8)     ,intent(out):: MinDelta
        real(8)     ,intent(out):: NewmarkDamp
        
        FollowerForce = Options%FollowerForce
        FollowerForceRig = Options%FollowerForceRig
        PrintInfo = Options%PrintInfo
        OutInBframe = Options%OutInBframe
        OutInaframe = Options%OutInaframe
        ElemProj = Options%ElemProj
        MaxIterations = Options%MaxIterations
        NumLoadSteps = Options%NumLoadSteps
        NumGauss = Options%NumGauss
        Solution = Options%Solution
        DeltaCurved = Options%DeltaCurved
        MinDelta = Options%MinDelta
        NewmarkDamp = Options%NewmarkDamp
        
        return
        
    end subroutine unpack_xbopts
    
    !---------------------------------------------------------------------------
    ! Extract sparse fields to separate variables
    ! Memory allocation is done before the call
    !---------------------------------------------------------------------------
    subroutine unpack_sparse(DimSparse,Vec1,IArray,JArray,Vec)
        
        integer     ,intent(in) :: DimSparse
        type(sparse),intent(in) :: Vec1(DimSparse)
        integer     ,intent(out):: IArray(DimSparse)
        integer     ,intent(out):: JArray(DimSparse)
        real(8)     ,intent(out):: Vec(DimSparse)
        
        integer:: k
        
        do k=1,DimSparse
            IArray(k) = Vec1(k)%i
            JArray(k) = Vec1(k)%j
            Vec(k)    = Vec1(k)%a
        end do
        
        return
        
    end subroutine unpack_sparse
    
    !---------------------------------------------------------------------------
    ! Convert a Fortran ordered array into a matrix of dimension (nx,ny)
    ! Vector elements are of type real(8)
    !---------------------------------------------------------------------------
    subroutine vec2mat_double(vec,mat,nx,ny)
        
        real(8),intent(in) :: vec(nx*ny)
        real(8),intent(out):: mat(nx,ny)
        integer,intent(in) :: nx,ny
        
        integer:: i,j
        
        mat = 0.d0
        do j=1,ny
            do i=1,nx
                mat(i,j)=vec(i+(j-1)*nx)
            end do
        end do
        
    end subroutine vec2mat_double
    
    !---------------------------------------------------------------------------
    ! Convert a Fortran ordered array into a 3d matrix of dimension (n1,n2,n3)
    ! Vector elements are of type real(8)
    !---------------------------------------------------------------------------
    subroutine vec2mat3d_double(Array,Mat,n1,n2,n3)
        
        integer,intent(in)   :: n1,n2,n3
        real(8),intent(in)   :: Array(n1*n2*n3)
        real(8),intent(inout):: Mat(n1,n2,n3)
        
        real(8),allocatable:: temp(:,:)
        integer:: i,i1,i2
        
        allocate(temp(n2,n3)); temp = 0.d0
        
        i1 = 0
        i2 = n1*(n2*n3-1)
        do i=1,n1
            i1 = i1+1
            i2 = i2+1
            call vec2mat_double(Array(i1:i2:n1), temp, n2, n3)
            Mat(i,:,:) = temp
        end do
        
        deallocate(temp)
        
        return
        
    end subroutine vec2mat3d_double
    
    !---------------------------------------------------------------------------
    ! Convert a Fortran ordered array into a matrix of dimension (nx,ny)
    ! Vector elements are of type int
    !---------------------------------------------------------------------------
    subroutine vec2mat_int(vec,mat,nx,ny)
        
        integer,intent(in) :: vec(nx*ny)
        integer,intent(out):: mat(nx,ny)
        integer,intent(in) :: nx,ny
        
        integer:: i,j
        
        mat=0
        do j=1,ny
            do i=1,nx
                mat(i,j)=vec(i+(j-1)*nx)
            end do
        end do
        
    end subroutine vec2mat_int
    
    !---------------------------------------------------------------------------
    ! Convert a Fortran ordered array into a matrix of dimension (nx,ny)
    ! Vector elements are of type real(8)
    !---------------------------------------------------------------------------
    subroutine mat2vec_double(mat,vec,nx,ny)
        
        real(8),intent(in) :: mat(nx,ny)
        real(8),intent(out):: vec(nx*ny)
        integer,intent(in) :: nx,ny
        
        integer:: i,j
        
        vec=0.d0
        do j=1,ny
            do i=1,nx
                vec(i+(j-1)*nx)=mat(i,j)
            end do
        end do
        
    end subroutine mat2vec_double
    
    !---------------------------------------------------------------------------
    ! Convert a 3d matrix of dimension (n1,n2,n3) into a Fortran ordered array
    ! Vector elements are of type real(8)
    !---------------------------------------------------------------------------
    subroutine mat2vec3d_double(Mat,Array,n1,n2,n3)
        
        integer,intent(in)   :: n1,n2,n3
        real(8),intent(in)   :: Mat(n1,n2,n3)
        real(8),intent(inout):: Array(n1*n2*n3)
        
        real(8),allocatable:: temp(:,:)
        integer:: i,i1,i2
        
        allocate(temp(n2,n3)); temp = 0.d0
        
        i1 = 0
        i2 = n1*(n2*n3-1)
        do i=1,n1
            i1 = i1+1
            i2 = i2+1
            temp = Mat(i,:,:)
            call mat2vec_double(temp, Array(i1:i2:n1), n2, n3)
        end do
        
        deallocate(temp)
        
        return
        
    end subroutine mat2vec3d_double
    
    !---------------------------------------------------------------------------
    ! Convert a Fortran ordered array into a matrix of dimension (nx,ny)
    ! Vector elements are of type int
    !---------------------------------------------------------------------------
    subroutine mat2vec_int(mat,vec,nx,ny)
        
        integer,intent(in) :: mat(nx,ny)
        integer,intent(out):: vec(nx*ny)
        integer,intent(in) :: nx,ny
        
        integer:: i,j
        
        vec=0
        do j=1,ny
            do i=1,nx
                vec(i+(j-1)*nx)=mat(i,j)
            end do
        end do
        
    end subroutine mat2vec_int
    
    !---------------------------------------------------------------------------
    ! From sparse type to full matrix
    !---------------------------------------------------------------------------
    subroutine sparse2full_rank(DimSprMat,SprMat,n1,n2,FulMat)
        
        integer     ,intent(in) :: DimSprMat
        type(sparse),intent(in) :: SprMat(DimSprMat)
        integer     ,intent(in) :: n1
        integer     ,intent(in) :: n2
        real(8)     ,intent(out):: FulMat(n1,n2)
        
        integer:: k
        
        FulMat = 0.d0
        do k=1,DimSprMat
            FulMat(SprMat(k)%i,SprMat(k)%j) = SprMat(k)%a
        end do
        
        return
        
    end subroutine sparse2full_rank
    
    !---------------------------------------------------------------------------
    ! Form the xbelem derived type
    !---------------------------------------------------------------------------
    subroutine do_xbelem_var(NumElems,Elem,NumNodes,MemNo,Conn,Master_Array,&
    &           Length,PreCurv,Psi,Vector,Mass_Array,Stiff_Array,&
    &           InvStiff_Array,RBMass_Array)
        
        integer     ,intent(in)   :: NumElems
        type(xbelem),intent(inout):: Elem(NumElems)
        integer     ,intent(in)   :: NumNodes(NumElems)
        integer     ,intent(in)   :: MemNo(NumElems)
        integer     ,intent(in)   :: Conn(MaxElNod*NumElems)
        integer     ,intent(in)   :: Master_Array(MaxElNod*NumElems*2)
        real(8)     ,intent(in)   :: Length(NumElems)
        real(8)     ,intent(in)   :: PreCurv(3*NumElems)
        real(8)     ,intent(in)   :: Psi(3*NumElems)
        real(8)     ,intent(in)   :: Vector(3*NumElems)
        real(8)     ,intent(in)   :: Mass_Array(6*NumElems*6)
        real(8)     ,intent(in)   :: Stiff_Array(6*NumElems*6)
        real(8)     ,intent(in)   :: InvStiff_Array(6*NumElems*6)
        real(8)     ,intent(in)   :: RBMass_Array(MaxElNod*NumElems*6*6)
        
        integer,allocatable:: Master(:,:)
        real(8),allocatable:: Mass(:,:)
        real(8),allocatable:: Stiff(:,:)
        real(8),allocatable:: InvStiff(:,:)
        real(8),allocatable:: RBMass(:,:,:)
        
        allocate(Master(MaxElNod*NumElems,2)); Master = 0
        allocate(Mass(6*NumElems,6)); Mass = 0.d0
        allocate(Stiff(6*NumElems,6)); Stiff = 0.d0
        allocate(InvStiff(6*NumElems,6)); InvStiff = 0.d0
        allocate(RBMass(MaxElNod*NumElems,6,6)); RBMass = 0.d0
        
        ! From one-dim array to multi-dim entity
        call vec2mat_int(Master_Array,Master,MaxElNod*NumElems,2)
        call vec2mat_double(Mass_Array,Mass,6*NumElems,6)
        call vec2mat_double(Stiff_Array,Stiff,6*NumElems,6)
        call vec2mat_double(InvStiff_Array,InvStiff,6*NumElems,6)
        call vec2mat3d_double(RBMass_Array,RBMass,MaxElNod*NumElems,6,6)
        
        ! Pack xbelem derived type
        call pack_xbelem(NumElems,Elem,NumNodes,MemNo,Conn,Master,Length,&
        &       PreCurv,Psi,Vector,Mass,Stiff,InvStiff,RBMass)
        
        deallocate(Master)
        deallocate(Mass)
        deallocate(Stiff)
        deallocate(InvStiff)
        deallocate(RBMass)
        
        return
        
    end subroutine do_xbelem_var
    
    !---------------------------------------------------------------------------
    ! Extract one-dim arrays from the xbelem derived type
    !---------------------------------------------------------------------------
    subroutine undo_xbelem_var(NumElems,Elem,NumNodes,MemNo,Conn,&
    &           Master_Array,Length,PreCurv,Psi,Vector,Mass_Array,&
    &           Stiff_Array,InvStiff_Array,RBMass_Array)
        
        integer     ,intent(in)   :: NumElems
        type(xbelem),intent(in)   :: Elem(NumElems)
        integer     ,intent(inout):: NumNodes(NumElems)
        integer     ,intent(inout):: MemNo(NumElems)
        integer     ,intent(inout):: Conn(MaxElNod*NumElems)
        integer     ,intent(inout):: Master_Array(MaxElNod*NumElems*2)
        real(8)     ,intent(inout):: Length(NumElems)
        real(8)     ,intent(inout):: PreCurv(3*NumElems)
        real(8)     ,intent(inout):: Psi(3*NumElems)
        real(8)     ,intent(inout):: Vector(3*NumElems)
        real(8)     ,intent(inout):: Mass_Array(6*NumElems*6)
        real(8)     ,intent(inout):: Stiff_Array(6*NumElems*6)
        real(8)     ,intent(inout):: InvStiff_Array(6*NumElems*6)
        real(8)     ,intent(inout):: RBMass_Array(MaxElNod*NumElems*6*6)
        
        integer,allocatable:: Master(:,:)
        real(8),allocatable:: Mass(:,:)
        real(8),allocatable:: Stiff(:,:)
        real(8),allocatable:: InvStiff(:,:)
        real(8),allocatable:: RBMass(:,:,:)
        
        allocate(Master(MaxElNod*NumElems,2)); Master = 0
        allocate(Mass(6*NumElems,6)); Mass = 0.d0
        allocate(Stiff(6*NumElems,6)); Stiff = 0.d0
        allocate(InvStiff(6*NumElems,6)); InvStiff = 0.d0
        allocate(RBMass(MaxElNod*NumElems,6,6)); RBMass = 0.d0
        
        ! Unpack xbelem derived type
        call unpack_xbelem(NumElems,Elem,NumNodes,MemNo,Conn,Master,Length,&
        &       PreCurv,Psi,Vector,Mass,Stiff,InvStiff,RBMass)
        
        ! From multi-dim entity to one-dim array
        call mat2vec_int(Master,Master_Array,MaxElNod*NumElems,2)
        call mat2vec_double(Mass,Mass_Array,6*NumElems,6)
        call mat2vec_double(Stiff,Stiff_Array,6*NumElems,6)
        call mat2vec_double(InvStiff,InvStiff_Array,6*NumElems,6)
        call mat2vec3d_double(RBMass,RBMass_Array,MaxElNod*NumElems,6,6)
        
        deallocate(Master)
        deallocate(Mass)
        deallocate(Stiff)
        deallocate(InvStiff)
        deallocate(RBMass)
        
        return
        
    end subroutine undo_xbelem_var
    
    !---------------------------------------------------------------------------
    ! Extract forces on unconstrained nodes and return one-dim array
    !---------------------------------------------------------------------------
    subroutine wrap_fem_m2v(N1,N2,Matrix_1d,N3,Vector,FilterIN)
        
        integer,intent(in) :: N1
        integer,intent(in) :: N2
        real(8),intent(in) :: Matrix_1d(N1*N2)
        integer,intent(in) :: N3
        real(8),intent(out):: Vector(N3)
        integer,intent(in) :: FilterIN(N1)
        
        real(8),allocatable:: Matrix(:,:)
        
        allocate(Matrix(N1,N2)); Matrix = 0.d0
        
        ! Convert PY data to F90 data
        call vec2mat_double(Matrix_1d,Matrix,N1,N2)
        
        ! Convert matrix to vector array
        Vector = fem_m2v(Matrix,N3,Filter=FilterIN)
        
        deallocate(Matrix)
        
        return
        
    end subroutine wrap_fem_m2v
    
    !---------------------------------------------------------------------------
    ! Extract forces on unconstrained nodes and return one-dim array
    ! No filter option
    ! Not able to get the optional argument to work properly in one single 
    ! routine
    !---------------------------------------------------------------------------
    subroutine wrap_fem_m2v_nofilter(N1,N2,Matrix_1d,N3,Vector)
        
        integer,intent(in) :: N1
        integer,intent(in) :: N2
        real(8),intent(in) :: Matrix_1d(N1*N2)
        integer,intent(in) :: N3
        real(8),intent(out):: Vector(N3)
        
        real(8),allocatable:: Matrix(:,:)
        
        allocate(Matrix(N1,N2)); Matrix = 0.d0
        
        ! Convert PY data to F90 data
        call vec2mat_double(Matrix_1d,Matrix,N1,N2)
        
        ! Convert matrix to vector array
        Vector = fem_m2v(Matrix,N3)
        
        deallocate(Matrix)
        
        return
        
    end subroutine wrap_fem_m2v_nofilter
    
    !---------------------------------------------------------------------------
    ! input_setup wrapper
    !---------------------------------------------------------------------------
    subroutine wrap_input_setup(NumElems,OutFile)
        
        integer          ,intent(in)   :: NumElems
        character(len=25),intent(inout):: OutFile
        
        type(xbopts):: Options
        
        ! Setup testcase
        call input_setup(NumElems,OutFile,Options)
        
        return
        
    end subroutine wrap_input_setup
    
    !---------------------------------------------------------------------------
    ! input_elem wrapper
    !---------------------------------------------------------------------------
    subroutine wrap_input_elem(NumElems,NumNodes_tot,NumNodes,MemNo,Conn,&
    &           Master_Array,Length,PreCurv,Psi,Vector,Mass_Array,&
    &           Stiff_Array,InvStiff_Array,RBMass_Array)
        
        integer,intent(in) :: NumElems
        integer,intent(out):: NumNodes_tot
        integer,intent(out):: NumNodes(NumElems)
        integer,intent(out):: MemNo(NumElems)
        integer,intent(out):: Conn(MaxElNod*NumElems)
        integer,intent(out):: Master_Array(MaxElNod*NumElems*2)
        real(8),intent(out):: Length(NumElems)
        real(8),intent(out):: PreCurv(3*NumElems)
        real(8),intent(out):: Psi(3*NumElems)
        real(8),intent(out):: Vector(3*NumElems)
        real(8),intent(out):: Mass_Array(6*NumElems*6)
        real(8),intent(out):: Stiff_Array(6*NumElems*6)
        real(8),intent(out):: InvStiff_Array(6*NumElems*6)
        real(8),intent(out):: RBMass_Array(MaxElNod*NumElems*6*6)
        
        type(xbelem),allocatable:: Elem(:)
        
        allocate(Elem(NumElems))
        
        ! Convert PY data to F90 data
        call do_xbelem_var(NumElems,Elem,NumNodes,MemNo,Conn,Master_Array,&
        &       Length,PreCurv,Psi,Vector,Mass_Array,Stiff_Array,&
        &       InvStiff_Array,RBMass_Array)
        
        ! Setup element properties
        call input_elem(NumElems,NumNodes_tot,Elem)
        
        ! Convert F90 data to PY data
        call undo_xbelem_var(NumElems,Elem,NumNodes,MemNo,Conn,Master_Array,&
        &       Length,PreCurv,Psi,Vector,Mass_Array,Stiff_Array,&
        &       InvStiff_Array,RBMass_Array)
        
        deallocate(Elem)
        
        return
        
    end subroutine wrap_input_elem
    
    !---------------------------------------------------------------------------
    ! input_node wrapper
    !---------------------------------------------------------------------------
    subroutine wrap_input_node(NumElems,NumNodes_tot,NumNodes,MemNo,Conn,&
    &           Master_Array,Length,PreCurv,Psi,Vector,Mass_Array,&
    &           Stiff_Array,InvStiff_Array,RBMass_Array,BoundConds,&
    &           PosIni_Array,ForceStatic_Array,PhiNodes)
        
        integer,intent(in)   :: NumElems
        integer,intent(in)   :: NumNodes_tot
        integer,intent(in)   :: NumNodes(NumElems)
        integer,intent(in)   :: MemNo(NumElems)
        integer,intent(in)   :: Conn(MaxElNod*NumElems)
        integer,intent(in)   :: Master_Array(MaxElNod*NumElems*2)
        real(8),intent(in)   :: Length(NumElems)
        real(8),intent(in)   :: PreCurv(3*NumElems)
        real(8),intent(in)   :: Psi(3*NumElems)
        real(8),intent(in)   :: Vector(3*NumElems)
        real(8),intent(in)   :: Mass_Array(6*NumElems*6)
        real(8),intent(in)   :: Stiff_Array(6*NumElems*6)
        real(8),intent(in)   :: InvStiff_Array(6*NumElems*6)
        real(8),intent(in)   :: RBMass_Array(MaxElNod*NumElems*6*6)
        integer,intent(inout):: BoundConds(NumNodes_tot)
        real(8),intent(inout):: PosIni_Array(NumNodes_tot*3)
        real(8),intent(inout):: ForceStatic_Array(NumNodes_tot*6)
        real(8),intent(inout):: PhiNodes(NumNodes_tot)
        
        type(xbelem),allocatable:: Elem(:)
        real(8)     ,allocatable:: PosIni(:,:)
        real(8)     ,allocatable:: ForceStatic(:,:)
        
        allocate(Elem(NumElems))
        allocate(PosIni(NumNodes_tot,3)); PosIni = 0.d0
        allocate(ForceStatic(NumNodes_tot,6)); ForceStatic = 0.d0
        
        ! Convert PY data to F90 data
        call do_xbelem_var(NumElems,Elem,NumNodes,MemNo,Conn,Master_Array,&
        &       Length,PreCurv,Psi,Vector,Mass_Array,Stiff_Array,&
        &       InvStiff_Array,RBMass_Array)

        call vec2mat_double(PosIni_Array,PosIni,NumNodes_tot,3)
        call vec2mat_double(ForceStatic_Array,ForceStatic,NumNodes_tot,6)
        
        ! Setup nodal properties
        call input_node(NumNodes_tot,Elem,BoundConds,PosIni,ForceStatic,PhiNodes)
        
        ! Convert F90 data to PY data
        call mat2vec_double(PosIni,PosIni_Array,NumNodes_tot,3)
        call mat2vec_double(ForceStatic,ForceStatic_Array,NumNodes_tot,6)
        
        deallocate(Elem)
        deallocate(PosIni)
        deallocate(ForceStatic)
        
        return
        
    end subroutine wrap_input_node
    
    !---------------------------------------------------------------------------
    ! Compute initial (undeformed) geometry
    !---------------------------------------------------------------------------
    subroutine wrap_xbeam_undef_geom(NumElems,NumNodes_tot,NumNodes,MemNo,&
    &           Conn,Master_Array,Length,PreCurv,Psi,Vector,Mass_Array,&
    &           Stiff_Array,InvStiff_Array,RBMass_Array,PosIni_Array,&
    &           PhiNodes,PsiIni_Array)
        
        integer,intent(in)   :: NumElems
        integer,intent(inout):: NumNodes_tot
        integer,intent(inout):: NumNodes(NumElems)
        integer,intent(inout):: MemNo(NumElems)
        integer,intent(inout):: Conn(MaxElNod*NumElems)
        integer,intent(inout):: Master_Array(MaxElNod*NumElems*2)
        real(8),intent(inout):: Length(NumElems)
        real(8),intent(inout):: PreCurv(3*NumElems)
        real(8),intent(inout):: Psi(3*NumElems)
        real(8),intent(inout):: Vector(3*NumElems)
        real(8),intent(inout):: Mass_Array(6*NumElems*6)
        real(8),intent(inout):: Stiff_Array(6*NumElems*6)
        real(8),intent(inout):: InvStiff_Array(6*NumElems*6)
        real(8),intent(inout):: RBMass_Array(MaxElNod*NumElems*6*6)
        real(8),intent(in)   :: PosIni_Array(NumNodes_tot*3)
        real(8),intent(in)   :: PhiNodes(NumNodes_tot)
        real(8),intent(inout):: PsiIni_Array(NumElems*MaxElNod*3)
        
        type(xbelem),allocatable:: Elem(:)
        real(8)     ,allocatable:: PosIni(:,:)
        real(8)     ,allocatable:: PsiIni(:,:,:)
        
        type(xbopts):: Options
        
        allocate(Elem(NumElems))
        allocate(PosIni(NumNodes_tot,3)); PosIni = 0.d0
        allocate(PsiIni(NumElems,MaxElNod,3)); PsiIni = 0.d0
        
        ! Convert PY data to F90 data
        call do_xbelem_var(NumElems,Elem,NumNodes,MemNo,Conn,Master_Array,&
        &       Length,PreCurv,Psi,Vector,Mass_Array,Stiff_Array,&
        &       InvStiff_Array,RBMass_Array)
        
        call vec2mat_double(PosIni_Array,PosIni,NumNodes_tot,3)
        call vec2mat3d_double(PsiIni_Array,PsiIni,NumElems,MaxElNod,3)
        
        ! Compute initial (undeformed) geometry
        call xbeam_undef_geom(Elem,PosIni,PhiNodes,PsiIni,Options)
        
        ! Convert F90 data to PY data
        call undo_xbelem_var(NumElems,Elem,NumNodes,MemNo,Conn,Master_Array,&
        &       Length,PreCurv,Psi,Vector,Mass_Array,Stiff_Array,&
        &       InvStiff_Array,RBMass_Array)
        
        call mat2vec3d_double(PsiIni,PsiIni_Array,NumElems,MaxElNod,3)
        
        deallocate(Elem)
        deallocate(PosIni)
        deallocate(PsiIni)
        
        return
        
    end subroutine wrap_xbeam_undef_geom
    
    !---------------------------------------------------------------------------
    ! Identify nodal degrees of freedom
    !---------------------------------------------------------------------------
    subroutine wrap_xbeam_undef_dofs(NumElems,NumNodes_tot,NumNodes,MemNo,&
    &           Conn,Master_Array,Length,PreCurv,Psi,Vector,Mass_Array,&
    &           Stiff_Array,InvStiff_Array,RBMass_Array,BoundConds,&
    &           Nod_Master,Nod_Vdof,Nod_Fdof,NumDof)
        
        integer,intent(in)   :: NumElems
        integer,intent(in)   :: NumNodes_tot
        integer,intent(in)   :: NumNodes(NumElems)
        integer,intent(in)   :: MemNo(NumElems)
        integer,intent(in)   :: Conn(MaxElNod*NumElems)
        integer,intent(in)   :: Master_Array(MaxElNod*NumElems*2)
        real(8),intent(in)   :: Length(NumElems)
        real(8),intent(in)   :: PreCurv(3*NumElems)
        real(8),intent(in)   :: Psi(3*NumElems)
        real(8),intent(in)   :: Vector(3*NumElems)
        real(8),intent(in)   :: Mass_Array(6*NumElems*6)
        real(8),intent(in)   :: Stiff_Array(6*NumElems*6)
        real(8),intent(in)   :: InvStiff_Array(6*NumElems*6)
        real(8),intent(in)   :: RBMass_Array(MaxElNod*NumElems*6*6)
        integer,intent(in)   :: BoundConds(NumNodes_tot)
        integer,intent(inout):: Nod_Master(2*NumNodes_tot)
        integer,intent(inout):: Nod_Vdof(NumNodes_tot)
        integer,intent(inout):: Nod_Fdof(NumNodes_tot)
        integer,intent(inout):: NumDof
        
        type(xbelem),allocatable:: Elem(:)
        type(xbnode),allocatable:: Node(:)
        
        allocate(Elem(NumElems))
        allocate(Node(NumNodes_tot))
        
        ! Convert PY data to F90 data
        call do_xbelem_var(NumElems,Elem,NumNodes,MemNo,Conn,Master_Array,&
        &       Length,PreCurv,Psi,Vector,Mass_Array,Stiff_Array,&
        &       InvStiff_Array,RBMass_Array)
        
        call pack_xbnode(NumNodes_tot,Node,Nod_Master,Nod_Vdof,Nod_Fdof)
        
        ! Identify nodal degrees of freedom
        call xbeam_undef_dofs(Elem,BoundConds,Node,NumDof)
        
        ! Convert F90 data to PY data
        call unpack_xbnode(NumNodes_tot,Node,Nod_Master,Nod_Vdof,Nod_Fdof)
        
        deallocate(Elem)
        deallocate(Node)
        
        return
        
    end subroutine wrap_xbeam_undef_dofs
    
    !---------------------------------------------------------------------------
    ! Assembly matrices for a nonlinear static problem
    !---------------------------------------------------------------------------
    subroutine wrap_cbeam3_asbly_static(NumElems,NumNodes_tot,&
    &           NumNodes,MemNo,Conn,Master_Array,Length,PreCurv,Psi,Vector,&
    &           Mass_Array,Stiff_Array,InvStiff_Array,RBMass_Array,&
    &           Nod_Master,Nod_Vdof,Nod_Fdof,&
    &           PosIni_Array,PsiIni_Array,PosDefor_Array,PsiDefor_Array,&
    &           ForceStatic_Array,DimMat,NumDof,&
    &           ks,Kglobal_Array,fs,Fglobal_Array,Qglobal)
        
        integer,intent(in)   :: NumElems
        integer,intent(in)   :: NumNodes_tot
        integer,intent(in)   :: NumNodes(NumElems)
        integer,intent(in)   :: MemNo(NumElems)
        integer,intent(in)   :: Conn(MaxElNod*NumElems)
        integer,intent(in)   :: Master_Array(MaxElNod*NumElems*2)
        real(8),intent(in)   :: Length(NumElems)
        real(8),intent(in)   :: PreCurv(3*NumElems)
        real(8),intent(in)   :: Psi(3*NumElems)
        real(8),intent(in)   :: Vector(3*NumElems)
        real(8),intent(in)   :: Mass_Array(6*NumElems*6)
        real(8),intent(in)   :: Stiff_Array(6*NumElems*6)
        real(8),intent(in)   :: InvStiff_Array(6*NumElems*6)
        real(8),intent(in)   :: RBMass_Array(MaxElNod*NumElems*6*6)
        integer,intent(in)   :: Nod_Master(2*NumNodes_tot)
        integer,intent(in)   :: Nod_Vdof(NumNodes_tot)
        integer,intent(in)   :: Nod_Fdof(NumNodes_tot)
        real(8),intent(in)   :: PosIni_Array(NumNodes_tot*3)
        real(8),intent(in)   :: PsiIni_Array(NumElems*MaxElNod*3)
        real(8),intent(in)   :: PosDefor_Array(NumNodes_tot*3)
        real(8),intent(in)   :: PsiDefor_Array(NumElems*MaxElNod*3)
        real(8),intent(in)   :: ForceStatic_Array(NumNodes_tot*6)
        integer,intent(in)   :: DimMat
        integer,intent(in)   :: NumDof
        integer,intent(inout):: ks
        real(8),intent(inout):: Kglobal_Array(NumDof*NumDof)
        integer,intent(inout):: fs
        real(8),intent(inout):: Fglobal_Array(NumDof*NumDof)
        real(8),intent(inout):: Qglobal(NumDof)
        
        type(xbelem),allocatable:: Elem(:)
        type(xbnode),allocatable:: Node(:)
        real(8)     ,allocatable:: Coords(:,:)
        real(8)     ,allocatable:: Psi0(:,:,:)
        real(8)     ,allocatable:: PosDefor(:,:)
        real(8)     ,allocatable:: PsiDefor(:,:,:)
        real(8)     ,allocatable:: AppForces(:,:)
        type(sparse),allocatable:: Kglobal(:)
        type(sparse),allocatable:: Fglobal(:)
        real(8)     ,allocatable:: Kglobal_Full(:,:)
        real(8)     ,allocatable:: Fglobal_Full(:,:)
        
        type(xbopts):: Options
        
        allocate(Elem(NumElems))
        allocate(Node(NumNodes_tot))
        allocate(Coords(NumNodes_tot,3)); Coords = 0.d0
        allocate(Psi0(NumElems,MaxElNod,3)); Psi0 = 0.d0
        allocate(PosDefor(NumNodes_tot,3)); PosDefor = 0.d0
        allocate(PsiDefor(NumElems,MaxElNod,3)); PsiDefor = 0.d0
        allocate(AppForces(NumNodes_tot,6)); AppForces = 0.d0
        allocate(Kglobal(DimMat*NumDof)); call sparse_zero(ks,Kglobal)
        allocate(Fglobal(DimMat*NumDof)); call sparse_zero(fs,Fglobal)
        allocate(Kglobal_Full(NumDof,NumDof)); Kglobal_Full = 0.d0
        allocate(Fglobal_Full(NumDof,NumDof)); Fglobal_Full = 0.d0
        
        ! Convert PY data to F90 data
        call do_xbelem_var(NumElems,Elem,NumNodes,MemNo,Conn,Master_Array,&
        &       Length,PreCurv,Psi,Vector,Mass_Array,Stiff_Array,&
        &       InvStiff_Array,RBMass_Array)
        
        call pack_xbnode(NumNodes_tot,Node,Nod_Master,Nod_Vdof,Nod_Fdof)
        
        call vec2mat_double(PosIni_Array,Coords,NumNodes_tot,3)
        call vec2mat3d_double(PsiIni_Array,Psi0,NumElems,MaxElNod,3)
        call vec2mat_double(PosDefor_Array,PosDefor,NumNodes_tot,3)
        call vec2mat3d_double(PsiDefor_Array,PsiDefor,NumElems,MaxElNod,3)
        call vec2mat_double(ForceStatic_Array,AppForces,NumNodes_tot,6)
        
        ! Assembly matrices for a nonlinear static problem
        call cbeam3_asbly_static(Elem,Node,Coords,Psi0,PosDefor,PsiDefor,&
        &       AppForces,ks,Kglobal,fs,Fglobal,Qglobal,Options)
        
        ! Convert F90 data to PY data
        call sparse2full_rank(ks,Kglobal,NumDof,NumDof,Kglobal_Full)
        call sparse2full_rank(fs,Fglobal,NumDof,NumDof,Fglobal_Full)
        
        call mat2vec_double(Kglobal_Full,Kglobal_Array,NumDof,NumDof)
        call mat2vec_double(Fglobal_Full,Fglobal_Array,NumDof,NumDof)
        
        deallocate(Elem)
        deallocate(Node)
        deallocate(Coords)
        deallocate(Psi0)
        deallocate(PosDefor)
        deallocate(PsiDefor)
        deallocate(AppForces)
        deallocate(Kglobal)
        deallocate(Fglobal)
        deallocate(Kglobal_Full)
        deallocate(Fglobal_Full)
        
        return
        
    end subroutine wrap_cbeam3_asbly_static
    
    !---------------------------------------------------------------------------
    ! Update global vectors 
    !---------------------------------------------------------------------------
    subroutine wrap_cbeam3_solv_update_static(NumElems,NumNodes_tot,NumNodes,&
    &           MemNo,Conn,Master_Array,Length,PreCurv,Psi,Vector,Mass_Array,&
    &           Stiff_Array,InvStiff_Array,RBMass_Array,Nod_Master,&
    &           Nod_Vdof,Nod_Fdof,NumDof,&
    &           DeltaX,PosIni_Array,PsiIni_Array,PosDefor_Array,PsiDefor_Array)
        
        integer,intent(in)   :: NumElems
        integer,intent(in)   :: NumNodes_tot
        integer,intent(in)   :: NumNodes(NumElems)
        integer,intent(in)   :: MemNo(NumElems)
        integer,intent(in)   :: Conn(MaxElNod*NumElems)
        integer,intent(in)   :: Master_Array(MaxElNod*NumElems*2)
        real(8),intent(in)   :: Length(NumElems)
        real(8),intent(in)   :: PreCurv(3*NumElems)
        real(8),intent(in)   :: Psi(3*NumElems)
        real(8),intent(in)   :: Vector(3*NumElems)
        real(8),intent(in)   :: Mass_Array(6*NumElems*6)
        real(8),intent(in)   :: Stiff_Array(6*NumElems*6)
        real(8),intent(in)   :: InvStiff_Array(6*NumElems*6)
        real(8),intent(in)   :: RBMass_Array(MaxElNod*NumElems*6*6)
        integer,intent(in)   :: Nod_Master(2*NumNodes_tot)
        integer,intent(in)   :: Nod_Vdof(NumNodes_tot)
        integer,intent(in)   :: Nod_Fdof(NumNodes_tot)
        integer,intent(in)   :: NumDof
        real(8),intent(in)   :: DeltaX(NumDof)
        real(8),intent(in)   :: PosIni_Array(NumNodes_tot*3)
        real(8),intent(in)   :: PsiIni_Array(NumElems*MaxElNod*3)
        real(8),intent(inout):: PosDefor_Array(NumNodes_tot*3)
        real(8),intent(inout):: PsiDefor_Array(NumElems*MaxElNod*3)
        
        type(xbelem),allocatable:: Elem(:)
        type(xbnode),allocatable:: Node(:)
        real(8)     ,allocatable:: Psi0(:,:,:)
        real(8)     ,allocatable:: PosDefor(:,:)
        real(8)     ,allocatable:: PsiDefor(:,:,:)
        real(8)     ,allocatable:: Coords(:,:)
        
        allocate(Elem(NumElems))
        allocate(Node(NumNodes_tot))
        allocate(Psi0(NumElems,MaxElNod,3)); Psi0 = 0.d0
        allocate(PosDefor(NumNodes_tot,3)); PosDefor = 0.d0
        allocate(PsiDefor(NumElems,MaxElNod,3)); PsiDefor = 0.d0
        allocate(Coords(NumNodes_tot,3)); Coords = 0.d0
        
        ! Convert PY data to F90 data
        call do_xbelem_var(NumElems,Elem,NumNodes,MemNo,Conn,Master_Array,&
        &       Length,PreCurv,Psi,Vector,Mass_Array,Stiff_Array,&
        &       InvStiff_Array,RBMass_Array)
        
        call pack_xbnode(NumNodes_tot,Node,Nod_Master,Nod_Vdof,Nod_Fdof)
        
        call vec2mat3d_double(PsiIni_Array,Psi0,NumElems,MaxElNod,3)
        call vec2mat_double(PosDefor_Array,PosDefor,NumNodes_tot,3)
        call vec2mat3d_double(PsiDefor_Array,PsiDefor,NumElems,MaxElNod,3)
        call vec2mat_double(PosIni_Array,Coords,NumNodes_tot,3)
        
        ! Update global vectors
        call cbeam3_solv_update_static(Elem,Node,Psi0,DeltaX,PosDefor,PsiDefor)
        
        ! Convert F90 data to PY data
        call mat2vec_double(PosDefor,PosDefor_Array,NumNodes_tot,3)
        call mat2vec3d_double(PsiDefor,PsiDefor_Array,NumElems,MaxElNod,3)
        
        deallocate(Elem)
        deallocate(Node)
        deallocate(Psi0)
        deallocate(PosDefor)
        deallocate(PsiDefor)
        deallocate(Coords)
        
        return
        
    end subroutine wrap_cbeam3_solv_update_static
    
    !---------------------------------------------------------------------------
    ! Dump out deformed configuration
    !---------------------------------------------------------------------------
    subroutine wrap_output_elems(NumElems,NumNodes_tot,NumNodes,MemNo,Conn,&
    &           Master_Array,Length,PreCurv,Psi,Vector,Mass_Array,&
    &           Stiff_Array,InvStiff_Array,RBMass_Array,PosDefor_Array,&
    &           PsiDefor_Array,OutFile,Ini0_Def1)
        
        integer,intent(in):: NumElems
        integer,intent(in):: NumNodes_tot
        integer,intent(in):: NumNodes(NumElems)
        integer,intent(in):: MemNo(NumElems)
        integer,intent(in):: Conn(MaxElNod*NumElems)
        integer,intent(in):: Master_Array(MaxElNod*NumElems*2)
        real(8),intent(in):: Length(NumElems)
        real(8),intent(in):: PreCurv(3*NumElems)
        real(8),intent(in):: Psi(3*NumElems)
        real(8),intent(in):: Vector(3*NumElems)
        real(8),intent(in):: Mass_Array(6*NumElems*6)
        real(8),intent(in):: Stiff_Array(6*NumElems*6)
        real(8),intent(in):: InvStiff_Array(6*NumElems*6)
        real(8),intent(in):: RBMass_Array(MaxElNod*NumElems*6*6)
        real(8),intent(in):: PosDefor_Array(NumNodes_tot*3)
        real(8),intent(in):: PsiDefor_Array(NumElems*MaxElNod*3)
        character(len=25),intent(in):: OutFile
        integer,intent(in):: Ini0_Def1
        
        type(xbelem),allocatable:: Elem(:)
        real(8)     ,allocatable:: PosDef(:,:)
        real(8)     ,allocatable:: PsiDef(:,:,:)
        real(8)     ,allocatable:: temp(:,:)
        
        allocate(Elem(NumElems))
        allocate(PosDef(NumNodes_tot,3)); PosDef = 0.d0
        allocate(PsiDef(NumElems,MaxElNod,3)); PsiDef = 0.d0
        allocate(temp(MaxElNod,3)); temp = 0.d0
        
        ! Convert PY data to F90 data
        call do_xbelem_var(NumElems,Elem,NumNodes,MemNo,Conn,Master_Array,&
        &                  Length,PreCurv,Psi,Vector,Mass_Array,Stiff_Array,&
        &                  InvStiff_Array,RBMass_Array)
        
        call vec2mat_double(PosDefor_Array,PosDef,NumNodes_tot,3)
        call vec2mat3d_double(PsiDefor_Array,PsiDef,NumElems,MaxElNod,3)
        
        ! Dump out deformed configuration
        ! Options: status='replace', access='append'
        if(Ini0_Def1.eq.0) then
            open(unit=11,file=OutFile(1:11)//'_und.txt',status='replace')
        else if(Ini0_Def1.eq.1) then
            open(unit=11,file=OutFile(1:11)//'_def.txt',status='replace')
        end if
        call output_elems(11,Elem,PosDef,PsiDef)
        close(11)
        
        deallocate(Elem)
        deallocate(PosDef)
        deallocate(PsiDef)
        deallocate(temp)
        
        return
        
    end subroutine wrap_output_elems
    
    !---------------------------------------------------------------------------
    ! Setup dynamic parameters
    !---------------------------------------------------------------------------
    subroutine wrap_input_dynsetup(NumSteps,t0,dt)
        
        integer,intent(out):: NumSteps
        real(8),intent(out):: t0
        real(8),intent(out):: dt
        type(xbopts)       :: Options
        
        call input_dynsetup(NumSteps,t0,dt,Options)
        
        return
        
    end subroutine wrap_input_dynsetup
    
    !---------------------------------------------------------------------------
    ! Define time-varying forcing terms
    !---------------------------------------------------------------------------
    subroutine wrap_input_dynforce(NumNodes,nrt,Time,ForceStatic,ForceDynAmp,&
    &           ForceTime)
        
        integer,intent(in) :: NumNodes
        integer,intent(in) :: nrt
        real(8),intent(in) :: Time(nrt)
        real(8),intent(in) :: ForceStatic(NumNodes,6)
        real(8),intent(out):: ForceDynAmp(NumNodes,6)
        real(8),intent(out):: ForceTime(nrt)
        
        call input_dynforce(NumNodes,Time,ForceStatic,ForceDynAmp,ForceTime)
        
        return
        
    end subroutine wrap_input_dynforce
    
    !---------------------------------------------------------------------------
    ! Setup time-varying forcing velocity
    !---------------------------------------------------------------------------
    subroutine wrap_input_forcedvel(NumNodes,nrt,Time,ForcedVel,ForcedVelDot)
        
        integer,intent(in) :: NumNodes
        integer,intent(in) :: nrt
        real(8),intent(in) :: Time(nrt)
        real(8),intent(out):: ForcedVel(nrt,6)
        real(8),intent(out):: ForcedVelDot(nrt,6)
        
        call input_forcedvel(NumNodes,Time,ForcedVel,ForcedVelDot)
        
        return
        
    end subroutine wrap_input_forcedvel
    
    !---------------------------------------------------------------------------
    ! Assembly matrices for dynamic problem
    !---------------------------------------------------------------------------
    subroutine wrap_cbeam3_asbly_dynamic(NumElems,NumNodes_tot,NumNodes,&
    &           MemNo,Conn,Master_Array,Length,PreCurv,Psi,Vector,Mass_Array,&
    &           Stiff_Array,InvStiff_Array,RBMass_Array,&
    &           Nod_Master,Nod_Vdof,Nod_Fdof,&
    &           Coords_Array,Psi0_Array,PosDefor_Array,PsiDefor_Array,&
    &           PosDeforDot_Array,PsiDeforDot_Array,&
    &           PosDeforDDot_Array,PsiDeforDDot_Array,&
    &           Force_Array,Vrel,VrelDot,NumDof,DimMat,&
    &           ms,Mglobal_Array,Mvel_Array,&
    &           cs,Cglobal_Array,Cvel_Array,&
    &           ks,Kglobal_Array,&
    &           fs,Fglobal_Array,Qglobal,Cao_Array)
        
        integer,intent(in) :: NumElems
        integer,intent(in) :: NumNodes_tot
        integer,intent(in) :: NumNodes(NumElems)
        integer,intent(in) :: MemNo(NumElems)
        integer,intent(in) :: Conn(MaxElNod*NumElems)
        integer,intent(in) :: Master_Array(MaxElNod*NumElems*2)
        real(8),intent(in) :: Length(NumElems)
        real(8),intent(in) :: PreCurv(3*NumElems)
        real(8),intent(in) :: Psi(3*NumElems)
        real(8),intent(in) :: Vector(3*NumElems)
        real(8),intent(in) :: Mass_Array(6*NumElems*6)
        real(8),intent(in) :: Stiff_Array(6*NumElems*6)
        real(8),intent(in) :: InvStiff_Array(6*NumElems*6)
        real(8),intent(in) :: RBMass_Array(MaxElNod*NumElems*6*6)
        integer,intent(in) :: Nod_Master(2*NumNodes_tot)
        integer,intent(in) :: Nod_Vdof(NumNodes_tot)
        integer,intent(in) :: Nod_Fdof(NumNodes_tot)
        real(8),intent(in) :: Coords_Array(NumNodes_tot*3)
        real(8),intent(in) :: Psi0_Array(NumElems*MaxElNod*3)
        real(8),intent(in) :: PosDefor_Array(NumNodes_tot*3)
        real(8),intent(in) :: PsiDefor_Array(NumElems*MaxElNod*3)
        real(8),intent(in) :: PosDeforDot_Array(NumNodes_tot*3)
        real(8),intent(in) :: PsiDeforDot_Array(NumElems*MaxElNod*3)
        real(8),intent(in) :: PosDeforDDot_Array(NumNodes_tot*3)
        real(8),intent(in) :: PsiDeforDDot_Array(NumElems*MaxElNod*3)
        real(8),intent(in) :: Force_Array(NumNodes_tot*6)
        real(8),intent(in) :: Vrel(6)
        real(8),intent(in) :: VrelDot(6)
        integer,intent(in) :: NumDof
        integer,intent(in) :: DimMat
        integer,intent(out):: ms
        real(8),intent(out):: Mglobal_Array(NumDof*NumDof)
        real(8),intent(out):: Mvel_Array(NumDof*6)
        integer,intent(out):: cs
        real(8),intent(out):: Cglobal_Array(NumDof*NumDof)
        real(8),intent(out):: Cvel_Array(NumDof*6)
        integer,intent(out):: ks
        real(8),intent(out):: Kglobal_Array(NumDof*NumDof)
        integer,intent(out):: fs
        real(8),intent(out):: Fglobal_Array(NumDof*NumDof)
        real(8),intent(out):: Qglobal(NumDof)
        !type(xbopts),intent(in):: Options ! DEPRECATED in PY framework
        real(8),intent(in) :: Cao_Array(3*3)
        
        type(xbelem),allocatable:: Elem(:)
        type(xbnode),allocatable:: Node(:)
        real(8)     ,allocatable:: Coords(:,:)
        real(8)     ,allocatable:: Psi0(:,:,:)
        real(8)     ,allocatable:: PosDefor(:,:)
        real(8)     ,allocatable:: PsiDefor(:,:,:)
        real(8)     ,allocatable:: PosDeforDot(:,:)
        real(8)     ,allocatable:: PsiDeforDot(:,:,:)
        real(8)     ,allocatable:: PosDeforDDot(:,:)
        real(8)     ,allocatable:: PsiDeforDDot(:,:,:)
        real(8)     ,allocatable:: Force(:,:)
        type(sparse),allocatable:: Mglobal(:)
        real(8)     ,allocatable:: Mglobal_Full(:,:)
        real(8)     ,allocatable:: Mvel(:,:)
        type(sparse),allocatable:: Cglobal(:)
        real(8)     ,allocatable:: Cglobal_Full(:,:)
        real(8)     ,allocatable:: Cvel(:,:)
        type(sparse),allocatable:: Kglobal(:)
        real(8)     ,allocatable:: Kglobal_Full(:,:)
        type(sparse),allocatable:: Fglobal(:)
        real(8)     ,allocatable:: Fglobal_Full(:,:)
        real(8)     ,allocatable:: Cao(:,:)
        
        type(xbopts):: Options
        
        allocate(Elem(NumElems))
        allocate(Node(NumNodes_tot))
        allocate(Coords(NumNodes_tot,3)); Coords = 0.d0
        allocate(Psi0(NumElems,MaxElNod,3)); Psi0 = 0.d0
        allocate(PosDefor(NumNodes_tot,3)); PosDefor = 0.d0
        allocate(PsiDefor(NumElems,MaxElNod,3)); PsiDefor = 0.d0
        allocate(PosDeforDot(NumNodes_tot,3)); PosDeforDot = 0.d0
        allocate(PsiDeforDot(NumElems,MaxElNod,3)); PsiDeforDot = 0.d0
        allocate(PosDeforDDot(NumNodes_tot,3)); PosDeforDDot = 0.d0
        allocate(PsiDeforDDot(NumElems,MaxElNod,3)); PsiDeforDDot = 0.d0
        allocate(Force(NumNodes_tot,6)); Force = 0.d0
        allocate(Mglobal(DimMat*NumDof)); call sparse_zero(ms,Mglobal)
        allocate(Mglobal_Full(NumDof,NumDof)); Mglobal_Full = 0.d0
        allocate(Mvel(NumDof,6)); Mvel = 0.d0
        allocate(Cglobal(DimMat*NumDof)); call sparse_zero(cs,Cglobal)
        allocate(Cglobal_Full(NumDof,NumDof)); Cglobal_Full = 0.d0
        allocate(Cvel(NumDof,6)); Cvel = 0.d0
        allocate(Kglobal(DimMat*NumDof)); call sparse_zero(ks,Kglobal)
        allocate(Kglobal_Full(NumDof,NumDof)); Kglobal_Full = 0.d0
        allocate(Fglobal(DimMat*NumDof)); call sparse_zero(fs,Fglobal)
        allocate(Fglobal_Full(NumDof,NumDof)); Fglobal_Full = 0.d0
        allocate(Cao(3,3)); Cao = 0.d0
        
        ! Convert PY data to F90 data
        call do_xbelem_var(NumElems,Elem,NumNodes,MemNo,Conn,Master_Array,&
        &       Length,PreCurv,Psi,Vector,Mass_Array,Stiff_Array,&
        &       InvStiff_Array,RBMass_Array)
        
        call pack_xbnode(NumNodes_tot,Node,Nod_Master,Nod_Vdof,Nod_Fdof)
        
        call vec2mat_double(Coords_Array,Coords,NumNodes_tot,3)
        call vec2mat3d_double(Psi0_Array,Psi0,NumElems,MaxElNod,3)
        call vec2mat_double(PosDefor_Array,PosDefor,NumNodes_tot,3)
        call vec2mat3d_double(PsiDefor_Array,PsiDefor,NumElems,MaxElNod,3)
        
        call vec2mat_double(PosDeforDot_Array,PosDeforDot,NumNodes_tot,3)
        call vec2mat3d_double(PsiDeforDot_Array,PsiDeforDot,NumElems,MaxElNod,3)
        call vec2mat_double(PosDeforDDot_Array,PosDeforDDot,NumNodes_tot,3)
        call vec2mat3d_double(PsiDeforDDot_Array,PsiDeforDDot,NumElems,MaxElNod,3)
        
        call vec2mat_double(Force_Array,Force,NumNodes_tot,6)
        call vec2mat_double(Cao_Array,Cao,3,3)
        
        ! Assembly matrices for a dynamic problem
        call cbeam3_asbly_dynamic(Elem,Node,Coords,Psi0,&
        &       PosDefor,PsiDefor,PosDeforDot,PsiDeforDot,&
        &       PosDeforDDot,PsiDeforDDot,Force,Vrel,VrelDot,&
        &       ms,Mglobal,Mvel,cs,Cglobal,Cvel,ks,Kglobal,fs,Fglobal,&
        &       Qglobal,Options,Cao)
        
        ! Convert F90 data to PY data
        call sparse2full_rank(ms,Mglobal,NumDof,NumDof,Mglobal_Full)
        call sparse2full_rank(cs,Cglobal,NumDof,NumDof,Cglobal_Full)
        call sparse2full_rank(ks,Kglobal,NumDof,NumDof,Kglobal_Full)
        call sparse2full_rank(fs,Fglobal,NumDof,NumDof,Fglobal_Full)
        
        call mat2vec_double(Mglobal_Full,Mglobal_Array,NumDof,NumDof)
        call mat2vec_double(Cglobal_Full,Cglobal_Array,NumDof,NumDof)
        call mat2vec_double(Kglobal_Full,Kglobal_Array,NumDof,NumDof)
        call mat2vec_double(Fglobal_Full,Fglobal_Array,NumDof,NumDof)
        
        call mat2vec_double(Mvel,Mvel_Array,NumDof,6)
        call mat2vec_double(Cvel,Cvel_Array,NumDof,6)
        
        deallocate(Elem)
        deallocate(Node)
        deallocate(Coords)
        deallocate(Psi0)
        deallocate(PosDefor)
        deallocate(PsiDefor)
        deallocate(PosDeforDot)
        deallocate(PsiDeforDot)
        deallocate(PosDeforDDot)
        deallocate(PsiDeforDDot)
        deallocate(Force)
        deallocate(Mglobal)
        deallocate(Mglobal_Full)
        deallocate(Mvel)
        deallocate(Cglobal)
        deallocate(Cglobal_Full)
        deallocate(Cvel)
        deallocate(Kglobal)
        deallocate(Kglobal_Full)
        deallocate(Fglobal)
        deallocate(Fglobal_Full)
        deallocate(Cao)
        
        return
        
    end subroutine wrap_cbeam3_asbly_dynamic
    
    !---------------------------------------------------------------------------
    ! From physical coordinates to state space vector
    !---------------------------------------------------------------------------
    subroutine wrap_cbeam3_solv_disp2state(&
    &           NumNodes_tot,NumDof,NumElems,Nod_Master,Nod_Vdof,Nod_Fdof,&
    &           Pos_Array,Psi_Array,PosDot_Array,PsiDot_Array,&
    &           X,dXdt)
        
        integer,intent(in) :: NumNodes_tot
        integer,intent(in) :: NumDof
        integer,intent(in) :: NumElems
        integer,intent(in) :: Nod_Master(2*NumNodes_tot)
        integer,intent(in) :: Nod_Vdof(NumNodes_tot)
        integer,intent(in) :: Nod_Fdof(NumNodes_tot)
        real(8),intent(in) :: Pos_Array(NumNodes_tot*3)
        real(8),intent(in) :: Psi_Array(NumElems*MaxElNod*3)
        real(8),intent(in) :: PosDot_Array(NumNodes_tot*3)
        real(8),intent(in) :: PsiDot_Array(NumElems*MaxElNod*3)
        real(8),intent(out):: X(NumDof)
        real(8),intent(out):: dXdt(NumDof)
        
        type(xbnode),allocatable:: Node(:)
        real(8)     ,allocatable:: Pos(:,:)
        real(8)     ,allocatable:: Psi(:,:,:)
        real(8)     ,allocatable:: PosDot(:,:)
        real(8)     ,allocatable:: PsiDot(:,:,:)
        
        allocate(Node(NumNodes_tot))
        allocate(Pos(NumNodes_tot,3)); Pos = 0.d0
        allocate(Psi(NumElems,MaxElNod,3)); Psi = 0.d0
        allocate(PosDot(NumNodes_tot,3)); PosDot = 0.d0
        allocate(PsiDot(NumElems,MaxElNod,3)); PsiDot = 0.d0
        
        ! Convert PY data to F90 data
        call pack_xbnode(NumNodes_tot,Node,Nod_Master,Nod_Vdof,Nod_Fdof)
        
        call vec2mat_double(Pos_Array,Pos,NumNodes_tot,3)
        call vec2mat3d_double(Psi_Array,Psi,NumElems,MaxElNod,3)
        call vec2mat_double(PosDot_Array,PosDot,NumNodes_tot,3)
        call vec2mat3d_double(PsiDot_Array,PsiDot,NumElems,MaxElNod,3)
        
        ! From physical coordinates to state space vector
        call cbeam3_solv_disp2state(Node,Pos,Psi,PosDot,PsiDot,X,dXdt)
        
        deallocate(Node)
        deallocate(Pos)
        deallocate(Psi)
        deallocate(PosDot)
        deallocate(PsiDot)
        
        return
        
    end subroutine wrap_cbeam3_solv_disp2state
    
    !---------------------------------------------------------------------------
    ! From state space vector to physical coordinates
    !---------------------------------------------------------------------------
    subroutine wrap_cbeam3_solv_state2disp(&
    &           NumElems,NumNodes_tot,NumNodes,&
    &           MemNo,Conn,Master_Array,Length,PreCurv,Psi,Vector,Mass_Array,&
    &           Stiff_Array,InvStiff_Array,RBMass_Array,&
    &           Nod_Master,Nod_Vdof,Nod_Fdof,&
    &           Coords_Array,Psi0_Array,&
    &           NumDof,X,dXdt,&
    &           Pos_Array,Psi_Array,PosDot_Array,PsiDot_Array)
        
        integer,intent(in) :: NumElems
        integer,intent(in) :: NumNodes_tot
        integer,intent(in) :: NumNodes(NumElems)
        integer,intent(in) :: MemNo(NumElems)
        integer,intent(in) :: Conn(MaxElNod*NumElems)
        integer,intent(in) :: Master_Array(MaxElNod*NumElems*2)
        real(8),intent(in) :: Length(NumElems)
        real(8),intent(in) :: PreCurv(3*NumElems)
        real(8),intent(in) :: Psi(3*NumElems)
        real(8),intent(in) :: Vector(3*NumElems)
        real(8),intent(in) :: Mass_Array(6*NumElems*6)
        real(8),intent(in) :: Stiff_Array(6*NumElems*6)
        real(8),intent(in) :: InvStiff_Array(6*NumElems*6)
        real(8),intent(in) :: RBMass_Array(MaxElNod*NumElems*6*6)
        integer,intent(in) :: Nod_Master(2*NumNodes_tot)
        integer,intent(in) :: Nod_Vdof(NumNodes_tot)
        integer,intent(in) :: Nod_Fdof(NumNodes_tot)
        real(8),intent(in) :: Coords_Array(NumNodes_tot*3)
        real(8),intent(in) :: Psi0_Array(NumElems*MaxElNod*3)
        integer,intent(in) :: NumDof
        real(8),intent(in) :: X(NumDof)
        real(8),intent(in) :: dXdt(NumDof)
        real(8),intent(out):: Pos_Array(NumNodes_tot*3)
        real(8),intent(out):: Psi_Array(NumElems*MaxElNod*3)
        real(8),intent(out):: PosDot_Array(NumNodes_tot*3)
        real(8),intent(out):: PsiDot_Array(NumElems*MaxElNod*3)
        
        type(xbelem),allocatable:: Elem(:)
        type(xbnode),allocatable:: Node(:)
        real(8)     ,allocatable:: Coords(:,:)
        real(8)     ,allocatable:: Psi0(:,:,:)
        real(8)     ,allocatable:: Pos(:,:)
        real(8)     ,allocatable:: Psi3d(:,:,:)
        real(8)     ,allocatable:: PosDot(:,:)
        real(8)     ,allocatable:: PsiDot(:,:,:)
        
        allocate(Elem(NumElems))
        allocate(Node(NumNodes_tot))
        allocate(Coords(NumNodes_tot,3)); Coords = 0.d0
        allocate(Psi0(NumElems,MaxElNod,3)); Psi0 = 0.d0
        allocate(Pos(NumNodes_tot,3)); Pos = 0.d0
        allocate(Psi3d(NumElems,MaxElNod,3)); Psi3d = 0.d0
        allocate(PosDot(NumNodes_tot,3)); PosDot = 0.d0
        allocate(PsiDot(NumElems,MaxElNod,3)); PsiDot = 0.d0
        
        ! Convert PY data to F90 data
        call do_xbelem_var(NumElems,Elem,NumNodes,MemNo,Conn,Master_Array,&
        &       Length,PreCurv,Psi,Vector,Mass_Array,Stiff_Array,&
        &       InvStiff_Array,RBMass_Array)
        
        call pack_xbnode(NumNodes_tot,Node,Nod_Master,Nod_Vdof,Nod_Fdof)
        
        call vec2mat_double(Coords_Array,Coords,NumNodes_tot,3)
        call vec2mat3d_double(Psi0_Array,Psi0,NumElems,MaxElNod,3)
        
        ! From state vector to physical coordinates
        call cbeam3_solv_state2disp(Elem,Node,Coords,Psi0,X,dXdt,Pos,Psi3d,PosDot,PsiDot)
        
        ! Convert F90 data to PY data
        call mat2vec_double(Pos,Pos_Array,NumNodes_tot,3)
        call mat2vec3d_double(Psi3d,Psi_Array,NumElems,MaxElNod,3)
        call mat2vec_double(PosDot,PosDot_Array,NumNodes_tot,3)
        call mat2vec3d_double(PsiDot,PsiDot_Array,NumElems,MaxElNod,3)
        
        deallocate(Elem)
        deallocate(Node)
        deallocate(Coords)
        deallocate(Psi0)
        deallocate(Pos)
        deallocate(Psi3d)
        deallocate(PosDot)
        deallocate(PsiDot)
        
        return
        
    end subroutine wrap_cbeam3_solv_state2disp
    
    !---------------------------------------------------------------------------
    ! Update results from linear dynamic solution
    !---------------------------------------------------------------------------
    subroutine wrap_cbeam3_solv_update_lindyn(&
    &           NumElems,NumNodes_tot,NumNodes,&
    &           MemNo,Conn,Master_Array,Length,PreCurv,Psi,Vector,Mass_Array,&
    &           Stiff_Array,InvStiff_Array,RBMass_Array,&
    &           Nod_Master,Nod_Vdof,Nod_Fdof,&
    &           PsiIni_Array,NumDof,DX,DXDt,PosDefor_Array,PsiDefor_Array,&
    &           PosDotDefor_Array,PsiDotDefor_Array)
        
        integer,intent(in)   :: NumElems
        integer,intent(in)   :: NumNodes_tot
        integer,intent(in)   :: NumNodes(NumElems)
        integer,intent(in)   :: MemNo(NumElems)
        integer,intent(in)   :: Conn(MaxElNod*NumElems)
        integer,intent(in)   :: Master_Array(MaxElNod*NumElems*2)
        real(8),intent(in)   :: Length(NumElems)
        real(8),intent(in)   :: PreCurv(3*NumElems)
        real(8),intent(in)   :: Psi(3*NumElems)
        real(8),intent(in)   :: Vector(3*NumElems)
        real(8),intent(in)   :: Mass_Array(6*NumElems*6)
        real(8),intent(in)   :: Stiff_Array(6*NumElems*6)
        real(8),intent(in)   :: InvStiff_Array(6*NumElems*6)
        real(8),intent(in)   :: RBMass_Array(MaxElNod*NumElems*6*6)
        integer,intent(in)   :: Nod_Master(2*NumNodes_tot)
        integer,intent(in)   :: Nod_Vdof(NumNodes_tot)
        integer,intent(in)   :: Nod_Fdof(NumNodes_tot)
        real(8),intent(in)   :: PsiIni_Array(NumElems*MaxElNod*3)
        integer,intent(in)   :: NumDof
        real(8),intent(in)   :: DX(NumDof)
        real(8),intent(in)   :: DXDt(NumDof)
        real(8),intent(inout):: PosDefor_Array(NumNodes_tot*3)
        real(8),intent(inout):: PsiDefor_Array(NumElems*MaxElNod*3)
        real(8),intent(inout):: PosDotDefor_Array(NumNodes_tot*3)
        real(8),intent(inout):: PsiDotDefor_Array(NumElems*MaxElNod*3)
        
        type(xbelem),allocatable:: Elem(:)
        type(xbnode),allocatable:: Node(:)
        real(8)     ,allocatable:: Psi0(:,:,:)
        real(8)     ,allocatable:: PosDefor(:,:)
        real(8)     ,allocatable:: PsiDefor(:,:,:)
        real(8)     ,allocatable:: PosDotDefor(:,:)
        real(8)     ,allocatable:: PsiDotDefor(:,:,:)
        
        allocate(Elem(NumElems))
        allocate(Node(NumNodes_tot))
        allocate(Psi0(NumElems,MaxElNod,3)); Psi0 = 0.d0
        allocate(PosDefor(NumNodes_tot,3)); PosDefor = 0.d0
        allocate(PsiDefor(NumElems,MaxElNod,3)); PsiDefor = 0.d0
        allocate(PosDotDefor(NumNodes_tot,3)); PosDotDefor = 0.d0
        allocate(PsiDotDefor(NumElems,MaxElNod,3)); PsiDotDefor = 0.d0
        
        ! Convert PY data to F90 data
        call do_xbelem_var(NumElems,Elem,NumNodes,MemNo,Conn,Master_Array,&
        &       Length,PreCurv,Psi,Vector,Mass_Array,Stiff_Array,&
        &       InvStiff_Array,RBMass_Array)
        
        call pack_xbnode(NumNodes_tot,Node,Nod_Master,Nod_Vdof,Nod_Fdof)
        
        call vec2mat3d_double(PsiIni_Array,Psi0,NumElems,MaxElNod,3)
        call vec2mat_double(PosDefor_Array,PosDefor,NumNodes_tot,3)
        call vec2mat3d_double(PsiDefor_Array,PsiDefor,NumElems,MaxElNod,3)
        call vec2mat_double(PosDotDefor_Array,PosDotDefor,NumNodes_tot,3)
        call vec2mat3d_double(PsiDotDefor_Array,PsiDotDefor,NumElems,MaxElNod,3)
        
        ! Update results from linear dynamic solution
        call cbeam3_solv_update_lindyn(Elem,Node,Psi0,DX,DXDt,PosDefor,&
        &       PsiDefor,PosDotDefor,PsiDotDefor)
        
        ! Convert F90 data to PY data
        call mat2vec_double(PosDefor,PosDefor_Array,NumNodes_tot,3)
        call mat2vec3d_double(PsiDefor,PsiDefor_Array,NumElems,MaxElNod,3)
        call mat2vec_double(PosDotDefor,PosDotDefor_Array,NumNodes_tot,3)
        call mat2vec3d_double(PsiDotDefor,PsiDotDefor_Array,NumElems,MaxElNod,3)
        
        deallocate(Elem)
        deallocate(Node)
        deallocate(Psi0)
        deallocate(PosDefor)
        deallocate(PsiDefor)
        deallocate(PosDotDefor)
        deallocate(PsiDotDefor)
        
        return 
        
    end subroutine wrap_cbeam3_solv_update_lindyn
    
    !---------------------------------------------------------------------------
    ! Wrapper to xbeam_asbly_orient
    !---------------------------------------------------------------------------
    subroutine wrap_xbeam_asbly_orient(&
    &           NumElems,NumNodes_tot,NumNodes,&
    &           MemNo,Conn,Master_Array,Length,PreCurv,Psi,Vector,Mass_Array,&
    &           Stiff_Array,InvStiff_Array,RBMass_Array,&
    &           Nod_Master,Nod_Vdof,Nod_Fdof,&
    &           PosDefor_Array,PsiDefor_Array,Vrel,Quat,&
    &           CQR_Array,CQQ_Array,NumDof,DimMat,fs,Frigid_Array,Cao_Array)
        
        integer,intent(in) :: NumElems
        integer,intent(in) :: NumNodes_tot
        integer,intent(in) :: NumNodes(NumElems)
        integer,intent(in) :: MemNo(NumElems)
        integer,intent(in) :: Conn(MaxElNod*NumElems)
        integer,intent(in) :: Master_Array(MaxElNod*NumElems*2)
        real(8),intent(in) :: Length(NumElems)
        real(8),intent(in) :: PreCurv(3*NumElems)
        real(8),intent(in) :: Psi(3*NumElems)
        real(8),intent(in) :: Vector(3*NumElems)
        real(8),intent(in) :: Mass_Array(6*NumElems*6)
        real(8),intent(in) :: Stiff_Array(6*NumElems*6)
        real(8),intent(in) :: InvStiff_Array(6*NumElems*6)
        real(8),intent(in) :: RBMass_Array(MaxElNod*NumElems*6*6)
        integer,intent(in) :: Nod_Master(2*NumNodes_tot)
        integer,intent(in) :: Nod_Vdof(NumNodes_tot)
        integer,intent(in) :: Nod_Fdof(NumNodes_tot)
        real(8),intent(in) :: PosDefor_Array(NumNodes_tot*3)
        real(8),intent(in) :: PsiDefor_Array(NumElems*MaxElNod*3)
        real(8),intent(in) :: Vrel(6)
        real(8),intent(in) :: Quat(4)
        real(8),intent(out):: CQR_Array(4*6)
        real(8),intent(out):: CQQ_Array(4*4)
        integer,intent(in) :: NumDof
        integer,intent(in) :: DimMat
        integer,intent(out):: fs
        real(8),intent(out):: Frigid_Array(6*(NumDof+6))
        real(8),intent(in) :: Cao_Array(3*3)
        
        type(xbelem),allocatable:: Elem(:)
        type(xbnode),allocatable:: Node(:)
        real(8)     ,allocatable:: PosDefor(:,:)
        real(8)     ,allocatable:: PsiDefor(:,:,:)
        real(8)     ,allocatable:: CQR(:,:)
        real(8)     ,allocatable:: CQQ(:,:)
        type(sparse),allocatable:: Frigid(:)
        real(8)     ,allocatable:: Frigid_Full(:,:)
        real(8)     ,allocatable:: Cao(:,:)
        
        type(xbopts):: Options
        
        allocate(Elem(NumElems))
        allocate(Node(NumNodes_tot))
        allocate(PosDefor(NumNodes_tot,3)); PosDefor = 0.d0
        allocate(PsiDefor(NumElems,MaxElNod,3)); PsiDefor = 0.d0
        allocate(CQR(4,6)); CQR = 0.d0
        allocate(CQQ(4,4)); CQQ = 0.d0
        allocate(Frigid(DimMat*NumDof)); call sparse_zero(fs,Frigid)
        allocate(Frigid_Full(6,NumDof+6)); Frigid_Full = 0.d0
        allocate(Cao(3,3)); Cao = 0.d0
        
        ! Convert PY data to F90 data
        call do_xbelem_var(NumElems,Elem,NumNodes,MemNo,Conn,Master_Array,&
        &       Length,PreCurv,Psi,Vector,Mass_Array,Stiff_Array,&
        &       InvStiff_Array,RBMass_Array)
        
        call pack_xbnode(NumNodes_tot,Node,Nod_Master,Nod_Vdof,Nod_Fdof)
        
        call vec2mat_double(PosDefor_Array,PosDefor,NumNodes_tot,3)
        call vec2mat3d_double(PsiDefor_Array,PsiDefor,NumElems,MaxElNod,3)
        call vec2mat_double(CQR_Array,CQR,4,6)
        call vec2mat_double(CQQ_Array,CQQ,4,4)
        call vec2mat_double(Cao_Array,Cao,3,3)
        
        ! Assembly rigid body matrices for changing orientation
        call xbeam_asbly_orient(Elem,Node,PosDefor,PsiDefor,Vrel,Quat,&
        &       CQR,CQQ,fs,Frigid,Options,Cao)
        
        ! Convert F90 data to PY data
        call mat2vec_double(CQR,CQR_Array,4,6)
        call mat2vec_double(CQQ,CQQ_Array,4,4)
        
        call sparse2full_rank(fs,Frigid,6,NumDof+6,Frigid_Full)
        call mat2vec_double(Frigid_Full,Frigid_Array,6,NumDof+6)
        
        deallocate(Elem)
        deallocate(Node)
        deallocate(PosDefor)
        deallocate(PsiDefor)
        deallocate(CQR)
        deallocate(CQQ)
        deallocate(Frigid)
        deallocate(Frigid_Full)
        deallocate(Cao)
        
        return
        
    end subroutine wrap_xbeam_asbly_orient
    
    !---------------------------------------------------------------------------
    ! Wrapper to out_outgrid
    !---------------------------------------------------------------------------
    subroutine wrap_xbeam_asbly_dynamic(&
    &           NumElems,NumNodes_tot,NumNodes,&
    &           MemNo,Conn,Master_Array,Length,PreCurv,Psi,Vector,Mass_Array,&
    &           Stiff_Array,InvStiff_Array,RBMass_Array,&
    &           Nod_Master,Nod_Vdof,Nod_Fdof,&
    &           Coords_Array,Psi0_Array,PosDefor_Array,PsiDefor_Array,&
    &           PosDeforDot_Array,PsiDeforDot_Array,&
    &           PosDeforDDot_Array,PsiDeforDDot_Array,Vrel,VrelDot,Quat,&
    &           NumDof,DimMat,ms,MRS_Array,MRR_Array,cs,CRS_Array,CRR_Array,&
    &           CQR_Array,CQQ_Array,ks,KRS_Array,fs,Frigid_Array,&
    &           Qrigid,Cao_Array)
        
        integer,intent(in) :: NumElems
        integer,intent(in) :: NumNodes_tot
        integer,intent(in) :: NumNodes(NumElems)
        integer,intent(in) :: MemNo(NumElems)
        integer,intent(in) :: Conn(MaxElNod*NumElems)
        integer,intent(in) :: Master_Array(MaxElNod*NumElems*2)
        real(8),intent(in) :: Length(NumElems)
        real(8),intent(in) :: PreCurv(3*NumElems)
        real(8),intent(in) :: Psi(3*NumElems)
        real(8),intent(in) :: Vector(3*NumElems)
        real(8),intent(in) :: Mass_Array(6*NumElems*6)
        real(8),intent(in) :: Stiff_Array(6*NumElems*6)
        real(8),intent(in) :: InvStiff_Array(6*NumElems*6)
        real(8),intent(in) :: RBMass_Array(MaxElNod*NumElems*6*6)
        integer,intent(in) :: Nod_Master(2*NumNodes_tot)
        integer,intent(in) :: Nod_Vdof(NumNodes_tot)
        integer,intent(in) :: Nod_Fdof(NumNodes_tot)
        real(8),intent(in) :: Coords_Array(NumNodes_tot*3)
        real(8),intent(in) :: Psi0_Array(NumElems*MaxElNod*3)
        real(8),intent(in) :: PosDefor_Array(NumNodes_tot*3)
        real(8),intent(in) :: PsiDefor_Array(NumElems*MaxElNod*3)
        real(8),intent(in) :: PosDeforDot_Array(NumNodes_tot*3)
        real(8),intent(in) :: PsiDeforDot_Array(NumElems*MaxElNod*3)
        real(8),intent(in) :: PosDeforDDot_Array(NumNodes_tot*3)
        real(8),intent(in) :: PsiDeforDDot_Array(NumElems*MaxElNod*3)
        real(8),intent(in) :: Vrel(6)
        real(8),intent(in) :: VrelDot(6)
        real(8),intent(in) :: Quat(4)
        integer,intent(in) :: NumDof
        integer,intent(in) :: DimMat
        integer,intent(out):: ms
        real(8),intent(out):: MRS_Array(6*NumDof)
        real(8),intent(out):: MRR_Array(6*6)
        integer,intent(out):: cs
        real(8),intent(out):: CRS_Array(6*NumDof)
        real(8),intent(out):: CRR_Array(6*6)
        real(8),intent(out):: CQR_Array(4*6)
        real(8),intent(out):: CQQ_Array(4*4)
        integer,intent(out):: ks
        real(8),intent(out):: KRS_Array(6*NumDof)
        integer,intent(out):: fs
        real(8),intent(out):: Frigid_Array(6*(NumDof+6))
        real(8),intent(out):: Qrigid(6)
        real(8),intent(in) :: Cao_Array(3*3)
        
        type(xbelem),allocatable:: Elem(:)
        type(xbnode),allocatable:: Node(:)
        real(8)     ,allocatable:: Coords(:,:)
        real(8)     ,allocatable:: Psi0(:,:,:)
        real(8)     ,allocatable:: PosDefor(:,:)
        real(8)     ,allocatable:: PsiDefor(:,:,:)
        real(8)     ,allocatable:: PosDeforDot(:,:)
        real(8)     ,allocatable:: PsiDeforDot(:,:,:)
        real(8)     ,allocatable:: PosDeforDDot(:,:)
        real(8)     ,allocatable:: PsiDeforDDot(:,:,:)
        type(sparse),allocatable:: MRS(:)
        real(8)     ,allocatable:: MRS_Full(:,:)
        real(8)     ,allocatable:: MRR(:,:)
        type(sparse),allocatable:: CRS(:)
        real(8)     ,allocatable:: CRS_Full(:,:)
        real(8)     ,allocatable:: CRR(:,:)
        real(8)     ,allocatable:: CQR(:,:)
        real(8)     ,allocatable:: CQQ(:,:)
        type(sparse),allocatable:: KRS(:)
        real(8)     ,allocatable:: KRS_Full(:,:)
        type(sparse),allocatable:: Frigid(:)
        real(8)     ,allocatable:: Frigid_Full(:,:)
        real(8)     ,allocatable:: Cao(:,:)
        
        type(xbopts):: Options
        
        allocate(Elem(NumElems))
        allocate(Node(NumNodes_tot))
        allocate(Coords(NumNodes_tot,3)); Coords = 0.d0
        allocate(Psi0(NumElems,MaxElNod,3)); Psi0 = 0.d0
        allocate(PosDefor(NumNodes_tot,3)); PosDefor = 0.d0
        allocate(PsiDefor(NumElems,MaxElNod,3)); PsiDefor = 0.d0
        allocate(PosDeforDot(NumNodes_tot,3)); PosDeforDot = 0.d0
        allocate(PsiDeforDot(NumElems,MaxElNod,3)); PsiDeforDot = 0.d0
        allocate(PosDeforDDot(NumNodes_tot,3)); PosDeforDDot = 0.d0
        allocate(PsiDeforDDot(NumElems,MaxElNod,3)); PsiDeforDDot = 0.d0
        allocate(MRS(DimMat*NumDof)); call sparse_zero(ms,MRS)
        allocate(MRS_Full(6,NumDof)); MRS_Full = 0.d0
        allocate(MRR(6,6)); MRR = 0.d0
        allocate(CRS(DimMat*NumDof)); call sparse_zero(cs,CRS)
        allocate(CRS_Full(6,NumDof)); CRS_Full = 0.d0
        allocate(CRR(6,6)); CRR = 0.d0
        allocate(CQR(4,6)); CQR = 0.d0
        allocate(CQQ(4,4)); CQQ = 0.d0
        allocate(KRS(DimMat*NumDof)); call sparse_zero(ks,KRS)
        allocate(KRS_Full(6,NumDof)); KRS_Full = 0.d0
        allocate(Frigid(DimMat*NumDof)); call sparse_zero(fs,Frigid)
        allocate(Frigid_Full(6,NumDof+6)); Frigid_Full = 0.d0
        allocate(Cao(3,3)); Cao = 0.d0
        
        ! Convert PY data to F90 data
        call do_xbelem_var(NumElems,Elem,NumNodes,MemNo,Conn,Master_Array,&
        &       Length,PreCurv,Psi,Vector,Mass_Array,Stiff_Array,&
        &       InvStiff_Array,RBMass_Array)
        
        call pack_xbnode(NumNodes_tot,Node,Nod_Master,Nod_Vdof,Nod_Fdof)
        
        call vec2mat_double(Coords_Array,Coords,NumNodes_tot,3)
        call vec2mat3d_double(Psi0_Array,Psi0,NumElems,MaxElNod,3)
        call vec2mat_double(PosDefor_Array,PosDefor,NumNodes_tot,3)
        call vec2mat3d_double(PsiDefor_Array,PsiDefor,NumElems,MaxElNod,3)
        call vec2mat_double(PosDeforDot_Array,PosDeforDot,NumNodes_tot,3)
        call vec2mat3d_double(PsiDeforDot_Array,PsiDeforDot,NumElems,MaxElNod,3)
        call vec2mat_double(PosDeforDDot_Array,PosDeforDDot,NumNodes_tot,3)
        call vec2mat3d_double(PsiDeforDDot_Array,PsiDeforDDot,NumElems,MaxElNod,3)
        call vec2mat_double(Cao_Array,Cao,3,3)
        
        ! Assembly rigid-body matrices for the dynamic problem
        call xbeam_asbly_dynamic(Elem,Node,Coords,Psi0,PosDefor,PsiDefor,&
        &       PosDeforDot,PsiDeforDot,PosDeforDDot,PsiDeforDDot,&
        &       Vrel,VrelDot,Quat,ms,MRS,MRR,cs,CRS,CRR,CQR,CQQ,ks,KRS,&
        &       fs,Frigid,Qrigid,Options,Cao)
        
        ! Convert F90 data to PY data
        call sparse2full_rank(ms,MRS,6,NumDof,MRS_Full)
        call mat2vec_double(MRS_Full,MRS_Array,6,NumDof)
        
        call mat2vec_double(MRR,MRR_Array,6,6)
        
        call sparse2full_rank(cs,CRS,6,NumDof,CRS_Full)
        call mat2vec_double(CRS_Full,CRS_Array,6,NumDof)
        
        call mat2vec_double(CRR,CRR_Array,6,6)
        call mat2vec_double(CQR,CQR_Array,4,6)
        call mat2vec_double(CQQ,CQQ_Array,4,4)
        
        call sparse2full_rank(ks,KRS,6,NumDof,KRS_Full)
        call mat2vec_double(KRS_Full,KRS_Array,6,NumDof)
        
        call sparse2full_rank(fs,Frigid,6,NumDof+6,Frigid_Full)
        call mat2vec_double(Frigid_Full,Frigid_Array,6,NumDof+6)
        
        deallocate(Elem)
        deallocate(Node)
        deallocate(Coords)
        deallocate(Psi0)
        deallocate(PosDefor)
        deallocate(PsiDefor)
        deallocate(PosDeforDot)
        deallocate(PsiDeforDot)
        deallocate(PosDeforDDot)
        deallocate(PsiDeforDDot)
        deallocate(MRS)
        deallocate(MRS_Full)
        deallocate(MRR)
        deallocate(CRS)
        deallocate(CRS_Full)
        deallocate(CRR)
        deallocate(CQR)
        deallocate(CQQ)
        deallocate(KRS)
        deallocate(KRS_Full)
        deallocate(Frigid)
        deallocate(Frigid_Full)
        deallocate(Cao)
        
        return
        
    end subroutine wrap_xbeam_asbly_dynamic
    
    !---------------------------------------------------------------------------
    ! Wrapper to fem_glob2loc_extract
    !---------------------------------------------------------------------------
    subroutine wrap_fem_glob2loc_extract(NumElems,NumNodes_tot,NumNodes,MemNo,&
    &           Conn,Master_Array,Length,PreCurv,Psi,Vector,Mass_Array,&
    &           Stiff_Array,InvStiff_Array,RBMass_Array,PosDefor_Array,&
    &           PsiDefor_Array,PosGlob_Array,NumNE_array)
        
        integer,intent(in) :: NumElems
        integer,intent(in) :: NumNodes_tot
        integer,intent(in) :: NumNodes(NumElems)
        integer,intent(in) :: MemNo(NumElems)
        integer,intent(in) :: Conn(MaxElNod*NumElems)
        integer,intent(in) :: Master_Array(MaxElNod*NumElems*2)
        real(8),intent(in) :: Length(NumElems)
        real(8),intent(in) :: PreCurv(3*NumElems)
        real(8),intent(in) :: Psi(3*NumElems)
        real(8),intent(in) :: Vector(3*NumElems)
        real(8),intent(in) :: Mass_Array(6*NumElems*6)
        real(8),intent(in) :: Stiff_Array(6*NumElems*6)
        real(8),intent(in) :: InvStiff_Array(6*NumElems*6)
        real(8),intent(in) :: RBMass_Array(MaxElNod*NumElems*6*6)
        real(8),intent(in) :: PosDefor_Array(NumNodes_tot*3)
        real(8),intent(in) :: PsiDefor_Array(NumElems*MaxElNod*3)
        real(8),intent(out):: PosGlob_Array(NumElems*MaxElNod*3)
        integer,intent(out):: NumNE_array(NumElems)
        
        type(xbelem),allocatable:: Elem(:)
        real(8)     ,allocatable:: PosDef(:,:)
        real(8)     ,allocatable:: PsiDef(:,:,:)
        real(8)     ,allocatable:: PosElem(:,:)
        real(8)     ,allocatable:: PosGlob(:,:)
        
        integer:: iElem,i1,i2,NumNE
        
        allocate(Elem(NumElems))
        allocate(PosDef(NumNodes_tot,3)); PosDef = 0.d0
        allocate(PsiDef(NumElems,MaxElNod,3)); PsiDef = 0.d0
        allocate(PosElem(MaxElNod,3)); PosElem = 0.d0
        allocate(PosGlob(NumElems*MaxElNod,3)); PosGlob = 0.d0
        
        ! Convert PY data to F90 data
        call do_xbelem_var(NumElems,Elem,NumNodes,MemNo,Conn,Master_Array,&
        &                  Length,PreCurv,Psi,Vector,Mass_Array,Stiff_Array,&
        &                  InvStiff_Array,RBMass_Array)
        
        call vec2mat_double(  PosDefor_Array,PosDef,NumNodes_tot,     3)
        call vec2mat3d_double(PsiDefor_Array,PsiDef,NumElems,MaxElNod,3)
        
        do iElem=1,NumElems
            call fem_glob2loc_extract(Elem(iElem)%Conn,PosDef,PosElem,NumNE)
            i1 = 1 + (iElem-1)*NumNE
            i2 = i1 + NumNE
            PosGlob(i1:i2,:) = PosElem
            NumNE_array(iElem) = NumNE
        end do
        
        call mat2vec_double(PosGlob,PosGlob_Array,NumElems*MaxElNod,3)

        deallocate(Elem)
        deallocate(PosDef)
        deallocate(PsiDef)
        deallocate(PosElem)
        deallocate(PosGlob)
        
        return
        
    end subroutine wrap_fem_glob2loc_extract
    
    !---------------------------------------------------------------------------
    ! Wrapper for cbeam3_solv_nlndyn under external forces
    !---------------------------------------------------------------------------

        subroutine wrap_cbeam3_solv_nlndyn(iOut,NumDof,NumSteps,Time,&
&                   NumElems, NumNodes, MemNo, Conn,        &!for do_xbelem_var
&                   Master_Array,                           &!for do_xbelem_var
&                   Length, PreCurv,                        &!for do_xbelem_var
&                   Psi, Vector, Mass_Array,                &!for do_xbelem_var
&                   Stiff_Array,                            &!for do_xbelem_var
&                   InvStiff_Array, RBMass_Array,           &!for do_xbelem_var
&                   NumNodes_tot, Master, Vdof, Fdof,       &!for pack_xbnode
&                   F0_Vec,Fa_Vec,Ftime,                            &
&                   Vrel_Vec, VrelDot_Vec, Coords_Vec, Psi0_Vec, PosDefor_Vec,  &
&                   PsiDefor_Vec, PosDotDefor_Vec, PsiDotDefor_Vec,     &
&                   PosPsiTime_Vec,VelocTime_Vec,DynOut_Vec,            &
&                   OutGrids,                               &
&                   FollowerForce, FollowerForceRig,        &!for pack_xbopts
&                   PrintInfo, OutInBframe, OutInaframe,    &!for pack_xbopts
&                   ElemProj, MaxIterations, NumLoadSteps,  &!for pack_xbopts
&                   NumGauss, Solution, DeltaCurved,        &!for pack_xbopts
&                   MinDelta, NewmarkDamp)                   !for pack_xbopts

        integer,    intent(in)  :: iOut                 ! Output file.
        integer,    intent(in)  :: NumDof               ! No. of independent DoF
        integer,    intent(in)  :: NumSteps         ! Number of timesteps
        real(8),    intent(in)  :: Time(NumSteps+1) ! Time steps.
        integer,    intent(in)  :: NumElems             !for do_xbelem_var
        integer,    intent(in)  :: NumNodes(NumElems)   !for do_xbelem_var+pkxbn !Number of nodes in each element
        integer,    intent(in)  :: MemNo(NumElems)      !for do_xbelem_var
        integer,    intent(in)  :: Conn(MaxElNod*NumElems)   !for do_xbelem_var
        integer,    intent(in)  :: Master_Array(MaxElNod*NumElems*2) !for do_xbe
        real(8),    intent(in)  :: Length(NumElems)     !for do_xbelem_var
        real(8),    intent(in)  :: PreCurv(3*NumElems)  !for do_xbelem_var
        real(8),    intent(in)  :: Psi(3*NumElems)      !for do_xbelem_var
        real(8),    intent(in)  :: Vector(3*NumElems)   !for do_xbelem_var
        real(8),    intent(in)  :: Mass_Array(6*NumElems*6) !for do_xbelem_var
        real(8),    intent(in)  :: Stiff_Array(6*NumElems*6)!for do_xbelem_var
        real(8),    intent(in)  :: InvStiff_Array(6*NumElems*6) !for do_xbelem_v
        real(8),    intent(in)  :: RBMass_Array(MaxElNod*NumElems*6*6)!for do_xb
        integer,    intent(in)  :: NumNodes_tot             !for pack_xbnode
        integer,    intent(in)  :: Master(2*NumNodes_tot)   !for pack_xbnode
        integer,    intent(in)  :: Vdof(NumNodes_tot)       !for pack_xbnode
        integer,    intent(in)  :: Fdof(NumNodes_tot)       !for pack_xbnode
        real(8),    intent(in)  :: F0_Vec(NumNodes_tot*6)   !Applied static nodal forces
        real(8),    intent(in)  :: Fa_Vec(NumNodes_tot*6)   !Amplitude of dynamic nodal forces
        real(8),    intent(in)  :: Ftime(NumSteps+1)    !Time history of applied forces
        real(8),    intent(in)  :: Vrel_Vec((NumSteps+1)*6)!Time history of vel of ref frame
        real(8),    intent(in)  :: VrelDot_Vec((NumSteps+1)*6)!Ti Hstry of accel of ref frame
        real(8),    intent(in)  :: Coords_Vec(NumNodes_tot*3)   !Undefrmd coords of grid points
        real(8),    intent(in)  :: Psi0_Vec(NumElems*MaxElNod*3)    !Undefrmd CRV of nodes in elems
        real(8),    intent(inout)   :: PosDefor_Vec(NumNodes_tot*3)   !Initial/final grid pts
        real(8),    intent(inout)   :: PsiDefor_Vec(NumElems*MaxElNod*3) !Init/Fnl CRVs
        real(8),    intent(inout)   :: PosDotDefor_Vec(NumNodes_tot*3)!d/dt of PosDefor
        real(8),    intent(inout)   :: PsiDotDefor_Vec(NumElems*MaxElNod*3) !d/dt of PsiDefor
        real(8),    intent(out)     :: PosPsiTime_Vec((NumSteps+1)*6)   !t-hstry of pos/crv at selected nodes.
        real(8),    intent(out)     :: VelocTime_Vec((NumSteps+1)*NumNodes_tot) !t-hstry of t-dervatvs at selected nodes.
        real(8),    intent(out)     :: DynOut_Vec(((NumSteps+1)*NumNodes_tot)*3)    !snapshot t-hstry of pos of all nodes to be appended
        logical,    intent(inout)   :: OutGrids(NumNodes_tot)   !Output grids
        logical,    intent(in) :: FollowerForce
        logical,    intent(in) :: FollowerForceRig
        logical,    intent(in) :: PrintInfo
        logical,    intent(in) :: OutInBframe
        logical,    intent(in) :: OutInaframe
        integer,    intent(in) :: ElemProj
        integer,    intent(in) :: MaxIterations
        integer,    intent(in) :: NumLoadSteps
        integer,    intent(in) :: NumGauss
        integer,    intent(in) :: Solution
        real(8),    intent(in) :: DeltaCurved
        real(8),    intent(in) :: MinDelta
        real(8),    intent(in) :: NewmarkDamp
        ! TODO: test wrapper.
        ! out and inout arguments must be then be converted from fortran matrix to vector after solve.
        ! TODO: write and test a wrapper for logicals and arrays of logicals
        ! F(t) = F0+Ftime(iStep)*Fa
        ! TODO: erors here are copy-pasted in nlnstatic

        ! Declare local variables
        type(xbelem),allocatable:: Elem(:) ! Initialise xbelem derived type
        type(xbnode),allocatable:: Node(:) ! Initialise xbnode derived type
        real(8)     ,allocatable:: F0(:,:)
        real(8)     ,allocatable:: Fa(:,:)
        real(8)     ,allocatable:: Vrel(:,:)
        real(8)     ,allocatable:: VrelDot(:,:)
        real(8)     ,allocatable:: Coords(:,:)
        real(8)     ,allocatable:: Psi0(:,:,:)
        real(8)     ,allocatable:: PosDefor(:,:)
        real(8)     ,allocatable:: PsiDefor(:,:,:)
        real(8)     ,allocatable:: PosDotDefor(:,:)
        real(8)     ,allocatable:: PsiDotDefor(:,:,:)
        real(8)     ,allocatable:: PosPsiTime(:,:)
        real(8)     ,allocatable:: VelocTime(:,:)
        real(8)     ,allocatable:: DynOut(:,:)
        ! create Options struct
        type(xbopts):: Options

        ! allocate memory for vectors of derived type
        allocate(Elem(NumElems))
        allocate(Node(NumNodes_tot))

        ! allocate memory for F90 Arrays
        allocate(F0(NumNodes_tot,6))
        allocate(Fa(NumNodes_tot,6))
        allocate(Vrel((NumSteps+1),6))
        allocate(VrelDot((NumSteps+1),6))
        allocate(Coords(NumNodes_tot,3))
        allocate(Psi0(NumElems,MaxElNod,3))
        allocate(PosDefor(NumNodes_tot,3))
        allocate(PsiDefor(NumElems,MaxElNod,3))
        allocate(PosDotDefor(NumNodes_tot,3))
        allocate(PsiDotDefor(NumElems,MaxElNod,3))
        allocate(PosPsiTime(NumSteps+1,6))
        allocate(VelocTime(NumSteps+1,NumNodes_tot))
        allocate(DynOut((NumSteps+1)*NumNodes_tot,3))

        ! Convert PY data to F90 data
        call do_xbelem_var(NumElems,Elem,NumNodes,MemNo,Conn,           &
&               Master_Array,Length,PreCurv,Psi,Vector,Mass_Array,      &
&               Stiff_Array,InvStiff_Array,RBMass_Array)

        ! Convert PY data to F90 data
        call pack_xbnode(NumNodes_tot,Node,Master,Vdof,Fdof)

        ! Convert PY data to F90 data
        call pack_xbopts(Options,FollowerForce,FollowerForceRig,            &
&                   PrintInfo,OutInBframe,OutInaframe,ElemProj,             &
&                   MaxIterations,NumLoadSteps,                             &
&                   NumGauss,Solution,DeltaCurved,MinDelta,NewmarkDamp)

        ! Convert vectors to multi-dim fortran arrays
        call vec2mat_double(F0_Vec,F0,NumNodes_tot,6)
        call vec2mat_double(Fa_Vec,Fa,NumNodes_tot,6)
        call vec2mat_double(Vrel_Vec,Vrel,size(Time),6)
        call vec2mat_double(VrelDot_Vec,VrelDot,size(Time),6)
        call vec2mat_double(Coords_Vec,Coords,NumNodes_tot,3)
        call vec2mat3d_double(Psi0_Vec,Psi0,NumElems,MaxElNod,3)
        call vec2mat_double(PosDefor_Vec,PosDefor,NumNodes_tot,3)
        call vec2mat3d_double(PsiDefor_Vec,PsiDefor,NumElems,MaxElNod,3)
        call vec2mat_double(PosDotDefor_Vec,PosDotDefor,NumNodes_tot,3)
        call vec2mat3d_double(PsiDotDefor_Vec,PsiDotDefor,NumElems,MaxElNod,3)
        call vec2mat_double(PosPsiTime_Vec,PosPsiTime,size(Time),6)
        call vec2mat_double(VelocTime_Vec,VelocTime,size(Time),NumNodes_tot)
        call vec2mat_double(DynOut_Vec,DynOut,size(Time)*NumNodes_tot,3)

        ! Call fortran solver
        call cbeam3_solv_nlndyn (iOut,NumDof,Time,Elem,Node,                &
                            F0,Fa,Ftime,                                    &
&                           Vrel,VrelDot,Coords,Psi0,PosDefor,PsiDefor,     &
&                           PosDotDefor,PsiDotDefor,PosPsiTime,VelocTime,   &
&                           DynOut,OutGrids,                                &
&                           Options)

        ! Convert multi-dim fortran arrays to vectors
        call mat2vec_double(PosDefor,PosDefor_Vec,NumNodes_tot,3)
        call mat2vec3d_double(PsiDefor,PsiDefor_Vec,NumElems,MaxElNod,3)
        call mat2vec_double(PosDotDefor,PosDotDefor_Vec,NumNodes_tot,3)
        call mat2vec3d_double(PsiDotDefor,PsiDotDefor_Vec,NumElems,MaxElNod,3)
        call mat2vec_double(PosPsiTime,PosPsiTime_Vec,size(Time),6)
        call mat2vec_double(VelocTime,VelocTime_Vec,size(Time),NumNodes_tot)
        call mat2vec3d_double(DynOut,DynOut_Vec,size(Time),NumNodes_tot,3)

        ! deallocate memory for vectors of derived types
        deallocate(Elem)
        deallocate(Node)

        ! deallocate memory for F90 Arrays
        deallocate(F0)
        deallocate(Fa)
        deallocate(Vrel)
        deallocate(VrelDot)
        deallocate(Coords)
        deallocate(Psi0)
        deallocate(PosDefor)
        deallocate(PsiDefor)
        deallocate(PosDotDefor)
        deallocate(PsiDotDefor)
        deallocate(PosPsiTime)
        deallocate(VelocTime)
        deallocate(DynOut)

    end subroutine wrap_cbeam3_solv_nlndyn


    !---------------------------------------------------------------------------
    ! Wrapper for cbeam3_solv_nlnstatic under external forces
    !---------------------------------------------------------------------------
    subroutine wrap_cbeam3_solv_nlnstatic (NumDof,&
&                   NumElems, NumNodes, MemNo, Conn,        &!for do_xbelem_var
&                   Master_Array,                           &!for do_xbelem_var
&                   Length, PreCurv,                        &!for do_xbelem_var
&                   Psi, Vector, Mass_Array,                &!for do_xbelem_var
&                   Stiff_Array,                            &!for do_xbelem_var
&                   InvStiff_Array, RBMass_Array,           &!for do_xbelem_var
&                   NumNodes_tot, Master, Vdof, Fdof,       &!for pack_xbnode
&                   AppForces_Vec,&
&                   Coords_Vec,Psi0_Vec,&
&                   PosDefor_Vec,PsiDefor_Vec,&
&                   FollowerForce, FollowerForceRig,        &!for pack_xbopts
&                   PrintInfo, OutInBframe, OutInaframe,    &!for pack_xbopts
&                   ElemProj, MaxIterations, NumLoadSteps,  &!for pack_xbopts
&                   NumGauss, Solution, DeltaCurved,        &!for pack_xbopts
&                   MinDelta, NewmarkDamp)                   !for pack_xbopts

        integer,    intent(in)  :: NumDof               ! No. of independent DoF
        integer,    intent(in)  :: NumElems             !for do_xbelem_var
        integer,    intent(in)  :: NumNodes(NumElems)   !for do_xbelem_var+pkxbn !Number of nodes in each element
        integer,    intent(in)  :: MemNo(NumElems)      !for do_xbelem_var
        integer,    intent(in)  :: Conn(MaxElNod*NumElems)   !for do_xbelem_var
        integer,    intent(in)  :: Master_Array(MaxElNod*NumElems*2) !for do_xbe
        real(8),    intent(in)  :: Length(NumElems)     !for do_xbelem_var
        real(8),    intent(in)  :: PreCurv(3*NumElems)  !for do_xbelem_var
        real(8),    intent(in)  :: Psi(3*NumElems)      !for do_xbelem_var
        real(8),    intent(in)  :: Vector(3*NumElems)   !for do_xbelem_var
        real(8),    intent(in)  :: Mass_Array(6*NumElems*6) !for do_xbelem_var
        real(8),    intent(in)  :: Stiff_Array(6*NumElems*6)!for do_xbelem_var
        real(8),    intent(in)  :: InvStiff_Array(6*NumElems*6) !for do_xbelem_v
        real(8),    intent(in)  :: RBMass_Array(MaxElNod*NumElems*6*6)!for do_xbe
        integer,    intent(in)  :: NumNodes_tot             !for pack_xbnode
        integer,    intent(in)  :: Master(2*NumNodes_tot)   !for pack_xbnode
        integer,    intent(in)  :: Vdof(NumNodes_tot)       !for pack_xbnode
        integer,    intent(in)  :: Fdof(NumNodes_tot)       !for pack_xbnode
        real(8),    intent(in)  :: AppForces_Vec(NumNodes_tot*6)!Applied static nodal forces
        real(8),    intent(in)  :: Coords_Vec(NumNodes_tot*3)   !Undefrmd coords of grid points
        real(8),    intent(in)  :: Psi0_Vec(NumElems*MaxElNod*3)    !Undefrmd CRV of nodes in elems
        real(8),    intent(inout)   :: PosDefor_Vec(NumNodes_tot*3)   !Initial/final grid pts
        real(8),    intent(inout)   :: PsiDefor_Vec(NumElems*MaxElNod*3) !Init/Fnl CRVs
        logical,    intent(in) :: FollowerForce
        logical,    intent(in) :: FollowerForceRig
        logical,    intent(in) :: PrintInfo
        logical,    intent(in) :: OutInBframe
        logical,    intent(in) :: OutInaframe
        integer,    intent(in) :: ElemProj
        integer,    intent(in) :: MaxIterations
        integer,    intent(in) :: NumLoadSteps
        integer,    intent(in) :: NumGauss
        integer,    intent(in) :: Solution
        real(8),    intent(in) :: DeltaCurved
        real(8),    intent(in) :: MinDelta
        real(8),    intent(in) :: NewmarkDamp

        ! Declare local variables
        type(xbelem),allocatable:: Elem(:) ! Initialise vec of xbelem derived type
        type(xbnode),allocatable:: Node(:) ! Initialise vec xbnode derived type
        real(8)     ,allocatable:: AppForces(:,:)
        real(8)     ,allocatable:: Coords(:,:)
        real(8)     ,allocatable:: Psi0(:,:,:)
        real(8)     ,allocatable:: PosDefor(:,:)
        real(8)     ,allocatable:: PsiDefor(:,:,:)

        ! create Options struct
        type(xbopts):: Options

        ! allocate memory for vectors of derived type
        allocate(Elem(NumElems))
        allocate(Node(NumNodes_tot))

        ! allocate memory for F90 Arrays
        allocate(AppForces(NumNodes_tot,6))
        allocate(Coords(NumNodes_tot,3))
        allocate(Psi0(NumElems,MaxElNod,3))
        allocate(PosDefor(NumNodes_tot,3))
        allocate(PsiDefor(NumElems,MaxElNod,3))

        ! Convert PY data to F90 data
        call do_xbelem_var(NumElems,Elem,NumNodes,MemNo,Conn,           &
&               Master_Array,Length,PreCurv,Psi,Vector,Mass_Array,      &
&               Stiff_Array,InvStiff_Array,RBMass_Array)

        ! Convert PY data to F90 data
        call pack_xbnode(NumNodes_tot,Node,Master,Vdof,Fdof)

        ! Convert PY data to F90 data
        call pack_xbopts(Options,FollowerForce,FollowerForceRig,            &
&                   PrintInfo,OutInBframe,OutInaframe,ElemProj,             &
&                   MaxIterations,NumLoadSteps,                             &
&                   NumGauss,Solution,DeltaCurved,MinDelta,NewmarkDamp)

        ! Convert vectors to multi-dim fortran arrays
        call vec2mat_double(AppForces_Vec,AppForces,NumNodes_tot,6)
        call vec2mat_double(Coords_Vec,Coords,NumNodes_tot,3)
        call vec2mat3d_double(Psi0_Vec,Psi0,NumElems,MaxElNod,3)
        call vec2mat_double(PosDefor_Vec,PosDefor,NumNodes_tot,3)
        call vec2mat3d_double(PsiDefor_Vec,PsiDefor,NumElems,MaxElNod,3)

        call cbeam3_solv_nlnstatic (NumDof,Elem,Node,AppForces,Coords,Psi0, &
&                                  PosDefor,PsiDefor,Options)

        ! Convert multi-dim fortran arrays to vectors
        call mat2vec_double(PosDefor,PosDefor_Vec,NumNodes_tot,3)
        call mat2vec3d_double(PsiDefor,PsiDefor_Vec,NumElems,MaxElNod,3)

        ! deallocate memory for vectors of derived types
        deallocate(Elem)
        deallocate(Node)

        ! deallocate memory for F90 Arrays
        deallocate(AppForces)
        deallocate(Coords)
        deallocate(Psi0)
        deallocate(PosDefor)
        deallocate(PsiDefor)

    end subroutine wrap_cbeam3_solv_nlnstatic

end module test

