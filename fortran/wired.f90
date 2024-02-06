! Wired Benchmarks - Fortran Implementation
! (c) 2024 ryan@freestatelabs


module wired
    use omp_lib

    implicit none 
    public cross, crosscols, crosscols2, subcols, multmv, cross2
    public createproblem, createproblem2
    public solve, solvecm, solve2, solve3

    contains

    ! Make the test problem with row-major ordering
    subroutine createproblem(Nn, Ns, nodes, sources)

        integer :: Nn, Ns, i
        real, dimension(:,:) :: nodes, sources 

        real :: zstart, zend, dz, z0, z1

        zstart = -1000.0 
        zend = 1000.0 
        dz = (zend - zstart)/Ns 
        z0 = zstart 
        z1 = dz 
        
        DO i = 1,Ns
            sources(1,i) = 0.0
            sources(2,i) = 0.0
            sources(3,i) = z0
            sources(4,i) = 0.0
            sources(5,i) = 0.0
            sources(6,i) = z1
            sources(7,i) = 5000.0
            z0 = z1 
            z1 = z1 + dz 
        END DO 
        DO i = 1,Nn
            nodes(1,i) = 3.0
            nodes(2,i) = 3.0
            nodes(3,i) = 0.0
        END DO 

    end subroutine createproblem 


    ! Create the problem with column-major ordering
    subroutine createproblem2(Nn, Ns, nodes, sources)

        integer :: Nn, Ns, i
        real, dimension(:,:) :: nodes, sources 

        real :: zstart, zend, dz, z0, z1

        zstart = -1000.0 
        zend = 1000.0 
        dz = (zend - zstart)/Ns 
        z0 = zstart 
        z1 = dz 
        
        DO i = 1,Ns
            sources(1,i) = 0.0
            sources(2,i) = 0.0
            sources(3,i) = z0
            sources(4,i) = 0.0
            sources(5,i) = 0.0
            sources(6,i) = z1
            sources(7,i) = 5000.0
            z0 = z1 
            z1 = z1 + dz 
        END DO 
        DO i = 1,Nn
            nodes(i,1) = 3.0
            nodes(i,2) = 3.0
            nodes(i,3) = 0.0
        END DO 

    end subroutine createproblem2 


    ! Compute c = a x b for 3-element vectors
    pure function cross(a, b) result(c)

        real, dimension(3) :: c
        real, dimension(3), intent(in) :: a, b
        
        c(1) = a(2) * b(3) - a(3) * b(2)
        c(2) = a(3) * b(1) - a(1) * b(3)
        c(3) = a(1) * b(2) - a(2) * b(1)

    end function cross


    ! Compute c = a x b for 3xN matrices, across the columns of a using vector b
    subroutine crosscols(c, a, b) 

        real, dimension(:,:) :: c
        real, dimension(:,:), intent(in) :: a 
        real, dimension(:), intent(in) :: b
        
        c(1,:) = a(2,:) * b(3) - a(3,:) * b(2)
        c(2,:) = a(3,:) * b(1) - a(1,:) * b(3)
        c(3,:) = a(1,:) * b(2) - a(2,:) * b(1)

    end subroutine crosscols


    ! Compute c = a x b for Nx3 matrices, across the rows of a using vector b
    subroutine crosscols2(c, a, b) 

        real, dimension(:,:) :: c
        real, dimension(:,:), intent(in) :: a 
        real, dimension(:), intent(in) :: b
        
        c(:,1) = a(:,2) * b(3) - a(:,3) * b(2)
        c(:,2) = a(:,3) * b(1) - a(:,1) * b(3)
        c(:,3) = a(:,1) * b(2) - a(:,2) * b(1)

    end subroutine crosscols2


    ! Subtract columns of b from vector a
    subroutine subcols(c, a, b) 

        real, dimension(:,:) :: c
        real, dimension(:), intent(in) :: a 
        real, dimension(:,:), intent(in) :: b
        integer i

        c(:,1) = -b(:,1) + a(1) 
        c(:,2) = -b(:,2) + a(2) 
        c(:,3) = -b(:,3) + a(3)

    end subroutine subcols


    ! Compute c = a x b by modifying vector c ("in-place") in an attempt to 
    ! reduce memory allocations.
    ! Does not appear to be any faster.
    subroutine cross2(c, a, b)

        real, dimension(3) :: c
        real, dimension(3), intent(in) :: a, b
        
        c(1) = a(2) * b(3) - a(3) * b(2)
        c(2) = a(3) * b(1) - a(1) * b(3)
        c(3) = a(1) * b(2) - a(2) * b(1)

    end subroutine cross2


    ! Take the norm of the columns of A and place them in b
    subroutine normcols(b, A) 

        real, dimension(:) :: b 
        real, dimension(:,:) :: A 

        b(:) = (A(:,1)**2 + A(:,2)**2 + A(:,3)**2)**0.5
    end subroutine normcols 


    ! Take the dot product of the columns of A with vector b and place them in c 
    subroutine dotcols(c, A, b) 

        real, dimension(:) :: c, b 
        real, dimension(:,:) :: A

        c(:) = A(:,1)*b(1) + A(:,2)*b(2) + A(:,3)*b(3) 

    end subroutine dotcols


    ! Multiply Matrix A row-wise by vector b and place in C
    subroutine multmv(C, A, b) 

        real, dimension(:,:) :: C, A 
        real, dimension(:) :: b 

        C(:,1) = A(:,1) * b(:) 
        C(:,2) = A(:,2) * b(:)
        C(:,3) = A(:,3) * b(:)

    end subroutine multmv


    ! Solve the Biot-Savart problem with opt-level 1: simple loops 
    ! This is the row-major method (slower)
    function solve(nodes, sources) result (bfield)

        real, dimension(:,:), intent(in) :: nodes, sources 
        real, dimension(:,:), allocatable :: bfield

        ! Instantiate intermediate values 
        real, dimension(3) :: a, b, c, cxa 
        real :: mu0, d
        integer :: Nn, Ns, i, j

        allocate(bfield(3,ubound(nodes,2)))
        bfield = 0.0    ! initialize all elements to zero
        mu0 = 1e-7

        Nn = size(nodes, 2) 
        Ns = size(sources, 2)

        do j = 1,Ns 
            d = mu0 * sources(7,j)

            do i = 1,Nn 
                a = sources(4:6,j) - sources(1:3,j) 
                b = sources(1:3,j) - nodes(:,i)
                c = sources(4:6,j) - nodes(:,i) 
                cxa = cross(c,a)
                cxa = cxa * d * (norm2(cxa)**(-2)) * (dot_product(a,c)/norm2(c) - dot_product(a,b)/norm2(b))
                bfield(:,i) = bfield(:,i) + cxa             
            end do  

        end do 

    end function solve 


    ! Solve the Biot-Savart problem with opt-level 1: simple loops in column-major
    ! ordering (faster)
    function solvecm(nodes, sources) result (bfield)

        real, dimension(:,:), intent(in) :: nodes, sources 
        real, dimension(:,:), allocatable :: bfield

        ! Instantiate intermediate values 
        real, dimension(3) :: a, b, c, cxa 
        real :: mu0, d
        integer :: Nn, Ns, i, j

        allocate(bfield(ubound(nodes,1),3))
        bfield = 0.0 
        mu0 = 1e-7

        Nn = size(nodes, 1) 
        Ns = size(sources, 2)

        do j = 1,Ns 
            d = mu0 * sources(7,j)

            do i = 1,Nn 
                a = sources(4:6,j) - sources(1:3,j) 
                b = sources(1:3,j) - nodes(i,:)
                c = sources(4:6,j) - nodes(i,:) 
                cxa = cross(c,a)
                cxa = cxa * d * (norm2(cxa)**(-2)) * (dot_product(a,c)/norm2(c) - dot_product(a,b)/norm2(b))
                bfield(i,:) = bfield(i,:) + cxa             
            end do  

        end do 

    end function solvecm


    ! Solve the biot-savart problem with opt-level 2: vectorization
    function solve2(nodes, sources) result (bfield)

        real, dimension(:,:), intent(in) :: nodes, sources 
        real, dimension(:,:), allocatable :: bfield

        ! Intermediate values 
        real, dimension(3) :: a
        real, dimension(:,:), allocatable :: b, c, cxa
        real, dimension(:), allocatable :: norm_cxa, dot_ac, norm_c, dot_ab, norm_b 
        real :: mu0, d
        integer :: Nn, Ns, i, j, thread_id

        Ns = size(sources,2) 
        Nn = size(nodes,1)

        allocate(bfield(Nn,3))
        allocate(b(Nn,3))
        allocate(c(Nn,3))
        allocate(cxa(Nn,3))
        allocate(norm_cxa(Nn))
        allocate(dot_ac(Nn))
        allocate(norm_c(Nn))
        allocate(dot_ab(Nn))
        allocate(norm_b(Nn))
        bfield = 0.0    ! initialize entire array to zero

        mu0 = 1e-7

        do j = 1,Ns 
            d = mu0 * sources(7,j)
            a = sources(4:6,j) - sources(1:3,j) 

            ! Subtract matrics from vector row by row 
            call subcols(b, sources(1:3,j), nodes)
            call subcols(c, sources(4:6,j), nodes)

            call crosscols2(cxa, c, a)
            call normcols(norm_cxa, cxa) 
            call dotcols(dot_ac, c, a) 
            call dotcols(dot_ab, b, a) 
            call normcols(norm_b, b) 
            call normcols(norm_c, c)
            call multmv(cxa,  cxa, d * (norm_cxa(:)**(-2)) * (dot_ac(:)/norm_c(:) - dot_ab(:)/norm_b(:)))
            
            bfield(:,:) = bfield(:,:) + cxa(:,:)             
        end do  

    end function solve2


    ! Partition a matrix with length `N` into smaller segments that can be processed
    ! on `Nt` number of threads. `it` is the current thread number. Returns 
    ! a 2-length vector corresponding to start/stop indices
    function partition(it, Nt, N) result(i) 

        integer, intent(in) :: it, Nt, N 
        integer :: rem, Nperthread
        integer,dimension(2) :: i

        Nperthread = N / Nt  ! Integer division returns only quotient 
        rem = mod(N, Nt)

        if (it == 1) then 
            i(1) = 1 
            i(2) = i(1) + Nperthread - 1
            ! print*,"it = ",it," spans ", i(1), i(2)
        else if (it == Nt) then 
            i(2) = N 
            i(1) = i(2) - Nperthread - rem + 1 
            ! print*,"it = ",it," spans ", i(1), i(2)
        else
            i(1) = (it-1)*Nperthread + 1
            i(2) = i(1) + Nperthread - 1 
            ! print*,"it = ",it," spans ", i(1), i(2)
        end if 

    end function partition 


    ! Solve using parallel processing via OpenMP
    function solve3(nodes, sources, Nt) result(bfield) 

        real, dimension(:,:), intent(in) :: nodes, sources 
        integer :: Nt
        real, dimension(:,:), allocatable :: bfield, bfield_partial
        integer :: Ns, Nn, it, thread_id 
        integer, dimension(2) :: ispan

        Ns = ubound(sources,2)
        Nn = ubound(nodes,1)

        allocate(bfield(Nn,3))
        allocate(bfield_partial(Nn,3))
        bfield = 0.0
        bfield_partial = 0.0

        if (Nt == 0) then
            Nt = OMP_GET_NUM_THREADS()
        end if 
        print*,"Solving with ", Nt, " threads."

        !$OMP parallel private(bfield_partial) shared(bfield)

        !$OMP do
        do it = 1,Nt
            ispan = partition(it, Nt, Ns) 
            bfield_partial = solve2(nodes, sources(:,ispan(1):ispan(2)))
        end do
        !$OMP end do

        !$OMP critical 
        bfield = bfield + bfield_partial 
        !$OMP end critical 

        !$OMP end parallel

    end function solve3 

end module wired