! Benchmark the Fortran implementation 
! (c) 2024 ryan@freestatelabs

program benchmark
    use wired
    use omp_lib

    INTEGER :: Nn, Ns, Nt, Nc, Ni, i
    REAL, ALLOCATABLE :: sources(:,:), nodes(:,:), bfield(:,:)
    REAL :: TSTART, TEND, CPU_TOTAL
    DOUBLE PRECISION :: OMP_START, OMP_END, WALL_TOTAL

    PRINT *, "Please enter the number of sources, nodes, iterations, and cores:"
    PRINT *, "Ns, Nn, Ni, Nc"
    READ(*,*) Ns, Nn, Ni, Nc
    print *, " "

    CALL OMP_SET_NUM_THREADS(Nc)

    PRINT *, "FWIRED program start..."
    !$OMP PARALLEL
    NT = OMP_GET_NUM_THREADS()
    !$OMP END PARALLEL
    print *, "  OMP thread count: ", NT
    print *, "  Sources:          ", Ns 
    print *, "  Nodes:            ", Nn
    print *, "  Iterations:       ", Ni

    ! Create the problem
    ALLOCATE(sources(7,Ns))
    ALLOCATE(nodes(Nn,3))
    call createproblem2(Nn, Ns, nodes, sources) 

    ! Do the calculation
    CALL cpu_time(TSTART)
    OMP_START = omp_get_wtime()

    do i = 1,Ni
        if (Nc > 1) then 
            bfield = solve3(nodes,sources, Nt)
        else if (Nc == 1) then 
            bfield = solve2(nodes, sources)
        else 
            bfield = solvecm(nodes,sources)
        end if

        if (i < Ni) then
            deallocate(bfield)
        end if

    end do 

    CALL cpu_time(TEND)
    OMP_END = OMP_GET_WTIME()
        
    CPU_TOTAL = (TEND - TSTART) / real(Ni)
    WALL_TOTAL = (OMP_END - OMP_START) / real(Ni)

    print *, " "
    print *, "Program complete."
    print *, "  CPU time  (s): ", CPU_TOTAL 
    print *, "  Wall time (s): ", WALL_TOTAL 
    print *, "  Ratio:         ", CPU_TOTAL/WALL_TOTAL
    print *, "  s/DOF:         ", (CPU_TOTAL/(Ns*Nn))

    print *, ""
    print *, "Bfield(:,1:10): "
    do i=1,10
        print *, "  ",bfield(i,:) 
    end do  

end program benchmark