""" Wired.jl benchmark: Python implementation 
    (c) 2024 ryan@freestatelabs 

    Run this file from the command line with the following arguments: 
    $ python3 benchmark.py NUM_NODES NUM_SOURCES NUM_ITERATIONS OPT_LEVEL NUM_THREADS

    OPT_LEVEL definitions:
    0 : Python
    1 : NumPy loops 
    2 : NumPy, vectorized
    3 : Numba 
    4 : Numba, parallel 
"""
from wired import * 
import time, sys

def runtest(N=1000, M=1000, Ni=10, opt_level=2, Nt=1): 
    set_num_threads(Nt)
    times = np.zeros((Ni,1))

    # Create the problem
    solver = solvers[opt_level]
    if opt_level > 0: 
        nodes, sources = create_problem2(N, M)
    else:
        nodes, sources = create_problem(N, M)
    
    if opt_level > 2:
        print("Precompiling...")
        t0 = time.perf_counter()
        _ = solver(nodes[0:2,:], sources[0:2,:])    # Precompile and throw out result
        t1 = time.perf_counter()
        print("Precompilation complete. Precompilation time: {0:.3E} s".format(t1-t0))

    for i in range(0, Ni): 

        t0 = time.perf_counter()
        B = solver(nodes, sources) 
        t1 = time.perf_counter()
        times[i] = t1 - t0

    print("First 10 rows of B: ")
    if opt_level > 0: 
        print(B[:,0:10])
    else:
        print(B[0:10])
    avg_time = np.mean(times)
    avg_time_per_DOF = avg_time/(N*M)

    print("Running tests with opt level {0}.".format(opt_level))
    print(" Iterations:   {0}".format(Ni))
    print(" Nodes:        {0}".format(N))
    print(" Sources:      {0}".format(M))
    print(" Threads:      {0}".format(Nt))
    print(" Avg time:     {0:.3E} s".format(avg_time))
    print(" Avg time/DOF: {0:.3E} s".format(avg_time_per_DOF))
    print("Tests complete.")


# Run the program 
if __name__ == "__main__":
    a = [int(x) for x in sys.argv[1:]]
    runtest(a[0], a[1], a[2], a[3], a[4])