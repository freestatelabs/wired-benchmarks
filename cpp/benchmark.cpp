#include <iostream>
#include <omp.h>
#include <wired.h>


void test_threads(int Nt)
{
    #pragma omp parallel for ordered
    for (int i=0; i<Nt; i++)
    {
        #pragma omp ordered
        {
            std::cout << "Hello from thread " << omp_get_thread_num(); << "\n";
        }
    }
}

void benchmark(int Nn, int Ns, int Ni, int Nt, bool verbose=false) 
{
    omp_set_num_threads(Nt);

    if (verbose)
    {
        test_threads(Nt);
    }
        
}