/*	Benchmark the Wired program in the C++ language
	(c) 2024 ryan@freestatelabs.com 
*/

#include <iostream>
#include <cstdlib>

#include "benchmark.h"


int main(int argc, char *argv[]) 
{
	if (argc == 5)
	{
		// Number of nodes, sources, iterations, and threads
		int Nn = atoi(argv[1]);
		int Ns = atoi(argv[2]);
		int Ni = atoi(argv[3]);
        int Nt = atoi(argv[4]);

        std::cout << "\nRunning benchmark with:\n";
        std::cout << "\tNodes:      " << Nn << "\n";
        std::cout << "\tSources:    " << Ns << "\n";
        std::cout << "\tIterations: " << Ni << "\n";
        std::cout << "\tThreads:    " << Nt << "\n";

        benchmark(Nn, Ns, Ni, Nt, true);
	}
	else
	{
		std::cout << "Error in input arguments. Please try again\n";
		std::cout << "\tOnly " << argc << " command line args provided.\n";
	}

	return 0;
}