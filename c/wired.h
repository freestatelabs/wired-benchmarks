// Wired.jl benchmarks - C implementation 
// (c) 2024 ryan@freestatelabs


// Convenience function for allocating a 2D array initialized to zero
double** zeros(int Nrows, int Ncols);

// Free a 2D array from memory
void free2darray(double** array, int Ncols);

// Create the sources for the test problem
double** createsources(int Ns);

// Create the output nodes for the test problem (observation points)
double** createnodes(int Nn);

// Solve the problem with simple (ideally auto-vectorized) loops
double** solve2(double** nodes, double** sources, int Nn, int Ns);

// Partition the problem for multi-threaded processing
void partition(int i[], int it, int Nt, int N);

// Solve the problem via multi-threaded processing
double** solve3(double** nodes, double** sources, int Nn, int Ns);

// Cross the columns of `A` with a 3-length vector `b` and place in `C`
//   `N` is the number of columns of `A`
void crosscols(double** C, double** A, double* b, int N);

// Take the vector norm of the columns of `A` and place in `b`
//   `N` is the number of columns of `A`
void normcols(double* b, double** A, int N);