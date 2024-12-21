// Wired.jl benchmarks - C implementation 
// (c) 2024 ryan@freestatelabs

#include <stdio.h>
#include <omp.h>
#include <math.h>
#include <stdlib.h>

// Convenience function to create a 2D array of zeros
double** zeros(int Nrows, int Ncols) {

    double** array; 
    array = calloc(Nrows, sizeof(double*));
    for (int i=0; i<Nrows; i++) {
        array[i] = calloc(Ncols, sizeof(double));
    }

    return array;
}


// Create an array of pointer arrays for nodes
double** createnodes(int Nn) {

    int M = 3;
    double** nodes = zeros(M, Nn);

    for (int i=0; i<Nn; i++) {
        nodes[0][i] = 3.0;
        nodes[1][i] = 3.0; 
        nodes[2][i] = 0.0;
    }

    return nodes;
}


// De-allocate a 2D array
void free2darray(double** array, int Ncols) {

    for (int i=0; i<Ncols; i++) {
        free(array[i]);
    }

    free(array);
}


// Norm of a 3xN matrix along direction N
void normcols(double* b, double** A, int N) {
    for (int i=0; i<N; i++) {
        b[i] = sqrt(pow(A[0][i],2) +pow(A[1][i],2) +pow(A[2][i],2));
    }
}


// Cross the columns of A with vector b
void crosscols(double** C, double** A, double* b, int N) {

    for (int i=0; i<N; i++) {
        C[0][i] = A[1][i] * b[2] - A[2][i] * b[1];
        C[1][i] = A[2][i] * b[0] - A[0][i] * b[2];
        C[2][i] = A[0][i] * b[1] - A[1][i] * b[0];
    }
}


// Return an array of pointer arrays for the sources
double** createsources(int Ns) {

    int M = 7; 
    double** sources = zeros(Ns, M);

    double zstart = -1000.0;
    double zend = 1000.0;
    double zstep = (zend-zstart)/Ns;
    double z = zstart;
    double Ic = 5000.0;

    for (int i=0; i<Ns; i++){
        sources[i][0] = 0.0;
        sources[i][1] = 0.0;
        sources[i][2] = z; 
        sources[i][3] = 0.0;
        sources[i][4] = 0.0;
        sources[i][5] = z + zstep; 
        sources[i][6] = Ic;
        z += zstep; 
    }

    return sources;
}


// Solve with optimization level 2: auto-vectorization using the compiler
double** solve2(double** nodes, double** sources, int Nn, int Ns) {

    double** bfield = zeros(3, Nn);

    double* a = calloc(Nn, sizeof(double)); 
    double** b = zeros(3,Nn);
    double** c = zeros(3,Nn);
    double** cxa = zeros(3,Nn);
    double* norm_cxa = calloc(Nn, sizeof(double)); 
    double* dot_ac = calloc(Nn, sizeof(double)); 
    double* dot_ab = calloc(Nn, sizeof(double)); 
    double* norm_c = calloc(Nn, sizeof(double)); 
    double* norm_b = calloc(Nn, sizeof(double)); 
    double* s = calloc(Nn, sizeof(double)); 
    double d;

    // Outer loop over sources
    for (int j=0; j<Ns; j++) {

        d = (1e-7) * sources[j][6];
        a[0] = sources[j][3] - sources[j][0];
        a[1] = sources[j][4] - sources[j][1];
        a[2] = sources[j][5] - sources[j][2];

        // Inner loop over nodes
        for (int i=0; i<Nn; i++) {
            b[0][i] = sources[j][0] - nodes[0][i];
            b[1][i] = sources[j][1] - nodes[1][i];
            b[2][i] = sources[j][2] - nodes[2][i];

            c[0][i] = sources[j][3] - nodes[0][i];
            c[1][i] = sources[j][4] - nodes[1][i];
            c[2][i] = sources[j][5] - nodes[2][i];
        }

        crosscols(cxa, c, a, Nn);
        normcols(norm_cxa, cxa, Nn);
        normcols(norm_b, b, Nn);
        normcols(norm_c, c, Nn);

        for (int i=0; i<Nn; i++) {
            dot_ac[i] = a[0]*c[0][i] + a[1]*c[1][i] + a[2]*c[2][i];
        }
            
        for (int i=0; i<Nn; i++) {
            dot_ab[i] = a[0]*b[0][i] + a[1]*b[1][i] + a[2]*b[2][i];
        }

        for (int i=0; i<Nn; i++) {
            s[i] = d * (dot_ac[i]/norm_c[i] - dot_ab[i]/norm_b[i]) * pow(norm_cxa[i],-2);
        }

        for (int i=0; i<Nn; i++) {
            bfield[0][i] += cxa[0][i] * s[i]; 
            bfield[1][i] += cxa[1][i] * s[i]; 
            bfield[2][i] += cxa[2][i] * s[i];
        }
    }    

    free(a); 
    free2darray(b, 3);
    free2darray(c, 3);
    free2darray(cxa, 3);
    free(norm_cxa);
    free(dot_ac);
    free(dot_ab);
    free(norm_b);
    free(norm_c);
    free(s);

    return bfield;
}


// Chunk a matrix of major length `N` for processing on `Nt` threads
//  `it` is the number of the current thread
// Modifies a 2-length array `i` of the indices for processing
void partition(int i[], int it, int Nt, int N){

    int rem = N % Nt; 
    int Nperthread = (N - rem) / Nt; 

    if (it == 0)  {
        i[0] = 0;
        i[1] = i[0] + Nperthread;
    }
    else if (it == Nt-1) {
        i[1] = N;
        i[0]= i[1] - Nperthread - rem;
    }
    else { 
        i[0] = it*Nperthread;
        i[1] = i[0] + Nperthread;
    }
}
         

// Solve with optimization level 3: parallel processing
double** solve3(double** nodes, double** sources, int Nn, int Ns) {

    double source[7];
    double** bfield = zeros(3, Nn);
    double** bfield_partial = zeros(3,Nn);
    double** sources_partial = zeros(Ns,7);
    int slice[2] = {0, 0};

    int i, j, k, it, row, Nt;
    #pragma omp parallel 
    {
        Nt = omp_get_num_threads();
    }
    int length_partition = 0;

    #pragma omp parallel private(i, j, k, it, row, slice, length_partition, bfield_partial, sources_partial) shared(bfield, sources) 
    {
        #pragma omp for ordered 
        for (it=0; it<Nt; it++) {
                partition(slice, it, Nt, Ns);
            }
    
        #pragma omp for
        for (it=0; it<Nt; it++) {
            partition(slice, it, Nt, Ns);
            length_partition = slice[1] - slice[0];

            // Create the partioned matrix 
            double** sources_partial = zeros(length_partition, 7);

            row = 0;
            for (j=slice[0]; j<slice[1]; j++) { 

                for (k=0; k<7; k++) {
                    sources_partial[row][k] = sources[j][k];
                }
                row++;
            }

            bfield_partial = solve2(nodes, sources_partial, Nn, length_partition); 
            free2darray(sources_partial, 7);

            #pragma omp critical 
            {
                for (i=0; i<3; i++) {
                    for (j=0; j<Nn; j++) {
                        bfield[i][j] += bfield_partial[i][j];
                    }
                }
            }
        }
    }
    
    return bfield;
}
