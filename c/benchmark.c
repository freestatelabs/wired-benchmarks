/*  Benchmark the Wired program in the C language
    (c) 2024 ryan@freestatelabs.com

    Run this program from the command line with four arguments: 
    $ benchmark.exe NUM_NODES NUM_SOURCES NUM_ITERATIONS NUM_THREADS
*/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>
#include "wired.h"


int main(int argc, char* argv[]) {

    int Nn = atoi(argv[1]);     
    int Ns = atoi(argv[2]);
    int Ni = atoi(argv[3]);
    int Nt = atoi(argv[4]);

    double tstart, tend, duration;

    printf("\nBiot-Savart benchmark program start: \n");
    printf("\t Nodes:   \t %d\n", Nn);
    printf("\t Sources: \t %d\n\n", Ns);
    printf("\t Threads: \t %d\n", Nt);

    omp_set_num_threads(Nt);

    printf("Testing threads:");
    #pragma omp parallel 
    #pragma omp for ordered
    for (int i=0; i<Nt; i++)
    {
        #pragma omp ordered
        printf(" %d", omp_get_thread_num());
    }
    printf("\n");

    double** nodes = createnodes(Nn);
    double** sources = createsources(Ns);
    double** bfield;

    tstart = omp_get_wtime();
    for (int i=0; i<Ni; i++) {
        if (Nt > 1) {
            bfield = solve3(nodes, sources, Nn, Ns);
        }
        else {
            bfield = solve2(nodes, sources, Nn, Ns);
        }

        if (i < Ni-1) {
            free2darray(bfield, 3);
        }
    }
    tend = omp_get_wtime();
    duration = (tend-tstart)/Ni;
    printf("Calculation duration: %f\n\n", duration);

    printf("Bfield: \n");
    for (int i = 0; i<10; i++) {
        for (int j = 0; j < 3; j++){
            printf("%f\t", bfield[j][i]);
        }
        printf("\n");
    }

    // printf("Nodes: \n");
    // for (int i = 0; i<10; i++) {
    //     for (int j = 0; j <3; j++){
    //         printf("%f\t", nodes[i][j]);
    //     }
    //     printf("\n");
    // }

    // printf("Sources: \n");
    // for (int i = 0; i<10; i++) {
    //     for (int j = 0; j <7; j++){
    //         printf("%f\t", sources[i][j]);
    //     }
    //     printf("\n");
    // }

    free2darray(bfield, 3);
    free2darray(nodes, 3); 
    free2darray(sources, Ns);
}