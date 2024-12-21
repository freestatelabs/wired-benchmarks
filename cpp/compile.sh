# g++ -O3 -fopenmp -march=native benchmark.cpp wired.cpp -o gcc.o -lm
icpx -O3 -qopenmp -march=native -ffast-math benchmark.cpp main.cpp -o icx.o -lm
# clang -O3 -fopenmp -march=znver5 -flto -fremap-arrays -zopt -ffast-math benchmark.cpp wired.cpp -o aocc.o -lm