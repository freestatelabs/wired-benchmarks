gcc -O3 -fopenmp -ffast-math -march=native benchmark.c wired.c -o gcc.o -lm
icx -O3 -qopenmp -march=native -ffast-math -qopt-report benchmark.c wired.c -o icx.o -lm
# clang -O3 -fopenmp -march=znver5 -mavx512f -mavx512dq -mavx512cd -mavx512bw -mavx512vl -ffast-math benchmark.c wired.c -o aocc.exe -lm
clang -O3 -fopenmp -march=znver5 -ffast-math -fgen-aor benchmark.c wired.c -o aocc.o -lm