gfortran -o3 -fopenmp -march=native -c wired.f90 
gfortran -o3 -fopenmp -march=native benchmark.f90 wired.f90 -o gcc.exe

ifx -o3 -qopenmp -march=native -c wired.f90 
ifx -o3 -qopenmp -march=native benchmark.f90 wired.f90 -o ifx.exe

flang -o3 -fopenmp -march=native -c wired.f90 
flang -o3 -fopenmp -march=native benchmark.f90 wired.f90 -o flang.exe