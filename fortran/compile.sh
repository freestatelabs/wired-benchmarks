gfortran -o3 -fopenmp -c wired.f90 
gfortran -o3 -fopenmp benchmark.f90 wired.f90 -o benchmark.exe