## Build

Serial:
gfortran -ffree-line-length-none biharmonic_wk_kernel.F90

GPU: 
pgf90 -O3 -Minfo=acc -ta=tesla,cuda7.5,pinned -acc biharmonic_wk_kernel.F90 

## Run
Serial:
./a.out

Titan:
export OMP_NUM_THREADS=16 
aprun -d 16 ./a.out
