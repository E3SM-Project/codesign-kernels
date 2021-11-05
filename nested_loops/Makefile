FC = 
FFLAGS =
LDFLAGS =

dummy:
	@echo "ERROR: Unknown target"

# On spock, both gnu and cray compilers use ftn wrapper from Cray PrgEnv
cray-cpu:
	make nested \
	"FC=ftn" \
	"FFLAGS= -O3" \
	"LDFLAGS= -O3"

cray-omp:
	make nested \
	"FC=ftn" \
	"FFLAGS=-DUSE_OMPOFFLOAD -homp -fopenmp -O3" \
	"LDFLAGS=-homp -fopenmp"

craygnu-omp:
	make nested \
	"FC=ftn" \
	"FFLAGS= -DUSE_OMPOFFLOAD -fopenmp -O3" \
	"LDFLAGS= -fopenmp"

xlf-omp:
	make nested \
	"FC=mpif90" \
	"FFLAGS=-DUSE_OMPOFFLOAD -qsmp -qoffload -O3" \
	"LDFLAGS=-qsmp -qoffload -lcudart -L$(CUDA_DIR)/lib64"

nvhpc-cpu:
	make nested \
	"FC=mpif90" \
	"FFLAGS=-g -O3" \
	"LDFLAGS=-O3"

nvhpc-acc:
	make nested \
	"FC=mpif90" \
	"FFLAGS=-DUSE_OPENACC -O3 -acc -Minfo=accel -ta=tesla:cc70,cc60,deepcopy,nollvm" \
	"LDFLAGS=-O3 -acc"

ifort:
	make nested \
	"FC=ifort" \
	"FFLAGS=-O3 -align array64byte -qopt-report=5" \
	"LDFLAGS=-O3"

gnu:
	make nested \
	"FC=gfortran" \
	"FFLAGS=-O3 -m64" \
	"LDFLAGS=-O3 -m64"

llvm:
	make nested \
	"FC=flang" \
	"FFLAGS=-O3" \
	"LDFLAGS=-O3"

nested: nested.o timerMod.o
	@echo "Linking"
	$(FC) $(LDFLAGS) -o nested nested.o timerMod.o

nested.o: nested.F90 timerMod.o
	@echo "Building nested source"
	$(FC) -c $(FFLAGS) nested.F90

timerMod.o: timerMod.f90
	@echo "Building timer module"
	$(FC) -c $(FFLAGS) timerMod.f90

clean:
	/bin/rm -f nested *.o *optrpt *.s *.mod *.MOD

