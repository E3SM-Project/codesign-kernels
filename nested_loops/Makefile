FC = 
FFLAGS =
LDFLAGS =

# if using C++/Kokkos, set
# * the EKAT installation root dir
EKAT = # no default
-include make.inc
CKE_FLAGS = -DUSE_CKE -I${EKAT}/include -I${EKAT}/include/kokkos
CKE_LIBS = -lstdc++ -L${EKAT}/lib -L${EKAT}/lib64 -lekat -lkokkoscore -ldl
CKE_OBJECTS = cke_mod.o cke.o cke_impl1.o cke_impl2.o

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
	"FC=mpif90" \
	"FFLAGS=-O3 -m64 -Wall -pedantic" \
	"LDFLAGS=-O3 -m64"

llvm:
	make nested \
	"FC=flang" \
	"FFLAGS=-O3" \
	"LDFLAGS=-O3"

gnu-cpu-cke:
	make nestedcke \
	"FC=mpif90" \
	"CXX=mpicxx" \
	"FFLAGS=-g -O3 -m64 -fopenmp -Wall -pedantic -DUSE_CKE -DF90_PACK_SIZE=4" \
	"CXXFLAGS=-std=c++14 -g -O3 -m64 -fopenmp -Wall -pedantic ${CKE_FLAGS} -DCKE_PACK_SIZE=4" \
	"LDFLAGS=-O3 -m64 -fopenmp ${CKE_LIBS}"

gnu-weaver-cke:
	make nestedcke \
	"FC=mpif90" \
	"CXX=mpicxx" \
	"FFLAGS=-g -O3 -m64 -fopenmp -Wall -pedantic -DUSE_CKE -DF90_PACK_SIZE=8" \
	"CXXFLAGS=-std=c++14 -g -O3 -m64 -fopenmp -Wall -pedantic \
		-arch=sm_70 --expt-extended-lambda ${CKE_FLAGS} -DCKE_PACK_SIZE=1" \
	"LDFLAGS=-O3 -m64 -fopenmp ${CKE_LIBS} -lcudadevrt -lcudart_static -lrt"

gnu-summit-cke:
	make nestedcke \
	"FC=mpif90" \
	"CXX=mpicxx" \
	"FFLAGS=-g -O3 -m64 -fopenmp -Wall -pedantic -DUSE_CKE -DF90_PACK_SIZE=8 -DNO_MPI" \
	"CXXFLAGS=-std=c++14 -g -O3 -m64 -fopenmp -Wall -pedantic \
		-arch=sm_70 --expt-extended-lambda ${CKE_FLAGS} -DCKE_PACK_SIZE=1" \
	"LDFLAGS=-O3 -m64 -fopenmp ${CKE_LIBS} -L${CUDAPATH}/lib64 -lcuda -lcudart"

intel-compy-cke:
	make nestedcke \
	"FC=mpiifort" \
	"CXX=mpiicpc" \
	"FFLAGS=-g -O3 -m64 -fopenmp -xcore-avx512 -DUSE_CKE -DNO_MPI -DF90_PACK_SIZE=16" \
	"CXXFLAGS=-std=c++14 -g -O3 -m64 -fopenmp -xcore-avx512 ${CKE_FLAGS} -DCKE_PACK_SIZE=8" \
	"LDFLAGS=-O3 -m64 -fopenmp ${CKE_LIBS}"

nested: nested.o timerMod.o
	@echo "Linking"
	$(FC) $(LDFLAGS) -o nested nested.o timerMod.o

nestedcke: nested.o timerMod.o $(CKE_OBJECTS)
	@echo "Linking"
	$(FC) -o nested nested.o timerMod.o $(CKE_OBJECTS) $(LDFLAGS)

nested.o: nested.F90 timerMod.o cke_mod.o
	@echo "Building nested source"
	$(FC) -c $(FFLAGS) nested.F90

timerMod.o: timerMod.f90
	@echo "Building timer module"
	$(FC) -c $(FFLAGS) timerMod.f90

cke_mod.o: cke_mod.F90
	@echo "Building cke module"
	$(FC) -c $(FFLAGS) cke_mod.F90

.cpp.o:
	$(CXX) $(CXXFLAGS) -c $< -o $@

cke.o: cke.hpp cke_impl.hpp
cke_impl1.o: cke.hpp cke_impl.hpp
cke_impl2.o: cke.hpp cke_impl.hpp

clean:
	/bin/rm -f nested *.o *optrpt *.s *.mod *.MOD
