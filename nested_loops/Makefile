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
	"FC=nvfortran" \
	"FFLAGS=-DUSE_OPENACC -DNO_MPI -O3 -acc -Minfo=accel -gpu=cc70,cc60,deepcopy" \
	"LDFLAGS=-O3 -acc"

nvhpc-omp:
	make nested \
	"FC=nvfortran" \
	"FFLAGS=-DUSE_OPENMP -DNO_MPI -O3 -mp=gpu -Minfo=accel -gpu=cc70,deepcopy" \
	"LDFLAGS=-O3 -mp=gpu"

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
	"FFLAGS=-O3 -DNO_MPI" \
	"LDFLAGS=-O3"

nested-yakl-amd:
	make nestedyakl \
	"FC=flang" \
	"CXX=hipcc" \
	"YARCH=HIP" \
	"CPPFLAGS=-DUSE_YAKL -DNO_MPI -P -traditional" \
	"FFLAGS=-g -O3 -m64 -fopenmp -Wall -pedantic -DUSE_YAKL -DNO_MPI -DF90_PACK_SIZE=4" \
	"CXXFLAGS=-g -O3 -m64 -fopenmp -Wall -pedantic -std=c++17 -DYAKL_ARCH_HIP -D__HIP_ROCclr__ -D__HIP_ARCH_GFX908__=1 --rocm-path=${ROCM_PATH} --offload-arch=gfx908 -x hip -I./YAKL/src" \
	"LDFLAGS=-O3 -m64 -L${HIP_PATH}/lib -L./YAKL/build -lyakl -lamdhip64 -lstdc++"

nested-yakl-amd-cpu:
	make nestedyakl \
	"FC=flang" \
	"CXX=clang++" \
	"YARCH=" \
	"CPPFLAGS=-DUSE_YAKL -DNO_MPI -P -traditional" \
	"FFLAGS=-g -O3 -m64 -fopenmp -Wall -pedantic -DUSE_YAKL -DNO_MPI -DF90_PACK_SIZE=4" \
	"CXXFLAGS=-g -O3 -m64 -fopenmp -Wall -pedantic -std=c++17 -I./YAKL/src" \
	"LDFLAGS=-O3 -m64 -L${HIP_PATH}/lib -L./YAKL/build -lyakl -lstdc++"

	
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


#
# Do not use these targets directly
#

nested: nested.o timerMod.o nested_vars.o
	@echo "Linking"
	$(FC) $(LDFLAGS) -o nested nested.o timerMod.o nested_vars.o

libYAKL:
	if [ -e ./YAKL/.git ]; then \
		(echo "Building YAKL..."); \
		(cd YAKL; git submodule update --init --recursive; mkdir -p build; cd build; FFLAGS="" LDFLAGS="" CXXFLAGS="" cmake -DYAKL_ARCH=${YARCH} ..; make -j 8 VERBOSE=ON; cd ../..) \
	else \
		(echo "Missing ./YAKL/.git, did you forget to 'git submodule update --init --recursive' ?"; exit 1) \
	fi
	
nestedyakl: libYAKL nested.o nested_yakl_c.o nested_vars.o
	@echo "Linking"
	$(FC) $(LDFLAGS) -o nested nested.o timerMod.o nested_vars.o nested_yakl_c.o

nestedcke: nested.o timerMod.o $(CKE_OBJECTS)
	@echo "Linking"
	$(FC) -o nested nested.o timerMod.o $(CKE_OBJECTS) $(LDFLAGS)

nested.o: nested.F90 timerMod.o cke_mod.o nested_vars.o
	@echo "Building nested source"
	$(FC) -c $(FFLAGS) nested.F90

timerMod.o: timerMod.f90
	@echo "Building timer module"
	$(FC) -c $(FFLAGS) timerMod.f90

nested_vars.o: nested_vars.F90
	@echo "Building nested_vars module"
	$(CPP) $(CPPFLAGS) $< > $*.f90
	$(FC) $(FFLAGS) -c $*.f90
	rm $*.f90

cke_mod.o: cke_mod.F90
	@echo "Building cke module"
	$(FC) -c $(FFLAGS) cke_mod.F90

.cpp.o:
	$(CXX) $(CXXFLAGS) -c $< -o $@

cke.o: cke.hpp cke_impl.hpp
cke_impl1.o: cke.hpp cke_impl.hpp
cke_impl2.o: cke.hpp cke_impl.hpp

clean:
	/bin/rm -rf nested *.o *optrpt *.s *.mod *.MOD ./YAKL/build/*
