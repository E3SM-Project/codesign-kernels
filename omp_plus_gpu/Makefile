.SUFFIXES: .cxx

tracers-yakl-amd:
	make tracers \
	"CXX=hipcc" \
	"YARCH=HIP" \
	"CPPFLAGS=-DUSE_YAKL -DNO_MPI -P -traditional" \
	"FFLAGS=-g -O3 -m64 -fopenmp -Wall -pedantic -DUSE_YAKL -DNO_MPI -DF90_PACK_SIZE=4" \
	"CXXFLAGS=-g -O3 -m64 -fopenmp -Wall -pedantic -std=c++17 -DYAKL_ARCH_HIP -D__HIP_ROCclr__ -D__HIP_ARCH_GFX90a__=1 --rocm-path=${ROCM_PATH} --offload-arch=gfx90a -x hip -I./YAKL/src -I./YAKL/external" \
	"LDFLAGS=-O3 -m64 -L${ROCM_PATH}/lib -L./YAKL/build -lyakl -lamdhip64 -fopenmp -lroctx64 -lrocfft "

tracers: libYAKL tracers.o
	$(CXX) $(LDFLAGS) -o tracers.x tracers.o

libYAKL:
	if [ -e ./YAKL/.git ]; then \
		(echo "Building YAKL..."); \
		(cd YAKL; git submodule update --init --recursive; mkdir -p build; cd build; FFLAGS="" LDFLAGS="" CXXFLAGS="" cmake -DYAKL_ARCH=${YARCH} ..; make -j 8 VERBOSE=ON; cd ../..) \
	else \
		(echo "Missing ./YAKL/.git, did you forget to 'git submodule update --init --recursive' ?"; exit 1) \
	fi

.cxx.o:
	$(CXX) $(CXXFLAGS) -c $< -o $@
