FC = 
FFLAGS =
LDFLAGS =

dummy:
	@echo "ERROR: Unknown target"

xlf-summit-omp-offload:
	make nested \
	"FC=xlf90_r" \
	"FFLAGS=-qsmp -qoffload -O3" \
	"LDFLAGS=-qsmp -qoffload -lcudart -L$(CUDA_DIR)/lib64"

pgi-cpu:
	make nested \
	"FC=pgf90" \
	"FFLAGS=-g -O3" \
	"LDFLAGS=-O3"

pgi-acc:
	make nested \
	"FC=pgf90" \
	"FFLAGS=-O3 -acc -Minfo=accel -ta=tesla:cc70,cc60,deepcopy,nollvm" \
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

nested.o: nested.f90 timerMod.o
	@echo "Building nested source"
	$(FC) -c $(FFLAGS) nested.f90

timerMod.o: timerMod.f90
	@echo "Building timer module"
	$(FC) -c $(FFLAGS) timerMod.f90

clean:
	/bin/rm nested *.o *optrpt

