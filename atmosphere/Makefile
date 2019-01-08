# Choose compiler when invoking make : make gnu=1 or pgiacc=1
ifdef gnu
	FC = gfortran
	FFLAGS = -O3 -ffree-line-length-none 
endif

ifdef intel
	FC = ifort
	FFLAGS += -O3 
endif	

ifdef pgi
	FC = pgf90
	FFLAGS += -fast 
endif	

ifdef pgiacc
	FC = pgf90
	FFLAGS += -fast -Minfo=acc -ta=tesla,cuda7.5,pinned -acc
endif	

ifdef cray
	FC = ftn
	FFLAGS += -O3
endif	

ifdef arm
	FC = armflang
	FFLAGS += -O3 -mcpu=thunderx2t99
endif	

ifdef armsve
	FC = armflang
	FFLAGS += -O3 -march=armv8-a+sve
endif	

LDFLAGS = 
# Name of executable
EXE = atm

.SUFFIXES: .c .f .f90 .F90 .o

# Source/Object Files
HEADERS  = $(wildcard *.h)
C_SRCS   = $(wildcard *.c)
F_SRCS   = $(wildcard *.F90)
C_OBJS   = $(patsubst %.c,%.o,$(wildcard *.c))
F_OBJS   = $(patsubst %.F90,%.o,$(wildcard *.F90))

OBJS =  $(F_OBJS)

all: $(EXE) 

$(EXE): $(HEADERS) $(OBJS)
	$(FC) $(FFLAGS) -o $@ $(OBJS) $(LDFLAGS)

%.o: %.f90
	$(FC) $(FFLAGS) -c $<

%.o: %.F90
	$(FC) $(FFLAGS) -c $<

# Utility targets
.PHONY: clean 

clean:
	rm -f $(EXE) *.o *.mod

