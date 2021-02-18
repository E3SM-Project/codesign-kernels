F90 := ifort
F90ACC := pgf90

nested: nested.f90
	cpp -P -traditional nested.f90 > nested-ref.f90
	cpp -P -traditional -DORIG nested.f90 > nested-orig.f90
	$(F90) -o nested-ref.x nested-ref.f90  -O3 -g -convert big_endian \
			-cpp -free -align array64byte -qopt-report=5 -xHost
	$(F90) -o nested-orig.x nested-orig.f90  -O3 -g -convert big_endian \
			-cpp -free -align array64byte -qopt-report=5 -xHost
	
nested-acc: nested.f90
	cpp -P -traditional -Uvector nested.f90 > nested-accref.f90
	cpp -P -traditional -Uvector -DORIG nested.f90 > nested-accorig.f90
	$(F90ACC) -o nested-accref.x nested-accref.f90  -O3 -g -acc -Minfo=accel \
		-ta=tesla:cc70,cc60,deepcopy,nollvm -byteswapio -Mfree
	$(F90ACC) -o nested-accorig.x nested-accorig.f90  -O3 -g -acc -Minfo=accel \
		-ta=tesla:cc70,cc60,deepcopy,nollvm -byteswapio -Mfree

clean:
	rm -f nested-ref.f90 nested-orig.f90 \
		nested-accref.f90 nested-accorig.f90 *.x *optrpt	
