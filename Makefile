# The compiler
FC = gfortran

# flags for debugging or for maximum performance, comment as necessary
#FCFLAGS = -g -fbounds-check
#FCFLAGS = -O2
# flags forall (e.g. look for system .mod files, required in gfortran)
#FCFLAGS += -I/usr/include

FCFLAGS = -I/usr/local/include/fgsl
# Debugging
DEBUG = -g -fcheck=all -Wall

#Below is equivalent to -I/usr/include/fgsl -lfgsl -lgsl -lgslcblas -lm
LDFLAGS = `pkg-config --cflags fgsl` 

#Libraries needed for linking, unused in the examples
LDLIBS = `pkg-config --libs fgsl`


# List of executables to be built within the package
PROGRAMS = mod_input.o mod_initial.o mod_output.o mod_magnetization.o mod_dynamics.o mod_read.o mod_self_absorption.o mod_lightcurves.o main.o

# "make" builds all
all: $(PROGRAMS)
	$(FC) $(LDFLAGS) $(DEBUG) $(PROGRAMS) -o torch $(LDLIBS)

#openfile.o:open.f90
#	$(FC) $(FCFLAGS) -c open.f90


#$conclusion
# ======================================================================
# And now the general rules, these should not require modification
# ======================================================================

# General rule for building prog from prog.o; $^ (GNU extension) is
# used in order to list additional object files on which the
# executable depends
%: %.o
	$(FC) $(FCFLAGS) -o $@ $^ $(LDFLAGS)

# General rules for building prog.o from prog.f90 or prog.F90; $< is
# used in order to list only the first prerequisite (the source file)
# and not the additional prerequisites such as module or include files
%.o: %.f90
	$(FC) $(FCFLAGS) -c $<

%.o: %.F90
	$(FC) $(FCFLAGS) -c $<

# Utility targets
.PHONY: clean veryclean

clean:
	rm -f *.o *.mod *.MOD

veryclean: clean
	rm -f *~ $(PROGRAMS)
	rm -f torch
	rm -r figs data

