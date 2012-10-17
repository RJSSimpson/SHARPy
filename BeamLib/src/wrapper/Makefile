
#-------------------------------------------------------------------------------
# Description: create a shared object from a collection of F90 source files
# 
# Usage: make [-f <filename>] [target]
# target:
#   clean       clean up object files and modules
#   'blank'     create a shared object
#
# Remarks:
#   -DNOLAPACK is used for conditional compilation
#   -cpp is used for preprocessing/conditional compilation
#-------------------------------------------------------------------------------

# Define dynamic library name
DYNLIB=../../bin/BeamLib.so

# Define flags
F90FLAG=-fPIC -Wextra -cpp -DNOLAPACK
LIBFLAG=-fPIC -shared -ffree-form

# Define path relative to cwd
PATHFOR=../fortran
PATHLIB=../lib

# Define source files
F90SRC=$(PATHFOR)/install/src/xbeam_shared.f90 \
       $(PATHLIB)/src/lib_sparse.f90 \
       $(PATHLIB)/src/lib_lu.f90 \
       $(PATHLIB)/src/lib_rot.f90 \
       $(PATHLIB)/src/lib_fem.f90 \
       $(PATHLIB)/src/lib_rotvect.f90 \
       $(PATHLIB)/src/lib_bgeom.f90 \
       $(PATHFOR)/install/src/xbeam_undef.f90 \
       $(PATHLIB)/src/lib_xbeam.f90 \
       $(PATHLIB)/src/lib_out.f90 \
       $(PATHLIB)/src/lib_cbeam3.f90 \
       $(PATHFOR)/install/src/cbeam3_asbly.f90 \
       $(PATHFOR)/install/src/cbeam3_solv.f90 \
       $(PATHLIB)/src/lib_out.f90 \
       $(PATHFOR)/install/src/xbeam_asbly.f90 \
       FxInput_PyFx.f90 \
       FxWrapper_PyFx.f90

# Define object files
F90OBJ=$(F90SRC:.f90=.o)

# Define targets
all: $(F90OBJ)
	gfortran -o $(DYNLIB) $^ $(LIBFLAG)
	rm -vf *.mod

%.o: %.f90
	gfortran -o $@ -c $< $(F90FLAG)

clean:
	rm -vf $(F90OBJ) *.mod $(DYNLIB)

