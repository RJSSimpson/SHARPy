#
#  Makefile for XBeam. HH & RPN 27/10/2011
#
#  Change TESTCASE to consider different test problems.
#
#  To link everything in C (using wrapC), only the compiler combination gcc and f95 works.
#
# Usage: 'make exe' to create the exe from sources (good for debugging)
#        'make lib' to update the library.
#        'make so' to create a shared library.
#        'make libexe' to create the exe from a library
#        'make wrapC' to create the exe from a Fortran library in C
#
###############################################################################
all: xbeam

exe: xbeam

clean:
	rm -f $(SRCS_ALL) $(SRC_WRAP) $(OBJDIR)/*.mod $(OBJDIR)/*.a

world: clean all

lib: xbeamlib
	echo 'Updating library'

so: xbeamso
	echo 'Updating shared library'

libexe: xbeamlib xbeamlibexe
	echo 'Creating xBeam from Library'

wrapC: xbeamlib xbeamC
	echo 'Creating xBeam in C from F90 Library'

################################################################################
## I/O SELECTION.

TESTCASE=rob
TESTDIR=../../main

ROOTDIR=../../..

SRCDIR=$(ROOTDIR)/fortran/install/src
OBJDIR=$(ROOTDIR)/fortran/install/linux/obj
WRAPDIR=$(ROOTDIR)/fortran/wrapper

LIBSSRCDIR=$(ROOTDIR)/lib/src
XBEAMLIBDIR=-L$(ROOTDIR)/fortran/install/linux/obj
LIBSDIR=-L/opt/lib

#######################################################################################################
##  COMPILER OPTIONS.

# *** Compilation with g95.
#F90COMP= g95
#F90OPTS= -O3 -fmod=$(OBJDIR) -ffree-form -Wobsolescent -Wuninitialized -Wunused-vars -ftrace=full -fzero

# *** Compilation with gfortran.
#CCCOMP=  gcc
#
F90COMP= gfortran
#F90OPTS= -fPIC -O3 -m64 -J $(OBJDIR)
F90OPTS= -g -m64 -cpp -fPIC -J $(OBJDIR)

# *** Compilation command with SUN Studio.
#CCCOMP= suncc
#
#F90COMP= sunf95
#F90OPTS= -C -g -moddir=$(OBJDIR)
#F90OPTS= -fast -moddir=$(OBJDIR)

LIBS= -llapack -lblas
XBEAMLIB= -lxbeam

################################################################################
## Source lists.

SRC_TOOLS=\
     $(OBJDIR)/lib_sparse.o \
     $(OBJDIR)/interface_lapack_v3.o \
     $(OBJDIR)/lib_out.o \
     $(OBJDIR)/lib_lu.o  \
     $(OBJDIR)/lib_fem.o \
     $(OBJDIR)/lib_rot.o \
     $(OBJDIR)/lib_rotvect.o \
     $(OBJDIR)/lib_bgeom.o   \
     $(OBJDIR)/lib_cbeam3.o  \
     $(OBJDIR)/lib_xbeam.o	 \
     $(OBJDIR)/lib_perturb.o

SRC_XBEAM=\
	 $(OBJDIR)/xbeam_shared.o	\
     $(OBJDIR)/xbeam_undef.o	\
     $(OBJDIR)/cbeam3_asbly.o	\
     $(OBJDIR)/cbeam3_solv.o	\
	 $(OBJDIR)/xbeam_asbly.o	\
	 $(OBJDIR)/xbeam_solv.o		\
	 $(OBJDIR)/xbeam_perturb.o	\
	 $(OBJDIR)/xbeam_fdiff.o

SRC_TEST=\
     $(OBJDIR)/input_$(TESTCASE).o	\
     $(OBJDIR)/main_$(TESTCASE).o

SRC_WRAP=\
     $(OBJDIR)/input_$(TESTCASE).o	\
	 $(OBJDIR)/xbeam_input_wrap.o	\
	 $(OBJDIR)/cbeam3_static_wrap.o

SRC_MAINLIB=\
     $(OBJDIR)/main.o

SRCS_MAIN= $(SRC_TEST)
SRCS_ALL=  $(SRC_TOOLS) $(SRC_XBEAM) $(SRC_TEST)
SRCS_LIB=  $(SRC_TOOLS) $(SRC_XBEAM) $(SRC_WRAP)

#########################################################################################
## LINK

xbeam: $(SRCS_ALL)
	$(F90COMP) $(F90OPTS) -o $(ROOTDIR)/fortran/bin/xbeam.exe $(SRCS_ALL) $(LIBS) $(LIBSDIR)

xbeamlib: $(SRCS_LIB)
	ar -ru $(ROOTDIR)/fortran/install/linux/obj/libxbeam.a $(SRCS_LIB)

xbeamlibexe: $(SRCS_MAIN)
	$(F90COMP) $(F90OPTS) -o $(ROOTDIR)/fortran/bin/xbeam.exe $(SRCS_MAIN) $(XBEAMLIB) $(LIBS) $(XBEAMLIBDIR) $(LIBSDIR)

xbeamso: $(SRCS_LIB)
	$(F90COMP) -shared -o libadd.so $(SRCS_LIB) $(LIBSDIR) $(LIBS)
	
xbeamC: $(SRC_MAINLIB)
	$(CCCOMP) -o $(ROOTDIR)/fortran/bin/xbeam.exe $(SRC_MAINLIB) $(XBEAMLIB) $(LIBS) $(XBEAMLIBDIR) $(LIBSDIR)

#########################################################################################
## TEST SOURCES

$(OBJDIR)/input_$(TESTCASE).o: $(TESTDIR)/input_$(TESTCASE).f90;
	$(F90COMP) -c $(F90OPTS) $(TESTDIR)/input_$(TESTCASE).f90 -o $(OBJDIR)/input_$(TESTCASE).o

$(OBJDIR)/main_$(TESTCASE).o: $(TESTDIR)/main_$(TESTCASE).f90;
	$(F90COMP) -c $(F90OPTS) $(TESTDIR)/main_$(TESTCASE).f90 -o $(OBJDIR)/main_$(TESTCASE).o

########################################################################################################
## WRAPPER SOURCES

$(OBJDIR)/main.o: $(WRAPDIR)/main.c;
	$(CCCOMP) -c  $(WRAPDIR)/main.c -o $(OBJDIR)/main.o

$(OBJDIR)/xbeam_input_wrap.o: $(WRAPDIR)/xbeam_input_wrap.f90;
	$(F90COMP) -c $(F90OPTS) $(WRAPDIR)/xbeam_input_wrap.f90 -o $(OBJDIR)/xbeam_input_wrap.o

$(OBJDIR)/cbeam3_static_wrap.o: $(WRAPDIR)/cbeam3_static_wrap.f90;
	$(F90COMP) -c $(F90OPTS) $(WRAPDIR)/cbeam3_static_wrap.f90 -o $(OBJDIR)/cbeam3_static_wrap.o

########################################################################################################
## XBEAM SOURCES

$(OBJDIR)/xbeam_shared.o: $(SRCDIR)/xbeam_shared.f90;
	$(F90COMP) -c $(F90OPTS) $(SRCDIR)/xbeam_shared.f90 -o $(OBJDIR)/xbeam_shared.o

$(OBJDIR)/xbeam_undef.o: $(SRCDIR)/xbeam_undef.f90;
	$(F90COMP) -c $(F90OPTS) $(SRCDIR)/xbeam_undef.f90 -o $(OBJDIR)/xbeam_undef.o

$(OBJDIR)/cbeam3_asbly.o: $(SRCDIR)/cbeam3_asbly.f90;
	$(F90COMP) -c $(F90OPTS) $(SRCDIR)/cbeam3_asbly.f90 -o $(OBJDIR)/cbeam3_asbly.o

$(OBJDIR)/cbeam3_solv.o: $(SRCDIR)/cbeam3_solv.f90;
	$(F90COMP) -c $(F90OPTS) $(SRCDIR)/cbeam3_solv.f90 -o $(OBJDIR)/cbeam3_solv.o

$(OBJDIR)/xbeam_asbly.o: $(SRCDIR)/xbeam_asbly.f90;
	$(F90COMP) -c $(F90OPTS) $(SRCDIR)/xbeam_asbly.f90 -o $(OBJDIR)/xbeam_asbly.o

$(OBJDIR)/xbeam_solv.o: $(SRCDIR)/xbeam_solv.f90;
	$(F90COMP) -c $(F90OPTS) $(SRCDIR)/xbeam_solv.f90 -o $(OBJDIR)/xbeam_solv.o

$(OBJDIR)/xbeam_perturb.o: $(SRCDIR)/xbeam_perturb.f90;
	$(F90COMP) -c $(F90OPTS) $(SRCDIR)/xbeam_perturb.f90 -o $(OBJDIR)/xbeam_perturb.o

$(OBJDIR)/xbeam_fdiff.o: $(SRCDIR)/xbeam_fdiff.f90;
	$(F90COMP) -c $(F90OPTS) $(SRCDIR)/xbeam_fdiff.f90 -o $(OBJDIR)/xbeam_fdiff.o

########################################################################################################
## LIBRARIES OF ROUTINES.

$(OBJDIR)/lib_rot.o: $(LIBSSRCDIR)/lib_rot.f90;
	$(F90COMP) -c $(F90OPTS) $(LIBSSRCDIR)/lib_rot.f90 -o $(OBJDIR)/lib_rot.o

$(OBJDIR)/lib_rodrigues.o: $(LIBSSRCDIR)/lib_rodrigues.f90;
	$(F90COMP) -c $(F90OPTS) $(LIBSSRCDIR)/lib_rodrigues.f90 -o $(OBJDIR)/lib_rodrigues.o

$(OBJDIR)/lib_rotvect.o: $(LIBSSRCDIR)/lib_rotvect.f90;
	$(F90COMP) -c $(F90OPTS) $(LIBSSRCDIR)/lib_rotvect.f90 -o $(OBJDIR)/lib_rotvect.o

$(OBJDIR)/lib_sparse.o: $(LIBSSRCDIR)/lib_sparse.f90;
	$(F90COMP) -c $(F90OPTS) $(LIBSSRCDIR)/lib_sparse.f90 -o $(OBJDIR)/lib_sparse.o

$(OBJDIR)/interface_lapack_v3.o: $(LIBSSRCDIR)/interface_lapack_v3.f90;
	$(F90COMP) -c $(F90OPTS) $(LIBSSRCDIR)/interface_lapack_v3.f90 -o $(OBJDIR)/interface_lapack_v3.o

$(OBJDIR)/lib_out.o: $(LIBSSRCDIR)/lib_out.f90;
	$(F90COMP) -c $(F90OPTS) $(LIBSSRCDIR)/lib_out.f90 -o $(OBJDIR)/lib_out.o

$(OBJDIR)/lib_lu.o: $(LIBSSRCDIR)/lib_lu.f90;
	$(F90COMP) -c $(F90OPTS) $(LIBSSRCDIR)/lib_lu.f90 -o $(OBJDIR)/lib_lu.o

$(OBJDIR)/lib_fem.o: $(LIBSSRCDIR)/lib_fem.f90;
	$(F90COMP) -c $(F90OPTS) $(LIBSSRCDIR)/lib_fem.f90 -o $(OBJDIR)/lib_fem.o

$(OBJDIR)/lib_bgeom.o: $(LIBSSRCDIR)/lib_bgeom.f90;
	$(F90COMP) -c $(F90OPTS) $(LIBSSRCDIR)/lib_bgeom.f90 -o $(OBJDIR)/lib_bgeom.o

$(OBJDIR)/lib_cbeam3.o: $(LIBSSRCDIR)/lib_cbeam3.f90;
	$(F90COMP) -c $(F90OPTS) $(LIBSSRCDIR)/lib_cbeam3.f90 -o $(OBJDIR)/lib_cbeam3.o

$(OBJDIR)/lib_xbeam.o: $(LIBSSRCDIR)/lib_xbeam.f90;
	$(F90COMP) -c $(F90OPTS) $(LIBSSRCDIR)/lib_xbeam.f90 -o $(OBJDIR)/lib_xbeam.o

$(OBJDIR)/lib_perturb.o: $(LIBSSRCDIR)/lib_perturb.f90;
	$(F90COMP) -c $(F90OPTS) $(LIBSSRCDIR)/lib_perturb.f90 -o $(OBJDIR)/lib_perturb.o

############################################################################################################
