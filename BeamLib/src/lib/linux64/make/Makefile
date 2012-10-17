#########################################################################################################
#
#  Makefile for xBeam low-level library HH & RPN 07Jan2011
#
#  Usage: 'make' or 'make lib' to update the library.
#
#########################################################################################################
all: xbeamlib

lib: xbeamlib
	echo 'Updating library'
	
clean:
	rm $(SRC_LIB)  $(OBJDIR)/*.mod

################################################################################
## I/O SELECTION.

SRCDIR=../../src
OBJDIR=../obj

########################################################################################################
##  COMPILER OPTIONS.
#
# *** Compilation with g95.
  F90COMP= g95
#  F90OPTS= -O3 -fmod=$(OBJDIR) -ffree-form -Wobsolescent -Wuninitialized -Wunused-vars -ftrace=full -fzero
#  F90OPTS= -O3 -fmod=$(OBJDIR)
  F90OPTS= -O3 -fmod=$(OBJDIR) -ffree-form -Wobsolescent -Wuninitialized -Wunused-vars -ftrace=full -fzero -fPIC
  LIBS= -llapack -lblas

# F90COMP = f95
# F90OPTS=  -O3 -m64 -moddir=$(OBJDIR)
# F90OPTS=  -g -m64 -moddir=$(OBJDIR)


# *** Compilation command with SUN Studio.
#
# F90COMP= sunf95
# F90OPTS= -C -g -moddir=$(OBJDIR)
# F90OPTS= -fast -moddir=$(OBJDIR)
# LIBS= -llapack -lblas

################################################################################
## Source lists.

SRC_LIB=\
     $(OBJDIR)/lib_sparse.o \
     $(OBJDIR)/interface_lapack_v3.o \
     $(OBJDIR)/lib_out.o \
     $(OBJDIR)/lib_lu.o \
     $(OBJDIR)/lib_fem.o \
     $(OBJDIR)/lib_rot.o \
     $(OBJDIR)/lib_rotvect.o \
     $(OBJDIR)/lib_bgeom.o   \
     $(OBJDIR)/lib_cbeam3.o  \
     $(OBJDIR)/lib_xbeam.o   \
     $(OBJDIR)/lib_perturb.o 
     
#########################################################################################
## LINK

xbeamlib: $(SRC_LIB)
	ar -ru ./libxbeam.a $(SRC_LIB)

########################################################################################################
## LIBRARIES OF ROUTINES.

$(OBJDIR)/lib_rot.o: $(SRCDIR)/lib_rot.f90;
	$(F90COMP) -c $(F90OPTS) $(SRCDIR)/lib_rot.f90 -o $(OBJDIR)/lib_rot.o

$(OBJDIR)/lib_rotvect.o: $(SRCDIR)/lib_rotvect.f90;
	$(F90COMP) -c $(F90OPTS) $(SRCDIR)/lib_rotvect.f90 -o $(OBJDIR)/lib_rotvect.o

$(OBJDIR)/lib_sparse.o: $(SRCDIR)/lib_sparse.f90;
	$(F90COMP) -c $(F90OPTS) $(SRCDIR)/lib_sparse.f90 -o $(OBJDIR)/lib_sparse.o

$(OBJDIR)/interface_lapack_v3.o: $(SRCDIR)/interface_lapack_v3.f90;
	$(F90COMP) -c $(F90OPTS) $(SRCDIR)/interface_lapack_v3.f90 -o $(OBJDIR)/interface_lapack_v3.o

$(OBJDIR)/lib_out.o: $(SRCDIR)/lib_out.f90;
	$(F90COMP) -c $(F90OPTS) $(SRCDIR)/lib_out.f90 -o $(OBJDIR)/lib_out.o

$(OBJDIR)/lib_lu.o: $(SRCDIR)/lib_lu.f90;
	$(F90COMP) -c $(F90OPTS) $(SRCDIR)/lib_lu.f90 -o $(OBJDIR)/lib_lu.o

$(OBJDIR)/lib_fem.o: $(SRCDIR)/lib_fem.f90;
	$(F90COMP) -c $(F90OPTS) $(SRCDIR)/lib_fem.f90 -o $(OBJDIR)/lib_fem.o

$(OBJDIR)/lib_bgeom.o: $(SRCDIR)/lib_bgeom.f90;
	$(F90COMP) -c $(F90OPTS) $(SRCDIR)/lib_bgeom.f90 -o $(OBJDIR)/lib_bgeom.o

$(OBJDIR)/lib_cbeam3.o: $(SRCDIR)/lib_cbeam3.f90;
	$(F90COMP) -c $(F90OPTS) $(SRCDIR)/lib_cbeam3.f90 -o $(OBJDIR)/lib_cbeam3.o

$(OBJDIR)/lib_xbeam.o: $(SRCDIR)/lib_xbeam.f90;
	$(F90COMP) -c $(F90OPTS) $(SRCDIR)/lib_xbeam.f90 -o $(OBJDIR)/lib_xbeam.o

$(OBJDIR)/lib_perturb.o: $(SRCDIR)/lib_perturb.f90;
	$(F90COMP) -c $(F90OPTS) $(SRCDIR)/lib_perturb.f90 -o $(OBJDIR)/lib_perturb.o

############################################################################################################
