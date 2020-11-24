# Makefile created by mkmf $Id: mkmf,v 18.0.18.4.2.1.2.3 2013/08/21 15:04:06 Seth.Underwood Exp $ 

SRCROOT = ../src.omp

include ../nvhpc.mk


.DEFAULT:
	-echo $@ does not exist.
all: nh_core.x
Main.o: ${SRCROOT}/Main.F90 nh_core.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	${SRCROOT}/Main.F90
nh_core.o: ${SRCROOT}/nh_core.F90
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	${SRCROOT}/nh_core.F90
./nh_core.F90: ${SRCROOT}/nh_core.F90
	cp ${SRCROOT}/nh_core.F90 .
./Main.F90: ${SRCROOT}/Main.F90
	cp /lustre/f1/unswept/Zhi.Liang/tests/gfdl_nh_core/src/./Main.F90 .
SRC = ${SRCROOT}/nh_core.F90 ${SRCROOT}/Main.F90
OBJ = nh_core.o Main.o
OFF = ${SRCROOT}/nh_core.F90 ${SRCROOT}/Main.F90
clean: neat
	-rm -f .nh_core.x.cppdefs $(OBJ) nh_core.x
neat:
	-rm -f $(TMPFILES)
localize: $(OFF)
	cp $(OFF) .
TAGS: $(SRC)
	etags $(SRC)
tags: $(SRC)
	ctags $(SRC)
nh_core.x: $(OBJ) 
	$(LD) $(OBJ) -o nh_core.x  $(LDFLAGS)
