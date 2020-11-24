# Makefile created by mkmf $Id: mkmf,v 18.0.18.4.2.1.2.3 2013/08/21 15:04:06 Seth.Underwood Exp $ 

SRCROOT = ../src.omp/

include ../nvhpc.mk


.DEFAULT:
	-echo $@ does not exist.
all: fv_mapz.x
Main.o: ${SRCROOT}/Main.F90 fv_mapz.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	${SRCROOT}/Main.F90
fv_mapz.o: ${SRCROOT}/fv_mapz.F90
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	${SRCROOT}/fv_mapz.F90
./fv_mapz.F90: ${SRCROOT}/fv_mapz.F90
	cp ${SRCROOT}/fv_mapz.F90 .
./Main.F90: ${SRCROOT}/Main.F90
	cp /lustre/f1/unswept/Zhi.Liang/tests/gfdl_fv_mapz/src/./Main.F90 .
SRC = ${SRCROOT}/fv_mapz.F90 ${SRCROOT}/Main.F90
OBJ = fv_mapz.o Main.o
OFF = ${SRCROOT}/fv_mapz.F90 ${SRCROOT}/Main.F90
clean: neat
	-rm -f .fv_mapz.x.cppdefs $(OBJ) fv_mapz.x
neat:
	-rm -f $(TMPFILES)
localize: $(OFF)
	cp $(OFF) .
TAGS: $(SRC)
	etags $(SRC)
tags: $(SRC)
	ctags $(SRC)
fv_mapz.x: $(OBJ) 
	$(LD) $(OBJ) -o fv_mapz.x  $(LDFLAGS)
