# Makefile created by mkmf $Id: mkmf,v 18.0.18.4.2.1.2.3 2013/08/21 15:04:06 Seth.Underwood Exp $ 

SRCROOT = ../src.acc/

include ../nvhpc.mk


.DEFAULT:
	-echo $@ does not exist.
all: fyppm.x
Main.o: ${SRCROOT}/Main.F90 fyppm.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	${SRCROOT}/Main.F90
fyppm.o: ${SRCROOT}/fyppm.F90
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	${SRCROOT}/fyppm.F90
./fyppm.F90: ${SRCROOT}/fyppm.F90
	cp ${SRCROOT}/fyppm.F90 .
./Main.F90: ${SRCROOT}/Main.F90
	cp /lustre/f1/unswept/Zhi.Liang/tests/gfdl_fyppm/src/./Main.F90 .
SRC = ${SRCROOT}/fyppm.F90 ${SRCROOT}/Main.F90
OBJ = fyppm.o Main.o
OFF = ${SRCROOT}/fyppm.F90 ${SRCROOT}/Main.F90
clean: neat
	-rm -f .fyppm.x.cppdefs $(OBJ) fyppm.x
neat:
	-rm -f $(TMPFILES)
localize: $(OFF)
	cp $(OFF) .
TAGS: $(SRC)
	etags $(SRC)
tags: $(SRC)
	ctags $(SRC)
fyppm.x: $(OBJ) 
	$(LD) $(OBJ) -o fyppm.x  $(LDFLAGS)
