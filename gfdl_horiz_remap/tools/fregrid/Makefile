# ------------------------------------------------------------------------------
# FMS/FRE Project: Makefile to build regridding executables
# ------------------------------------------------------------------------------
# Copyright (C) NOAA Geophysical Fluid Dynamics Laboratory, 2009-2010, 2017
# Designed and written by V. Balaji, Amy Langenhorst, Aleksey Yakovlev and
# Seth Underwood
# ------------------------------------------------------------------------------

# MPICC and CC default settings
MPICC    := nvc
CC       := nvc

# Get the site specific settings, CC and MPICC may be overridden
include ../../nvhpc.mk

# Default Directories
SRCDIR := .
PREFIX := .
PREFIX_EXEC = $(PREFIX)/$(SITE)/bin
PREFIX_SHARE = $(PREFIX)/share
PREFIX_DOCS = $(PREFIX_SHARE)
PREFIX_SHARE_EXEC = $(PREFIX_SHARE)/bin

VPATH = $(SRCDIR):.

CFLAGS   := -O3 -g -traceback $(CFLAGS_SITE)
CFLAGS_O2:= -O2 -g -traceback $(CFLAGS_SITE)
INCLUDES := -I${NETCDF_HOME}/include -I$(SRCDIR)/ -I$(SRCDIR)/../../shared
CLIBS     := $(LDFLAGS) -L${NETCDF_HOME}/lib64 -L${HDF5_HOME}/lib -lnetcdf -lhdf5_hl -lhdf5 -lz $(CLIBS_SITE) $(STATIC)

TARGETS  := fregrid
TARGETS_PARALLEL := fregrid_parallel
ifndef NOPARALLEL
TARGETS += $(TARGETS_PARALLEL)
endif

SOURCES  := fregrid.c bilinear_interp.c conserve_interp.c fregrid_util.c gradient_c2l.c

OBJECTS  := $(SOURCES:c=o)

SHARED_OBJS = ../../shared/mpp.o ../../shared/mpp_domain.o  \
	      ../../shared/mpp_io.o ../../shared/tool_util.o \
	      ../../shared/create_xgrid.o ../../shared/interp.o  \
	      ../../shared/mosaic_util.o ../../shared/read_mosaic.o \
              ../../shared/affinity.o

#SHARED_HDRS = ../../shared/constant.h  \
#	       ../../shared/mosaic_util.h ../../shared/mpp.h \
#	       ../../shared/read_mosaic.h ../../shared/create_xgrid.h \
#	       ../../shared/interp.h ../../shared/mpp_domain.h \
#	       ../../shared/mpp_io.h ../../shared/tool_util.h

HEADERS = ./bilinear_interp.h ./conserve_interp.h ./fregrid_util.h ./globals.h ./gradient_c2l.h

all: $(TARGETS)

fregrid: $(OBJECTS)
	$(CC) -o $@ $^ $(SHARED_OBJS) $(CLIBS)

fregrid_parallel: $(OBJECTS)
	$(MPICC) -o $@ $^ $(SHARED_OBJS) $(CLIBS)

conserve_interp.o: conserve_interp.c $(HEADERS)
	$(CC) $(CFLAGS_O2) $(INCLUDES) -c $<

bilinear_interp.o: bilinear_interp.c $(HEADERS)
	$(CC) $(CFLAGS_O2) $(INCLUDES) -c $<

gradient_c2l.o: gradient_c2l.c $(HEADERS)
	$(CC) $(CFLAGS) $(INCLUDES) -c $<

fregrid_util.o: fregrid_util.c $(HEADERS)
	$(CC) $(CFLAGS_O2) $(INCLUDES) -c $<

%.o: %.c
	$(CC) $(CFLAGS) $(INCLUDES) -c $<

install: all
	install -m 755 -d $(PREFIX)
	install -m 755 -d $(PREFIX_EXEC)
	install -m 755 $(TARGETS) $(PREFIX_EXEC)

docs:

install-docs: docs

clean:
	-rm -f *.o $(TARGETS)
