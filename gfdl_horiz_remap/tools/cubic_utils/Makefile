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

TARGETS  := make_remap_file

SOURCES  := make_remap_file.c

OBJECTS  := $(SOURCES:c=o)

SHARED_OBJS = ../../shared/mpp.o ../../shared/mpp_domain.o  \
	      ../../shared/mpp_io.o ../../shared/tool_util.o \
	      ../../shared/create_xgrid.o ../../shared/interp.o  \
	      ../../shared/mosaic_util.o ../../shared/read_mosaic.o \
	      ../../shared/affinity.o

#SHARED_HDRS = ../../shared/constant.h \
#	       ../../shared/mosaic_util.h ../../shared/mpp.h \
#	       ../../shared/read_mosaic.h ../../shared/create_xgrid.h \
#	       ../../shared/interp.h ../../shared/mpp_domain.h \
#	       ../../shared/mpp_io.h ../../shared/tool_util.h

all: $(TARGETS)

make_remap_file: $(OBJECTS)
	$(CC) -o $@ $^ $(SHARED_OBJS) $(CLIBS)

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
