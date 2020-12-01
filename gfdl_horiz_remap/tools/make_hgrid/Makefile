# ------------------------------------------------------------------------------
# FMS/FRE Project: Makefile to build make_hgrid
# ------------------------------------------------------------------------------
# Copyright (C) NOAA Geophysical Fluid Dynamics Laboratory, 2009-2011, 2017
# This program is distributed under the terms of the GNU General Public
# License. See the file COPYING contained in this directory
# ------------------------------------------------------------------------------
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

CFLAGS   := -O3 -traceback $(CFLAGS_SITE)
CFLAGS_O2:= -O2 -traceback $(CFLAGS_SITE)
INCLUDES := -I${NETCDF_HOME}/include -I$(SRCDIR)/ -I$(SRCDIR)/../../shared
CLIBS     := $(LDFLAGS) $(shell nc-config --libs) $(CLIBS_SITE) $(STATIC)

TARGETS  := make_hgrid.x
TARGETS_PARALLEL := make_hgrid_parallel
ifndef NOPARALLEL
TARGETS += $(TARGETS_PARALLEL)
endif

SOURCES  := make_hgrid.c create_conformal_cubic_grid.c create_gnomonic_cubic_grid.c
SOURCES  += create_grid_from_file.c create_lonlat_grid.c

OBJECTS  := $(SOURCES:c=o)

SHARED_OBJS = ../../shared/mpp.o ../../shared/mpp_domain.o  \
	      ../../shared/mpp_io.o ../../shared/tool_util.o \
	      ../../shared/create_xgrid.o ../../shared/interp.o  \
	      ../../shared/mosaic_util.o ../../shared/read_mosaic.o \
              ../../shared/affinity.o

#SHARED_HDRS = ../../shared/constant.h \
#	  ../../shared/mosaic_util.h ../../shared/mpp.h \
#	  ../../shared/read_mosaic.h ../../shared/create_xgrid.h \
#	  ../../shared/interp.h ../../shared/mpp_domain.h \
#	  ../../shared/mpp_io.h ../../shared/tool_util.h

HEADERS = ./create_hgrid.h

all: $(TARGETS)

make_hgrid.x: $(OBJECTS)
	$(CC) -o $@ $^ $(SHARED_OBJS) $(CLIBS)

make_hgrid_parallel: $(OBJECTS)
	$(MPICC) -o $@ $^ $(SHARED_OBJS) $(CLIBS)

make_hgrid.o: make_hgrid.c
	$(CC) $(CFLAGS) $(INCLUDES) -c $<

create_gnomonic_cubic_grid.o: create_gnomonic_cubic_grid.c
	$(CC) $(CFLAGS) $(INCLUDES) -c $<

create_conformal_cubic_grid.o: create_conformal_cubic_grid.c
	$(CC) $(CFLAGS) $(INCLUDES) -c $<

create_lonlat_grid.o: create_lonlat_grid.c
	$(CC) $(CFLAGS) $(INCLUDES) -c $<

create_grid_from_file.o: create_grid_from_file.c
	$(CC) $(CFLAGS) $(INCLUDES) -c $<

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
