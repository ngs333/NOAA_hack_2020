# ------------------------------------------------------------------------------
# FMS/FRE Project: Makefile to build make_solo_mosaic
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
CLIBS     := -L${NETCDF_HOME}/lib64 -L${HDF5_HOME}/lib -lnetcdf -lhdf5_hl -lhdf5 -lz $(CLIBS_SITE) $(STATIC)

TARGETS  := make_solo_mosaic

SOURCES  := make_solo_mosaic.c get_contact.c

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

HEADERS =  ./get_contact.h

all: $(TARGETS)

make_solo_mosaic: $(OBJECTS)
	$(CC) -o $@ $^ $(SHARED_OBJS) $(CLIBS)

make_solo_mosaic.o: make_solo_mosaic.c $(HEADERS)
	$(CC) $(CFLAGS) $(INCLUDES) -c $<

get_contact.o: get_contact.c $(HEADERS)
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
