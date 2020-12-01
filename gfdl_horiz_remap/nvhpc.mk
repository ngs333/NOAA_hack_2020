# Make macros for the ncrc4 site

MPICC := nvc

CFLAGS_SITE := -tp haswell
FFLAGS_SITE := -tp haswell

CLIBS_SITE :=
FLIBS_SITE :=

NETCDF_HOME := $(shell nc-config --prefix)
HDF5_HOME := $$HDF5_ROOT

NOPARALLEL := t
