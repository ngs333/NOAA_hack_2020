# Make macros for the ncrc4 site

MPICC := nvc

CFLAGS_SITE := -msse2
FFLAGS_SITE := -msse2

CLIBS_SITE :=
FLIBS_SITE :=

NETCDF_HOME := $$NETCDF_ROOT
HDF5_HOME := $$HDF5_ROOT

NOPARALLEL := t
