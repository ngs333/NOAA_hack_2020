# Make macros for the ncrc4 site

MPICC := nvc

#CFLAGS_SITE := -tp haswell -Minfo -Mbounds
#CLIBS_SITE := 

CFLAGS_SITE := -tp haswell -acc -ta=nvidia:managed -Minfo -Mbounds -g  -Mnoinline
CLIBS_SITE := -Mcuda -ta=nvidia

NETCDF_HOME := $(shell nc-config --prefix)

NOPARALLEL := t
