COM_CFLAGS = -std=c99 -O2 -W -Wall

bin_PROGRAMS = qcd_mass qcd_pcac qcd_prop qcd_pseudosc qcd_splitsum qcd_phyfit\
qcd_error

qcd_mass_SOURCES     = mass.c qcd_arg_parse.h qcd_arg_parse.c
qcd_mass_CFLAGS      = $(COM_CFLAGS)

qcd_prop_SOURCES     = prop.c qcd_arg_parse.h qcd_arg_parse.c
qcd_prop_CFLAGS      = $(COM_CFLAGS)

qcd_pseudosc_SOURCES = pseudosc.c qcd_arg_parse.h qcd_arg_parse.c
qcd_pseudosc_CFLAGS  = $(COM_CFLAGS)

qcd_pcac_SOURCES     = pcac.c qcd_arg_parse.h qcd_arg_parse.c
qcd_pcac_CFLAGS      = $(COM_CFLAGS)

qcd_splitsum_SOURCES = splitsum.c qcd_arg_parse.h qcd_arg_parse.c
qcd_splitsum_CFLAGS  = $(COM_CFLAGS)

qcd_phyfit_SOURCES   = phyfit.c                             \
                       parameters.h parameters.c            \
                       data_loader.h data_loader.c          \
                       models.h models.c                    \
                       output.h output.c                    \
                       qed_fvol_tabfunc.h qed_fvol_tabfunc.c
qcd_phyfit_CFLAGS    = $(COM_CFLAGS)

qcd_error_SOURCES    = error.c
qcd_error_CFLAGS     = $(COM_CFLAGS)

ACLOCAL_AMFLAGS      = -I .buildutils/m4
