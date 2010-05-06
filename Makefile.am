COM_CFLAGS = -O2 -W -Wall

bin_PROGRAMS = qcd_prop qcd_pcac

qcd_prop_SOURCES = prop.c qcd_arg_parse.h qcd_arg_parse.c
qcd_prop_CFLAGS = $(COM_CFLAGS)

qcd_pcac_SOURCES = pcac.c qcd_arg_parse.h qcd_arg_parse.c
qcd_pcac_CFLAGS = $(COM_CFLAGS)

ACLOCAL_AMFLAGS = -I .buildutils/m4
