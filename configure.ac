# Process this file with autoconf to produce a configure script.

# Initialization
AC_PREREQ([2.64])
AC_INIT([QCDutils],[2.0rc1],[antonin.portelli@me.com],[qcdutils])
AC_CONFIG_AUX_DIR([.buildutils])
AC_CONFIG_SRCDIR([prop.c])
AC_CONFIG_MACRO_DIR([.buildutils/m4])
AM_INIT_AUTOMAKE([foreign -Wall -Werror])
AM_CONFIG_HEADER([config.h])

# Check for OpenMP
AM_CFLAGS="$AM_CFLAGS $OPENMP_CFLAGS"
CFLAGS="$AM_CFLAGS $CFLAGS"

# Checks for programs.
AC_PROG_CC
AM_PROG_CC_C_O

# Configure parameters non-standard install prefix of libs
AC_ARG_WITH([latan],
    [AS_HELP_STRING([--with-latan=prefix],
	[try this for a non-standard install prefix of the LatAnalyze library])],
    [AM_CFLAGS="$AM_CFLAGS -I$with_latan/include"]
	[AM_LDFLAGS="$AM_LDFLAGS -L$with_latan/lib"],
	[]
)
AC_ARG_WITH([argtable],
    [AS_HELP_STRING([--with-argtable=prefix],
	[try this for a non-standard install prefix of the argtable library])],
    [AM_CFLAGS="$AM_CFLAGS -I$with_argtable/include"]
	[AM_LDFLAGS="$AM_LDFLAGS -L$with_argtable/lib"],
	[]
)

LX_FIND_MPI

CFLAGS="$AM_CFLAGS $MPI_CFLAGS $CFLAGS"
LDFLAGS="$AM_LDFLAGS $MPI_CLDFLAGS $LDFLAGS"

# Checks for libraries.
AC_CHECK_LIB([m],[cos],[],[AC_MSG_ERROR([libm library not found])])
AC_CHECK_LIB([latan],[latan_get_version],[],[AC_MSG_ERROR([LatAnalyze library not found])])
AC_CHECK_LIB([argtable2],[arg_parse],[],[AC_MSG_ERROR([argtable library not found])])

LDFLAGS="$LDFLAGS $ARGTABLE_LDFLAGS"

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([string.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_TYPE_SIZE_T

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
