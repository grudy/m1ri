#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT(m1ri,1)
AC_CONFIG_SRCDIR([m1ri/config.h.in])
AC_CONFIG_HEADERS([m1ri/config.h])

# Checks for programs.
AC_PROG_AWK
AC_PROG_CC_C99
AC_PROG_INSTALL
AC_PROG_LN_S
AC_PROG_MAKE_SET

AM_INIT_AUTOMAKE([subdir-objects])
AC_CONFIG_MACRO_DIR([m4])
LT_INIT
# Checks for libraries.
# FIXME: Replace `main' with a function in `-lm':
AC_CHECK_LIB([m], [log])

# Checks for header files.
AC_CHECK_HEADERS([stdint.h stdlib.h string.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_CHECK_HEADER_STDBOOL
AC_C_INLINE
AC_TYPE_SIZE_T
# OpenMP support
AC_ARG_ENABLE([openmp],
        AS_HELP_STRING( [--enable-openmp],[add support for OpenMP multicore support.]))

AS_IF([test "x$enable_openmp" = "xyes"], [
   AC_OPENMP()
])
AC_SUBST(OPENMP_CFLAGS)
if test -n "$OPENMP_CFLAGS"; then
  M1RI_HAVE_OPENMP=1
else
  M1RI_HAVE_OPENMP=0
fi
AC_SUBST(M1RI_HAVE_OPENMP)

# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_REALLOC

RELEASE="AC_PACKAGE_VERSION"
AC_SUBST(RELEASE)

AC_CONFIG_FILES([Makefile m1ri/m1ri_config.h m1ri.pc])
AC_OUTPUT

