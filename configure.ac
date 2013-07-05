#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([libm1ri], [0.01])
AC_CANONICAL_SYSTEM
AC_CONFIG_MACRO_DIR([m4])
AM_INIT_AUTOMAKE
AC_PROG_LIBTOOL
LT_INIT
AC_PROG_INSTALL
AC_PROG_CC_C99()
AC_CONFIG_SRCDIR([m1ri/m1ri.c])

AC_CONFIG_HEADERS([config.h])
# Checks for programs.

# Checks for libraries.


# Checks for header files.
AC_CHECK_HEADERS([stdint.h stdlib.h string.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_C_INLINE
AC_TYPE_SIZE_T

# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_REALLOC

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
