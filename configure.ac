#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT(m1ri,1)
AC_CONFIG_SRCDIR([m1ri/config.h.in])
AC_CONFIG_HEADERS([m1ri/config.h])

# Checks for programs.
AC_PROG_AWK
AC_PROG_CC_C99
AC_PROG_LIBTOOL
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
   AX_OPENMP()
])
AC_SUBST(OPENMP_CFLAGS)
if test -n "$OPENMP_CFLAGS"; then
  M1RI_HAVE_OPENMP=1
else
  M1RI_HAVE_OPENMP=0
fi
AC_SUBST(m1RI_HAVE_OPENMP)




# SSE2 support
AC_ARG_ENABLE([sse2],
        AS_HELP_STRING([--disable-sse2], [don't use SSE2 instruction set.]),
	, [if test "$m1ri_wrapword" = "yes"; then enable_sse2="no"; else enable_sse2="yes"; fi])

AS_IF([test "x$enable_sse2" != "xno"], [
   if test "$m1ri_wrapword" = "yes"; then
      AC_MSG_ERROR([SSE2 cannot be supported when wrapping word in a C++ class.])
   fi
   case $host_cpu in i[[3456]]86*|x86_64*)
          AX_CPU_VENDOR()
          if test "x$ax_cv_cpu_vendor" = "xIntel"; then
             AX_EXT() # SSE2 is slower on the Opteron
          fi
   esac
])
if test x"$ax_cv_have_sse2_ext" = x"yes"; then
  m1ri_HAVE_SSE2=1
else
  m1ri_HAVE_SSE2=0
fi
# PNG

have_libpng="no"
AC_ARG_ENABLE([png],
   [AC_HELP_STRING([--disable-png], [disable PNG support @<:@default=enabled@:>@])],
   [
    if test "x${enableval}" = "xyes" ; then
       want_png="yes"
    else
       want_png="no"
    fi
   ],
   [want_png="yes"])

AC_MSG_CHECKING([whether to build with PNG support])
AC_MSG_RESULT([${want_png}])

if test "x${want_png}" = "xyes" ; then
   PKG_CHECK_MODULES([PNG], [libpng],
      [have_libpng="yes"; LIBPNG_LIBADD=`pkg-config --libs libpng`],
      [have_libpng="no"])
   if ! test "x${have_libpng}" = "xyes" ; then
      AC_CHECK_LIB([png],
         [png_create_write_struct],
         [have_libpng="yes"; LIBPNG_LIBADD="-lpng"],
         [AC_CHECK_LIB([png14],
            [png_create_write_struct],
            [have_libpng="yes"; LIBPNG_LIBADD="-lpng14"],
            [AC_CHECK_LIB([png12],
               [png_create_write_struct],
               [have_libpng="yes"; LIBPNG_LIBADD="-lpng12"],
               [AC_CHECK_LIB([png10],
                  [png_create_write_struct],
                  [have_libpng="yes"; LIBPNG_LIBADD="-lpng10"],
                  [have_libpng="no"])
               ])
            ])
        ])
   fi
   if test "x${have_libpng}" = "xno" ; then
      AC_MSG_WARN([Can not find a usuable PNG library. Make sure that CPPFLAGS and LDFLAGS are correctly set.])
   fi
fi

if test "x${have_libpng}" = "xyes" ; then
   M1RI_HAVE_LIBPNG=1
   AC_SUBST(M1RI_HAVE_LIBPNG)
   AC_SUBST(LIBPNG_LIBADD)
else
   M1RI_HAVE_LIBPNG=0
   AC_SUBST(M1RI_HAVE_LIBPNG)
fi



AC_SUBST(m1ri_HAVE_SSE2)

# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_REALLOC

RELEASE="AC_PACKAGE_VERSION"
AC_SUBST(RELEASE)

AC_CONFIG_FILES([Makefile m1ri.pc])

AC_OUTPUT

