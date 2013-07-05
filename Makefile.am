AUTOMAKE_OPTIONS = gnu
ACLOCAL_AMFLAGS = -I m4

AM_CFLAGS=${SIMD_CFLAGS} ${OPENMP_CFLAGS} ${DEBUG_FLAGS}

lib_LTLIBRARIES = libm1ri.la

libm1ri_la_SOURCES =   m1ri_3dt.c\
m1riarith.c\
m1ri_classical.c\
m1ri_strassen.c\
m1ri_cubes.c\
m1ri_small.c\
m1ri_combine.c\
m7d.c\
m5d.c\
m1ri_test.c


pkgincludesubdir = $(includedir)/m1ri
pkgincludesub_HEADERS = m5d.h\
m1ri_3dt.h\
m7d.c.h\
m1ri_combine.h\
m1ri_small.h\
m1ri_cubes.h\
m1ri_strassen.h\
m1riwrappers.h\
m1ri_classical.h\
m1riarith.h


