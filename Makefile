AUTOMAKE_OPTIONS = gnu
ACLOCAL_AMFLAGS = -I m4

AM_CFLAGS=${SIMD_CFLAGS} ${OPENMP_CFLAGS} ${DEBUG_FLAGS}

lib_LTLIBRARIES = libm1ri.la

libm1ri_la_SOURCES = \
	m1ri_3dt.h \



pkgincludesubdir = $(includedir)/m1ri
pkgincludesub_HEADERS = m	 m1riwrappers.h \
	 m1rielarith.h \
	 m1ristrassen.h \ 
	 m1ri_small.h \
 		m1ri_classical.h \
	 m7d.h \
	 m5d.h \
	 m1ri_cubes.h \
	 m1ri_combine.h \



