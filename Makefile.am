AUTOMAKE_OPTIONS = gnu
ACLOCAL_AMFLAGS = -I m4
AM_CFLAGS=${SIMD_CFLAGS} ${OPENMP_CFLAGS} ${DEBUG_FLAGS}
lib_LTLIBRARIES = libm1ri.la

libm1ri_la_SOURCES = m1ri/m3d.c\
	m1ri/m1ri_strassen.c\
	m1ri/m7d.c\
	m1ri/m5d.c\
	m1ri/m1ri_io.c \
	m1ri/misc.c\
	m1ri/m1ri_permutations.c\
	m1ri/m1ri_decom.c
	
BUILT_SOURCES = m1ri/m1ri_config.h

pkgincludesubdir = $(includedir)/m1ri
pkgincludesub_HEADERS =  m1ri/m1ri.h
	m1ri/m5d.h\
	m1ri/m3d.h\
	m1ri/m7d.h\
	m1ri/m1ri_strassen.h\
	m1ri/m1riwrappers.h\
	m1ri/m1ri_io.h \
	m1ri/misc.h\
	m1ri/m1ri_permutations.h\
	m1ri/m1ri_decom.h

nodist_pkgincludesub_HEADERS = m1ri/m1ri_config.h

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = m1ri.pc

libm1ri_la_LDFLAGS = -release 0.0.$(RELEASE) -no-undefined
libm1ri_la_LIBADD = $(LIBPNG_LIBADD)



check_PROGRAMS =  test_arithmetic test_hadamard test_create_m3d test_m3d test_m5d test_m7d 

test_arithmetic_SOURCES=tests/test_arithmetic.c 
test_arithmetic_LDFLAGS=-lm1ri -lm
test_arithmetic__CFLAGS=$(AM_CFLAGS)

test_hadamard_SOURCES=tests/test_hadamard.c 
test_hadamard_LDFLAGS=-lm1ri -lm
test_hadamard__CFLAGS=$(AM_CFLAGS)


test_m3d_SOURCES=tests/test_m3d.c 
test_m3d_LDFLAGS=-lm1ri -lm
test_m3d__CFLAGS=$(AM_CFLAGS)

test_m5d_SOURCES=tests/test_m5d.c 
test_m5d_LDFLAGS=-lm1ri -lm
test_m5d__CFLAGS=$(AM_CFLAGS)

test_m7d_SOURCES=tests/test_m7d.c 
test_m7d_LDFLAGS=-lm1ri -lm
test_m7d__CFLAGS=$(AM_CFLAGS)

TESTS =  test_arithmetic test_hadamard test_create_m3d test_m3d test_m5d test_m7d 


SUBDIRS = . bench

clean-local:
	(cd tests; make clean; cd ..)
