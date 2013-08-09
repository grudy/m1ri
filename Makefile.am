AUTOMAKE_OPTIONS = gnu
ACLOCAL_AMFLAGS = -I m4
AM_CFLAGS=${SIMD_CFLAGS} ${OPENMP_CFLAGS} ${DEBUG_FLAGS}
lib_LTLIBRARIES = libm1ri.la

libm1ri_la_SOURCES = m1ri/m1ri_3dt.c\
	m1ri/m1riarith.c\
	m1ri/m1ri_classical.c\
	m1ri/m1ri_strassen.c\
	m1ri/m1ri_cubes.c\
	m1ri/m1ri_small.c\
	m1ri/m1ri_combine.c\
	m1ri/m7d.c\
	m1ri/m5d.c\
	m1ri/m1ri_io.c \
	m1ri/misc.c

pkgincludesubdir = $(includedir)/m1ri
pkgincludesub_HEADERS = m1ri/m5d.h\
	m1ri/m1ri_3dt.h\
	m1ri/m7d.h\
	m1ri/m1ri_combine.h\
	m1ri/m1ri_small.h\
	m1ri/m1ri_cubes.h\
	m1ri/m1ri_strassen.h\
	m1ri/m1riwrappers.h\
	m1ri/m1riarith.h\
	m1ri/m1ri_classical.h\
	m1ri/m1ri_io.h \
	m1ri/misc.h

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = m1ri.pc

libm1ri_la_LDFLAGS = -release 0.0.$(RELEASE) -no-undefined

check_PROGRAMS =  hadamard m3d_create_test m3d_functions m5d_tests m7d_tests 
hadamard_SOURCES=tests/hadamard.c 
hadamard_LDFLAGS=-lm1ri -lm
hadamard__CFLAGS=$(AM_CFLAGS)

m3d_create_test_SOURCES=tests/m3d_create_test.c 
m3d_create_test_LDFLAGS=-lm1ri -lm
m3d_create_test__CFLAGS=$(AM_CFLAGS)

m3d_functions_SOURCES=tests/m3d_functions.c 
m3d_functions_LDFLAGS=-lm1ri -lm
m3d_functions__CFLAGS=$(AM_CFLAGS)

m5d_tests_SOURCES=tests/m5d_tests.c 
m5d_tests_LDFLAGS=-lm1ri -lm
m5d_tests__CFLAGS=$(AM_CFLAGS)

m7d_tests_SOURCES=tests/m7d_tests.c 
m7d_tests_LDFLAGS=-lm1ri -lm
m7d_tests__CFLAGS=$(AM_CFLAGS)

TESTS =  hadamard m3d_create_test m3d_functions m5d_tests m7d_tests 

clean-local:
