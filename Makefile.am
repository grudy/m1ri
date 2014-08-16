AUTOMAKE_OPTIONS = gnu subdir-objects
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

BUILT_SOURCES = m1ri/config.h

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
	m1ri/m1ri_decom.h \
	m1ri/config.h



pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = m1ri.pc

libm1ri_la_LDFLAGS = -release 0.0.$(RELEASE) -no-undefined
libm1ri_la_LIBADD = $(LIBPNG_LIBADD)

check_PROGRAMS =  test_arithmetic test_hadamard m3d_mul_benchmark test_m3d test_m5d test_m7d

test_arithmetic_SOURCES=tests/test_arithmetic.c
test_arithmetic_LDFLAGS=-lm1ri -lm
test_arithmetic_CFLAGS=$(AM_CFLAGS)

test_hadamard_SOURCES=tests/test_hadamard.c
test_hadamard_LDFLAGS=-lm1ri -lm
test_hadamard_CFLAGS=$(AM_CFLAGS)

test_m3d_SOURCES=tests/test_m3d.c
test_m3d_LDFLAGS=-lm1ri -lm
test_m3d_CFLAGS=$(AM_CFLAGS)

test_m5d_SOURCES=tests/test_m5d.c
test_m5d_LDFLAGS=-lm1ri -lm
test_m5d_CFLAGS=$(AM_CFLAGS)

test_m7d_SOURCES=tests/test_m7d.c
test_m7d_LDFLAGS=-lm1ri -lm
test_m7d_CFLAGS=$(AM_CFLAGS)

TESTS =  test_arithmetic test_hadamard test_m3d test_m5d test_m7d 


bin_PROGRAMS = m3d_mul_benchmark

m3d_mul_benchmark_SOURCES=bench/m3d_mul_benchmark.c
m3d_mul_benchmark_CFLAGS=$(AM_CFLAGS)
m3d_mul_benchmark_LDADD=-lm1ri -lm
m3d_mul_benchmark_LDFLAGS=-no-install


# benchmarketing

SUBDIRS = .

clean-local:
	(cd tests; make clean; cd ..)
