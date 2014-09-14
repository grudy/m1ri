CC = gcc -std=cc99
DEFINES =
TOPSRCDIR = ..
CFLAGS = -I$(TOPSRCDIR) -D_XOPEN_SOURCE=600 -g -O2 $(DEFINES) 
STAGEDIR := $(realpath -s $(TOPSRCDIR)/.libs)
LDFLAGS = -L$(STAGEDIR) -Wl,-rpath,$(STAGEDIR) -lm1ri -lm
PAPI_FLAGS =  
PAPI_CFLAGS = 
DEBUG = -ggdb

TESTS = \
	test_hadamard.c	\
	test_create_m3d.c	\
	test_m5d.c	\
	m5d_tests.c	\
	test_m7d.c
    test_arithmetic.c\
    test_permutations.c\
    test_decom.c

pkgincludesubdir = $(includedir)/m1ri
M1RI_HEADERS = ../src/m1ri/m5d.h\
	../src/m1ri/m3d.h\
	../src/m7d.h\
	../src/m1ri_combine.h\
	 ../src/m1ri_small.h\
	../src/m1ri_cubes.h\
	../src/m1ri_strassen.h\
	.../src/m1riwrappers.h\
	../src//m1ri_classical.h\
	../src/m1ri_io.h\
	../src/m1ri_permutations.h\
	../src/m1ri_decom.h

distclean: clean
