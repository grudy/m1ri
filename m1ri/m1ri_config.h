#ifndef M1RI_M1RI_CONFIG_H
#define M1RI_M1RI_CONFIG_H

// Defines determined during configuration of M1RI.
#define __M1RI_HAVE_OPENMP		@M1RI_HAVE_OPENMP@

#define __M1RI_CC                       "gcc -std=gnu99"
#define __M1RI_CFLAGS                   "@SIMD_CFLAGS@ -fopenmp -g -O2"
#define __M1RI_SIMD_CFLAGS              "@SIMD_CFLAGS@"
#define __M1RI_OPENMP_CFLAGS            "-fopenmp"

#endif // M1RI_M1RI_CONFIG_H
