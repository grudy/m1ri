#ifndef M1RIPROJECT_M1RIWRAPPERS_H
#define M1RIPROJECT_M1RIWRAPPERS_H
/** 
 
Function wrappers
TOMAS J. BOOTHBY AND ROBERT W. BRADSHAW "BITSLICING AND THE METHOD OF FOUR
RUSSIANS OVER LARGER FINITE FIELDS"

Copyright 2013 William Andrew Alumbaugh <williamandrewalumbaugh@gmail.com>
 
 This program is free software; you can redistribute it and/or
 modify it under the terms of the GNU General Public License
 as published by the Free Software Foundation; either version 2
 of the License, or (at your option) any later version.
 
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 
 m1riwrappers.h
 */

#ifdef HAVE_CONFIG_H

#include <m1ri/config.h>
#include <m1ri/config.h>
#endif

#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <stdbool.h>
#include <stdio.h>
#include <mm_malloc.h>
#include <assert.h>

#include <m1ri/misc.h>

#define M1RI_FN(a, b, c, d) ((a)^(b))&((c)^(d)) /* for finding R[0]# (the first half of the value representingthe sum of vectory and vectorx, vectorr) */
#define M1RI_ST(a, b , c) ((a)^(b)^(c)) /* performing the (S= x[0] XOR y[1] XOR [x1]) and (T = x[1] XOR Y[0] XOR Y[1]) operations of addition */
#define M1RI_DN(a, n) ((a)/(n)) + ((1) && (a%n))/* division by n rounded up */
#define M1RI_MAX(a,b)	((a > b)?	a: b)
#define M1RI_MIN(a, b)	((a > b)? b: a)

typedef u_int64_t vec;

#define M1RI_RADIX 64
/* static const int M1RI_RADIX = 64;	*/
/* ======= */



static	const u_int64_t leftbit	= (1ULL)<<(M1RI_RADIX-1);
static	const u_int64_t rightbit = 1;



static u_int64_t const ibits = 0x8040201008040201;
typedef int rci_t;
typedef int wi_t;
/* typedef unsigned int vbit; */


/*
	Computes the size for padding
*/
static inline u_int32_t powerof2(u_int32_t v)
{
	v--;
	v |= v >> 1;
	v |= v >> 2;
	v |= v >> 4;
	v |= v >> 8;
	v |= v >> 16;
	v++; 
	return v; 
}

static inline void * m1ri_malloc(size_t size) {
	void * allocate = malloc(size);
	if(allocate == NULL)
	{
		m1ri_die("Out of memory, exiting\n");
	}	
	return allocate;
}

/*
	Wrapper for calloc
*/
static inline void * m1ri_calloc(size_t nobj, size_t size) 
{
	void * allocate = calloc(nobj, size);
	if(allocate == NULL)
	{
		m1ri_die("Out of memory, exiting\n");
	}	
	return allocate;
}

/** 
	Wrapper for realloc	
 */
static inline void * m1ri_realloc(void * val, size_t size) 
{
	
	void * reallocate =	realloc(val, size);
	if(reallocate == NULL)
	{
		m1ri_die("Out of memory, exiting\n");
	}
		
	return reallocate;
}

/** 
 \Releases a value into the wilderness
 \param val pointer to be freed
*/

static inline void m1ri_free(void * val) {
free(val);	
}

/** 
 \brief For testing if windowed
 */
static u_int8_t const iswindowed = 0x1; 

static u_int8_t const notwindowed = 0x2;


/** 
 \brief Wrapper for rand
 \
 \Made for working on 64 bit variables "vecs" 
 */
static inline u_int64_t	m1ri_rand() 
{

	assert(RAND_MAX >= ((1ULL<<25)-1));
	u_int64_t randomword = random();
	randomword ^= (u_int64_t)random() << 25;
	randomword ^= (u_int64_t)random() << 50;
	return randomword;
}

static inline void m1ri_sort( const void *ptr, size_t count, size_t size, int (*comp)(const void *, const void *))
{
	qsort(	&ptr,	count, size,	comp); 
}

#endif
