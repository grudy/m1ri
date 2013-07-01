/*
 
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

#ifndef M1RIPROJECT_M1RIWRAPPERS_H
#define M1RIPROJECT_M1RIWRAPPERS_H
#define FN(a, b, c, d) ((a)^(b))&((c)^(d)) //for finding R[0]# (the first half of the value representingthe sum of vectory and vectorx, vectorr)
#define ST(a, b , c) ((a)^(b)^(c)) //performing the (S= x[0] XOR y[1] XOR [x1]) and (T = x[1] XOR Y[0] XOR Y[1]) operations of addition
#define RU64(a) ((a/64) + ((1) && (a%64)))//division by 64 rounded up
#define DN(a, n) ((a/n) + ((1) && (a%n)))//division by n rounded up


#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>

u_int64_t const leftbit = 0x8000000000000000;
u_int64_t const rightbit = 0x1;
typedef  u_int64_t vec ;


u_int64_t const  ibits = 0x8040201008040201;



typedef int rci_t ;

typedef int wi_t ;


typedef unsigned int vbit;

typedef u_int64_t m1ri_word;

static inline void * m1ri_malloc(size_t size)
{
    void * allocate = malloc(size);
    return  allocate;
    
}


static inline void * m1ri_calloc(size_t nobj, size_t size)
{
    void * call = calloc(nobj, size);
    return call;

}




/*
 
 */


static inline void * m1ri_realloc(void * val, size_t size)
{
    void  * reallocate =  realloc(val, size);
    return reallocate;
}


/*
 
 */

static inline void m1ri_free(void * val)
{
    free(val);
    
}



/*
 For testing if windowed
 */

static u_int8_t const iswindowed = 0x1; 
static u_int8_t const notwindowed = 0x2;




static inline int  m1ri_rand()
{
     
    int randomword = rand();
    return randomword;
}


u_int64_t static const bc_l[64] =  { 0x8000000000000000,    0xc000000000000000,  0xe000000000000000,    0xf000000000000000,  0xf800000000000000,
    0xfc00000000000000,  0xfe00000000000000,    0xff00000000000000,  0xff80000000000000,    0xffc0000000000000,  0xffe0000000000000,    0xfff0000000000000,
    0xfff8000000000000,    0xfffc000000000000,  0xfffe000000000000,    0xffff000000000000,  0xffff800000000000,    0xffffc00000000000,  0xffffe00000000000,
    0xfffff00000000000,  0xfffff80000000000,    0xfffffc0000000000,  0xfffffe0000000000,    0xffffff0000000000,  0xffffff8000000000,    0xffffffc000000000,
    0xffffffe000000000,    0xfffffff000000000,  0xfffffff800000000,    0xfffffffc00000000,  0xfffffffe00000000,    0xffffffff00000000,  0xffffffff80000000,
    0xffffffffc0000000,  0xffffffffe0000000,    0xfffffffff0000000,  0xfffffffff8000000,    0xfffffffffc000000,  0xfffffffffe000000,    0xffffffffff000000,
    0xffffffffff800000,    0xffffffffffc00000,  0xffffffffffe00000,    0xfffffffffff00000,  0xfffffffffff80000,    0xfffffffffffc0000,  0xfffffffffffe0000,
    0xffffffffffff0000,  0xffffffffffff8000,    0xffffffffffffc000,  0xffffffffffffe000,    0xfffffffffffff000,  0xfffffffffffff800,    0xfffffffffffffc00,
    0xfffffffffffffe00,    0xffffffffffffff00,  0xffffffffffffff80,    0xffffffffffffffc0,  0xffffffffffffffe0,    0xfffffffffffffff0,  0xfffffffffffffff8,
    0xfffffffffffffffc,  0xfffffffffffffffe,    0xffffffffffffffff};

u_int64_t static const bc_r[64] = {  0x1,    0x3,  0x7,    0xf,  0x1f,    0x3f,  0x7f,    0xff,  0x1ff,    0x3ff,  0x7ff,    0xfff,  0x1fff,
    0x3fff,  0x7fff,    0xffff,  0x1ffff,    0x3ffff,  0x7ffff,    0xfffff,  0x1fffff,    0x3fffff,  0x7fffff,    0xffffff,  0x1ffffff,    0x3ffffff,
    0x7ffffff,    0xfffffff,  0x1fffffff,    0x3fffffff,  0x7fffffff,    0xffffffff,  0x1ffffffff,    0x3ffffffff,  0x7ffffffff,    0xfffffffff,
    0x1fffffffff,    0x3fffffffff,  0x7fffffffff,    0xffffffffff,  0x1ffffffffff,    0x3ffffffffff,  0x7ffffffffff,    0xfffffffffff,  0x1fffffffffff,
    0x3fffffffffff,  0x7fffffffffff,    0xffffffffffff,  0x1ffffffffffff,    0x3ffffffffffff,  0x7ffffffffffff,    0xfffffffffffff,  0x1fffffffffffff,
    0x3fffffffffffff,  0x7fffffffffffff,    0xffffffffffffff,  0x1ffffffffffffff,    0x3ffffffffffffff,  0x7ffffffffffffff,    0xfffffffffffffff,
    0x1fffffffffffffff,    0x3fffffffffffffff,  0x7fffffffffffffff,    0xffffffffffffffff};


static inline void m1ri_swap_vec(vec *a, vec *b)
{
    vec temp;
    temp = *a;
    *a = *b;
    *b = temp;

}



static inline void m1ri_sort( const void *ptr, size_t count, size_t size, int (*comp)(const void *, const void *))
               {
                    qsort(  &ptr,  count, size,  comp);
                         
               
               }
               
               
               
               
#endif
