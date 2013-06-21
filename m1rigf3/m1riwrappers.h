/*
Function wrapper
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
#define RU64(a) ((a/64) + ((1) & (a%64)))//division by 64 rounded up



#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>


typedef  u_int64_t vec ;


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


static inline void * m1ri_realloc(void * val, size_t size)
{
    void  * reallocate =  realloc(val, size);
    return reallocate;
}

static inline void m1ri_free(void * val)
{
    free(val);
    
}


static inline int  m1ri_rand()
{
     
    int randomword = rand();
    return randomword;
}


#endif
