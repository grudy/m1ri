/*
 Matrix Represenations and basic operations
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
 
 m5d.h
 */

#ifndef M1RIPORJECT_M5D_H
#define M1RIPROJECT_M5D_H

#include <stdlib.h>
#include "m1riwrappers.h"







/********************************************
 Creates  a union of 192 bits
 ********************************************/

typedef struct {
    
    vec units;
    
    vec middle;
    
    vec sign;
    
} vfd;




/*
 GF(5) Matrix structure
 
 */

typedef struct {
    
    rci_t nrows; //< number of rows
    
    rci_t ncols; //< number of columns
    
    wi_t width; //< the number vfd's needed to hold columns
    
    
    vfd * block;  //< block containing the data contiguous in memory
    
    vfd ** rows;  // < pointers to rows of the matrix
    
    
    u_int8_t flags;
    
    
    
    
    
    
    
} m5d_t;


/*
 Matrix Windows
 ______________
 
 
 |   [A0 | A1]
 A = | --------
 |   |A2 | A3]
 
 
 
 */
typedef  struct{
    m5d_t  a0;
    m5d_t  a1;
    m5d_t  a2;
    m5d_t a3;
    
    
}m5d_windows;


/*
 Read n bits from a sign portion of an element
 x = rows
 y = columns
 M = Matrix read from
 */

vec m5d_rm_bits(m5d_t *M, rci_t  x, rci_t  y, int  n) {
    
    
    
    
    
    wi_t  block = (y  ) / 64;
    
    int  spill = (y  % 64) + n - 64;
    
    vec bits;
    
    bits = (spill <= 0) ? M->rows[x][block].middle << -spill : (M->rows[x][block + 1].sign << (64 - spill)) | (M->rows[x][block].middle >> spill);
    
    
    return bits;
    
    
}


vec m5d_rs_bits(m5d_t *M, rci_t  x, rci_t  y, int  n) {
    
    
    
    
    
    wi_t  block = (y  ) / 64;
    
    int  spill = (y  % 64) + n - 64;
    
    vec bits;
    
    bits = (spill <= 0) ? M->rows[x][block].sign << -spill : (M->rows[x][block + 1].sign << (64 - spill)) | (M->rows[x][block].sign >> spill);
    
    
    return bits;
    
    
}
/*
 Read n bits from units
 x = rows
 y = columns
 M = Matrix read from
 */

vec m5d_ru_bits(m5d_t *M, rci_t  x, rci_t  y, int  n) {
    
    
    
    
    
    wi_t  block = (y  ) / 64;
    
    int  spill = (y  % 64) + n - 64;
    
    vec bits;
    
    bits = (spill <= 0) ? M->rows[x][block].units << -spill : (M->rows[x][block + 1].units<< (64 - spill)) | (M->rows[x][block].units>> spill);
    
    
    
    
    
    return bits;
    
    
}





/*
 Read n elements
 x = rows
 y = columns
 M = Matrix read from
 */

vfd m5d_read_elems(m5d_t *M, rci_t  x, rci_t  y, int  n) {
    
    
    
    
    wi_t  block = (y  ) / 64;
    
    int  spill = (y  % 64) + n - 64;
    
    vfd elem;
    
    elem.units = (spill <= 0) ? M->rows[x][block].units << -spill : ((M->rows[x][(block + 1)].units<< (64 - spill)) | (M->rows[x][block].units >> spill));
    
    elem.sign = (spill <= 0) ?  (M->rows[x][block].sign << -spill) : (M->rows[x][block + 1].sign << (64 - spill)) | (M->rows[x][block].sign>> spill);
    
    elem.middle = (spill <= 0) ?  (M->rows[x][block].middle << -spill) : (M->rows[x][block + 1].middle << (64 - spill)) | (M->rows[x][block].middle>> spill);
    
    elem.middle = (elem.middle >> (64 - n));
    
    elem.units = (elem.units >> (64 - n));
    
    elem.sign = (elem.sign >> (64 - n));
    
    
    
    return elem;
    
    
}





/*
 Swap rows in a matrix;
 */
void * m5d_rowswap (m5d_t * M, rci_t row_a, rci_t  row_b)
{
    
    
    if((M->nrows >= (row_a ) && (M->nrows >= row_b)))
    {
        vfd * temp =  m1ri_malloc(M->width * sizeof(vfd));
        temp =  M->rows[row_a -1];
        M->rows[row_a -1] = M->rows[row_b -1];
        M->rows[row_b -1] =  temp;
        
        
        
    }
    
    
    {
        
        
    }
    return 0;
}


/*
 
 */


//unfinished
void *  m5d_write_elem( m5d_t * M,rci_t x, rci_t y, vec s, vec u )
{
    
    
    
    wi_t  block = (y  ) / 64;
    
    int   spill = (y  % 64) - 63;
    
    
    
    s = ~(s == 0);
    u = ~(u == 0);
    
    
    M->rows[x][block].units  = (u == 0) ? (~(rightbit << -spill) &  (M->rows[x][block].units))  : ((u << (64 - spill)) | (M->rows[x][block].units));
    
    M->rows[x][block].sign  = (s == 0) ? (~(rightbit << -spill) &  (M->rows[x][block].units))  : ((u << (64 - spill)) | (M->rows[x][block].units));
    M->rows[x][block].middle  = (s == 0) ? (~(rightbit << -spill) &  (M->rows[x][block].units))  : ((u << (64 - spill)) | (M->rows[x][block].units));
    
    return 0;
    
    
}


/*
 
 */



vfd  * m5d_block_allocate(vfd * block, rci_t  nrows,  wi_t  width)
{
    
    
    block  = m1ri_malloc(nrows * width * sizeof(vfd) );
    
    
    
    return block;
    
    
    
}

/*
 
 */




vfd ** m5d_row_alloc(vfd * block, vfd ** rows, wi_t width, rci_t nrows)
{
    
    
    
    
    rows = m1ri_malloc( nrows * width * sizeof(vfd *));
    
    
    for (int i = 0; i <  nrows;  i++ )
    {
        rows[i]  = (block + (i * width));
        
        
    };
    
    return rows;
}

/*
 
 */

void * m5d_create( m5d_t * a, rci_t nrows, rci_t ncols)
{
    
    
    a->ncols = ncols;
    a->nrows = nrows;
    a->width =  RU64(ncols);
    a->block = m5d_block_allocate(a->block,  a->nrows,    a->width);
    a->rows  = m5d_row_alloc(a->block, a->rows, a->width, a->nrows);
    a->flags = iswindowed;
    
    return a;
    
}

/*
 
 */

vfd * m5d_rand(m5d_t * a)
{
    
    for(int i = 0; i < (a->nrows * a->width); i++)
    {
        
        a->block[i].sign = m1ri_rand();
        
        
        a->block[i].middle = m1ri_rand();
        
        a->block[i].units = m1ri_rand();
        
        
        
        
    }
    return a->block;
}


/*
 Make an Identity Matrix
 a = Identity matrix
 n = matrix size (row length and column width)
 
 
 */


m5d_t  m5d_identity_set(m5d_t * a)

{
    if (a->nrows == a->ncols)
    {
        
        
        
        for(int i = 0; i < (a->nrows/64); i++ )
        {
            
            a->rows[i][i].units = ibits;
            
        }
        
        
    }
    return *a;
}

/*
 
 */


m5d_t   m5d_identity(m5d_t  *a, rci_t n)
{
    a = m5d_create(a, n, n);
    *a = m5d_identity_set(a);
    
    return *a;
    
    
}


/*
 
 */





m5d_windows m5d_windows_create(m5d_t *a)
{
    m5d_windows b;
    int demi= DN((a->width * a->nrows * a->ncols), 4 );
    
    b.a0.block =  m1ri_malloc(demi);
    b.a1.block =  m1ri_malloc(demi);
    b.a2.block =  m1ri_malloc(demi);
    b.a3.block =  m1ri_malloc(demi);
    {
        
        b.a0.block = a->block;
        b.a1.block  = a->block + demi;
        b.a2.block = a->block + (2 * demi);
        b.a3.block = a->block + ( 3 * demi);
        
        
        
        
        
        
        
        
    }
    b.a0.flags = notwindowed;
    b.a1.flags = notwindowed;
    b.a2.flags = notwindowed;
    b.a3.flags = notwindowed;
    
    return b;
    
}






/*
 
 Releases a m5d_t into the wilderness.
 */



void m5d_free( m5d_t *  tofree)
{
    
    
    m1ri_free(tofree->rows);
    m1ri_free(tofree->block);
    
}


void addgf5(vfd * r, vfd * x, vfd * y)

{
    
    vec s;
    r->sign = x->sign ^ y->sign;
    r->middle = ~r->sign^ x->middle ^ y->middle;
    r->units = ~r->middle^ x->units ^ y->units;
    s = ~r->units;
    r->sign = s ^ x->sign ^ y->sign;
    r->middle = ~r->sign^ x->middle ^ y->middle;
    r->units = ~r->middle^ x->units ^ y->units;
    
    
}



vfd addgf5r(vfd  *x, vfd *y)
{
    vfd  r;
    r.sign = x->sign ^ y->sign;
    r.middle = ~r.sign^ x->middle ^ y->middle;
    r.units = ~r.middle^ x->units ^ y->units;
    r.sign = ~r.units ^ x->sign ^ y->sign;
    r.middle = ~r.sign^ x->middle ^ y->middle;
    r.units = ~r.middle^ x->units ^ y->units;
    
    return r;
}

void subgf5( vfd *r, vfd *x, vfd *y)               //multiply matrix x by by matrix y.   The product is matrix r.

{
    
    
    
}



vfd subgf5r(vfd x, vfd y)               //multiply matrix x by by matrix y.   The product is matrix r.

{
    vfd r;
    return r;
}




/********************************************
 matrix r = (direct sum matrix r + matrix x)
 ********************************************/
void iaddgf5(vfd *r,vfd *x)
{
    
    
    
    
    
    
    
}

void isubgf5(vfd *r,vfd *x)  //matrix r = (matrix r - matrix x)
{
    
    
    
}



void  m5d_mul( vfd *r, vfd *x, vfd *y)             //multiply matrix x by y assinging the output to r
{
    
}




//return the value of the matrix multiplied


vfd m5d_mul_i(vfd x, vfd y)
{
   
    vfd r;
    r.units = x.units & y.units;
    r.sign  = (y.sign ^ x.sign) & (r.units);
    
    return r;
    
}






/* * * * * * * * * * * * * * * * * * * *
 Subtract a 1 Megabyte Matrix from another
 1 megabyte Matrix
 * * * * * * * * * * * * * * * * * * * * */



void sub_64gf5(vfd *R, vfd *A, vfd *B)
{
    int i;
    for (i = 0; i < (sizeof(vec)); i++ )
    {
        R[i] = subgf5r(A[i], B[i]);
    }
}


void add_64gf5(vfd *R, vfd *A, vfd *B)
{
    int i;
    for (i = 0; i < (sizeof(vec)); i++ )
    {
        R[i] = addgf5r(&A[i], &B[i]);
    }
    
}














#endif
