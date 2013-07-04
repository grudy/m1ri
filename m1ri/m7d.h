//
//  m7d.h
//  m1riproject
//  Created by grudy on 6/21/13.
//  Copyright (c) 2013 William Alumbaugh. All rights reserved.
//

#ifndef M1RIPROJECT_M7D_H
#define M1RIPROJECT_M7D_H

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
 
 m7d.h
 */


#include <stdlib.h>
#include "m1riwrappers.h"







/********************************************
 Creates  a union of 192 bits
 ********************************************/

typedef struct {
    
    vec units;
    
    vec middle;
    
    vec sign;
    
} vtri;




/*
 GF(7) Matrix structure
 
 */

typedef struct {
    
    rci_t nrows; //< number of rows
    
    rci_t ncols; //< number of columns
    
    wi_t width; //< the number of vtri's needed to hold columns
    
    
    vtri * block;  //< block containing the data contiguous in memory
    
    vtri ** rows;  // < pointers to rows of the matrix
    
    
    u_int8_t flags;
    
    
    
    
    

    
} m7d_t;


/*
 Matrix Windows
 ______________
 

    | [A0 | A1]
A = | --------
    | |A2 | A3]
 
 
 
 */
typedef  struct{
    m7d_t  a0;
    m7d_t  a1;
    m7d_t  a2;
    m7d_t a3;
    
    
}m7d_windows;


/*
 Read n bits from a sign portion of an element
 x = rows
 y = columns
 M = Matrix read from
 */

vec m7d_rm_bits(m7d_t *M, rci_t  x, rci_t  y, int  n) {
    
    
    
    
    
    wi_t  block = (y  ) / 64;
    
    int  spill = (y  % 64) + n - 64;
    
    vec bits;
    
    bits = (spill <= 0) ? M->rows[x][block].middle << -spill : (M->rows[x][block + 1].sign << (64 - spill)) | (M->rows[x][block].middle >> spill);
    
    
    return bits;
    
    
}


vec m7d_rs_bits(m7d_t *M, rci_t  x, rci_t  y, int  n) {
    
    
    
    
    
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

vec m7d_ru_bits(m7d_t *M, rci_t  x, rci_t  y, int  n) {
    
    
    
    
    
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

vtri m7d_read_elems(m7d_t *M, rci_t  x, rci_t  y, int  n) {
    
    
    
    
    wi_t  block = (y  ) / 64;
    
    int  spill = (y  % 64) + n - 64;
    
    vtri elem;
    
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
void * m7d_rowswap (m7d_t * M, rci_t row_a, rci_t  row_b)
{
    
    
    if((M->nrows >= (row_a ) && (M->nrows >= row_b)))
    {
        vtri * temp =  m1ri_malloc(M->width * sizeof(vtri));
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
void *  m7d_write_elem( m7d_t * M,rci_t x, rci_t y, vec s, vec u )
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



vtri  * m7d_block_allocate(vtri * block, rci_t  nrows,  wi_t  width)
{
    
    
    block  = m1ri_malloc(nrows * width * sizeof(vtri) );
    
    
    
    return block;
    
    
    
}

/*
 
 */




vtri ** m7d_row_alloc(vtri * block, vtri ** rows, wi_t width, rci_t nrows)
{
    
    
    
    
    rows = m1ri_malloc( nrows * width * sizeof(vtri *));
    
    
    for (int i = 0; i <  nrows;  i++ )
    {
        rows[i]  = (block + (i * width));
        
        
    };
    
    return rows;
}

/*
 
 */

m7d_t m7d_create( m7d_t * a, rci_t nrows, rci_t ncols)
{
    
    
    a->ncols = ncols;
    a->nrows = nrows;
    a->width =  RU64(ncols);
    a->block = m7d_block_allocate(a->block,  a->nrows,    a->width);
    a->rows  = m7d_row_alloc(a->block, a->rows, a->width, a->nrows);
    a->flags = iswindowed;
    
    return *a;
    
}

/*
 
 */

vtri * m7d_rand(m7d_t * a)
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


m7d_t  m7d_identity_set(m7d_t * a)

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


m7d_t   m7d_identity(m7d_t  *a, rci_t n)
{
    *a = m7d_create(a, n, n);
    *a = m7d_identity_set(a);
    
    return *a;
    
    
}


/*
 
 */





m7d_windows m7d_windows_create(m7d_t *c)
{
    m7d_windows b;
    int demi= DN((c->width * c->nrows * c->ncols), 4 );
    
    b.a0.block =  m1ri_malloc(demi);
    b.a1.block =  m1ri_malloc(demi);
    b.a2.block =  m1ri_malloc(demi);
    b.a3.block =  m1ri_malloc(demi);
    {
        
        b.a0.block = c->block;
        b.a1.block  = c->block + demi;
        b.a2.block = c->block + (2 * demi);
        b.a3.block = c->block + ( 3 * demi);
        
        
        
        
        
        
        
        
    }
    
    
    return b;
    
}






/*
 
 Releases a m7d_t into the wilderness.
 */



void m7d_free( m7d_t *  tofree)
{
    
    
    m1ri_free(tofree->rows);
    m1ri_free(tofree->block);
    
}


void addgf7(vtri * r, vtri * x, vtri * y)

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



vtri addgf7r(vtri  *x, vtri *y)
{
    vtri  r;
    r.sign = x->sign ^ y->sign;
    r.middle = ~r.sign^ x->middle ^ y->middle;
    r.units = ~r.middle^ x->units ^ y->units;
    r.sign = ~r.units ^ x->sign ^ y->sign;
    r.middle = ~r.sign^ x->middle ^ y->middle;
    r.units = ~r.middle^ x->units ^ y->units;
    
    return r; 
}

void subgf7( vtri *r, vtri *x, vtri *y)               //multiply matrix x by by matrix y.   The product is matrix r.

{

    
    
}



vtri subgf7r(vtri x, vtri y)               //multiply matrix x by by matrix y.   The product is matrix r.

{
    vtri r;
    return r;
}




/********************************************
 matrix r = (direct sum matrix r + matrix x)
 ********************************************/
void iaddgf7(vtri *r,vtri *x)
{
    
    
    
    
    
    
    
}

void isubgf7(vtri *r,vtri *x)  //vector  r = (vector r - vector x)
{
   
    
    
}



void  m7d_mul( vtri *r, vtri *x, vtri *y)             //multiply vector x by y assinging the output to r
{
    
}




//return the value of the vector multiplied


vtri m7d_mul_i(vtri x, vtri y)
{
    
    vtri r;
    r.units = x.units & y.units;
    r.sign  = (y.sign ^ x.sign) & (r.units);
    
    return r;
    
}




m7d_t m7d_transpose(m7d_t * a)
{
    m7d_t b = m7d_create(a, a->ncols , a->nrows);
    int i, j;
    
    for(i = 0; i < a->nrows; i++)
        
        
        for(j = 0; j < a->ncols; j++)
        {
            a->rows[i][j] = b.rows[j][i];
            
            
            
        }
    
    
    return b;
    
    
}


/* * * * * * * * * * * * * * * * * * * *
 Subtract a 1 kilobyte Matrix from another
 1 kilobyte Matrix
 * * * * * * * * * * * * * * * * * * * * */



void sub_64gf7(vtri *R, vtri *A, vtri *B)
{
    int i;
    for (i = 0; i < (sizeof(vec)); i++ )
    {
        R[i] = subgf7r(A[i], B[i]);
    }
}


void add_64gf7(vtri *R, vtri *A, vtri *B)
{
    int i;
    for (i = 0; i < (sizeof(vec)); i++ )
    {
        R[i] = addgf7r(&A[i], &B[i]);
    }
    
}







#endif
