/*  m1riproject
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
#ifndef M1RIPROJECT_M7D_H
#define M1RIPROJECT_M7D_H
#include <stdlib.h>
#include <m1ri/m1riwrappers.h>







/********************************************
 Creates  a structure of 192 bits
 ********************************************/

typedef struct {
    
    vec units;
    vec middle;
    vec sign;
    
    
    
} vtri;




/*
 GF(7) Matrix structure
 GF7

000 = 0    
100  = 1
010 = 2
110 = 3
001 = 4
101 = 5
011 = 6

 */

typedef struct {
    
    rci_t nrows; //< number of rows
    
    rci_t ncols; //< number of columns
    
    wi_t width; //< the number of vbg's needed to hold columns
    
    vtri * block;  //< block containing the data contiguous in memory
    
    vtri ** rows;  // < pointers to rows of the matrix
    
    vec  svtri;   //Identifies first vbg used in row
    
    // wi_t rowstride;  //vbg's in block to traverse to  get to first
    
    u_int32_t  lblock; //  first block pointed to in a window
    u_int32_t fcol;  //column offset of first block
    u_int8_t flags;    //IsWindowed, NotWindowed

    
    
} m7d_t;


/*
 Matrix Windows
 ______________
  
 | [A0 | A1]
 A = | --------
 | |A2 | A3]

 */

//Read 'middle' bits
vec m7d_rm_bits(m7d_t *M, rci_t  x, rci_t  y, int  n) ;


// Read n 'sign' bits
vec m7d_rs_bits(m7d_t *M, rci_t  x, rci_t  y, int  n);
// Read  n 'unit' bits
vec m7d_ru_bits(m7d_t *M, rci_t  x, rci_t  y, int  n);
// Read   n elements
vtri m7d_read_elems(m7d_t *M, rci_t  x, rci_t  y, int  n);

void  m7d_rowswap (m7d_t * M, rci_t row_a, rci_t  row_b);

/*
 
 */
void m7d_copypadding(m7d_t  * , m7d_t const * );
void m7d_putpadding(m7d_t  * , m7d_t const * );
void add_64_m7d(vtri ** , vtri **, vtri **);
//unfinished
void *  m7d_write_elem( m7d_t * M,rci_t x, rci_t y, vec s,  vec m , vec u );

/*
 
 */


vtri  * m7d_block_allocate(vtri * block, rci_t  nrows,  wi_t  width);
/*
 
 */

vtri ** m7d_row_alloc(vtri * block, vtri ** rows, wi_t width, rci_t nrows);
/*
 
 */

m7d_t m7d_create( m7d_t * a, rci_t nrows, rci_t ncols);
/*
 
 */
void reduce_vtri( vtri * );
m7d_t m7d_rand(m7d_t * a);
/*
 Make an Identity Matrix
 a = Identity matrix
 n = matrix size (row length and column width)
 
 
 */


m7d_t  m7d_identity_set(m7d_t * a);
/*
 
 */


m7d_t   m7d_identity(m7d_t  *a, rci_t n);

/*
 
 */



/*
 Releases a m7d_t into the wilderness.
 */



int m7d_equal(m7d_t const *, m7d_t const *);
void m7d_free( m7d_t *  );
void vtri_negate( vtri * );


/********************************************
 matrix r = (direct sum matrix r + matrix x)
 ********************************************/
void add_vtri(vtri *, vtri * , vtri *);

void iadd_vtri(vtri  *, vtri *);

void isub_m7d(vtri  *, vtri *);

//Scalar  multiplication
vtri m7d_mul_2(vtri);
vtri m7d_mul_3(vtri);
vtri m7d_mul_4(vtri);
vtri m7d_mul_5(vtri);
vtri m7d_mul_6(vtri);


// negate  r0, r1, r2 ← a0, a1, a2


/********************************************
 matrix r = (difference matrix r - matrix x)                //Or x will the function be  r  = x- r???
 ********************************************/

/*
	GF(7) Addition on a single M1RI word.
*/
void m7d_add_r(m7d_t *, m7d_t *, m7d_t *);




m7d_t   m7d_window(m7d_t *, rci_t , rci_t , rci_t , rci_t );
void   m7d_window_create(m7d_t *, m7d_t * , rci_t , rci_t , rci_t , rci_t);



#endif
