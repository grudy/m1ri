/** 
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
#include <m1ri/m1riwrappers.h>

#ifndef _VFDM5_
#define _VFDM5_
/** 
The Representation for GF5
(Boothby's Additive Method) 
000 0
001 1
011 2
101 3
111 4

*/

/** *******************************************
 Creates  a structure for GF(5) Matrices
 ********************************************/

typedef struct {
    
    vec units;
    
    vec middle;
    
    vec sign;
    
} vfd;



/** 

This is unused because of the current representation 
static inline  vfd fold5( vec s3, vec s2, vec s1, vec s0)
{
    vfd r;
    vec t = s2 | s1;
    r.units = s0 ^ t;
    r.middle = (r.units & s0) ^ (s3 ^ s1);
    r.sign = (t ^ s2 ) | (r.middle & s3  );
    return r;
}
*/

/** 
 GF(5) Matrix structure
 
 */

typedef struct m5d_t{
    
    rci_t nrows; //< number of rows
    
    rci_t ncols; //< number of columns
    
    wi_t width; //< the number vfd's needed to hold columns

    vfd * block;  //< block containing the data contiguous in memory
    
    vfd ** rows;  // < pointers to rows of the matrix
    u_int32_t  fblock; //  first block pointed to in a window
    u_int32_t fcol;  //column offset of first block
 
    u_int8_t flags;
	vec  svfd;   //Identifies first vfd used in row
 
    u_int32_t  lblock; //  first block pointed to in a window
    
    
    
} m5d_t;

 
 
typedef struct
{
	/** 
	
	 */
    m5d_t * block;
    /** 
    (slicesize ^ 2) * 64
    */
    m5d_t ** row;
    // (slicesize ^ 2) * 64
    wi_t slicesize;
    //width in slices horizaontally per row
    wi_t width;   
    rci_t nrows;
    rci_t ncols;
    
}m5_slice;




#endif


void m5d_copypadding(m5d_t  * , m5d_t const * );
void m5d_putpadding(m5d_t  * , m5d_t const * );

/** 
 

/** 
 Matrix Windows
 ______________
 
 
 |   [A0 | A1]
 A = | --------
 |   |A2 | A3]
 
 
 
 */


vec m5d_rm_bits(m5d_t *, rci_t  , rci_t  , int  ) ;
vec m5d_rs_bits(m5d_t *, rci_t  , rci_t  , int  );
vec m5d_ru_bits(m5d_t *, rci_t  , rci_t  , int  );
vfd m5d_read_elems(m5d_t *M, rci_t  x, rci_t  y, int  );
void * m5d_rowswap (m5d_t * , rci_t , rci_t  );
/** 
 
*/
//unfinished
void *  m5d_write_elem( m5d_t * ,rci_t , rci_t , vec , vec , vec);


/** 
 
 */



vfd  * m5d_block_allocate(vfd * block, rci_t  nrows,  wi_t  width);
/** 
 
 */




vfd ** m5d_row_alloc(vfd * block, vfd ** rows, wi_t width, rci_t nrows);
/** 
 
 */

m5d_t m5d_create( m5d_t * , rci_t nrows, rci_t ncols);
/** 
 
 */

vfd * m5d_rand(m5d_t * );
/** 
 Make an Identity Matrix
 a = Identity matrix
 n = matrix size (row length and column width)
*/


m5d_t  m5d_identity_set(m5d_t * );

/** 
		 
*/
m5d_t   m5d_identity(m5d_t  *, rci_t );

/** 
 
 */


/** 
 
 Releases a m5d_t into the wilderness.
 */

void m5d_free( m5d_t *  );



void vfd_sub( vfd *, vfd *, vfd *);               //subtract vector x by by vector y.   The product is vector r.

/** *******************************************
 matrix r = (direct sum matrix r + matrix x)
 ********************************************/
void iadd_vfd(vfd *,vfd *);

void m5d_add_r(m5d_t *, m5d_t *, m5d_t *);
void m5d_sub(m5d_t *, m5d_t *, m5d_t *);
void m5d_sub_d(m5d_t  * a , m5d_t * b);
void add_vfd(vfd *, vfd * , vfd *);
void sub_vfd(vfd *, vfd * , vfd *);
void m5d_add2(vfd * , vfd * , vfd * );

void m5d_add2_i(vfd * , vfd * );
int m5d_equal(m5d_t const *, m5d_t const *);
void m5d_add_64(vfd **, vfd **  , vfd ** );
m5d_t   m5d_window(m5d_t *, rci_t , rci_t , rci_t , rci_t );
void   m5d_window_create(m5d_t *, m5d_t * , rci_t , rci_t , rci_t , rci_t );

void m5d_sub_64(m5d_t * c ,m5d_t  * a , m5d_t * b);


        
void m5d_sub_i(vfd *,vfd *);
vfd m5d_mul2(vfd);
vfd m5d_mul3(vfd);
vfd m5d_mul4(vfd);




/* These Functions are expanded upon Toms 
 implementation of the method of four RUSSIANS
 */
void *  m5d_combine3(vfd *, vfd * );
void m5d_combine4(vfd *, vfd * );
void m5d_combine5(vfd *, vfd * );
void m5d_combine6(vfd *, vfd * );
void m5d_combine7(vfd *, vfd * );
void m5d_combine8(vfd *, vfd *);

/* ****************************************************************************
								GF(5) base cases for multiplication
**************************************************************************** */

void m5d_mul_64(vfd **, vfd **, vfd **);
void m5d_mul_32(vfd *, vfd *, vfd *);
void m5d_mul_16(vfd *, vfd *, vfd *);
void m5d_mul_8(vfd *, vfd *, vfd *);
void m5d_mul_4(vfd *R, vfd *A, vfd *B);
vfd * m5d_transpose_vfd(vfd  **, vfd **);




/*	These functions work with lots of partitions*/
 
/*Pointers to submatrices*/


void  m5d_slices(m5_slice *  , m5d_t * , wi_t );
void  m5d_quarter(m5_slice *  , m5d_t * );
m5d_t m5d_transpose_sliced(m5d_t * );
m5d_t  * m5_blockslice_allocate(m5d_t * , rci_t  ,  wi_t  );
m5d_t ** m5_rowslice_allocate(m5d_t * , m5d_t ** , wi_t , rci_t );
m5d_t * m5d_hadamard(m5d_t const * , m5d_t const *  );
void m5d_copy(m5d_t *, m5d_t const *);
void * m5d_colswap_capped_row(m5d_t *, rci_t , rci_t, rci_t );

#endif
