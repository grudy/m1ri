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
#include <m1ri/m1riwrappers.h>







/*
 Matrix Windows
 ______________
 
 
 |   [A0 | A1]
 A = | --------
 |   |A2 | A3]
 
 
 
 */


vec m5d_rm_bits(m5d_t *M, rci_t  x, rci_t  y, int  n) ;
vec m5d_rs_bits(m5d_t *M, rci_t  x, rci_t  y, int  n);
vec m5d_ru_bits(m5d_t *M, rci_t  x, rci_t  y, int  n);
vfd m5d_read_elems(m5d_t *M, rci_t  x, rci_t  y, int  n);
void * m5d_rowswap (m5d_t * M, rci_t row_a, rci_t  row_b);
/*
 
*/
//unfinished
void *  m5d_write_elem( m5d_t * ,rci_t , rci_t , vec , vec , vec);


/*
 
 */



vfd  * m5d_block_allocate(vfd * block, rci_t  nrows,  wi_t  width);
/*
 
 */




vfd ** m5d_row_alloc(vfd * block, vfd ** rows, wi_t width, rci_t nrows);
/*
 
 */

m5d_t m5d_create( m5d_t * a, rci_t nrows, rci_t ncols);
/*
 
 */

vfd * m5d_rand(m5d_t * a);
/*
 Make an Identity Matrix
 a = Identity matrix
 n = matrix size (row length and column width)
 
 
 */


m5d_t  m5d_identity_set(m5d_t * a);
/*
 
 */


m5d_t   m5d_identity(m5d_t  *a, rci_t n);

/*
 
 */


/*
 
 Releases a m5d_t into the wilderness.
 */



void m5d_free( m5d_t *  );

void addgf5(vfd *, vfd * , vfd *);



void m5d_add_r(m5d_t  *, m5d_t *, m5d_t *);


void m5d_sub( vfd *, vfd *, vfd *);               //subtract vector x by by vector y.   The product is vector r.





vfd m5d_sub_r(vfd , vfd );          

vfd m5d_mul2(vfd);


vfd m5d_mul3(vfd);


vfd m5d_mul4(vfd);



/********************************************
 matrix r = (direct sum matrix r + matrix x)
 ********************************************/
void iaddgf5(vfd *,vfd *);

void isubgf5(vfd *,vfd *);




int m5d_equal(m5d_t const *, m5d_t const *);



m5d_t   m5d_window(m5d_t *, rci_t , rci_t , rci_t , rci_t );
void   m5d_window_create(m5d_t *, m5d_t * , rci_t , rci_t , rci_t , rci_t );

void sub_64gf5(vfd *, vfd *, vfd *);

void add_m1riw_gff5(vfd *, vfd *, vfd *);











#endif
