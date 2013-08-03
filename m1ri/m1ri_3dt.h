 
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
 
 m1ri_3dt.h
 */

#ifndef M1RIGF3_M1RIGF3_H
#define M1RIGF3_M1RIGF3_H

#include <stdlib.h>
#include <m1ri/m1riwrappers.h>






/********************************************
Creates  a struct of 128 bits
********************************************/

typedef struct {
    
    vec units;
    vec sign;
} vbg;




/*
    GF(3) Matrix structure
 
*/

typedef struct {
    
    rci_t nrows; //< number of rows
    
    rci_t ncols; //< number of columns
    
    wi_t width; //< the number of vbg's needed to hold columns
    
    vbg * block;  //< block containing the data contiguous in memory
    
    vbg ** rows;  // < pointers to rows of the matrix
    
    vec  svbg;   //Identifies first vbg used in row
    u_int64_t a;
    // wi_t rowstride;  //vbg's in block to traverse to  get to first
    
    u_int32_t  lblock; //  first block pointed to in a window
    u_int32_t fcol;  //column offset of first block
    u_int8_t flags;    //IsWindowed, NotWindowed    
    
    
    
    
} m3d_t;


/*
 Read n bits from a s portion of an element
 x = rows
 y = columns
 M = Matrix read from
 */
m3d_t transpose(m3d_t  const *);
vec m3d_rs_bits(m3d_t const *, rci_t  , rci_t  , int  );

/*
 Read n bits from units
 x = rows
 y = columns
 M = Matrix read from
 */

vec m3d_ru_bits(m3d_t const  *, rci_t  , rci_t  , int  );




/*
 Read n elements
 x = rows
 y = columns
 M = Matrix read from 
*/

vbg m3d_read_elems(m3d_t const *, rci_t  , rci_t  , int  );
m3d_t m3d_transposewin(m3d_t  const * );


/*
Swap rows in a matrix;
*/

void * m3d_rowswap (m3d_t  * , rci_t , rci_t );



/*
Naive column swapping
*/



void * m3d_colswap(m3d_t *, rci_t , rci_t );




 void   m3d_write_elem( m3d_t * ,rci_t , rci_t , vec , vec  );

static inline void  * m3d_block_allocate(vbg * block, rci_t  nrows,  wi_t  width)
{
    block  = m1ri_calloc(nrows * width ,  sizeof(vbg) );
    return block;
}

/*
 
 */


static inline vbg ** m3d_row_alloc(vbg * block, vbg ** rows, wi_t width, rci_t nrows)
{
	int i;
    rows = m1ri_malloc( nrows * width * sizeof(vbg *));
    for ( i = 0; i <  nrows;  i++ )
    {
        rows[i]  = block + i * width;
    };
    return rows;
}


m3d_t m3d_create( m3d_t *  , rci_t nrows, rci_t );



m3d_t m3d_rand(m3d_t * );


/*
 Make an Identity Matrix
 a = Identity matrix 
 n = matrix size (row length and column width)

 
*/
m3d_t    m3d_identity_set(m3d_t * );


m3d_t   m3d_identity(m3d_t  *, rci_t );


/*
 windows in m1ri_word rows * m1ri_word column incriments
 stvbg = the vbg/width offset from the base matrix
 strow = row offset in increments of 64
 sizecol  = cols * 64
 sizerow  = rows * 64
 */
m3d_t  m3d_window(m3d_t  *, rci_t , rci_t , rci_t , rci_t );


/*
 Same as m3d_window but the second argument is made into the window
 */

 void   m3d_window_create(m3d_t *, m3d_t * , rci_t , rci_t , rci_t , rci_t );

/*
 Concat b on the end of a, the result is c
 
 
 
   [a] [b] ----->  [a b]   ===  C

*/
m3d_t m3d_concat(m3d_t * , m3d_t * , m3d_t * );


/*
    Stacks a on b, resulting matrix is c
    [a]
         ===  C
    [b]
 
*/
m3d_t m3d_stack(m3d_t * ,  m3d_t * , m3d_t * );

    
  
/*
 
 Releases a m3d_t into the wilderness.
 */






int m3d_equal(m3d_t const  *, m3d_t const  *);





void m3d_free( m3d_t *  );




#endif
