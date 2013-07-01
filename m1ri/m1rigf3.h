 
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
 
 m1rigf3.h
 */

#ifndef M1RIGF3_M1RIGF3_H
#define M1RIGF3_M1RIGF3_H

#include <stdlib.h>
#include "m1riwrappers.h"







/********************************************
Creates  a union of 128 bits
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
    
    wi_t width; //< the number of m1ri_words needed to hold columns
    
    
    vbg * block;  //< block containing the data contiguous in memory
    
    vbg ** rows;  // < pointers to rows of the matrix
    
    
} m3d_t;
/*
 Read n bits from a s portion of an element
 x = rows
 y = columns
 M = Matrix read from
 */

vec m3d_rs_bits(m3d_t *M, rci_t  x, rci_t  y, int  n) {
    
    
    
    
    
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

vec m3d_ru_bits(m3d_t *M, rci_t  x, rci_t  y, int  n) {
    
    
    
    
    
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

vbg m3d_read_elems(m3d_t *M, rci_t  x, rci_t  y, int  n) {
   
    
    
  
    wi_t  block = (y  ) / 64;
    
    int  spill = (y  % 64) + n - 64;
    
    vbg elem;
    
    elem.units = (spill <= 0) ? M->rows[x][block].units << -spill : ((M->rows[x][(block + 1)].units<< (64 - spill)) | (M->rows[x][block].units >> spill));
    
    elem.sign = (spill <= 0) ?  (M->rows[x][block].sign << -spill) : (M->rows[x][block + 1].sign << (64 - spill)) | (M->rows[x][block].sign>> spill);
    
    elem.units = (elem.units >> (64 - n));
   
    elem.sign = (elem.units >> (64 - n));
    
    return elem;
    
    
}


vbg *  m3d_write_elem( m3d_t * M,rci_t x, rci_t y, vec s, vec u )
{
    
    
    
    wi_t  block = (y  ) / 64;
    
    int   spill = (y  % 64) - 63;
    
    vbg * elem = 0;
    
     s = ~(s == 0);
     u = ~(u == 0);
    
    
    M->rows[x][block].units  = (u == 0) ? (~(u << -spill) &  (M->rows[x][block].units))  : ((u << (64 - spill)) | (M->rows[x][block].units));
    
    M->rows[x][block].sign  = (s == 0) ? (~(s << -spill) &  (M->rows[x][block].units))  : ((u << (64 - spill)) | (M->rows[x][block].units));
    
    return elem;
    
    
}



vbg  * m3d_block_allocate(vbg * block, rci_t  nrows,  wi_t  width)
{
    
    
    block  = m1ri_malloc(nrows * width * sizeof(vbg) );
    

    
    return block;
    
    
    
}







vbg ** m3d_row_alloc(vbg * block, vbg ** rows, wi_t width, rci_t nrows)
{
    
    
    
   
    rows = m1ri_malloc( nrows * width * sizeof(vbg *));

    
    for (int i = 0; i <  nrows;  i++ )
    {
                rows[i]  = block + (i * width);
        
    
    };
        
    return rows;
}

m3d_t * m3d_create( m3d_t * a, rci_t nrows, rci_t ncols)
{
    
    
    a->ncols = ncols;
    a->nrows = nrows;
    a->width =  RU64(ncols);
    a->block  = m3d_block_allocate(a->block,  a->nrows,    a->width);
    a->rows = m3d_row_alloc(a->block, a->rows, a->width, a->nrows);
    
    
    return a;
    
}



vbg * m3d_rand(m3d_t * a)
{
   
    for(int i = 0; i < (a->nrows * a->width); i++)
    {
        
        a->block[i].sign = m1ri_rand();
        
        
        
        
        a->block[i].units = m1ri_rand();
    
        
        
        
    }
    return a->block;
}


/*
 
    Releases a m3d_t into the wilderness.
*/
void m3d_free( m3d_t *  tofree)
{
    
    
    m1ri_free(tofree->rows);
    m1ri_free(tofree->block);
   
}








#endif
