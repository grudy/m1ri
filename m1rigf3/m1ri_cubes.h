/*  Cube Form 
 //  
 //  m1riproject
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
 
 //  Copyright (c) 2013 William Alumbaugh. All rights reserved.
 
 m1ri_cubes.h
 */


#ifndef m1riproject_m1ri_cubes_h
#define m1riproject_m1ri_cubes_h
#include "m1riwrappers.h"
#include "m1ri_3dt.h"


typedef struct {
    
    rci_t nrows; //< number of rows
    
    rci_t ncols; //< number of columns
    
    wi_t width; //< the number of vbg's needed to hold columns
    
    
    
    vbg * block;  //< block containing the data contiguous in memory
    
    vbg ** rows;  // < pointers to rows of the matrix
    

    
    u_int8_t flags;
    
    
    
    
    
} m3d_cbe;

vbg  * m3d_cube_block(vbg * block, rci_t  nrows,  wi_t  width)
{
    
    
    block  = m1ri_malloc(nrows * width * sizeof(vbg) );
    
    
    
    return block;
    
    
    
}

/*
 
 */




vbg ** m3d_cube_rows(vbg * block, vbg ** rows, wi_t width, rci_t nrows)
{
    
    
    
    
    rows = m1ri_malloc( nrows * width * sizeof(vbg *));
    
    
    for (int i = 0; i <  nrows;  i++ )
    {
        rows[i]  = (block + (i * width));
        
        
    };
    
    return rows;
}



m3d_cbe m3dcbe_create( m3d_cbe * a, rci_t nrows, rci_t ncols)
{
    
    
    a->ncols = ncols;
    a->nrows = nrows;
    a->width =  RU64(ncols);
    a->block = m3d_block_allocate(a->block,  a->nrows,    a->width);
    a->rows  = m3d_row_alloc(a->block, a->rows, a->width, a->nrows);
    a->flags = iswindowed;
    
    return *a;
    
}
/*
m3d_cbe m3d_cbe_form(m3d_t *a)
{
   
   

    
    for(int i = 0; i < a->nrows; i = i + 64)
    {
        

    
    }
        
        
        
        
}

*/

#endif
