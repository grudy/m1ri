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




#include "m1ri_3dt.h"
#include "m5d.h"


typedef struct
{
 
    
    m3d_t * block;
    m3d_t ** row;
    wi_t slicesize;// (slicesize ^ 2) * 64
    wi_t width;   ///width in slices horizaontally per row
    rci_t nrows;
    rci_t ncols;
    
}m3_slice;



typedef struct
{
    
    
    m5d_t * block;
    m5d_t ** row;
    wi_t slicesize;// (slicesize ^ 2) * 64
    wi_t width;   ///width in slices horizaontally per row
    rci_t nrows;
    rci_t ncols;
    
}m5_slice;


 void m3d_swap(m3d_t *  , m3d_t * );


vbg * m3d_transpose_vbg(vbg  **a , vbg **b);
 m3d_t  * m3_blockslice_allocate(m3d_t * block, rci_t  nrows,  wi_t  width);

 m3d_t ** m3_rowslice_allocate(m3d_t * block, m3d_t ** rows, wi_t width, rci_t nrows);

m3d_t  m3d_cubes(m3d_t *, m3d_t   * , rci_t );
void  m3d_slices(m3_slice *  , m3d_t * , wi_t );
//A direct transpose, using no windows
void  m3d_transpose(m3d_t * restrict, m3d_t * restrict);
m3d_t m3d_transpose_sliced(m3d_t * );
m5d_t  * m5_blockslice_allocate(m5d_t * , rci_t  ,  wi_t  );

m5d_t ** m5_rowslice_allocate(m5d_t * , m5d_t ** , wi_t , rci_t );

#endif
