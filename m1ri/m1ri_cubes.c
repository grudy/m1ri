
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
 
 m1ri_cubes.h
 */

#include "m1ri_cubes.h"


static inline m3d_t  * m3_blockslice_allocate(m3d_t * block, rci_t  nrows,  wi_t  width)
{
    block  = m1ri_calloc(nrows *  width  ,   sizeof(m3d_t  ) );
    return block;
}

 static inline m3d_t ** m3_rowslice_allocate(m3d_t * block, m3d_t ** rows, wi_t width, rci_t nrows)
{
	int i;
    rows = m1ri_calloc( nrows , width  * sizeof(m3d_t *));
    for ( i = 0; i <  nrows;  i++ )
    {
        rows[i]  = block + (i * width);
    };
    return rows;
}



void  m3d_slices(m3_slice *  c, m3d_t * a, wi_t slicesize)
{
    wi_t l,  r,  colroundeddown;
    int  i,  f,extrarows ,  extracols;
    extracols = a->width%slicesize;
    colroundeddown = a->width/slicesize;
    c->width = a->width;
    extrarows = a->nrows%(  slicesize);
    c->nrows = M1RI_DN(a->nrows, (M1RI_RADIX * slicesize));
    c->ncols = M1RI_DN(a->width, slicesize);
    l = a->nrows / (M1RI_RADIX * slicesize);
    l = l * slicesize;
    c->slicesize = slicesize;
 	c->block = m3_blockslice_allocate(c->block,  c->nrows,   c->width);
    c->row = m3_rowslice_allocate(c->block,  c->row,   c->width, c->nrows);
    r = 0 ;
     
    for ( i = 0; i <  l;  i = i + slicesize)
    {       
    	for( f = 0; f <colroundeddown ; f++)
        {
        	m3d_window_create(a, &c->row[r][f],i , (f * slicesize), slicesize, slicesize);
        }
        
        if(extracols > 0)
        {
        	m3d_window_create(a, &c->row[r][f],i , (f * slicesize), slicesize, extracols);
		}
        r++;
        
    }
    
    if(extrarows >0 )
    {
		for( f = 0; f <colroundeddown ; f++)
        {
           m3d_window_create(a, &c->row[r][f],i , (f * slicesize), extrarows, slicesize);
        }

    	if(extracols > 0)
    	{
           m3d_window_create(a, &c->row[r][f],i , (f * slicesize), extrarows, extracols);
    	}
    }
}


static inline vbg *  m3d_transpose_vbg(vbg  **a, vbg  **b  )
{
    int i, x;
    vbg temp;
    for(i = 0; i <64; i ++)
    {
        for(x = 0; x < 64; x ++)
        {
          temp.units =  (a[x][0].units & (leftbit >> i) );
          temp.sign =  (a[x][0].sign & (leftbit >> i) );
          b[i][0].units = (temp.units) ?  b[i][0].units | (leftbit >> x) : b[i][0].units ;
          b[i][0].sign = (temp.sign) ? b[i][0].sign | (leftbit >> x) : b[i][0].sign  ;
        }
    }
    return *b;
}

m3d_t m3d_transpose_sliced(m3d_t * a)
{
    int x, y;
    m3d_t c;
    c = m3d_create(&c, a->ncols, a->nrows);
    m3_slice * b, *d;
    d = malloc(sizeof(m3_slice));
    b = malloc(sizeof(m3_slice));
    m3d_slices(b, a, 1);
    m3d_slices(d, &c, 1);
    for (x = 0; x < b->nrows; x++) {
        for (y = 0; y < b->ncols; y ++) {
         m3d_transpose_vbg(b->row[x][y].rows, d->row[y][x].rows);
            
        }
    }
    return c;
  
}

void m3d_quarter(m3_slice * c , m3d_t * a)
{
	
	//int arows, acols;
	c->block = m3_blockslice_allocate(c->block,  2,   2);
    c->row = m3_rowslice_allocate(c->block,  c->row,   2, 2);
    m3d_window_create(a, &c->row[0][0], 0, 0 , a->nrows/128, a->ncols/128);
	m3d_window_create(a, &c->row[0][1], 0, a->ncols/128 , a->nrows/128, a->ncols/128);   
    m3d_window_create(a, &c->row[1][0], a->nrows/128, 0 , a->nrows/128, a->ncols/128);
	m3d_window_create(a, &c->row[1][1], a->nrows/128,a->ncols/128,  a->nrows/128, a->ncols/128);
    
}
m5d_t  * m5_blockslice_allocate(m5d_t * block, rci_t  nrows,  wi_t  width)
{
    block  = m1ri_calloc(nrows * width ,  sizeof(m5d_t  ) );
    return block;
}

m5d_t ** m5_rowslice_allocate(m5d_t * block, m5d_t ** rows, wi_t width, rci_t nrows)
{
	int i;
    rows = m1ri_malloc( nrows * width * sizeof(m5d_t *));
    for ( i = 0; i <  nrows;  i++ )
    {
        rows[i]  = block + (i * width);
    };
    return rows;
}
void  m5d_slices(m5_slice *  c, m5d_t * a, wi_t slicesize)
{
    wi_t l,  r,  colroundeddown;
    int  i,  f,extrarows ,  extracols;
    extracols = a->width%slicesize;
    colroundeddown = a->width/slicesize;
    c->width = a->width;
    extrarows = a->nrows%(  slicesize);
    c->nrows = M1RI_DN(a->nrows, (M1RI_RADIX * slicesize));
    c->ncols = M1RI_DN(a->width, slicesize);
    l = a->nrows / (M1RI_RADIX * slicesize);
    l = l * slicesize;
    c->slicesize = slicesize;
 	c->block = m5_blockslice_allocate(c->block,  c->nrows,   c->width);
    c->row = m5_rowslice_allocate(c->block,  c->row,   c->width, c->nrows);
    r = 0 ;
     
    for ( i = 0; i <  l;  i = i + slicesize)
    {       
    	for( f = 0; f <colroundeddown ; f++)
        {
        	m5d_window_create(a, &c->row[r][f],i , (f * slicesize), slicesize, slicesize);
        }
        
        if(extracols > 0)
        {
        	m5d_window_create(a, &c->row[r][f],i , (f * slicesize), slicesize, extracols);
		}
        r++;
        
    }
    
    if(extrarows >0 )
    {
		for( f = 0; f <colroundeddown ; f++)
        {
           m5d_window_create(a, &c->row[r][f],i , (f * slicesize), extrarows, slicesize);
        }

    	if(extracols > 0)
    	{
           m5d_window_create(a, &c->row[r][f],i , (f * slicesize), extrarows, extracols);
    	}
    }
}

vfd *  m5d_transpose_vfd(vfd  **a, vfd  **b  )
{
    int i, x;
    vfd temp;
    for(i = 0; i <64; i ++)
    {
      for(x = 0; x < 64; x ++)
       {
       
        temp.units =  (a[x][0].units & (leftbit >> i) );
        temp.sign =  (a[x][0].sign & (leftbit >> i) );
        temp.middle =  (a[x][0].middle & (leftbit >> i) );    
        b[i][0].units = (temp.units) ?  b[i][0].units | (leftbit >> x) : b[i][0].units ;
        b[i][0].sign = (temp.sign) ? b[i][0].sign | (leftbit >> x) : b[i][0].sign  ;
        b[i][0].middle = (temp.middle) ? b[i][0].middle | (leftbit >> x) : b[i][0].middle;    
        }

    }
    
    return *b;
}

m5d_t m5d_transpose_sliced(m5d_t * a)
{
    int x, y;
    m5d_t c;
    c = m5d_create(&c, a->ncols, a->nrows);
    m5_slice * b, *d;
    d = malloc(sizeof(m5_slice));
    b = malloc(sizeof(m5_slice));
    m5d_slices(b, a, 1);
    m5d_slices(d, &c, 1);
    for (x = 0; x < b->nrows; x++) {
        for (y = 0; y < b->ncols; y ++) {
			m5d_transpose_vfd(b->row[x][y].rows, d->row[y][x].rows);  
        }
    }
    return c;
}
void m5d_quarter(m5_slice * c , m5d_t * a)
{

	//int arows, acols;
	c->block = m5_blockslice_allocate(c->block,  2,   2);
    c->row = m5_rowslice_allocate(c->block,  c->row,   2, 2);
    m5d_window_create(a, &c->row[0][0], 0, 0 , a->nrows/128, a->ncols/128);
	m5d_window_create(a, &c->row[0][1], 0, a->ncols/128 , a->nrows/128, a->ncols/128);   
    m5d_window_create(a, &c->row[1][0], a->nrows/128, 0 , a->nrows/128, a->ncols/128);
	m5d_window_create(a, &c->row[1][1], a->nrows/128,a->ncols/128,  a->nrows/128, a->ncols/128);
    
}

/**

*/
m7d_t  * m7_blockslice_allocate(m7d_t * block, rci_t  nrows,  wi_t  width)
{
    block  = m1ri_calloc(nrows * width ,  sizeof(m7d_t  ) );
    return block;
}

m7d_t ** m7_rowslice_allocate(m7d_t * block, m7d_t ** rows, wi_t width, rci_t nrows)
{
	int i;
    rows = m1ri_malloc( nrows * width * sizeof(m7d_t *));
    for ( i = 0; i <  nrows;  i++ )
    {
        rows[i]  = block + (i * width);
    };
    return rows;
}



vtri *  m7d_transpose_vtri(vtri  **a, vtri  **b  )
{
    int i, x;
    vtri temp;
    for(i = 0; i <64; i ++)
    {
      for(x = 0; x < 64; x ++)
      {
            
        temp.units =  (a[x][0].units & (leftbit >> i) );
        temp.sign =  (a[x][0].sign & (leftbit >> i) );
        temp.middle =  (a[x][0].middle & (leftbit >> i) );
        b[i][0].units = (temp.units) ?  b[i][0].units | (leftbit >> x) : b[i][0].units ;
        b[i][0].sign = (temp.sign) ? b[i][0].sign | (leftbit >> x) : b[i][0].sign  ;
        b[i][0].middle = (temp.middle) ? b[i][0].middle | (leftbit >> x) : b[i][0].middle;
    
      }
    }
    
    return *b;
}

void  m7d_slices(m7_slice *  c, m7d_t * a, wi_t slicesize)
{
    wi_t l, z, r,  colroundeddown;
    int  i,  f,extrarows ,  extracols;
    extracols = a->width%slicesize;
    colroundeddown = a->width/slicesize;
    c->width = M1RI_DN(a->width, slicesize);
    extrarows = a->nrows%(M1RI_RADIX * slicesize);
    c->nrows = M1RI_DN(a->nrows, (M1RI_RADIX * slicesize));
    c->ncols = M1RI_DN(a->width, slicesize);
    l = a->nrows / (M1RI_RADIX * slicesize);
    l = l * slicesize;
    c->slicesize = slicesize;
    c->block = m7_blockslice_allocate(c->block,  c->nrows,   c->width);
    c->row = m7_rowslice_allocate(c->block,  c->row,   c->width, c->nrows);
    z = 0;
    r = 0 ;
    
    for ( i = 0; i <  l;  i = i + slicesize)
    {
   
        for( f = 0; f <colroundeddown ; f++)
        {
            c->block[z] =  m7d_window(a, ( i ) , (f * slicesize), slicesize, slicesize);
            z++;    
        }
        
        c->block[z] =  m7d_window(a, ( i ) , (f * slicesize),slicesize, extracols);
        z++;
        c->row[r] =  c->block + (c->ncols * i );
        r++;
        
    }

    for( f = 0; f <colroundeddown ; f++)
    {
        c->block[z] =  m7d_window(a, ( i ) , (f * slicesize),extrarows, slicesize);
        z++;
    }

    c->block[z] =  m7d_window(a, ( i ) , (f * slicesize),extrarows, extracols);
	c->row[r] =   c->block + (c->ncols * i );
   
}

m7d_t m7d_transpose_sliced(m7d_t * a)
{
    int x, y;
    m7d_t c;
    c = m7d_create(&c, a->ncols, a->nrows);
    m7_slice * b, *d;
    d = malloc(sizeof(m7_slice));
    b = malloc(sizeof(m7_slice));
    m7d_slices(b, a, 1);
    m7d_slices(d, &c, 1);
    for (x = 0; x < b->nrows; x++)
     {
        for (y = 0; y < b->ncols; y ++)
         {
         
         	m7d_transpose_vtri(b->row[x][y].rows, d->row[y][x].rows);  
         	 
        }
    }
    
    return c;
}

void m7d_quarter(m7_slice * c , m7d_t * a)
{
	//int arows, acols;
	c->block = m7_blockslice_allocate(c->block,  2,   2);
    c->row = m7_rowslice_allocate(c->block,  c->row,   2, 2);
    m7d_window_create(a, &c->row[0][0], 0, 0 , a->nrows/128, a->ncols/128);
	m7d_window_create(a, &c->row[0][1], 0, a->ncols/128 , a->nrows/128, a->ncols/128);   
    m7d_window_create(a, &c->row[1][0], a->nrows/128, 0 , a->nrows/128, a->ncols/128);
	m7d_window_create(a, &c->row[1][1], a->nrows/128,a->ncols/128,  a->nrows/128, a->ncols/128);
    
}
