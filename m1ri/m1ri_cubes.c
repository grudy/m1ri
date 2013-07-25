
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
 
 m1ri_cubes.h
 */

#include "m1ri_cubes.h"

//This copies a matrix as a contingous matrix of in (slicesize ^2) * 64 slices 
m3d_t  m3d_cubes(m3d_t * c, m3d_t  *a, rci_t slicesize )
{

          u_int64_t l, i, f, x, y, z, lf, r, extracols, extrarows, colroundeddown;
       // int l, i, f, x, y, z, lf, r, extracols, extrarows, colroundeddown;
        extracols = a->width%slicesize;
        colroundeddown = a->width/slicesize;
        c->ncols = DN(a->width, slicesize);
        extrarows = a->nrows%(m1ri_word * slicesize);
        c->nrows = DN(a->nrows, (m1ri_word * slicesize));
        c->width =  a->width;
    l = a->nrows / (m1ri_word * slicesize);
   l  = l  * m1ri_word * slicesize;
    
        c->block = m3d_block_allocate(a->block,  a->nrows,   a->width);
    z = 0;
    r = 0 ;
    
      c->rows = m1ri_malloc( a->nrows * a->width * sizeof(vbg *));
        for ( i = 0; i <  l;  i = i + (slicesize* m1ri_word))
        {
            for( f = 0; f <colroundeddown ; f++)
            {   
                 lf = f * slicesize; 
                for( x = 0; x <(slicesize  * m1ri_word ) ; x++)
                {
                    for (y = 0 ; y < slicesize; y++) {
                       
                       c->block[z] = a->rows[i+ x][lf  + y];
                        z++;
                        
                    }

                }
   
                
            }
            
            for( x = 0; x <(slicesize  * m1ri_word ) ; x++)
            {
                for (y = 0 ; y < extracols; y++) {
                    
                  c->block[z] = a->rows[i+ x][(f * slicesize) + y];
                    z++;
                    
                }
                
                
            }
            
        
        c->rows[r]  = c->block + (i * slicesize );
            r++;
            
        }
    
    
    for( f = 0; f <colroundeddown ; f++)
    {
        for( x = 0; x <(extrarows   ) ; x++)
        {
            for (y = 0 ; y <  slicesize; y++) {
                
                c->block[z] = a->rows[i+ x][(f * slicesize) + y];
                z++;
                
            }
            
            
            
        }
    
        
    }

    for( x = 0; x <(extrarows  ) ; x++)
    {
        for (y = 0 ; y < extracols; y++) {
            
            c->block[z] = a->rows[i+ x][(f * slicesize) + y];
            z++;
            
        }
        
        
    }
    
       /// r->rows  = m3d_row_alloc(a->block, a->rows, a->width, a->nrows);
        c->flags = notwindowed;
        c->lblock = a->ncols%64;
        c->fcol = 0;
        c->svbg = 0;
        return *c;
      c->rows[r]  = c->block + (i * slicesize );  
    

    
    
    
}





 m3d_t  * m3_blockslice_allocate(m3d_t * block, rci_t  nrows,  wi_t  width)
{
        block  = m1ri_calloc(nrows *  width  ,   sizeof(m3d_t  ) );
    return block;
}

 m3d_t ** m3_rowslice_allocate(m3d_t * block, m3d_t ** rows, wi_t width, rci_t nrows)
{
    rows = m1ri_calloc( nrows , width  * sizeof(m3d_t *));
    for (int i = 0; i <  nrows;  i++ )
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
    c->nrows = DN(a->nrows, (m1ri_word * slicesize));
    c->ncols = DN(a->width, slicesize);
    //printf("extracols: %d    extrarows %d", extracols, extrarows);
    
    l = a->nrows / (m1ri_word * slicesize);
    l = l * slicesize;
    //printf("%d",  l);
    //colpassed = (slicesize * m1ri_word) - 1;
    c->slicesize = slicesize;
 
    c->block = m3_blockslice_allocate(c->block,  c->nrows,   c->width);
    c->row = m3_rowslice_allocate(c->block,  c->row,   c->width, c->nrows);
    //z = 0;
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
         //c->row[ r][f]    = m3d_window(a, i , (f * slicesize), slicesize, extracols);
            
     
        }
            
            
       
        r++;
        
    }
    
    
    if(extrarows >0 )
    {
   for( f = 0; f <colroundeddown ; f++)
        {
           m3d_window_create(a, &c->row[r][f],i , (f * slicesize), extrarows, slicesize);
       //c->row[r][f] =  m3d_window(a, i , (f * slicesize),extrarows, slicesize);
   
        }
        
        
    
    
    if(extracols > 0)
    {
   
           m3d_window_create(a, &c->row[r][f],i , (f * slicesize), extrarows, extracols);
        //c->row[r][f] =  m3d_window(a,i  , (f * slicesize),extrarows, extracols);
  
    }
    
    
    

    }
    
    
    

    
}

// 
  void m3d_transpose(m3d_t  * restrict a, m3d_t *  restrict b)
{
   
    *b = m3d_create(b, a->ncols, a->nrows);
   
    //int  nrows;
    int x,  z, i, l, mod;

    //nrows = DN(a->nrows, (m1ri_word ));
   
    
    vbg temp;
    for(i = 0; i < a->nrows; i++)
    {
        l = i/64;
        mod  = i%64;
       
            for(z = 0; z < a->width; z ++)
            {
                for(x = 0; x < m1ri_word; x ++)
                {
                  
                    temp.units =  (a->rows[i][z].units & (leftbit >> x) );
                    temp.sign =  (a->rows[i][z].sign & (leftbit >> x) );
                    
                    b->rows[(z * 64) + x][l].units = (temp.units) ?   b->rows[(z * 64) + x][l].units| (leftbit >> mod) :  b->rows[(z * 64) + x][l].units;
                    b->rows[(z * 64) + x][l].sign = (temp.sign) ?b->rows[(z * 64) + x][l].sign | (leftbit >> mod) : b->rows[(z * 64) + x][l].sign  ;
                                      
                    
                    
                }
                
                
            }
            
            
            
        }
    
 
   // return *b;
  //  m1ri_free(b);

}
vbg *  m3d_transpose_vbg(vbg  **a, vbg  **b  )
{
    //m3d_t *b = m1ri_malloc(sizeof(m3d_t));
    
    //m3d_create(b, 64, 64);
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
            // b->rows[i][0].units = (temp.units) ? : ;
            
            
            
            
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
    //m1ri_free(temp);
}

m5d_t  * m5_blockslice_allocate(m5d_t * block, rci_t  nrows,  wi_t  width)
{
    block  = m1ri_calloc(nrows * width ,  sizeof(m5d_t  ) );
    return block;
}

m5d_t ** m5_rowslice_allocate(m5d_t * block, m5d_t ** rows, wi_t width, rci_t nrows)
{
    rows = m1ri_malloc( nrows * width * sizeof(m5d_t *));
    for (int i = 0; i <  nrows;  i++ )
    {
        rows[i]  = block + i * width;
    };
    return rows;
}


void  m5d_slices(m5_slice *  c, m5d_t * a, wi_t slicesize)
{
    wi_t l, z, r,  colroundeddown;
    int  i,  f,extrarows ,  extracols;
    extracols = a->width%slicesize;
    colroundeddown = a->width/slicesize;
    c->width = DN(a->width, slicesize);
    extrarows = a->nrows%(m1ri_word * slicesize);
    c->nrows = DN(a->nrows, (m1ri_word * slicesize));
    c->ncols = DN(a->width, slicesize);
    l = a->nrows / (m1ri_word * slicesize);
    l = l * slicesize;
    //colpassed = (slicesize * m1ri_word) - 1;
    c->slicesize = slicesize;
    c->block = m5_blockslice_allocate(c->block,  c->nrows,   c->width);
    c->row = m5_rowslice_allocate(c->block,  c->row,   c->width, c->nrows);
    z = 0;
    r = 0 ;
    
    
    for ( i = 0; i <  l;  i = i + slicesize)
    {
   
        for( f = 0; f <colroundeddown ; f++)
        {
            
            c->block[z] =  m5d_window(a, ( i ) , (f * slicesize), slicesize, slicesize);
            // c->block[z] = a->rows[i+ x][lf  + y];
            z++;
            
            
        }
        
        
        c->block[z] =  m5d_window(a, ( i ) , (f * slicesize),slicesize, extracols);
        z++;
        

        c->row[r] =  c->block + (c->ncols * i );
        //   c->rows[r]  = c->block + (i * slicesize );
        r++;
        
    }
    

    
    for( f = 0; f <colroundeddown ; f++)
    {
        c->block[z] =  m5d_window(a, ( i ) , (f * slicesize),extrarows, slicesize);
        //  c->block[z] = a->rows[i+ x][(f * slicesize) + y];
        z++;
        
    }
    
    
    c->block[z] =  m5d_window(a, ( i ) , (f * slicesize),extrarows, extracols);
    
  
    
    c->row[r] =   c->block + (c->ncols * i );
    //   c->rows[r]  = c->block + (i * slicesize );
    

    
}





m5d_t  m5d_cubes(m5d_t * c, m5d_t  *a, rci_t slicesize )
{
    
    int l, i, f, x, y, z, lf, r, extracols, extrarows, colroundeddown;
    extracols = a->width%slicesize;
    colroundeddown = a->width/slicesize;
    c->ncols = DN(a->width, slicesize);
    extrarows = a->nrows%(m1ri_word * slicesize);
    c->nrows = DN(a->nrows, (m1ri_word * slicesize));
    c->width =  a->width;
    l = a->nrows / (m1ri_word * slicesize);
    l  = l  * m1ri_word * slicesize;
    
    c->block = m5d_block_allocate(a->block,  a->nrows,   a->width);
    //c->rows  = m5d_row_alloc(c->block, <#vfd **rows#>, <#wi_t width#>, <#rci_t nrows#>)
    z = 0;
    r = 0 ;
    f = 0;
    c->rows = m1ri_malloc( a->nrows * a->width * sizeof(vbg *));
    for ( i = 0; i <  l;  i = i + (slicesize* m1ri_word))
    {
        for( f = 0; f <colroundeddown ; f++)
        {
            
            for( x = 0; x <(slicesize  * m1ri_word ) ; x++)
            {
                for (y = 0 ; y < slicesize; y++) {
                    lf = f * slicesize;
                    c->block[z] = a->rows[i+ x][lf  + y];
                    z++;
                    
                }
                
            }
            
            
        }
        
        for( x = 0; x <(slicesize  * m1ri_word ) ; x++)
        {
            for (y = 0 ; y < extracols; y++) {
                
                c->block[z] = a->rows[i+ x][(f * slicesize) + y];
                z++;
                
            }
            
            
        }
        
        
        c->rows[r]  = c->block + (i * slicesize );
        r++;
        
    }
    
    
    
    
    for( x = 0; x <(extrarows   ) ; x++)
    {
        for (y = 0 ; y <  slicesize; y++) {
            
            c->block[z] = a->rows[i+ x][(f * slicesize) + y];
            z++;
            
        }
        
        
        
    }
    
    
    
    
    
    
    
    /// r->rows  = m5d_row_alloc(a->block, a->rows, a->width, a->nrows);
    c->flags = notwindowed;

    c->fcol = 0;
   
    return *c;
    
    
    
    
    
    
}




