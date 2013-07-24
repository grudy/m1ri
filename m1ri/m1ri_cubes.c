
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



inline  void m3d_swap(m3d_t *  a, m3d_t *  b)
{
    vbg temp[64];
    temp[0] = a->rows[0][0]; temp[1] = a->rows[1][0];
    temp[2] = a->rows[2][0]; temp[3] = a->rows[3][0];
    temp[4] = a->rows[4][0]; temp[5] = a->rows[5][0];
    temp[6] = a->rows[6][0]; temp[7] = a->rows[7][0];
    temp[8] = a->rows[8][0]; temp[9] = a->rows[9][0];
    temp[10] = a->rows[10][0]; temp[11] = a->rows[11][0];
    temp[12] = a->rows[12][0]; temp[13] = a->rows[13][0];
    temp[14] = a->rows[14][0]; temp[15] = a->rows[15][0];
    temp[16] = a->rows[16][0]; temp[17] = a->rows[17][0];
    temp[18] = a->rows[18][0]; temp[19] = a->rows[19][0];
    temp[20] = a->rows[20][0]; temp[21] = a->rows[21][0];
    temp[22] = a->rows[22][0]; temp[23] = a->rows[23][0];
    temp[24] = a->rows[24][0]; temp[25] = a->rows[25][0];
    temp[26] = a->rows[26][0]; temp[27] = a->rows[27][0];
    temp[28] = a->rows[28][0]; temp[29] = a->rows[29][0];
    temp[30] = a->rows[30][0]; temp[31] = a->rows[31][0];
    temp[32] = a->rows[32][0]; temp[33] = a->rows[33][0];
    temp[34] = a->rows[34][0]; temp[35] = a->rows[35][0];
    temp[36] = a->rows[36][0]; temp[37] = a->rows[37][0];
    temp[38] = a->rows[38][0]; temp[39] = a->rows[39][0];
    temp[40] = a->rows[40][0]; temp[41] = a->rows[41][0];
    temp[42] = a->rows[42][0]; temp[43] = a->rows[43][0];
    temp[44] = a->rows[44][0]; temp[45] = a->rows[45][0];
    temp[46] = a->rows[46][0]; temp[47] = a->rows[47][0];
    temp[48] = a->rows[48][0]; temp[49] = a->rows[49][0];
    temp[50] = a->rows[50][0]; temp[51] = a->rows[51][0];
    temp[52] = a->rows[52][0]; temp[53] = a->rows[53][0];
    temp[54] = a->rows[54][0]; temp[55] = a->rows[55][0];
    temp[56] = a->rows[56][0]; temp[57] = a->rows[57][0];
    temp[58] = a->rows[58][0]; temp[59] = a->rows[59][0];
    temp[60] = a->rows[60][0]; temp[61] = a->rows[61][0];
    temp[62] = a->rows[62][0]; temp[63] = a->rows[63][0];
    a->rows[0][0] = b->rows[0][0]; a->rows[1][0] = b->rows[1][0];
    a->rows[2][0] = b->rows[2][0]; a->rows[3][0] = b->rows[3][0];
    a->rows[4][0] = b->rows[4][0]; a->rows[5][0] = b->rows[5][0];
    a->rows[6][0] = b->rows[6][0]; a->rows[7][0] = b->rows[7][0];
    a->rows[8][0] = b->rows[8][0]; a->rows[9][0] = b->rows[9][0];
    a->rows[10][0] = b->rows[10][0]; a->rows[11][0] = b->rows[11][0];
    a->rows[12][0] = b->rows[12][0]; a->rows[13][0] = b->rows[13][0];
    a->rows[14][0] = b->rows[14][0]; a->rows[15][0] = b->rows[15][0];
    a->rows[16][0] = b->rows[16][0]; a->rows[17][0] = b->rows[17][0];
    a->rows[18][0] = b->rows[18][0]; a->rows[19][0] = b->rows[19][0];
    a->rows[20][0] = b->rows[20][0]; a->rows[21][0] = b->rows[21][0];
    a->rows[22][0] = b->rows[22][0]; a->rows[23][0] = b->rows[23][0];
    a->rows[24][0] = b->rows[24][0]; a->rows[25][0] = b->rows[25][0];
    a->rows[26][0] = b->rows[26][0]; a->rows[27][0] = b->rows[27][0];
    a->rows[28][0] = b->rows[28][0]; a->rows[29][0] = b->rows[29][0];
    a->rows[30][0] = b->rows[30][0]; a->rows[31][0] = b->rows[31][0];
    a->rows[32][0] = b->rows[32][0]; a->rows[33][0] = b->rows[33][0];
    a->rows[34][0] = b->rows[34][0]; a->rows[35][0] = b->rows[35][0];
    a->rows[36][0] = b->rows[36][0]; a->rows[37][0] = b->rows[37][0];
    a->rows[38][0] = b->rows[38][0]; a->rows[39][0] = b->rows[39][0];
    a->rows[40][0] = b->rows[40][0]; a->rows[41][0] = b->rows[41][0];
    a->rows[42][0] = b->rows[42][0]; a->rows[43][0] = b->rows[43][0];
    a->rows[44][0] = b->rows[44][0]; a->rows[45][0] = b->rows[45][0];
    a->rows[46][0] = b->rows[46][0]; a->rows[47][0] = b->rows[47][0];
    a->rows[48][0] = b->rows[48][0]; a->rows[49][0] = b->rows[49][0];
    a->rows[50][0] = b->rows[50][0]; a->rows[51][0] = b->rows[51][0];
    a->rows[52][0] = b->rows[52][0]; a->rows[53][0] = b->rows[53][0];
    a->rows[54][0] = b->rows[54][0]; a->rows[55][0] = b->rows[55][0];
    a->rows[56][0] = b->rows[56][0]; a->rows[57][0] = b->rows[57][0];
    a->rows[58][0] = b->rows[58][0]; a->rows[59][0] = b->rows[59][0];
    a->rows[60][0] = b->rows[60][0]; a->rows[61][0] = b->rows[61][0];
    a->rows[62][0] = b->rows[62][0]; a->rows[63][0] = b->rows[63][0];
    b->rows[0][0] = temp[0]; b->rows[1][0] = temp[1];
    b->rows[2][0] = temp[2]; b->rows[3][0] = temp[3];
    b->rows[4][0] = temp[4]; b->rows[5][0] = temp[5];
    b->rows[6][0] = temp[6]; b->rows[7][0] = temp[7];
    b->rows[8][0] = temp[8]; b->rows[9][0] = temp[9];
    b->rows[10][0] = temp[10]; b->rows[11][0] = temp[11];
    b->rows[12][0] = temp[12]; b->rows[13][0] = temp[13];
    b->rows[14][0] = temp[14]; b->rows[15][0] = temp[15];
    b->rows[16][0] = temp[16]; b->rows[17][0] = temp[17];
    b->rows[18][0] = temp[18]; b->rows[19][0] = temp[19];
    b->rows[20][0] = temp[20]; b->rows[21][0] = temp[21];
    b->rows[22][0] = temp[22]; b->rows[23][0] = temp[23];
    b->rows[24][0] = temp[24]; b->rows[25][0] = temp[25];
    b->rows[26][0] = temp[26]; b->rows[27][0] = temp[27];
    b->rows[28][0] = temp[28]; b->rows[29][0] = temp[29];
    b->rows[30][0] = temp[30]; b->rows[31][0] = temp[31];
    b->rows[32][0] = temp[32]; b->rows[33][0] = temp[33];
    b->rows[34][0] = temp[34]; b->rows[35][0] = temp[35];
    b->rows[36][0] = temp[36]; b->rows[37][0] = temp[37];
    b->rows[38][0] = temp[38]; b->rows[39][0] = temp[39];
    b->rows[40][0] = temp[40]; b->rows[41][0] = temp[41];
    b->rows[42][0] = temp[42]; b->rows[43][0] = temp[43];
    b->rows[44][0] = temp[44]; b->rows[45][0] = temp[45];
    b->rows[46][0] = temp[46]; b->rows[47][0] = temp[47];
    b->rows[48][0] = temp[48]; b->rows[49][0] = temp[49];
    b->rows[50][0] = temp[50]; b->rows[51][0] = temp[51];
    b->rows[52][0] = temp[52]; b->rows[53][0] = temp[53]; 
    b->rows[54][0] = temp[54]; b->rows[55][0] = temp[55]; 
    b->rows[56][0] = temp[56]; b->rows[57][0] = temp[57]; 
    b->rows[58][0] = temp[58]; b->rows[59][0] = temp[59]; 
    b->rows[60][0] = temp[60]; b->rows[61][0] = temp[61]; 
    b->rows[62][0] = temp[62]; b->rows[63][0] = temp[63]; 

   


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
   
    int  nrows;
    int x,  z, i, l, mod;

    nrows = DN(a->nrows, (m1ri_word ));
   
    
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




