
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
 
 m1ri_3dt.c
 */

#include "m1ri_3dt.h"
#include "m1riarith.h"
#include "m1ri_classical.h"

void * m3d_rowswap (m3d_t * M, rci_t row_a, rci_t  row_b)
{
    vbg  temp ;
    
    if((M->nrows >= (row_a ) && (M->nrows >= row_b)))
    {
        
        temp =  *M->rows[row_a -1];
        M->rows[row_a -1] = M->rows[row_b -1];
        M->rows[row_b -1] =  &temp;
        
        
        
    }
    
    
    
    return 0;
}


vec m3d_rs_bits(m3d_t const *M, rci_t  x, rci_t  y, int  n) {
    
    
    
    
    
    wi_t  block = (y  ) / m1ri_word;
    
    int  spill = (y  % m1ri_word) + n - m1ri_word;
    
    vec bits;
    
    bits  = (n == 0) ? (~(leftbit >> spill) &  (M->rows[x][block].sign))  : ((leftbit >> spill) | (M->rows[x][block].units));
    
    
    return bits;
    
    
}


vec m3d_ru_bits(m3d_t const *M, rci_t  x, rci_t  y, int  n) {
    
    
    
    
    
    wi_t  block = (y  ) / m1ri_word;
    
    int  spill = (y  % m1ri_word) + n - m1ri_word;
    
    vec bits;
    
    
    
    bits  = (n == 0) ? (~(leftbit >> spill) &  (M->rows[x][block].units))  : ((leftbit >> spill) | (M->rows[x][block].units));
    
    
    
    
    
    
    return bits;
    
    
}

vbg m3d_read_elems(m3d_t const *M, rci_t  x, rci_t  y, int  n) {
    
    
    
    
    wi_t  block = (y  ) / m1ri_word;
    
    int  spill = (y  % m1ri_word) + n - m1ri_word;
    
    vbg elem;
    
    elem.units = (spill <= 0) ? M->rows[x][block].units << -spill : ((M->rows[x][(block + 1)].units<< (64 - spill)) | (M->rows[x][block].units >> spill));
    
    elem.sign = (spill <= 0) ?  (M->rows[x][block].sign << -spill) : (M->rows[x][block + 1].sign << (64 - spill)) | (M->rows[x][block].sign>> spill);
    
    elem.units = (elem.units >> (m1ri_word - n));
    
    elem.sign = (elem.sign >> (m1ri_word - n));
    
    
    
    return elem;
    
    
}


void * m3d_colswap(m3d_t *M, rci_t col_a, rci_t col_b)
{
    if((M->ncols >= (col_a ) && (M->nrows >= col_b)))
    {
        
        vec block_a = col_a/m1ri_word;
        vec block_b = col_b/m1ri_word;
        vec dif_a = col_a%m1ri_word;
        vec dif_b = col_b%m1ri_word;
        vec a_place =  leftbit >>  dif_a ;
        vec b_place =  leftbit >> dif_b ;
        vbg tempa;
        vbg tempb;
        
        
        for(int i = 0; i > M->nrows; i++)
        {
            
            
            
            tempa.sign = (M->rows[i][block_a].sign) & a_place;
            tempa.units = (M->rows[i][block_a].units) & a_place;
            tempb.units = (M->rows[i][block_b].units) & b_place;
            tempb.sign = (M->rows[i][block_b].sign) & b_place;
            M->rows[i][block_b].units = (tempa.sign == 0)? (~(leftbit >> dif_b) &  (M->rows[i][block_b].units))  : ((leftbit >> dif_b) | (M->rows[i][block_b].units));
            M->rows[i][block_b].sign = (tempa.sign == 0)? (~(leftbit >> dif_b) &  (M->rows[i][block_b].sign))  : ((leftbit >> dif_b) | (M->rows[i][block_b].sign));
            
            M->rows[i][block_a].units = (tempb.sign == 0)? (~(leftbit >> dif_a) &  (M->rows[i][block_a].sign))  : ((leftbit >> dif_a) | (M->rows[i][block_a].sign));
            M->rows[i][block_a].sign = (tempb.units == 0)? (~(leftbit >> dif_a) &  (M->rows[i][block_a].units))  : ((leftbit >> dif_a) | (M->rows[i][block_a].units));
        }
    }
    
    return M;
}
void  m3d_write_elem( m3d_t * M,rci_t x, rci_t y, vec s, vec u )
{
    wi_t  block = (y  ) / m1ri_word;
    int   spill =  (y  % m1ri_word) ;
    M->rows[x][block].units  = (u == 0) ? (~(leftbit >> spill) &  (M->rows[x][block].units))  : ((leftbit >> spill) | (M->rows[x][block].units));
    M->rows[x][block].sign  = (s == 0) ? (~(leftbit  >> spill) &  (M->rows[x][block].sign))  : ((leftbit  >> spill) | (M->rows[x][block].sign));
    
}


/*
 
 */



vbg  * m3d_block_allocate(vbg * block, rci_t  nrows,  wi_t  width)
{
    block  = m1ri_calloc(nrows * width ,  sizeof(vbg) );
    return block;
}

/*
 
 */


vbg ** m3d_row_alloc(vbg * block, vbg ** rows, wi_t width, rci_t nrows)
{
    rows = m1ri_malloc( nrows * width * sizeof(vbg *));
    for (int i = 0; i <  nrows;  i++ )
    {
        rows[i]  = block + i * width;
    };
    return rows;
}

/*
 
 */

m3d_t m3d_create( m3d_t * a, rci_t nrows, rci_t ncols)
{
    
    
    a->ncols = ncols;
    a->nrows = nrows;
    a->width =  RU64(ncols);
    a->block = m3d_block_allocate(a->block,  a->nrows,    a->width);
    a->rows  = m3d_row_alloc(a->block, a->rows, a->width, a->nrows);
    a->flags = notwindowed;
    a->lblock = ncols%64;
    a->fcol = 0;
    a->svbg = 0;
    return *a;
    
}

 


m3d_t  m3d_rand(m3d_t * a)
{
    
    for(int i = 0; i < (a->nrows * a->width); i++)
    {
        
       a->block[i].units = m1ri_rand();
        a->block[i].units = m1ri_rand();
    }
    
    return *a;
    
}


/*
 Make an Identity Matrix
 a = Identity matrix
 n = matrix size (row length and column width)
 
 */


m3d_t    m3d_identity_set(m3d_t * a)

{
    if(a->ncols == a->nrows)
    {
        int k,i,  j,l;
        for( i  = 1; i < (a->width  ) ; ++i)
        {
            l =  ((i - 1) * m1ri_word);
            j = i - 1;
            for ( k = 0 ; k < m1ri_word; k++)
            {
                
                a->rows[l][j].units = lbit[k];
                l++;
                
            }
            
        }
        
        if((a->ncols%m1ri_word) != 0)
        {
            l = a->ncols %m1ri_word;
            k = ((a->width -1) * m1ri_word);
            l = m1ri_word - l;
            for(i = 0; i < (m1ri_word - l); i++)
            {
                
                a->rows[k + i][a->width-1].units = lbit[i];
            }
            
        }
        if ((a->ncols%64) == 0)
        {
            
            l = (a->width - 1) * 64;
            for(i  = 0; i < 64; i++)
                
            {
                a->rows[l][a->width -1].units = lbit[i];
                l++;
                
            }
            
        }
        
        
    }
    return *a;
}




m3d_t   m3d_identity(m3d_t  *a, rci_t n)
{
    *a = m3d_create(a, n, n);
    *a = m3d_identity_set(a);
    return *a;
}




/*
 windows in 64 rows * 64 column incriments
 stvbg = the vbg or width offset from the base matrix
 strow = row offset in increments of 64
 sizecol  = cols * 64
 sizerow  = rows * 64
 */

m3d_t   m3d_window(m3d_t *c, rci_t strow, rci_t stvbg, rci_t sizerows, rci_t sizecols)
{
    
   
    m3d_t  submatrix;
    
    if((strow + sizerows) > c->width)
    {
    
    
        return submatrix;
    }
    
    
    
  
    if((stvbg + sizecols) > c->width)
    {
        
        
        return submatrix;
    }
    
        
    
    submatrix.nrows =   m1ri_word * sizecols; 
    submatrix.ncols =  m1ri_word * sizecols;
    submatrix.flags = iswindowed;
    submatrix.width =  sizecols;
    submatrix.block = &c->block[(stvbg * stvbg)];
    submatrix.rows = m1ri_calloc(m1ri_word * sizerows ,  sizecols * sizeof(vbg *));
    submatrix.lblock = ( (sizerows +  strow)  ==  c->width)? c->lblock:  0;
    submatrix.fcol   = 0;
    submatrix.svbg = stvbg;
   

    
    int i;
    
    for(  i =   strow; i < (strow + (m1ri_word * sizerows)) ; i++)
    {
     
       submatrix.rows[i - strow] = c->rows[i];
        
    }
    
    return submatrix;
    
}

void   m3d_window_create(m3d_t *c, m3d_t * submatrix, rci_t strow, rci_t stvbg, rci_t sizerows, rci_t sizecols)
{
  
    
    submatrix->nrows =   m1ri_word * sizecols;
    submatrix->ncols =  m1ri_word * sizecols;
    submatrix->flags = iswindowed;
    submatrix->width = 1 * sizecols;
    submatrix->block = &c->block[(stvbg * stvbg)];
    submatrix->rows = m1ri_calloc(m1ri_word * sizerows ,  sizecols * sizeof(vbg *));
    submatrix->lblock = ( (sizerows +  strow)  ==  c->width)? c->lblock:  0;
    submatrix->fcol   = 0;
    submatrix->svbg = stvbg;
    
    
    int i;
    
    for(  i =   strow; i < strow + (m1ri_word * sizerows) ; i++)
    {
        
        submatrix->rows[i - strow] = c->rows[i];
        
    }
    
    
    
}


/*
 Concat b on the end of a, the result is c
 
 
 
 [a] [b] ----->  [a b]   ===  C
 
 */
m3d_t m3d_concat(m3d_t * c, m3d_t * a, m3d_t * b)
{
    if (a->nrows != b->nrows)
    {
        /*if concat hath failed*/
        return *c;
        
    }
    
    if(a->nrows == b->nrows)
    {
        *c = m3d_create(c, a->nrows , a->ncols + b->ncols);
        int x, y;
        x =  0;
        
        
        while (x < a->nrows) {
            
            
            for(y = 0; y < a->width; y++)
            {
                c->rows[x] = a->rows[x];
            }
            
            
            for(y = a->width; y  < c->width; y++)
            {
                
                
                
            }
            x++;
        }
        
        
    }
    return  *c;
}


/*
 Stacks a on b, resulting matrix is c
 [a]
 ===  C
 [b]
 
 */
m3d_t m3d_stack(m3d_t * c,  m3d_t * a, m3d_t * b)
{
    if (a->ncols != b->ncols)
    {
        /*If a stacked matrix cannot be created*/
        return *c;
        
    }
    
    if(a->ncols == b->ncols)
    {
        *c = m3d_create(c, (a->nrows + b->nrows), a->ncols);
        int x =  0;
        while (x < a->nrows) {
            c->rows[x] = a->rows[x];
            x++;
        }
        while(x < (a->nrows + b->nrows)){
            
            c->rows[x] = a->rows[x - a->nrows];
            
            x++;
        }
        
    }
    
    return *c;
    
    
}





/*
 
 Releases a m3d_t into the wilderness.
 */




int m3d_equal(m3d_t const *a, m3d_t const *b)
{
    //  for a->nrows
    if ((a->nrows != b->nrows)    || ( a->ncols != b->ncols)  )
    {
        return 0;
    }
    int i, j;
    for( i = 0; i < a->nrows; i++)
    {
        
        for(j = 0; j < b->width; j++)
        {
            if((a->rows[i][j].sign != b->rows[i][j].sign) || (a->rows[i][j].sign != b->rows[i][j].sign))
            {
                return 0;
            }
            
        }
    }
    return 1;
}

void m3d_free( m3d_t *  tofree)
{
    
    
    m1ri_free(tofree->rows);
    m1ri_free(tofree->block);
    
}


