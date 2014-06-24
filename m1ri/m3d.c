
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
 
  Matrix Represenations and basic operations over GF(3)
 m3d.c
 */

#include "m3d.h"
#include "m1riarith.h"


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


vec m3d_rs_bits(m3d_t const *M, rci_t  x, rci_t  y, int  n) 
{
    
    wi_t  block = (y  ) / M1RI_RADIX;
    int  spill = (y  % M1RI_RADIX) + n - M1RI_RADIX;
    vec bits;
    bits  = (n == 0) ? (~(leftbit >> spill) &  (M->rows[x][block].sign))  : ((leftbit >> spill) | (M->rows[x][block].units));
    return bits;
}


vec m3d_ru_bits(m3d_t const *M, rci_t  x, rci_t  y, int  n) 
{
	wi_t  block = (y  ) / M1RI_RADIX;
	int  spill = (y  % M1RI_RADIX) + n - M1RI_RADIX;
    vec bits;
    bits  = (n == 0) ? (~(leftbit >> spill) &  (M->rows[x][block].units))  : ((leftbit >> spill) | (M->rows[x][block].units));
    return bits;   
}

vbg m3d_read_elems(m3d_t const *M, rci_t  x, rci_t  y, int  n) 
{
    wi_t  block = (y  ) / M1RI_RADIX;
    int  spill = (y  % M1RI_RADIX) + n - M1RI_RADIX;
    vbg elem;
    
    elem.units = (spill <= 0) ? M->rows[x][block].units << -spill : ((M->rows[x][(block + 1)].units<< (64 - spill)) | (M->rows[x][block].units >> spill));
    
    elem.sign = (spill <= 0) ?  (M->rows[x][block].sign << -spill) : (M->rows[x][block + 1].sign << (64 - spill)) | (M->rows[x][block].sign>> spill);
    
    elem.units = (elem.units >> (M1RI_RADIX - n));
    
    elem.sign = (elem.sign >> (M1RI_RADIX - n));
    
    
    
    return elem;
    
    
}


void  m3d_colswap(m3d_t *M, rci_t col_a, rci_t col_b)
{
    if((M->ncols >= (col_a ) && (M->ncols >= col_b)))
    {
        int i;
        vec block_a, block_b, dif_a, dif_b, a_place, b_place; 
        vbg tempa, tempb;
         block_a = (col_a-1)/M1RI_RADIX;
         block_b = (col_b-1)/M1RI_RADIX;
         dif_a = col_a%M1RI_RADIX;
         dif_b = col_b%M1RI_RADIX;
         a_place =  leftbit >>  dif_a ;
         b_place =  leftbit >> dif_b ;
        if(block_a == block_b)
        { 

              
          for( i = 0; i < M->nrows; i++)
          {
		     
		  
		           tempa.units  = (b_place  & M->rows[i][block_b].units) ? (a_place  ): 0;
		     tempb.units  = (a_place  & M->rows[i][block_a].units) ? (b_place  ): 0; 
		       M->rows[i][block_a].units  = (tempa.units)  ? M->rows[i][block_a].units  | tempa.units :   M->rows[i][block_a].units  & ~a_place; 
		       M->rows[i][block_b].units  = (tempb.units)  ? M->rows[i][block_a].units  | tempb.units :   M->rows[i][block_b].units  & ~b_place;  
		         tempa.sign  = (b_place  & M->rows[i][block_b].sign) ? (a_place  ): 0;
		     tempb.sign  = (a_place  & M->rows[i][block_a].sign) ? (b_place  ): 0; 
		       M->rows[i][block_a].sign  = (tempa.sign)  ? (M->rows[i][block_a].sign  | tempa.sign) :   M->rows[i][block_a].sign  & ~a_place; 
		       M->rows[i][block_b].sign  = (tempb.sign)  ? (M->rows[i][block_a].sign  | tempb.sign) :   M->rows[i][block_b].sign  & ~b_place; 
		     
		       

		       
          }
    
        }
        
      
        
        
    }
    
}

void m3d_colswap_capped_row(m3d_t *M, rci_t col_a, rci_t col_b, rci_t start_row)
{
    if((M->ncols >= (col_a ) && (M->ncols >= col_b)))
    {
        int i;
        vec block_a, block_b, dif_a, dif_b, a_place, b_place; 
        vbg tempa, tempb;
         block_a = (col_a-1)/M1RI_RADIX;
         block_b = (col_b-1)/M1RI_RADIX;
         dif_a = col_a%M1RI_RADIX;
         dif_b = col_b%M1RI_RADIX;
         a_place =  leftbit >>  dif_a ;
         b_place =  leftbit >> dif_b ;
        if(block_a == block_b)
        { 

              
          for( i = start_row; i < M->nrows; i++)
          {
		     
		  
		           tempa.units  = (b_place  & M->rows[i][block_b].units) ? (a_place  ): 0;
		     tempb.units  = (a_place  & M->rows[i][block_a].units) ? (b_place  ): 0; 
		       M->rows[i][block_a].units  = (tempa.units)  ? M->rows[i][block_a].units  | tempa.units :   M->rows[i][block_a].units  & ~a_place; 
		       M->rows[i][block_b].units  = (tempb.units)  ? M->rows[i][block_a].units  | tempb.units :   M->rows[i][block_b].units  & ~b_place;  
		         tempa.sign  = (b_place  & M->rows[i][block_b].sign) ? (a_place  ): 0;
		     tempb.sign  = (a_place  & M->rows[i][block_a].sign) ? (b_place  ): 0; 
		       M->rows[i][block_a].sign  = (tempa.sign)  ? (M->rows[i][block_a].sign  | tempa.sign) :   M->rows[i][block_a].sign  & ~a_place; 
		       M->rows[i][block_b].sign  = (tempb.sign)  ? (M->rows[i][block_a].sign  | tempb.sign) :   M->rows[i][block_b].sign  & ~b_place; 
		     
		       

		       
          }
    
        }
        
      
        
        
    }
    


}


void  m3d_write_elem( m3d_t * M,rci_t x, rci_t y, vec s, vec u )
{
    wi_t  block = (y  ) / M1RI_RADIX;
    int   spill =  (y  % M1RI_RADIX) ;
    M->rows[x][block].units  = (u == 0) ? (~(leftbit >> spill) &  (M->rows[x][block].units))  : ((leftbit >> spill) | (M->rows[x][block].units));
    M->rows[x][block].sign  = (s == 0) ? (~(leftbit  >> spill) &  (M->rows[x][block].sign))  : ((leftbit  >> spill) | (M->rows[x][block].sign));
    
}

m3d_t m3d_create( m3d_t *  a, rci_t nrows, rci_t ncols)
{
    
    a->ncols = ncols;
    a->nrows = nrows;
    a->width =  M1RI_DN(ncols, M1RI_RADIX);
    a->block =  NULL;
    a->rows =  NULL;
    a->block = m3d_block_allocate(a->block,  a->nrows,    a->width);
    a->rows  = m3d_row_alloc(a->block, a->rows, a->width, a->nrows);
    a->flags = 0;
    a->lblock = ncols%64;
    a->fcol = 0;
    a->svbg = 0;
    return *a;
    
}
m3d_t  m3d_rand(m3d_t * a)
{
    int i,  z;
    for(i = 0; i < (a->nrows); i++)
    {
        for( z = 0; z  < (a->width); z++)
            {  
       			a->rows[i][z].sign = m1ri_rand();
       			a->rows[i][z].units = m1ri_rand();
       			a->rows[i][z].sign =  a->rows[i][z].sign & a->rows[i][z].units;
            
            
            }
    
    }
    
    return *a;    
}


/** 
 Make an Identity Matrix
 a = Identity matrix
 n = matrix size (row length and column width)
 
 */
m3d_t m3d_transposewin(m3d_t const *a )
{
    m3d_t *b = m1ri_malloc(sizeof(m3d_t));
    m3d_create(b, a->nrows, a->ncols);
	int i, x;
	vbg temp;
for(i = 0; i < a->nrows; i ++)
{
	for(x = 0; x < a->ncols; x ++)
	{
	
		temp.units =  (a->rows[x][0].units & (leftbit >> i) );
		temp.sign =  (a->rows[x][0].sign & (leftbit >> i) );
		b->rows[i][0].units = (temp.units) ?  b->rows[i][0].units | (leftbit >> x) : b->rows[i][0].units ;
        b->rows[i][0].sign = (temp.sign) ? b->rows[i][0].sign | (leftbit >> x) : b->rows[i][0].sign  ;      
	}
	

    }
                                                
        return *b;
                                                
}

                                                                                        
m3d_t m3d_identity_set(m3d_t * a)

{
    if(a->ncols == a->nrows)
    {
        int k,i, j,l;
        for( i = 1; i < (a->width ) ; ++i)
        {
            l = ((i - 1) * M1RI_RADIX);
            j = i - 1;
            for ( k = 0 ; k < M1RI_RADIX; k++)
            {
                
              a->rows[l][j].units = (leftbit)>>k;
                
                l++;
                
            }
            
        }
        
        if((a->ncols%M1RI_RADIX) != 0)
        {
            l = a->ncols %M1RI_RADIX;
            k = ((a->width -1) * M1RI_RADIX);
            l = M1RI_RADIX - l;
            for(i = 0; i < (M1RI_RADIX - l); i++)
            {
                
              a->rows[k + i][a->width-1].units = (leftbit)>>i;
            }
            
        }
        if ((a->ncols%64) == 0)
        {
            
            l = (a->width - 1) * 64;
            for(i = 0; i < 64; i++)
                
            {
              a->rows[l][a->width -1].units = (leftbit)>>i;
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

/** 
 windows in M1RI_RADIX rows * M1RI_RADIX column incriments
 stvbg = the vbg or width offset from the base matrix
 strow = row offset in increments of 64
 sizecol  = cols * 64
 sizerow  = rows * 64
 */

void   m3d_window(m3d_t *c,m3d_t * submatrix, rci_t strow, rci_t stvbg, rci_t sizerows, rci_t sizecols)
{

    
    /** c->width should not be compared twice */
    if((strow + sizerows) > c->width)
    {    
        return ;
    }
    
    if((stvbg + sizecols) > c->width)
    {
        return ;
    }
    int i, f;
	f = strow * M1RI_RADIX;
    submatrix->nrows =   M1RI_RADIX * sizerows;
    submatrix->ncols =  M1RI_RADIX * sizecols;
    submatrix->flags = iswindowed;
    submatrix->width =  sizecols;
    submatrix->rows = m1ri_calloc(M1RI_RADIX * sizerows ,  sizecols * sizeof(vbg *));
    submatrix->lblock = ( (sizerows +  strow)  ==  c->width)? c->lblock:  0;
    submatrix->fcol   = 0;
    submatrix->svbg = stvbg;
    
    
    
    for(  i =   f; i < (f + (M1RI_RADIX * sizerows)) ; i++)
    {
        submatrix->rows[i - f] = c->rows[i] + stvbg;   
    }
    
    
}
static inline void   m3d_window_create(m3d_t *c, m3d_t * submatrix, rci_t strow, rci_t stvbg, rci_t sizerows, rci_t sizecols)
{
     /** c->width should not be compared twice */
    
    if((strow + sizerows) > c->width)
    {   
        return;
    }
    
    if((stvbg + sizecols) > c->width)
    {
        return;    
    }
    int f = strow << 6;
    
    submatrix->nrows =   M1RI_RADIX * sizerows;
    submatrix->ncols =  M1RI_RADIX * sizecols;
    submatrix->flags = iswindowed;
    submatrix->width =  sizecols;
    submatrix->block =      m1ri_calloc(sizecols * sizerows, sizeof(m3d_t));
    submatrix->block = &c->block[(strow * stvbg)];
    submatrix->rows = m1ri_calloc(M1RI_RADIX * sizerows ,  sizecols * sizeof(vbg *));
    submatrix->lblock = ( (sizerows +  strow)  ==  c->width)? c->lblock:  0;
    submatrix->fcol   = 0;
    submatrix->svbg = stvbg;
   
    for( int i =   f; i < (f + (M1RI_RADIX * sizerows)) ; i++)
    {    
        submatrix->rows[i - f] = c->rows[i] + stvbg;
    }
    
}

/** 
 Concat b on the end of a, the result is c 
 [a] [b] ----->  [a b]   ===  C
 */
 /* This function still needs work*/
 
 void m3d_concat(m3d_t * c, m3d_t * a, m3d_t * b)
{
    if (a->nrows != b->nrows)
    {
        /* if concat hath failed*/
        return;  
    }
    
    if(a->nrows == b->nrows)
    {
    	m3d_create(c, a->nrows , a->ncols + b->ncols);
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
    return ;
}

/** 
 Stacks a on b, resulting matrix is c
 [a]
 ===  C
 [b]
 */
 
m3d_t m3d_stack(m3d_t * c,  m3d_t * a, m3d_t * b)
{
    if (a->ncols != b->ncols)
    {
        /* If a stacked matrix cannot be created*/
        return *c;
    }
    
    if(a->ncols == b->ncols)
    {
        *c = m3d_create(c, (a->nrows + b->nrows), a->ncols);
        int x =  0;
        while (x < a->nrows)
        {
            c->rows[x] = a->rows[x];
            x++;
        }
        while(x < (a->nrows + b->nrows))
        {
        
            c->rows[x] = a->rows[x - a->nrows];    
            x++;
        }
        
    }
    
    return *c;
}

void m3d_copypadding(m3d_t  * r, m3d_t  const * x)
{
		int i, s;
        for( i = 0; i < x->nrows; i++)
        {
        	 for( s = 0; s < x->width; s++)
        	 {
            r->rows[i][s] = x->rows[i][s];
            }
            
        }
        
}


void m3d_putpadding(m3d_t  * r, m3d_t  const * x)
{
		int i, s;
        for( i = 0; i < r->nrows; i++)
        {
        	for( s = 0; s < r->width; s++)
        	{
            	r->rows[i][s] = x->rows[i][s];
            }          
        }
	
}

/**
  Checks if an m3d_t is equal to another.
*/
int m3d_equal(m3d_t const *a, m3d_t const *b)
{
    if ((a->nrows != b->nrows)    || ( a->ncols != b->ncols)  )
    {
        return 0;
    }
    int i, j;
    for( i = 0; i < a->nrows; i++)
    {
        
        for(j = 0; j < b->width; j++)
        {
            if((a->rows[i][j].sign != b->rows[i][j].sign) || (a->rows[i][j].units != b->rows[i][j].units))
            {
                printf("row [%d][%d] not equal \n", i, j);
                return 0;
            }
            
        }
    }
    return 1;
}

/** 
 Releases a m3d_t into the wilderness.
 */

void m3d_free( m3d_t *  tofree)
{ 		
    m1ri_free(tofree->rows);
    m1ri_free(tofree->block);
}




/**
   Allocates blocks for an m3d_slice
*/
static inline m3d_t  * m3_blockslice_allocate(m3d_t * block, rci_t  nrows,  wi_t  width)
{
    block  = m1ri_calloc(nrows *  width  ,   sizeof(m3d_t  ) );
    return block;
}



/**
   Allocates rows for an m3d_slice
   
   
*/
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
/**
  Allocates a m3_slice to consist of four equally sized windows
*/
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

void  m3d_transpose(m3d_t   * a)
{

   
  int x, y;
     m3d_t * c;
    m3d_create(c, a->ncols, a->nrows);
    m3_slice * b, *d;
    d = malloc(sizeof(m3_slice));
    b = malloc(sizeof(m3_slice));
    m3d_slices(b, a, 1);
    m3d_slices(d, c, 1);
    for (x = 0; x < b->nrows; x++) {
        for (y = 0; y < b->ncols; y ++) {
         m3d_transpose_vbg(b->row[x][y].rows, d->row[y][x].rows);
            
        }
    }


   
}

void m3d_copy(m3d_t * a, m3d_t const *b)
{
  m3d_create(a, b->ncols, b->nrows);
  for(int i = 0; i < a->nrows; i++)
  {
    for(int j = 0; j < b->ncols; j++)
    {
    
      a->rows[i][j] = b->rows[i][j];
 
    }
    
     a->lblock = b->lblock; //  first block pointed to in a window
     a->fcol = b->fcol;  ///column offset of first block
     a->flags = b->flags;
  
  }
  
  


}




