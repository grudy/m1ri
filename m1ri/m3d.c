
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

m3d_t * m3d_create( rci_t nrows, rci_t ncols)
{
    m3d_t * a;
    a = m1ri_malloc(sizeof(m3d_t));
    a->ncols = ncols;
    a->nrows = nrows;
    a->width =  M1RI_DN(ncols, M1RI_RADIX);
    a->block =  NULL;
    a->rows =  NULL;
    a->block = m3d_block_allocate(a->nrows,    a->width);
    a->rows  = m3d_row_alloc(a->block, a->width, a->nrows);
    a->flags = 0;
    a->lblock = ncols%64;
    a->fcol = 0;
    a->svbg = 0;
    return a;
    
}
void  m3d_rand(m3d_t * a)
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
    
   
}



m3d_t * m3d_transposewin(const m3d_t  *a )
{
    m3d_t *b = m1ri_malloc(sizeof(m3d_t));
    b = m3d_create(a->nrows, a->ncols);
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
                                                
        return b;
                                                
}

                                                                                        
static inline void m3d_set_one_ui(m3d_t * a)

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

}



static inline void m3d_set_two_ui(m3d_t * a)

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
              a->rows[l][j].sign = (leftbit)>>k;
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
              a->rows[k + i][a->width-1].sign = (leftbit)>>i;
            }
            
        }
        if ((a->ncols%64) == 0)
        {
            
            l = (a->width - 1) * 64;
            for(i = 0; i < 64; i++)
                
            {
              a->rows[l][a->width -1].units = (leftbit)>>i;
              a->rows[l][a->width -1].sign = (leftbit)>>i;
                l++;
                
            }
            
        }
        
        
    }

}


void m3d_set_ui(m3d_t * a, unsigned int scalar)
{
	unsigned int value = scalar%3;
	switch(value)
	{
		case 0: 
				m3d_set_zero;
			break ;
		case 1:
				m3d_set_one_ui(a);;
			break;
  		case 2: m3d_set_two_ui(a);
  			break;
		
  		
  	}	
  


}

m3d_t  * m3d_identity(m3d_t  *a, rci_t n)
{
    a = m3d_create( n, n);
    m3d_set_ui(a, 1);
    return a;
}

/** 
 \brief windows in M1RI_RADIX rows * M1RI_RADIX column incriments
 \param stvbg = the vbg or width offset from the base matrix
 \param strow = row offset in increments of 64
 \param sizecol  = cols * 64
 \param sizerow  = rows * 64
 */
m3d_t *    m3d_init_window(const m3d_t *c, const rci_t strow, const rci_t stvbg,const rci_t sizerows,const rci_t sizecols)
{

    m3d_t * submatrix = m1ri_malloc(sizeof(m3d_t));
    /** c->width should not be compared twice */
    if((strow + sizerows) > c->width)
    {    
        return  0;
    }
    
    if((stvbg + sizecols) > c->ncols)
    {
        return  0;
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
    
    
   
    for(  i =   f; i < (f + submatrix->nrows) ; i++)
    {
        submatrix->rows[i - f] = c->rows[i] + stvbg;   
    }
    return submatrix;
    
}

/** 
 Concat b on the end of a, the result is c 
 [a] [b] ----->  [a b]   ===  C
 */
 
 m3d_t *  m3d_concat(m3d_t * c,const   m3d_t * a, const m3d_t * b)
{
    if(a->nrows != b->nrows)
    {
    	m1ri_die("m3d_concat: bad arguments to concat!\n");
    
    }
    if(c == NULL)
    { 
      c = m3d_create( (a->ncols + b->ncols),  a->nrows );
    }
    else if(c->nrows != a->nrows || c->ncols != (a->nrows + b->nrows))
    {
    	m1ri_die("m3d_concat: c has wrong dimensions!\n");
    
    }
	
    	
	c = m3d_create( a->nrows ,a->ncols + b->ncols);
	int x, y;
	x =  0;
    
    while (x < a->nrows) 
    {
    	
    	for(y = 0; y < a->width; y++)
        {
            c->rows[x] = a->rows[x];
        }
            
        for(y = a->width; y  < c->width; y++)
        {
                
        }
        x++;
            
    }
        
    
    return c;
}

/** 
 Stacks a on b, resulting matrix is c
 [a]
 ===  C
 [b]
 */
 
m3d_t *  m3d_stack(m3d_t * c,const   m3d_t * a,const m3d_t * b)
{
    
    if(a->ncols != b->ncols)
    {
    	m1ri_die("m3d_stack: bad arguments to stack!\n");
    
    }
    if(c == NULL)
    { 
      c = m3d_create( a->ncols,  (a->nrows + b->nrows));
    }
    else if(c->nrows != a->nrows || c->ncols != (a->nrows + b->nrows))
    {
    	m1ri_die("m3d_stack: c has wrong dimensions!\n");
    
    }
    int x =  0;
    int y = 0;
    while (x < a->nrows)
    {
            c->rows[x] = a->rows[x];
            x++;
    }
    while(x < (a->nrows + b->nrows))
    {
        
        	c->rows[x] = b->rows[y];
            y++;    
            x++;
    }
        
        
        
    
    
    return c;
}



m3d_t * m3d_copy_cutoff(m3d_t  * a, m3d_t  const * b)
{
	int i, s;
    for( i = 0; i < b->nrows; i++)
    {
    	for( s = 0; s < b->width; s++)
        {
            a->rows[i][s] = b->rows[i][s];
        }          
    }
    return a;
	
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
                //printf("row [%d][%d] not equal \n", i, j);
                return 0;
            }
            
        }
    }
    return 1;
}

/** 
 Releases a m3d_t into the wilderness.
 */

void m3d_free( m3d_t *  a)
{ 		
    m1ri_free(a->rows);
    m1ri_free(a->block);
    a == NULL;
    
}




/**
   Allocates blocks for an m3d_slice
*/
static inline m3d_t  * m3_blockslice_allocate(rci_t  nrows,  wi_t  width)
{
    m3d_t * block  = m1ri_calloc(nrows *  width  ,   sizeof(m3d_t  *) );
    return block;
}



/**
   Allocates rows for an m3d_slice
   
   
*/
static inline m3d_t ** m3_rowslice_allocate(m3d_t * block,  wi_t width, rci_t nrows)
{
	
	int i;
    m3d_t ** rows = m1ri_calloc( nrows , width  * sizeof(m3d_t **));
    for ( i = 0; i <  nrows;  i++ )
    {
        rows[i]  = block + (i * width);
    };
    return rows;
}



void  m3d_slices(m3_slice *  c, const m3d_t * a, wi_t slicesize)
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
 	c->block = m3_blockslice_allocate(c->nrows,   c->width);
    c->row = m3_rowslice_allocate(c->block , c->width, c->nrows);
    r = 0 ;
     
    for ( i = 0; i <  l;  i = i + slicesize)
    {       
    	for( f = 0; f <colroundeddown ; f++)
        {
        	c->row[(r * f) + f] = m3d_init_window(a,i , (f * slicesize), slicesize, slicesize);
        }
        
        if(extracols > 0)
        {
        	c->row[(r * f) + f] = m3d_init_window(a,i , (f * slicesize), slicesize, extracols);
		}
        r++;
        
   	}
    
    if(extrarows >0 )
    {
		for( f = 0; f <colroundeddown ; f++)
        {
           c->row[(r * f) + f] = m3d_init_window(a, i , (f * slicesize), extrarows, slicesize);
        }

    	if(extracols > 0)
    	{
           c->row[(r * f) + f] = m3d_init_window(a, i , (f * slicesize), extrarows, extracols);
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

m3d_t * m3d_transpose_sliced(m3d_t * a)
{
    int x, y;
    m3d_t * c;
    c = m3d_create(a->ncols, a->nrows);
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
    return c;
  
}
/**
  Allocates a m3_slice to consist of four equally sized windows
*/
#include "m1ri_io.h"
m3_slice * m3d_quarter(const m3d_t * a)
{
	 m3_slice * c = m1ri_malloc(sizeof(m3_slice));
	 
	 
	 c->block = m3_blockslice_allocate(  2,   2);
     
     
     c->row = m1ri_malloc(4 * sizeof(m3d_t *));
     c->row[0] = m3d_init_window(a,  0, 0 , a->nrows/128, a->ncols/128);
	 c->row[1] = m3d_init_window(a, 0, a->ncols/128 , a->nrows/128, a->ncols/128);   
     c->row[2] = m3d_init_window(a, a->nrows/128, 0 , a->nrows/128, a->ncols/128);
	 c->row[3] = m3d_init_window(a, a->nrows/128,a->ncols/128,  a->nrows/128, a->ncols/128);
	 

    return c;
}

void  m3d_transpose(m3d_t   * a)
{

   
  	int x, y;
    m3d_t * c;
    c = m3d_create( a->ncols, a->nrows);
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

m3d_t *  m3d_copy(m3d_t * a, m3d_t const *b)
{
  if(a == NULL)
  {	
  	a = m3d_create( b->ncols, b->nrows);
  }
  
  if((a->ncols < b->ncols) || (a->nrows < a->nrows));
  {
  	//m1ri_die("m3d_copy: Provided return matrix has wrong dimensions.\n");
  
  }
  
  
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
  
  
  return a;

}


static inline void add_vbg(vbg *   r, vbg const *   x, vbg const * y)

{
    /* r->units = (x->units ^ y->sign) & (x->sign ^ y->units); // ///r0 ← (x0 ⊕y->1)∧(x1 ⊕y->0); */
    /* r->sign = (M1RI_ST(x->units, y->sign, x->sign ) | M1RI_ST(x->sign, y->units, y->sign)); //// r1 ← s XOR t. */



 	r->units = y->sign ^ x->units;
    r->sign = y->units ^ x->sign;
    r->sign = r->sign & r->units;
    r->units = r->units ^ x->sign;
    r->units = (y->units ^ x->units) | r->units;




}

inline void sub_m3d( vbg *r, vbg const *x, vbg const *y)
{
   /*  r->units = ((x->units^y->units) | (x->sign^y->sign)); */
   /*  r->sign = (((x->units^y->units)^x->sign)&(y->units ^ x->sign)); */
    r->sign = y->units ^ x->units;
    r->units = y->sign ^ x->sign;
    r->units = r->units | r->sign;
    r->sign = r->sign ^ y->sign;
    r->sign = (y->units ^ x->sign) & r->sign;   
   
   
   
   
}




void vbg_negation(vbg *r)
{
     r->sign = r->sign ^ r->units;
}





static inline void iadd_vbg(vbg *r,vbg  const *  x)
{
    vec t;
    t = x->units ^ r->sign;
    r->sign = x->units ^ r->units;
    r->units = x->units ^ r->units;
    r->sign = r->sign & t;
    t = t ^ x->sign;
    r->units = t | r->units;
}

static inline void isub_m3d(vbg  *r,vbg  *x)
{
    vec t;
    r->units = x->units ^ r->units;
    t  = r->units | r->sign;
    t = t ^ x->sign;
    r->sign = x->units ^ r->sign;
    r->sign = r->sign & t;
    r->units = t | r->units;     
}


void  vbg_mul( vbg *r, vbg  *x, vbg  *y)            {
    r->units = y->units ^ x->units ;
    r->sign = (y->sign ^ x->sign) & (r->units);
    
}

vbg vbg_mul_i(vbg const x, vbg const y)
{
    vbg r;
    r.units = x.units & y.units;
    r.sign  = (y.sign ^ x.sign) & (r.units);
    
    return r;   
}

m3d_t * m3d_hadamard(m3d_t * c, m3d_t const *  a, m3d_t const *  b)
{
    
    if (c == NULL)
	{
		c = m3d_create(a->nrows, b->ncols);

	} 
	else if( (c->nrows != a->nrows || c->ncols != b->ncols)) 
	{
		m1ri_die("m3d_hadamard: Provided return matrix has wrong dimensions.\n");	
	}
    if((a->nrows != b->nrows) || ( b->ncols != a->ncols))
    {
       
      m1ri_die("m3d_hadamard: Input Matrices must have same dimension.\n");
    }

    int i, j;
    if(a->ncols < 256)
    { 
    	for( i = 0; i < a->nrows; i++)
    	{
    		for(j = 0; j < (a->width ); j++)
    		{	  
        		c->rows[i][j] = vbg_mul_i(a->rows[i][j], b->rows[i][j]);
    		}  
    	}

    }
        
        
   
    
    return c;
}








m3d_t *  m3d_sub(m3d_t * r,   const  m3d_t  *x, const m3d_t  *y)
{
  int n , i;
  if (r == NULL)
	{
	r = m3d_create(x->nrows, y->ncols);

	} 
	else if( (r->nrows != x->nrows || (r->ncols != y->ncols))) 
	{
		m1ri_die("m3d_sub: Provided return matrix has wrong dimensions.\n");
    	
	
	}
    if((x->nrows != y->nrows) || ( y->ncols != x->ncols))
    {
       
      m1ri_die("m3d_sub: Input Matrices must have same dimension.\n");
    }
    
    for(i = 0; i < x->nrows; i++)
    {
    	for(n = 0; n < x->width; n++)
        {
		  sub_m3d(r->rows[i] + n, x->rows[i] + n , y->rows[i] + n );
        }
    }
   
   return r;
}




m3d_t  * m3d_add(m3d_t *c, const m3d_t  *a,const m3d_t  *b)
{
    
    if (c == NULL)
	{
	c = m3d_create(a->nrows, b->ncols);

	} 
	else if( (c->nrows != a->nrows || c->ncols != b->ncols)) 
	{
		m1ri_die("m3d_add: Provided return matrix has wrong dimensions.\n");
    	
	
	}
    if((a->nrows != b->nrows) || ( b->ncols != a->ncols))
    {
       
      m1ri_die("m3d_add: Input Matrices must have same dimension.\n");
    }
     c = m3d_create( a->nrows , b->ncols);
        int i, j;
        for( i = 0; i < a->nrows; i++)
        {
            for(j = 0; j < (a->width ); j++)
            {
            	add_vbg(&c->rows[i][j], &a->rows[i][j], &b->rows[i][j]);
        	}   
        }
    return c;  
}





m3d_t * m3d_submatrix(m3d_t *S, const m3d_t *M, const rci_t lowr, const rci_t lowc, const rci_t highr, const rci_t highc)
{
  	rci_t s_rows =  highr -  lowr ;
  	rci_t s_cols =  highc - lowc ;
  	if( S != NULL)
	{
		if((S->nrows < s_rows) || (S->ncols < s_cols))
		{
			m1ri_die("m3d_submatrix: S has too small of dimensions");
    	}
  
	}
  	else
  	{

  		S = m3d_create(s_rows, s_cols);
 	}
 	
 	vec temp_mask_l = ((leftbit >> ((lowc%M1RI_RADIX) -1)) - 1);
	vec temp_mask_r = ~temp_mask_l;
 	rci_t s_width = lowc/M1RI_RADIX;
 	if(!(lowc % M1RI_RADIX))
 	{
 		
 		if(s_cols/64 != 0)
 		{
 		
 			for(int i = 0; i < s_rows; i++)
 			{
				
						for(int m = 0; m <= s_width; m++)
						{
								S->rows[i][m].units	= M->rows[i][m].units;
 								S->rows[i][m].sign	= M->rows[i][m].sign;
 						}				
 			}
 		}
 		
 		if(s_cols%64)
 		{
 	
				   vec temp_mask = ~((leftbit >> ((s_cols%64) -1)) - 1);
 		   vbg temp;
 		   for(int i = 0; i < s_rows; i++)
 		   {
 					 
 				temp.units = M->rows[lowr + i][s_width + lowc / M1RI_RADIX].units & temp_mask;
 				temp.sign  =  M->rows[lowr + i][s_width + lowc / M1RI_RADIX].sign & temp_mask;
       			S->rows[i][s_cols / M1RI_RADIX].units = temp.units;
       			S->rows[i][s_cols / M1RI_RADIX].sign = temp.sign;

 		   } 
 		
 		}
 		
 		
 			
	}
	else 
	{
		
		vec temp_mask_l = ((leftbit >> ((lowc%M1RI_RADIX) -1)) - 1);
		vec temp_mask_r = ~temp_mask_l;
		//vbg tempr;
		vbg temp;
		rci_t j;
	
		for(int i = 0; i < s_rows; i++)
		{
 	  
			for(j=0; j+M1RI_RADIX <=s_cols; j+= M1RI_RADIX)
			{
		   		temp.units = (M->rows[s_rows + i][(lowc + j) / M1RI_RADIX].units & temp_mask_l) << (lowc%M1RI_RADIX);
 			   	temp.sign  =  (M->rows[s_rows + i][(lowc + j)/ M1RI_RADIX].sign & temp_mask_l) << (lowc%M1RI_RADIX);
       	   		S->rows[i][j].units = temp.units;
       	   		S->rows[i][j].sign = temp.sign;
       	   		S->rows[i][j].units |= (M->rows[s_rows + i][((lowc + j ) / M1RI_RADIX) + 1].units & temp_mask_r) >> (M1RI_RADIX - (lowc%M1RI_RADIX));
       	   		S->rows[i][j].sign  |= (M->rows[s_rows + i][((lowc + j) / M1RI_RADIX) + 1].sign & temp_mask_r) >> (M1RI_RADIX - (lowc%M1RI_RADIX));
				
		
			}
        	
 			S->rows[i][j/M1RI_RADIX].units &= temp_mask_l;
 			S->rows[i][j/M1RI_RADIX].sign  &= temp_mask_l;
      		S->rows[i][j/M1RI_RADIX].units |= m3d_ru_bits(M, lowr+i, lowc+j, s_cols - j) & temp_mask_r;
 		    S->rows[i][j/M1RI_RADIX].sign |= m3d_rs_bits(M, lowr+i, lowc+j, s_cols - j) & temp_mask_r;
 	   
		}	 
	}
  return S;	
}

int m3d_is_zero(const m3d_t *A)
{
	if(A == NULL)
	{
	
		m1ri_die("m3d_is_zero: A cannot be null!\n");

	
	}
	int i, j;

	for(i = 0; i < (A->width ); i++)
   	{
        for(j = 0; j < (A->width ); j++)
        {
            if((A->rows[i][j].units != 0) || (A->rows[i][j].sign != 0));
            return 0;
            
        }   
    
    }
  
  
  
  return 1;

}

m3d_t *m3d_mul_scalar(m3d_t *C, const long a, const m3d_t *B);


/*
  Matrix Multiplication Lookup tables

*/
void *  m3d_combine3(vbg *table, vbg *input )
{
    vbg t, a, b, c;
    t.sign = t.units = 0;
    a = input[0];
    b = input[1];
    c = input[2];
    table[0] = t;
    table[1] = a;
    table[2] = b;
    table[4] = c;

    add_vbg(&t, &a, &b);
    table[3] = t;
    iadd_vbg(&t, &c);
    table[7] = t;
    isub_m3d(&t, &a);
    table[6] = t;
    add_vbg((table + 5), &a , &b);

    return 0;
    
}


void m3d_combine4(vbg *table, vbg *input )
{
    vbg t, a, b, c , d;
    t.sign = t.units = 0;
    a = input[0];
    b = input[1];
    c = input[2];
    d = input[3];
    
    table[0] = t;
    table[1] = a;
    table[2] = b;
    table[4] = c;
    table[8] = d;
    
    add_vbg(&t, &c, &d);
    
    table[12] = t;
    
    add_vbg(&t,&b,&c);
    table[6] = t;
    iadd_vbg(&t,&d);
    table[14] = t;
    isub_m3d(&t,&c);
    table[10] = t;
    
    add_vbg(&t,&b,&c);
    table[3] = t;
    iadd_vbg(&t, &d);
    
    
    
    table[11] = t;
    iadd_vbg(&t, &c);
    table[15] = t;
    isub_m3d(&t, &d);
    table[7] = t;
    isub_m3d(&t, &b);
    table[5] = t;
    iadd_vbg(&t, &d);
    table[13] = t;
    isub_m3d(&t, &c);
    table[9] = t;
    
    
}


void m3d_combine5(vbg *table, vbg *input )
{
	int i;
    vbg e, *t4;

    m3d_combine4(table, input);
    e = input[4];
    t4 = table + 16;
    table[16] = e;
    
    for ( i = 1; i < 16 ; i ++ ) {
        add_vbg(t4 + i, table + i, &e);
    }
}


void m3d_combine6(vbg *table, vbg *input )

{
    vbg f, *t5;
    int i;
    m3d_combine5(table, input);
    f = input[5];
    t5 = (table + 32);
    table [32] = f;
    for (i = 1; i < 32; i++)
        add_vbg((t5 + i), (table + i), &f);
    
}

 void m3d_combine7(vbg *table, vbg *input )

{
    
    vbg g, *t6;
    int i;
    m3d_combine6(table, input);
    g = input[6];
    t6 = (table+64);
    table[64] = g;
    for (i = 1; i < 64; i = i +1) {
        add_vbg((t6 + i), (table + i), &g );
    }
 
}


void m3d_combine8(vbg *table, vbg *input)

{
    vbg h, *t7;
    int i;
    
    m3d_combine7(table, input);
    h = input[7];
    t7 = (table+128);
    table[128] = h;
    
    for (i = 1; i < 128; i++)
        add_vbg((t7 + i), (table+i), &h);
}


void m3d_mul_64(vbg *R, vbg *  A, vbg *  B)
{
    int i;
    vbg t1, t2, r1, r2, a;
    vec v1, v2;
    
    vbg  tables6[9][64];
    vbg tables5[2][32];
    
	for (i = 0; i < 9; i ++)
    {
        m3d_combine6(tables6[i], B + (6*i));
    }
   
    for (i = 0; i < 2; i ++)
    {
        m3d_combine5(tables5[i], B + (54 + (5 * i)));
    }

    for (i = 0; i < 64; i ++  )/* i from 0 <= i < 64 */
    {
        a = A[i];
        v2 = a.sign;
    
        v1 = (a.units ^ v2);		
        r1 = tables6[0][v1&63];
        v1 >>= 6;
        r2 = tables6[0][v2&63];
        v2 >>= 6;
        t1 = tables6[1][v1&63]; iadd_vbg(&r1, &t1);v1 >>= 6;
        t2 = tables6[1][v2&63]; iadd_vbg(&r2, &t2); v2 >>= 6;
        t1 = tables6[2][v1&63]; iadd_vbg(&r1, &t1); v1 >>= 6;
        t2 = tables6[2][v2&63]; iadd_vbg(&r2, &t2); v2 >>= 6;
        t1 = tables6[3][v1&63]; iadd_vbg(&r1, &t1); v1 >>= 6;
        t2 = tables6[3][v2&63]; iadd_vbg(&r2, &t2); v2 >>= 6;
        t1 = tables6[4][v1&63]; iadd_vbg(&r1, &t1); v1 >>= 6;
        t2 = tables6[4][v2&63]; iadd_vbg(&r2, &t2); v2 >>= 6;
        t1 = tables6[5][v1&63]; iadd_vbg(&r1, &t1); v1 >>= 6;
        t2 = tables6[5][v2&63]; iadd_vbg(&r2, &t2); v2 >>= 6;
        t1 = tables6[6][v1&63]; iadd_vbg(&r1, &t1); v1 >>= 6;
        t2 = tables6[6][v2&63]; iadd_vbg(&r2, &t2); v2 >>= 6;
        t1 = tables6[7][v1&63]; iadd_vbg(&r1, &t1); v1 >>= 6;
        t2 = tables6[7][v2&63]; iadd_vbg(&r2, &t2); v2 >>= 6;
        t1 = tables6[8][v1&63]; iadd_vbg(&r1, &t1); v1 >>= 6;
        t2 = tables6[8][v2&63]; iadd_vbg(&r2, &t2); v2 >>= 6;
        t1 = tables5[0][v1&31]; iadd_vbg(&r1, &t1); v1 >>= 5;
        t2 = tables5[0][v2&31]; iadd_vbg(&r2, &t2); v2 >>= 5;
        t1 = tables5[1][v1&31]; iadd_vbg(&r1, &t1);
        t2 = tables5[1][v2&31]; iadd_vbg(&r2, &t2);
        
        isub_m3d(&r1, &r2);
        R[i] = r1;
       /*  */
    }
    
}

/* 32 * 64,2048 bit, 256 byte matrix(slice) multiplication */
void mul_32_m3d(vbg *R, vbg *A, vbg *B)
{
    long i;
    vbg t1, t2, r1, r2, a;
    long v1, v2;
    
    vbg tables5[4][32];
    vbg tables4[3][16];
    for (i = 1; i < 4; i ++)
        
        m3d_combine5(tables5[i], B + 0 + 5*i);
    for (i = 0; i < 3; i++)
        m3d_combine4(tables4[i], B + 20 + 4*i);
    
    for (i = 0;i < 32; i++)
    {
        
        a = A[i];
        v2 = a.sign;
        v1 = a.units ^ v2;
        t1 = tables5[0][v1&31]; v1 >>= 5;
        t2 = tables5[0][v2&31]; v2 >>= 5;
        t1 = tables5[1][v1&31]; iadd_vbg(&r1, &t1); v1 >>= 5;
        t2 = tables5[1][v2&31]; iadd_vbg(&r2, &t2); v2 >>= 5;
        t1 = tables5[2][v1&31]; iadd_vbg(&r1, &t1); v1 >>= 5;
        t2 = tables5[2][v2&31]; iadd_vbg(&r2, &t2); v2 >>= 5;
        t1 = tables5[3][v1&31]; iadd_vbg(&r1, &t1); v1 >>= 5;
        t2 = tables5[3][v2&31]; iadd_vbg(&r2, &t2); v2 >>= 5;
        t1 = tables4[0][v1&15]; iadd_vbg(&r1, &t1); v1 >>= 4;
        t2 = tables4[0][v2&15]; iadd_vbg(&r2, &t2); v2 >>= 4;
        t1 = tables4[1][v1&15]; iadd_vbg(&r1, &t1); v1 >>= 4;
        t2 = tables4[1][v2&15]; iadd_vbg(&r2, &t2); v2 >>= 4;
        t1 = tables4[2][v1&15]; iadd_vbg(&r1, &t1);
        t2 = tables4[2][v2&15]; iadd_vbg(&r2, &t2);
        
        isub_m3d(&r1, &r2);
        R[i] = r1;
    }
    
}

/* 16 * 64,1024 bit, 128 byte matrix(slice) multiplication */
void mul_16_m3d(vbg *R, vbg *A, vbg *B)
{
    long i;
    vbg t1, t2, r1, r2, a;
    long v1, v2;
    
    vbg tables4[4][16];
    for (i = 0; i < 4; i++)
        m3d_combine4(tables4[i], B + (4*i));
    for (i = 0;  i < 16; i++)
    {
        a = A[i];
        v2 = a.sign;
        v1 = a.units ^ v2;
        r1 = tables4[0][v1&15]; v1 >>= 4;
        r2 = tables4[0][v2&15]; v2 >>= 4;
        t1 = tables4[1][v1&15]; iadd_vbg(&r1, &t1); v1 >>= 4;
        t2 = tables4[1][v2&15]; iadd_vbg(&r2, &t2); v2 >>= 4;
        t1 = tables4[2][v1&15]; iadd_vbg(&r1, &t1); v1 >>= 4;
        t2 = tables4[2][v2&15]; iadd_vbg(&r2, &t2); v2 >>= 4;
        t1 = tables4[3][v1&15]; iadd_vbg(&r1, &t1);
        t2 = tables4[3][v2&15]; iadd_vbg(&r2, &t2);
    
        isub_m3d(&r1, &r2);
        R[i] = r1;
    }
}

/* 8 * 64,512 bit, m1ri_word byte matrix(slice) multiplication */
void mul_8_m3d(vbg *R, vbg *A, vbg *B)

{
    int i;
    vbg t1, t2, r1, r2, a;
    vec v1, v2;
    
    vbg tables4[2][16];
    for (i = 0; i < 2; i++)
        m3d_combine4(tables4[i], B + (4*i));
    for (i = 0; i < 8; i++)
    {
        a = A[i];
    v2 = a.sign;
    v1 = a.units ^ v2;
    r1 = tables4[0][v1&15]; v1 >>= 4;
    r2 = tables4[0][v2&15]; v2 >>= 4;
    t1 = tables4[1][v1&15]; iadd_vbg(&r1, &t1);
    t2 = tables4[1][v2&15]; iadd_vbg(&r2, &t2);
    
    isub_m3d(&r1, &r2);
    R[i] = r1;
    }
}






/* 4 * 64,256 bit, 32 byte matrix(slice) multiplication */
void mul_4_m3d(vbg *R, vbg *A, vbg *B)
{
    int i;
    vbg r1, r2, a;
    vec v1, v2;
    
    vbg table4[16];
    for (i = 0; i < 1; i++)
        m3d_combine4(table4, B + (4*i));
    for(i = 0; i < 4; i++)
    {
        a = A[i];
        v2 = a.sign;
        v1 = a.units ^ v2;
        r1 = table4[v1&15];
        r2 = table4[v2&15];
        
        isub_m3d(&r1, &r2);
        R[i] = r1;
    }
    
}



inline void m3d_mul_zero(m3d_t * a)
{
	
	int i, j;
	for(i = 0; i < a->nrows; i++)
	{
		for(j = 0; j < a->ncols; j++)
		{
		  a->rows[i][j].units = 0;
		  a->rows[i][j].sign  = 0;
		  
		
		}
	
	}

}


inline void m3d_mul_two(m3d_t * a, const m3d_t * b)
{
	int i, j;
	for(i = 0; i < a->nrows; i++)
	{
		for(j = 0; j < a->ncols; j++)
		{
		  a->rows[i][j].units = ~(b->rows[i][j].units) && b->rows[i][j].sign ;
		
		  
		
		}
	
	}


}


m3d_t *m3d_mul_scalar(m3d_t *C, const long a, const m3d_t *B)
{
	

    if(C == NULL)
    { 
      C = m3d_create( B->ncols, B->nrows);
    }
    
    else if(C->nrows != B->nrows || C->ncols !=  B->ncols)
    {
    	m1ri_die("m3d_mul_scalar: C has wrong dimensions!\n");
    
    }
    
	long m = a%3;
	switch(m)
	{
		case 0: m3d_mul_zero(C);
			break ;
  		case 2: m3d_mul_two(C, B);
			break ;
  	}
  
  return C;
}




void m3d_add_row(m3d_t *A, rci_t ar, const m3d_t *B, rci_t br, rci_t start_col)
{
  


}





