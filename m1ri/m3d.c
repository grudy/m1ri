
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
    bits  = (n == 0) ? (~(rightbit << spill) &  (M->rows[x][block].sign))  : ((rightbit << spill) | (M->rows[x][block].units));
    return bits;
}


vec m3d_ru_bits(m3d_t const *M, rci_t  x, rci_t  y, int  n) 
{
	wi_t  block = (y  ) / M1RI_RADIX;
	int  spill = (y  % M1RI_RADIX) + n - M1RI_RADIX;
    vec bits;
    bits  = (n == 0) ? (~(rightbit << spill) &  (M->rows[x][block].units))  : ((rightbit << spill) | (M->rows[x][block].units));
    return bits;   
}

vbg m3d_read_elems(m3d_t const *M, rci_t  x, rci_t  y, int  n) 
{
    wi_t  block = (y  ) / M1RI_RADIX;
    int  spill = (y  % M1RI_RADIX) + n - M1RI_RADIX;
    vbg elem;
    if(spill <= 0)
    {
    
    	elem.sign =  M->rows[x][block ].sign <<  - spill; 
    	elem.units =  M->rows[x][block].units << -spill  ;

    }
    else
    {
    	elem.sign = (M->rows[x][block + 1].sign << (M1RI_RADIX - spill)) | (M->rows[x][block].sign>> spill);
    	elem.units = (M->rows[x][block + 1].units << (M1RI_RADIX - spill)) | (M->rows[x][block].units>> spill);
    }
    
    
    
    elem.sign = (elem.sign >> (M1RI_RADIX - n));
    
    elem.units = (elem.units >> (M1RI_RADIX - n));

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
         a_place =  rightbit << dif_a ;
         b_place =  rightbit <<  dif_b ;
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
         a_place =  rightbit <<  dif_a ;
         b_place =  rightbit <<  dif_b ;
       

              
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


void  m3d_write_elem( m3d_t * M,rci_t x, rci_t y, vec s, vec u )
{
    wi_t  block = (y  ) / M1RI_RADIX;
    int   spill =  (y  % M1RI_RADIX) ;
    M->rows[x][block].units  = (u == 0) ? (~(rightbit <<  spill) &  (M->rows[x][block].units))  : ((rightbit <<  spill) | (M->rows[x][block].units));
    M->rows[x][block].sign  = (s == 0) ? (~(rightbit << spill) &  (M->rows[x][block].sign))  : ((rightbit << spill) | (M->rows[x][block].sign));
    
}

void m3d_add_i(m3d_t * x, m3d_t *y) 
{
	 int i, j;
        for( i = 0; i < x->nrows; i++)
        {
            for(j = 0; j < (x->width ); j++)
            {
            	m3d_inc(&x->rows[i][j], &y->rows[i][j]);
        	}   
        }
	

}

void m3d_sub_i(m3d_t * x, m3d_t *y) 
{
	 int i, j;
        for( i = 0; i < x->nrows; i++)
        {
            for(j = 0; j < (x->width ); j++)
            {
            	m3d_dec(&x->rows[i][j], &y->rows[i][j]);
        	}   
        }
	

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
    a->flags = notwindowed;
    a->lblock = ncols%64;
    a->fcol = 0;
    a->svbg = 0;
    return a;
    
}
void  m3d_rand(m3d_t * a)
{
    int i,  z;
    rci_t cutoff = a->ncols% 64;
    if(cutoff)
    {
    
    
    	vec mask_rand = (rightbit  << (cutoff )) - 1;
    	//mask_rand = ~mask_rand;
    	
    	for(i = 0; i < (a->nrows); i++)
   	 	{
        	for( z = 0; z  < (a->width - 1 ); z++)
            {  
       			a->rows[i][z].sign = m1ri_rand();
       			a->rows[i][z].units = m1ri_rand();
       			a->rows[i][z].sign =  a->rows[i][z].sign & a->rows[i][z].units;
            
            
            }
            	
    			a->rows[i][z].sign = m1ri_rand();
       			a->rows[i][z].units = m1ri_rand();
       			a->rows[i][z].sign =  a->rows[i][z].sign & a->rows[i][z].units;
       			a->rows[i][z].sign = a->rows[i][z].sign & mask_rand;
       			a->rows[i][z].units = a->rows[i][z].units & mask_rand;
            
   	 	}
    
    
    
    }
    
    else
    {
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
                
              a->rows[l][j].units = (rightbit)<<k;
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
                
              a->rows[k + i][a->width-1].units = (rightbit)<<i;
            }
            
        }
        if ((a->ncols%64) == 0)
        {
            
            l = (a->width - 1) * 64;
            for(i = 0; i < 64; i++)
                
            {
              a->rows[l][a->width -1].units = (rightbit)<<i;
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
                
              a->rows[l][j].units = (rightbit )<<k;
              a->rows[l][j].sign = (rightbit)<<k;
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
                
              a->rows[k + i][a->width-1].units = (rightbit)<<i;
              a->rows[k + i][a->width-1].sign = (rightbit)<<i;
            }
            
        }
        if ((a->ncols%64) == 0)
        {
            
            l = (a->width - 1) * 64;
            for(i = 0; i < 64; i++)
                
            {
              a->rows[l][a->width -1].units = (rightbit)<<i;
              a->rows[l][a->width -1].sign = (rightbit)<<i;
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

 m3d_t *    m3d_init_window_unshackled(const m3d_t *c, const rci_t strow, const rci_t stvbg,const rci_t sizerows,const rci_t sizecols)
{

	

	m3d_t * submatrix = m1ri_malloc(sizeof(m3d_t));
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
   // if(c == NULL)
    //{ 
      //c = m3d_create( (a->ncols + b->ncols),  a->nrows );
  //  }
    //else if(c->nrows != a->nrows || c->ncols != (a->nrows + b->nrows))
    //{
    //	m1ri_die("m3d_concat: c has wrong dimensions!\n");
    
    //}
	
    	
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
    for( i = 0; i < a->nrows; i++)
    {
    	for( s = 0; s < a->width; s++)
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
    if(a->flags == notwindowed)
    {
    
		m1ri_free(a->block); 	
	}
    
    m1ri_free(a);
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
    m3d_t ** rows = m1ri_calloc( nrows , sizeof(m3d_t **));
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


/**
  Allocates a m3_slice to consist of four equally sized windows
*/
m3_slice * m3d_quarter( m3d_t  const * a)
{
	 m3_slice * c = m1ri_malloc(sizeof(m3_slice));
	 
	 
	
     c->row = m1ri_calloc( 4 , sizeof(m3d_t **));

     
     
     c->row[0] = m3d_init_window(a,  0, 0 , a->nrows/128, a->ncols/128);
	 c->row[1] = m3d_init_window(a, 0, a->ncols/128 , a->nrows/128, a->ncols/128);   
     c->row[2] = m3d_init_window(a, a->nrows/128, 0 , a->nrows/128, a->ncols/128);
	 c->row[3] = m3d_init_window(a, a->nrows/128,a->ncols/128,  a->nrows/128, a->ncols/128);
	 
	 
	 

    return c;
}





m3d_t *  m3d_copy(m3d_t * a, m3d_t const *b)
{
  if(a == NULL)
  {	
  	a = m3d_create( b->nrows, b->ncols);
  }
  
  if((a->ncols < b->ncols) || (a->nrows < b->nrows))
  {
  	m1ri_die("m3d_copy: Provided return matrix has wrong dimensions.\n");
  
  }
  
  
  for(int i = 0; i < b->nrows; i++)
  {
    for(int j = 0; j < b->width; j++)
    {
    
		a->rows[i][j].units = b->rows[i][j].units;
		a->rows[i][j].sign = b->rows[i][j].sign;
 
    }
    
  
  }
  
  
     a->lblock = b->lblock; //  first block pointed to in a window
     a->fcol = b->fcol;  ///column offset of first block
     a->flags = b->flags;
  
  
  return a;

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
    //if(a->ncols < 256)
    //{ 
    
    for( i = 0; i < a->nrows; i++)
    {
    	for(j = 0; j < (a->width ); j++)
    	{	  
    		c->rows[i][j] = vbg_mul_i(a->rows[i][j], b->rows[i][j]);
    	}  
    }

	//}
        
        
   
    
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







m3d_t  * m3d_add(m3d_t *c,  m3d_t const *a, m3d_t const *b)
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



/*

	\brief incremental subtraction, but where the subtrahend is changed
	
*/


 

static inline void sub_m3d_r(vbg  const *r,vbg   *x)
{
	vbg m;
	m.units = x->units;
	m.sign =  x->sign;
    vec t;
    x->units = m.units ^ r->units;
     t  = x->units | r->sign;
    t = t ^ m.sign;
    x->sign = m.units ^ r->sign;
    x->sign = x->sign & t;
    x->units = t | x->units;  
         
}


void m3d_sub_r(m3d_t   *x , m3d_t   const *r)
{
  
	int n , i;
	for(i = 0; i < x->nrows; i++)
    {
    	for(n = 0; n < x->width; n++)
        {
        
       
		  sub_m3d_r(r->rows[i] + n,  x->rows[i] + n );
        }
    }
   



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
 	
 	vec temp_mask_l = ((rightbit << ((lowc%M1RI_RADIX) -1)) - 1);
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
 	
				   vec temp_mask = ~((rightbit << ((s_cols%64) -1)) - 1);
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
		
		vec temp_mask_l = ((rightbit <<((lowc%M1RI_RADIX) -1)) - 1);
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





static inline void m3d_combine4(vbg *table, vbg  ** const input) 
{
    vbg  t,   a,  b,  c,  d;
    t.sign = t.units = 0;
    a = *input[0];
    b = *input[1];
    c = *input[2];
    d = *input[3];

    table[0] = t;
    table[1] = a;
    table[2] = b;
    table[4] = c;
    table[8] = d;

    add_vbg(&t,&c,&d);
    table[12] = t;
    
    add_vbg(&t,&b,&c);
    table[6] = t;
    m3d_inc(&t,&d);
    table[14] = t;
    m3d_dec(&t,&c);
    table[10] = t;

    add_vbg(&t,&a,&b);
    table[3] = t;
    m3d_inc(&t,&d);
    
    table[11] = t;
    m3d_inc(&t,&c);
    table[15] = t;
    m3d_dec(&t,&d);
    
    table[7] = t;
    m3d_dec(&t,&b);
    table[5] = t;
    m3d_inc(&t,&d);
    table[13] = t;
    m3d_dec(&t,&c);
    table[9] = t;
}

static inline void m3d_combine5(vbg *table, vbg  **  const input) 
{
    vbg *e, *t4;
    int i;

    m3d_combine4(table, input);
    e = input[4];
    t4 = table+16;
    table[16] = *e;
    
    for(i=1;i<16;i++)
        add_vbg(t4 + i, table+i, e);
}
    
static inline void m3d_combine6(vbg *table, vbg  ** const input) 
{
    vbg *e, *t4;
    vbg * f, *t5;
    int i;
    
    m3d_combine4(table, input);
    e = input[4];
    t4 = table+16;
    table[16] = *e;

    f = input[5];
    t5 = table+32;
    table[32] = *f;
    
    for(i=1;i<16;i++)
        add_vbg(t4 + i, table+i, e);

    for(i=1;i<32;i++)
        add_vbg(t5 + i, table+i, f);

}


void m3d_mul_64(vbg **R, vbg  ** const A, vbg   ** const B)
{
   
    vbg t, r, * a;
    vec v;
    int i;

    vbg tables6[4][64];
    vbg tables5[8][32];

    for(i=0;i<4;i++)
    {
        m3d_combine6(tables6[i], B + 6*i);
    }
    for(i=0;i<8;i++)
    {
        m3d_combine5(tables5[i], B + 24 + (5*i));
	}
    for(i=0;i<64;i++) 
    {
    
        a = A[i];
        v = a->sign;// finds values equal to one
        r = tables6[0][v&63];                 v >>= 6;
        t = tables6[1][v&63]; m3d_inc(&r, &t); v >>= 6;
        t = tables6[2][v&63]; m3d_inc(&r, &t); v >>= 6;
        t = tables6[3][v&63]; m3d_inc(&r, &t); v >>= 6;
        
        t = tables5[0][v&31]; m3d_inc(&r, &t); v >>= 5;
        t = tables5[1][v&31]; m3d_inc(&r, &t); v >>= 5;
        t = tables5[2][v&31]; m3d_inc(&r, &t); v >>= 5;
        t = tables5[3][v&31]; m3d_inc(&r, &t); v >>= 5;
        t = tables5[4][v&31]; m3d_inc(&r, &t); v >>= 5;
        t = tables5[5][v&31]; m3d_inc(&r, &t); v >>= 5;
        t = tables5[6][v&31]; m3d_inc(&r, &t); v >>= 5;
        t = tables5[7][v&31]; m3d_inc(&r, &t);
		// The ones and twos are added above
		r.sign^=r.units;  //this multiplies r by two
		v = a->units ^ a->sign;	//this finds values that are equal to one

        t = tables6[0][v&63]; m3d_inc(&r, &t); v >>= 6;
        t = tables6[1][v&63]; m3d_inc(&r, &t); v >>= 6;
        t = tables6[2][v&63]; m3d_inc(&r, &t); v >>= 6;
        t = tables6[3][v&63]; m3d_inc(&r, &t); v >>= 6;
       
        t = tables5[0][v&31]; m3d_inc(&r, &t); v >>= 5;
        t = tables5[1][v&31]; m3d_inc(&r, &t); v >>= 5;
        t = tables5[2][v&31]; m3d_inc(&r, &t); v >>= 5;
        t = tables5[3][v&31]; m3d_inc(&r, &t); v >>= 5;
        t = tables5[4][v&31]; m3d_inc(&r, &t); v >>= 5;
        t = tables5[5][v&31]; m3d_inc(&r, &t); v >>= 5;
        t = tables5[6][v&31]; m3d_inc(&r, &t); v >>= 5;
    	t = tables5[7][v&31]; m3d_inc(&r, &t);
    	
    	
    	
    	
    	
    	
		R[i][0] = r;
        
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





