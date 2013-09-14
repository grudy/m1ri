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
 
 m7d.c
 */

#include "m7d.h"


vec m7d_rm_bits(m7d_t *M, rci_t  x, rci_t  y, int  n) {
    
    
    wi_t  block = (y  ) / M1RI_RADIX;
    int  spill = (y  % M1RI_RADIX) + n - M1RI_RADIX;
    vec bits;
    bits = (spill <= 0) ? M->rows[x][block].middle << -spill : (M->rows[x][block + 1].sign << (M1RI_RADIX - spill)) | (M->rows[x][block].middle >> spill);
    return bits;
    
    
}


vec m7d_rs_bits(m7d_t *M, rci_t  x, rci_t  y, int  n) {
    
    
    
    
    
    wi_t  block = (y  ) / M1RI_RADIX;
    int  spill = (y  % M1RI_RADIX) + n - M1RI_RADIX;
    vec bits;
    bits = (spill <= 0) ? M->rows[x][block].sign << -spill : (M->rows[x][block + 1].sign << (M1RI_RADIX - spill)) | (M->rows[x][block].sign >> spill);
    return bits;
    
    
}
/*
 Read n bits from units
 x = rows
 y = columns
 M = Matrix read from
 */

vec m7d_ru_bits(m7d_t *M, rci_t  x, rci_t  y, int  n) {

    wi_t  block = (y  ) / M1RI_RADIX;
    int  spill = (y  % M1RI_RADIX) + n - M1RI_RADIX;
    vec bits;
    bits = (spill <= 0) ? M->rows[x][block].units << -spill : (M->rows[x][block + 1].units<< (M1RI_RADIX - spill)) | (M->rows[x][block].units>> spill);
    return bits;
    
    
}





/*
 Read n elements
 x = rows
 y = columns
 M = Matrix read from
 */

vtri m7d_read_elems(m7d_t *M, rci_t  x, rci_t  y, int  n) 
{
    wi_t  block = (y  ) / M1RI_RADIX;
    int  spill = (y  % M1RI_RADIX) + n - M1RI_RADIX;
    vtri elem;
    elem.units = (spill <= 0) ? M->rows[x][block].units << -spill : ((M->rows[x][(block + 1)].units<< (M1RI_RADIX - spill)) | (M->rows[x][block].units >> spill));
    elem.sign = (spill <= 0) ?  (M->rows[x][block].sign << -spill) : (M->rows[x][block + 1].sign << (M1RI_RADIX - spill)) | (M->rows[x][block].sign>> spill);
    elem.middle = (spill <= 0) ?  (M->rows[x][block].middle << -spill) : (M->rows[x][block + 1].middle << (M1RI_RADIX - spill)) | (M->rows[x][block].middle>> spill);
    elem.middle = (elem.middle >> (M1RI_RADIX - n));
    elem.units = (elem.units >> (M1RI_RADIX - n));
    elem.sign = (elem.sign >> (M1RI_RADIX - n));
    
     
    return elem;
    
    
}





/*
 Swap rows in a matrix;
 */
void m7d_rowswap (m7d_t * M, rci_t row_a, rci_t  row_b)
{
    
    
    if((M->nrows >= (row_a ) && (M->nrows >= row_b)))
    {
        vtri  temp ;
        temp =  *M->rows[row_a -1];
        M->rows[row_a -1] = M->rows[row_b -1];
        M->rows[row_b -1] =  &temp;
        
    }
    
    else
    {
        
        return;
    }
    
}


/*
 
 */

m7d_t    m7d_identity_set(m7d_t * a)

{
    if(a->ncols == a->nrows)
    {
        int k,i,  j,l;
        for( i  = 1; i < (a->width  ) ; ++i)
        {
            l =  ((i - 1) * M1RI_RADIX);
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
        if ((a->ncols%M1RI_RADIX) == 0)
        {
            
            l = (a->width - 1) * M1RI_RADIX;
            for(i  = 0; i < M1RI_RADIX; i++)
                
            {
                a->rows[l][a->width -1].units = (leftbit)>>i;
                l++;
                
            }
            
        }
        
        
    }
    return *a;
}




m7d_t   m7d_identity(m7d_t  *a, rci_t n)
{
    *a = m7d_create(a, n, n);
    *a = m7d_identity_set(a);
    return *a;
}



//unfinished
void *  m7d_write_elem( m7d_t * M,rci_t x, rci_t y, vec s, vec m,  vec u )
{
    wi_t  block = (y  ) / M1RI_RADIX;
    int   spill =  (y  % M1RI_RADIX) ;
    M->rows[x][block].units  = (u == 0) ? (~(leftbit >> spill) &  (M->rows[x][block].units))  : ((leftbit >> spill) | (M->rows[x][block].units));
    M->rows[x][block].middle  = (m == 0) ? (~(leftbit  >> spill) &  (M->rows[x][block].middle))  : ((leftbit  >> spill) | (M->rows[x][block].middle));
    M->rows[x][block].sign  = (s == 0) ? (~(leftbit  >> spill) &  (M->rows[x][block].sign))  : ((leftbit  >> spill) | (M->rows[x][block].sign));
    return 0;

}


/*
 
 */



vtri  * m7d_block_allocate(vtri * block, rci_t  nrows,  wi_t  width)
{
    
    block  = m1ri_malloc(nrows * width * sizeof(vtri) );
    return block;
}

/*
 
 */




vtri ** m7d_row_alloc(vtri * block, vtri ** rows, wi_t width, rci_t nrows)
{
	int i;
    rows = m1ri_malloc( nrows * width * sizeof(vtri *));
    for ( i = 0; i <  nrows;  i++ )
    {
        rows[i]  = (block + (i * width));
    };
    return rows;
    
}

/*
 
 */

m7d_t m7d_create( m7d_t * a, rci_t nrows, rci_t ncols)
{
    a->ncols = ncols;
    a->nrows = nrows;
    a->width = M1RI_DN(ncols, M1RI_RADIX);
    a->block = m7d_block_allocate(a->block,  a->nrows,    a->width);
    a->rows  = m7d_row_alloc(a->block, a->rows, a->width, a->nrows);
    a->flags = iswindowed;
    return *a;
    
}


m7d_t   m7d_window(m7d_t *c, rci_t strow, rci_t stvtri, rci_t sizerows, rci_t sizecols)
{
    
    
    m7d_t  submatrix;
    
    if((strow + sizerows) > c->width)
    {
        
        
        return submatrix;
    }
    
    
    
    
    if((stvtri + sizecols) > c->width)
    {
        
        
        return submatrix;
    }
    
    
    
    submatrix.nrows =   M1RI_RADIX * sizecols;
    submatrix.ncols =  M1RI_RADIX * sizecols;
    submatrix.flags = iswindowed;
    submatrix.width =  sizecols;
    submatrix.block = &c->block[(stvtri * stvtri)];
    submatrix.rows = m1ri_calloc(M1RI_RADIX * sizerows ,  sizecols * sizeof(vtri *));
    submatrix.lblock = ( (sizerows +  strow)  ==  c->width)? c->lblock:  0;
    submatrix.fcol   = 0;
    submatrix.svtri = stvtri;
    
    
    
    int i;
    
    for(  i =   strow; i < (strow + (M1RI_RADIX * sizerows)) ; i++)
    {
        
        submatrix.rows[i - strow] = c->rows[i];
        
    }
    
    return submatrix;
    
}

void   m7d_window_create(m7d_t *c, m7d_t * submatrix, rci_t strow, rci_t stvtri, rci_t sizerows, rci_t sizecols)
{
    
    
    submatrix->nrows =   M1RI_RADIX * sizecols;
    submatrix->ncols =  M1RI_RADIX * sizecols;
    submatrix->flags = iswindowed;
    submatrix->width = 1 * sizecols;
    submatrix->block = &c->block[(stvtri * stvtri)];
    submatrix->rows = m1ri_calloc(M1RI_RADIX * sizerows ,  sizecols * sizeof(vtri *));
    submatrix->lblock = ( (sizerows +  strow)  ==  c->width)? c->lblock:  0;
    submatrix->fcol   = 0;
    submatrix->svtri = stvtri;
    
    
    int i;
    
    for(  i =   strow; i < strow + (M1RI_RADIX * sizerows) ; i++)
    {
        
        submatrix->rows[i - strow] = c->rows[i];
        
    }
    
    
    
}

/*
 
 */
vtri * m7d_allone(m7d_t * a)
{
    int i;
    for( i = 0; i < (a->nrows * a->width); i++)
    {
        
        a->block[i].sign = 0xffffffffffffffff;
        
        a->block[i].middle = 0xffffffffffffffff;
        
        a->block[i].units = 0xffffffffffffffff;
        
    
        
    }
    return a->block;
}

void vtri_negate(vtri * a)
{
	a->units = !a->units;
	a->middle = !a->middle;
	a->sign = !a->sign;


}

m7d_t  m7d_rand(m7d_t * a)
{
    int i;
    for( i = 0; i < (a->nrows * a->width); i++)
    {
        
        a->block[i].sign = m1ri_rand();
        a->block[i].middle = m1ri_rand() ;
        a->block[i].units = m1ri_rand();
   
    }
    return *a;
}



/*
 
 Releases a m7d_t into the wilderness.
 */



void m7d_free( m7d_t *  tofree)
{
    m1ri_free(tofree->rows);
    m1ri_free(tofree->block);
   
}


void add_vtri(vtri * r, vtri * x, vtri * y)

{
    
    vec s;
    vec t;
    
    r->sign = x->units ^ y->units;

    s = (x->units & y->units);
    r->middle = s^ x->middle ^ y->middle;
    t = (((s) & (x->middle | y->middle)) | (x->middle & y->middle) );
    r->sign = x->sign ^ y->sign ^ t;
    s = x->sign | y->sign | t;
    
    t = (r->units & s );
    r->units = s ^ r->units;
    
   r->middle = r->middle ^ t ;
    
   r->sign = r->sign ^ ( t & r->middle);
    
}

void m7d_sub_i(vtri  *r, vtri *y)
{
	/*
	Subtraction function
	*/

}


void iadd_vtri(vtri  *x, vtri *y)
{
    vtri  r;
    vec s;
    vec t;

    r.sign = x->units ^ y->units;
    s = (x->units & y->units);
    r.middle = s^ x->middle ^ y->middle;
    t = (((s) & (x->middle | y->middle)) | (x->middle & y->middle) );
    r.sign = x->sign ^ y->sign ^ t;
    s = x->sign | y->sign | t;
    
    t = (r.units & s );
    x->units  = s ^ r.units;
    x->middle  = x->middle ^ t ;
    x->sign  = x->sign ^ (  t & x->middle);
    /*Cleanup*/
      

}


void reduce_vtri( vtri * a)
{
    vtri b = *a ;
    a->units  = b.units ^ (b.units  | b.sign | b.middle) ;//  )
    a->middle  = b.middle ^ (b.units  | b.sign | b.middle) ;
    a->sign  = b.units ^ (b.units  | b.sign | b.middle) ;
}





vtri m7d_mul_2(vtri a)
{
    vec temp;
    temp = a.units;
    a.units = a.middle;
    a.middle = a.sign;
    a.sign = temp;
    return a;
    
}
vtri m7d_mul_3(vtri a)
{
  
    vec z = a.units | a.middle | a.sign;
    vec temp = a.units;
    a.units =  a.sign ^ z;
    a.sign =   a.middle  ^ z;
    a.middle = temp ^ z;
   
    return a;
    
}
vtri m7d_mul_4(vtri a)
{
    vec temp = a.units;
    a.units = a.sign;
    a.sign = a.middle;
    a.middle = temp;
    
    return a;
    
    
}
vtri m7d_mul_5(vtri a)
{
    vec z = a.units| a.middle | a.sign;
    vec temp = a.units;
    a.units = a.middle ^ z;
    a.middle = a.sign ^ z;
    a.sign = z ^ temp;
    return a;
    
}
vtri m7d_mul_6(vtri a)
{
    vec z = a.units | a.middle | a.sign;
	a.units =  a.units ^ z;
    a.sign =   a.sign  ^ z;
    a.middle = a.middle ^ z;
    
    return a;
     
}
void m7d_add_4r(vtri *x, vtri * y)
{
 	vtri  r;
    vec s;
    vec t;
	r.sign = x->units ^ y->units;
    s = (x->units & y->units);
    r.middle = s^ x->middle ^ y->middle;
    t = (((s) & (x->middle | y->middle)) | (x->middle & y->middle) );
    r.sign = x->sign ^ y->sign ^ t;
    s = x->sign | y->sign | t;
    t = (r.units & s );
    x->units  = s ^ r.units;
    x->middle  = x->middle ^ t ;
    x->sign  = x->sign ^ (  t & x->middle);
    
    //Optimize the multiplication later
    r.units = x->units;
    x->units = x->sign;
    x->sign = x->middle;
    x->middle = r.units;
}





/*************
sub_m7dr is a placeholder for now
**********************/
vtri sub_m7dr(vtri const x, vtri const y)
{
    vtri r;
    r.units = ((x.units^y.units) | (x.sign^y.sign));
    r.sign = (((x.units^y.units)^x.sign)&(y.units ^ x.sign));
    
    return r;
}

void m7d_sub_64(vtri **R, vtri  **A, vtri  **B)
{
    int i;
    for (i= 0; i < M1RI_RADIX; i++ )
    {
        R[i][0] = sub_m7dr(A[i][0], B[i][0]);
    }
    
}

int m7d_equal(m7d_t const *a, m7d_t const *b)
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
            if((a->rows[i][j].sign != b->rows[i][j].sign) || (a->rows[i][j].units != b->rows[i][j].units ) || (a->rows[i][j].middle != b->rows[i][j].middle ))
            {
                return 0;
            }
            
        }
    }
    return 1;
}


void m7d_copypadding(m7d_t  * r, m7d_t  const * x)
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

void m7d_putpadding(m7d_t  * r, m7d_t  const * x)
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

void m7d_add_64(vtri **R, vtri   **A, vtri  **B)
{
    int i;
    for (i = 0; i < M1RI_RADIX; i++ )
    {
    	add_vtri(&R[i][0], &A[i][0], &B[i][0]);
       // R[i][0] = add_m7dr(A[i][0], B[i][0]);
    }

}
void m7d_add_r(m7d_t *c, m7d_t *a, m7d_t *b)
{

    if((a->nrows == b->nrows) && ( b->ncols == a->ncols))
    {
        int i, j;
     
        for( i = 0; i < a->nrows; i++)
        {
            for(j = 0; j < (a->width ); j++)
            {
                add_vtri(&c->rows[i][j], &a->rows[i][j], &b->rows[i][j]);    
            }
        }
        
    }
      
}


