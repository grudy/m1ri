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
 
 m5d.c
 */

#include "m5d.h"

/** 
 Read n bits from a sign portion of an element
 x = rows
 y = columns
 M = Matrix read from
 */

vec m5d_rm_bits(m5d_t *M, rci_t  x, rci_t  y, int  n) {
    
	wi_t  block = (y  ) / M1RI_RADIX;
    int  spill = (y  % M1RI_RADIX) + n - M1RI_RADIX;
    vec bits;
    bits = (spill <= 0) ? M->rows[x][block].middle << -spill : (M->rows[x][block + 1].sign << (M1RI_RADIX - spill)) | (M->rows[x][block].middle << spill);
    return bits;
      
}


vec m5d_rs_bits(m5d_t *M, rci_t  x, rci_t  y, int  n) {
    
	wi_t  block = (y  ) / M1RI_RADIX;
    int  spill = (y  % M1RI_RADIX) + n - M1RI_RADIX;
    vec bits;
    bits = (spill <= 0) ? M->rows[x][block].sign << -spill : (M->rows[x][block + 1].sign << (M1RI_RADIX - spill)) | (M->rows[x][block].sign << spill);

    return bits;
}
/** 
 Read n bits from units
 x = rows
 y = columns
 M = Matrix read from
 */

vec m5d_ru_bits(m5d_t *M, rci_t  x, rci_t  y, int  n) {
    
    wi_t  block = (y  ) / M1RI_RADIX;
    int  spill = (y  % M1RI_RADIX) + n - M1RI_RADIX;
    vec bits;
    bits = (spill <= 0) ? M->rows[x][block].units << -spill : (M->rows[x][block + 1].units<< (M1RI_RADIX - spill)) | (M->rows[x][block].units << spill);
  
    return bits;
}


/** 
 Read n elements
 x = rows
 y = columns
 M = Matrix read from
 */

vfd m5d_read_elems(m5d_t *M, rci_t  x, rci_t  y, int  n) {
    
    
    wi_t  block = (y  ) / M1RI_RADIX;
    int  spill = (y  % M1RI_RADIX) + n - M1RI_RADIX;
    vfd elem;
    
    elem.units = (spill <= 0) ? M->rows[x][block].units << -spill : ((M->rows[x][(block + 1)].units<< (M1RI_RADIX - spill)) | (M->rows[x][block].units >> spill));
    elem.sign = (spill <= 0) ?  (M->rows[x][block].sign << -spill) : (M->rows[x][block + 1].sign << (M1RI_RADIX - spill)) | (M->rows[x][block].sign>> spill);
    elem.middle = (spill <= 0) ?  (M->rows[x][block].middle << -spill) : (M->rows[x][block + 1].middle << (M1RI_RADIX - spill)) | (M->rows[x][block].middle>> spill);
    elem.middle = (elem.middle >> (M1RI_RADIX - n));
    elem.units = (elem.units >> (M1RI_RADIX - n));
    elem.sign = (elem.sign >> (M1RI_RADIX - n));

    return elem;
    
    
}





/** 
 Swap rows in a matrix;  Finish the else statement here
 */
void * m5d_rowswap (m5d_t * M, rci_t row_a, rci_t  row_b)
{
    
    
    if((M->nrows >= (row_a ) && (M->nrows >= row_b)))
    {
        vfd temp;
        temp =  *M->rows[row_a -1];
        *M->rows[row_a -1] = *M->rows[row_b -1];
        *M->rows[row_b -1] =  temp;
      
        
    }

    return 0;
}


/** 
 
 */


/* unfinished */
void *  m5d_write_elem( m5d_t * M,rci_t x, rci_t y, vec s, vec m , vec u)
{
      
    wi_t  block = (y  ) / M1RI_RADIX;
    int   spill = (y  % M1RI_RADIX) - (M1RI_RADIX -1);
	s = ~(s == 0);
    u = ~(u == 0);
	M->rows[x][block].units  = (u == 0) ? (~(rightbit << -spill) &  (M->rows[x][block].units))  : ((u << (M1RI_RADIX - spill)) | (M->rows[x][block].units));
    M->rows[x][block].sign  = (s == 0) ? (~(rightbit << -spill) &  (M->rows[x][block].units))  : ((u << (M1RI_RADIX - spill)) | (M->rows[x][block].units));
    M->rows[x][block].middle  = (m == 0) ? (~(rightbit << -spill) &  (M->rows[x][block].units))  : ((u << (M1RI_RADIX - spill)) | (M->rows[x][block].units));
    
    return 0;
}


/** 
 
 */





/** 
 
 */

vfd ** m5d_row_alloc(vfd * block, vfd ** rows, wi_t width, rci_t nrows)
{
 
    int i;
    rows = m1ri_calloc( nrows ,  sizeof(vfd *));
    for ( i = 0; i <  nrows;  i++ )
    {
        rows[i]  = (block + (i * width));
        
        
    };
    
    return rows;
}



m5d_t  * m5d_create( rci_t nrows, rci_t ncols)
{
    
    m5d_t * a = m1ri_malloc(sizeof(m5d_t));
    a->ncols = ncols;
    a->nrows = nrows;
    a->width = M1RI_DN(ncols, M1RI_RADIX);
    a->block = m5d_block_allocate(a->block,  a->nrows,    a->width);
    a->rows  = m5d_row_alloc(a->block, a->rows, a->width, a->nrows);
    a->flags = notwindowed;
    a->fcol = 0;
    return a;
    
}



m5d_t *    m5d_init_window(const m5d_t *c,const rci_t strow, const rci_t stvfd, const rci_t sizerows, const  rci_t sizecols)
{
   
    m5d_t * submatrix = m1ri_malloc(sizeof(m5d_t));
    /** c->width should not be compared twice */
    if((strow + sizerows) > c->width)
    {    
        return  0;
    }
    
    if((stvfd + sizecols) > c->ncols)
    {
        return  0;
    }
    int i, f;
	f = strow * M1RI_RADIX;
    submatrix->nrows =   M1RI_RADIX * sizerows;
    submatrix->ncols =  M1RI_RADIX * sizecols;
    submatrix->flags = iswindowed;
    submatrix->width =  sizecols;
    submatrix->rows = m1ri_calloc(M1RI_RADIX * sizerows ,  sizecols * sizeof(vfd *));
    submatrix->lblock = ( (sizerows +  strow)  ==  c->width)? c->lblock:  0;
    submatrix->fcol   = 0;
    submatrix->svfd = stvfd;
    
    
    
    for(  i =   f; i < (f + (M1RI_RADIX * sizerows)) ; i++)
    {
        submatrix->rows[i - f] = c->rows[i] + stvfd;   
    }
    return submatrix;
}

vfd * m5d_rand(m5d_t * a)
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
       			a->rows[i][z].middle = m1ri_rand();
       			
       			a->rows[i][z].sign =  a->rows[i][z].sign | a->rows[i][z].units | a->rows[i][z].middle;
            
            
            }
				a->rows[i][z].sign = m1ri_rand();
       			a->rows[i][z].units = m1ri_rand();
       			a->rows[i][z].middle = m1ri_rand();
       			
       			a->rows[i][z].sign =  a->rows[i][z].sign | a->rows[i][z].units | a->rows[i][z].middle;
            
       			a->rows[i][z].sign = a->rows[i][z].sign & mask_rand;
       			a->rows[i][z].middle  = a->rows[i][z].middle & mask_rand;

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
       			a->rows[i][z].middle = m1ri_rand();
       			
       			a->rows[i][z].sign =  a->rows[i][z].sign | a->rows[i][z].units | a->rows[i][z].middle;
            
            
            }
    
   	 	}
 	}   
}



void   m5d_set_ui(m5d_t * a, rci_t value)

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
                a->rows[k + i][a->width-1].units = 0;
                a->rows[l][a->width -1].middle = 0;
                a->rows[l][a->width -1].sign = (rightbit)<<i;
            }
        }
        
        if ((a->ncols%M1RI_RADIX) == 0)
        {
            
            l = (a->width - 1) * M1RI_RADIX;
            for(i  = 0; i < M1RI_RADIX; i++)
                
            {
                a->rows[l][a->width -1].units = 0;
                a->rows[l][a->width -1].middle = 0;
                a->rows[l][a->width -1].sign = (rightbit)<<i;
                l++;    
            }
            
        }
        
    }
    
}


m5d_t *  m5d_identity(m5d_t * a, rci_t n)
{
	
    a = m5d_create( n, n);
    m5d_set_ui(a, 1);
    
    return a;   
}


void m5d_free( m5d_t *  a)
{ 		
    m1ri_free(a->rows);
    if(a->flags == notwindowed)
    {
	  m1ri_free(a->block); 	
	}
    
    m1ri_free(a);
    a == NULL;
    
}

void add_vfd(vfd * r, vfd * a, vfd * b)
{
    vec c, d, e, f, g, h, i, j, k, l, m, n ,o, p, q;
    c = b->units ^ a->units;
    d = b->middle ^ a->middle;
    e = b->sign ^ a->sign;
    f = d & c;
    g = f | b->middle;
    h = f ^ a->sign;
    i = h | e;
	r->sign = i; /** */
    j = i ^ b->units;
    k = j ^ a->middle;
    l = k | c;
    m = l ^ e;
    n = m ^ g;
	r->middle = n; /** */
    o = m | d;
    p = o ^ c;
    q = p^n;
    r->units = q;
    
        /** */

	/** 
	more optimized? 
    vec c, d, e, f,   j,  m;
    c = b->units ^ a->units;
    d = b->middle ^ a->middle;
    e = b->sign ^ a->sign;
    f = d & c;
    r->sign = (f ^ a->sign) | e;
    j = r->sign ^ b->units;
    m = ((j ^ a->middle)| c) ^ e;
	r->middle  = m ^ (f | b->middle);
    r->units = ((m | d) ^ c)^(r->middle)) ;
     = q; */
    /** 
	*/
    
   /** * def add(a,b):
    c = b[0] ^ a[0]
    d = b[1] ^ a[1]
    e = b[2] ^ a[2]
    f = d & c
    g = f | b[1]
    h = f ^ a[2]
    i = h | e
    j = i ^ b[0]
    k = j ^ a[1]
    l = k | c
    m = l ^ e
    n = m ^ g
    o = m | d
    p = o ^ c
    q = p ^ n
    return q,n,i
    */
    
}

void m5d_add2(vfd * r, vfd * a, vfd * b)
{
    
    vec c, d, e, f, g, h, i, j, k, l, m, n ,o, p, q;
    c = b->units ^ a->middle;
    d = b->middle ^ a->units;
    e = c & a->middle;
    f = d ^ b->sign;
    g = d ^ a->sign;
    h = g | f;
    i = h | c;
	r->sign = i; /** */
    j = h | e;
    k = h ^ a->sign;
    l = k | b->middle;
    m = l & j;
    n = m ^ c;
    r->middle = n;
    o = n | c;
    p = o ^ f;
    q = p ^ e;
	r->units = q;
    /** 
	def add2(a,b):
    c = b[0] ^ a[1]
    d = b[1] ^ a[0]
    e = c & a[1]
    f = d ^ b[2]
    g = d ^ a[2]
    h = g | f
    i = h | c
    j = h | e
    k = h ^ a[2]
    l = k | b[1]
    m = l & j
    n = m ^ c
    o = n | c
    p = o ^ f
    q = p ^ e
    return q,n,i
	*/
}

void m5d_add2_i(vfd * a, vfd * b)
{
	vec c, d, e, f, g, h, i, j, k, l, m, o, p, q;
    c = b->units ^ a->middle;
    d = b->middle ^ a->units;
    e = c & a->middle;
    f = d ^ b->sign;
    g = d ^ a->sign;
    h = g | f;
    i = h | c;
    j = h | e;
    k = h ^ a->sign;
    l = k | b->middle;
    m = l & j;
    a->middle = m ^ c;
    o = a->middle  | c;
    p = o ^ f;
    q = p ^ e;
	a->units = q;
	a->sign = i;
    /** 
	def add2(a,b):
    c = b[0] ^ a[1]
    d = b[1] ^ a[0]
    e = c & a[1]
    f = d ^ b[2]
    g = d ^ a[2]
    h = g | f
    i = h | c
    j = h | e
    k = h ^ a[2]
    l = k | b[1]
    m = l & j
    n = m ^ c
    o = n | c
    p = o ^ f
    q = p ^ e
    return q,n,i
	*/

}

void vfd_sub( vfd *r, vfd *a, vfd *b)
{
    vec c, d, e, f, g, h, i, j, k, l, m, n ,o, p, q;
    c = b->units ^ a->units;
    d = b->middle ^ a->middle;
    e = c ^ b->sign;
    f = e | b->units;
    g = f ^ a->sign;
    h = g | d;
    r->sign = h | c;
    j = r->sign ^ c;
    k = j & a->sign;
    l = k | b->middle;
    m = l ^ g;
    r->middle = m ^ a->middle;
    o = m | d;
    p = o ^ c;
    r->units = p ^ a->middle;

    
    /** 
def sub(a,b):
    c = b[0] ^ a[0]
    d = b[1] ^ a[1]
    e = c ^ b[2]
    f = e | b[0]
    g = f ^ a[2]
    h = g | d
    i = h | c  last
    j = i ^ c
    k = j & a[2]
    l = k | b[1]
    m = l ^ g
    n = m ^ a[1]
    o = m | d
    p = o ^ c
    q = p ^ a[1]
    return q,n,i
    */

}


 /* matrix r = (matrix r - matrix x) */
void fb_i(vfd *r,vfd *x) 
{
    vec c, d, e, f, g, h, i, j, k, l, m, n ,o, p, q;
    c = x->units ^ r->units;
    d = x->middle ^ r->middle;
    e = c ^ x->sign;
    f = e ^ x->middle;
    g = f ^ r->sign;
    h = g | d;
    i = h | c;
    r->sign = i;
    j = i ^ c;
    k = j & r->sign;
    l = k | x->middle;
    m = l ^ g;
    n = m ^ r->middle;
    r->middle = n;
    o = m | d;
    p = o ^ c;
    q = p ^ r->middle;
    r->units = q;
       
}

/** 
	Scalar Multiplication
*/
void m5d_mul2(vfd * a)
{
    vec t = a->units;
	a->units = a->middle;
	a->middle = ~t & a->sign;
    
    
}
/** 
	Scalar Multiplication
*/
void m5d_mul3(vfd * a)
{
    vec t = a->units;
    a->units =  (~a->middle) & a->sign;
    a->middle = t;
    
    
}
/** 
	Scalar Multiplication
*/
void m5d_mul4(vfd * a)
{
  
    a->units = (~a->units) & a->sign;
    a->middle = (~a->middle) & a->sign;
    

}
 /* matrix r = (matrix r - matrix x) */
void m5d_dec(vfd *r,vfd *x) 
{
    vec c, d, e, f, g, h,  i, j, k, l, m, o, p, q;
    c = x->units ^ r->units;
    d = x->middle ^ r->middle;
    e = c ^ x->sign;
    f = e | x->units;
    g = f ^ r->sign;
    h = g | d;
    
	i = r->sign;
    r->sign =  h | c;
    j = r->sign ^ c;
    k = j & i;
    l = k | x->middle;
    m = l ^ g;
	q = r->middle;
    r->middle = m ^ r->middle;
    o = m | d;
    p = o ^ c;
    r->units = p ^ q;
    
       
}

m5d_t * m5d_add(m5d_t * c,  const m5d_t  *a,const m5d_t  *b)
{

 	if (c == NULL)
	{
		c = m5d_create(a->nrows, b->ncols);

	} 
	else if( (c->nrows != a->nrows || c->ncols != b->ncols)) 
	{
		m1ri_die("m5d_add: Provided return matrix has wrong dimensions.\n");
    	
	
	}	         

	if((a->nrows != b->nrows) || ( b->ncols != a->ncols))
    {
       
      m1ri_die("m5d_add: Input Matrices must have same dimension.\n");
    }
    
    
    int i, j;
	for( i = 0; i < a->nrows; i++)
        {
            for(j = 0; j < (a->width ); j++)
            {   
        
               add_vfd(&c->rows[i][j], &a->rows[i][j], &b->rows[i][j]);
               
            }
        }
        
    
    return c;
    
}
/**
  Checks if an m5d_t is equal to another. 
*/
int m5d_equal(m5d_t const *a, m5d_t const *b)
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

void m5d_add_64(vfd **R, vfd   **A, vfd  **B)
{
    
    for (int i = 0; i < M1RI_RADIX; i++ )
    {
		add_vfd(R[i], A[i], B[i]);
    }

}




m5d_t * m5d_copy_cutoff(m5d_t  * a, m5d_t  const * b)
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
void m5d_sub_64(vfd **c ,vfd **  a , vfd **b)
{
    int i;
    for (i= 0; i < M1RI_RADIX; i++ )
    {
    
        vfd_sub( c[i], a[i], b[i]);
    }
}



m5d_t * m5d_sub(m5d_t * r, const m5d_t  * a ,const  m5d_t * b)
{
  int i, j;
  if (r == NULL)
	{
		r = m5d_create(a->nrows, b->ncols);

	} 
	else if( (r->nrows != a->nrows || r->ncols != b->ncols)) 
	{
		m1ri_die("m5d_sub: Provided return matrix has wrong dimensions.\n");
    	
	
	}
    if((a->nrows != b->nrows) || ( a->ncols != b->ncols))
    {
       
      m1ri_die("m5d_sub: Input Matrices must have same dimension.\n");
    }

    for( i = 0; i < a->nrows; i++)
    {
        for(j = 0; j < (a->width ); j++)
        {   
                vfd_sub(&r->rows[i][j], &a->rows[i][j], &b->rows[i][j]);
        }
    }
        

    return r;
}

void m5d_sub_d(m5d_t  * a , m5d_t * b)
{
    if((a->nrows == b->nrows) && ( b->ncols == a->ncols))
    {
        int i, j;
        for( i = 0; i < a->nrows; i++)
        {
            for(j = 0; j < (a->width ); j++)
            {   
                m5d_dec(&a->rows[i][j], &b->rows[i][j]);
            }
        }
        
    }


}




static inline void m5d_combine4(vfd *table, vfd  ** const input) 
{
    vfd  t,   a,  b,  c,  d;
    t.sign = t.units =  t.middle = 0;
    a = *input[0];
    b = *input[1];
    c = *input[2];
    d = *input[3];

    table[0] = t;
    table[1] = a;
    table[2] = b;
    table[4] = c;
    table[8] = d;

    add_vfd(&t,&c,&d);
    table[12] = t;
    
    add_vfd(&t,&b,&c);
    table[6] = t;
    m5d_inc(&t,&d);
    table[14] = t;
    m5d_dec(&t,&c);
    table[10] = t;

    add_vfd(&t,&a,&b);
    table[3] = t;
    m5d_inc(&t,&d);
    
    table[11] = t;
    m5d_inc(&t,&c);
    table[15] = t;
    m5d_dec(&t,&d);
    
    table[7] = t;
    m5d_dec(&t,&b);
    table[5] = t;
    m5d_inc(&t,&d);
    table[13] = t;
	m5d_dec(&t,&c);
    table[9] = t;
}

static inline void m5d_combine5(vfd *table, vfd  **  const input) 
{
    vfd *e, *t4;
    int i;

    m5d_combine4(table, input);
    e = input[4];
    t4 = table+16;
    table[16] = *e;
    
    for(i=1;i<16;i++)
        add_vfd(t4 + i, table+i, e);
}
    
static inline void m5d_combine6(vfd *table, vfd  ** const input) 
{
    vfd *e, *t4;
    vfd * f, *t5;
    int i;
    
    m5d_combine4(table, input);
    e = input[4];
    t4 = table+16;
    table[16] = *e;

    f = input[5];
    t5 = table+32;
    table[32] = *f;
    
    for(i=1;i<16;i++)
        add_vfd(t4 + i, table+i, e);

    for(i=1;i<32;i++)
        add_vfd(t5 + i, table+i, f);

}


void m5d_mul_64(vfd **R, vfd  ** const A, vfd   ** const B)
{
   
    vfd t, r, r2, * a;
    vec v;
    int i;

    vfd tables6[4][64];
    vfd tables5[8][32];

    for(i=0;i<4;i++)
    {
        m5d_combine6(tables6[i], B + 6*i);
    }
    for(i=0;i<8;i++)
    {
        m5d_combine5(tables5[i], B + 24 + (5*i));
	}
    for(i=0;i<64;i++) 
    {
    /*
    	001 1	
		011 2
		101 3
		111 4
    
    */
    
    
    	/*   first part is ones */
        a = A[i];
        
		
		v = ~(a->units) &  ~(a->middle) & a->sign;

		
        
        
        r = tables6[0][v&63];                 v >>= 6;
        t = tables6[1][v&63]; m5d_inc(&r, &t); v >>= 6;
        t = tables6[2][v&63]; m5d_inc(&r, &t); v >>= 6;
        t = tables6[3][v&63]; m5d_inc(&r, &t); v >>= 6;
        
        t = tables5[0][v&31]; m5d_inc(&r, &t); v >>= 5;
        t = tables5[1][v&31]; m5d_inc(&r, &t); v >>= 5;
        t = tables5[2][v&31]; m5d_inc(&r, &t); v >>= 5;
        t = tables5[3][v&31]; m5d_inc(&r, &t); v >>= 5;
        t = tables5[4][v&31]; m5d_inc(&r, &t); v >>= 5;
        t = tables5[5][v&31]; m5d_inc(&r, &t); v >>= 5;
        t = tables5[6][v&31]; m5d_inc(&r, &t); v >>= 5;
        t = tables5[7][v&31]; m5d_inc(&r, &t);

      
      
      
      	v =  (~a->units)   & (a->middle) & a->sign;
		
		
	
        r2 = tables6[0][v&63];                 v >>= 6;
        t = tables6[1][v&63]; m5d_inc(&r2, &t); v >>= 6;
        t = tables6[2][v&63]; m5d_inc(&r2, &t); v >>= 6;
        t = tables6[3][v&63]; m5d_inc(&r2, &t); v >>= 6;
        
        t = tables5[0][v&31]; m5d_inc(&r2, &t); v >>= 5;
        t = tables5[1][v&31]; m5d_inc(&r2, &t); v >>= 5;
        t = tables5[2][v&31]; m5d_inc(&r2, &t); v >>= 5;
        t = tables5[3][v&31]; m5d_inc(&r2, &t); v >>= 5;
        t = tables5[4][v&31]; m5d_inc(&r2, &t); v >>= 5;
        t = tables5[5][v&31]; m5d_inc(&r2, &t); v >>= 5;
        t = tables5[6][v&31]; m5d_inc(&r2, &t); v >>= 5;
        t = tables5[7][v&31]; m5d_inc(&r2, &t);

		

		m5d_mul2(&r2);
		v = (a->units) &   (~a->middle) & a->sign;
		m5d_inc(&r, &r2);

		
		r2 = tables6[0][v&63];                 v >>= 6;
        t = tables6[1][v&63]; m5d_inc(&r2, &t); v >>= 6;
        t = tables6[2][v&63]; m5d_inc(&r2, &t); v >>= 6;
        t = tables6[3][v&63]; m5d_inc(&r2, &t); v >>= 6;
        
        t = tables5[0][v&31]; m5d_inc(&r2, &t); v >>= 5;
        t = tables5[1][v&31]; m5d_inc(&r2, &t); v >>= 5;
        t = tables5[2][v&31]; m5d_inc(&r2, &t); v >>= 5;
        t = tables5[3][v&31]; m5d_inc(&r2, &t); v >>= 5;
        t = tables5[4][v&31]; m5d_inc(&r2, &t); v >>= 5;
        t = tables5[5][v&31]; m5d_inc(&r2, &t); v >>= 5;
        t = tables5[6][v&31]; m5d_inc(&r2, &t); v >>= 5;
        t = tables5[7][v&31]; m5d_inc(&r2, &t);

    	
    		

		v  = a->units   & a->middle & a->sign;
		m5d_mul3(&r2);
		m5d_inc(&r, &r2);



    	r2 = tables6[0][v&63];                 v >>= 6;
        t = tables6[1][v&63]; m5d_inc(&r2, &t); v >>= 6;
        t = tables6[2][v&63]; m5d_inc(&r2, &t); v >>= 6;
        t = tables6[3][v&63]; m5d_inc(&r2, &t); v >>= 6;
        
        t = tables5[0][v&31]; m5d_inc(&r2, &t); v >>= 5;
        t = tables5[1][v&31]; m5d_inc(&r2, &t); v >>= 5;
        t = tables5[2][v&31]; m5d_inc(&r2, &t); v >>= 5;
        t = tables5[3][v&31]; m5d_inc(&r2, &t); v >>= 5;
        t = tables5[4][v&31]; m5d_inc(&r2, &t); v >>= 5;
        t = tables5[5][v&31]; m5d_inc(&r2, &t); v >>= 5;
        t = tables5[6][v&31]; m5d_inc(&r2, &t); v >>= 5;
        t = tables5[7][v&31]; m5d_inc(&r2, &t);

    	m5d_mul4(&r2);
    	
    	m5d_inc(&r, &r2);
    	
		R[i][0] = r;
        
    }
    
}



 
m5d_t  * m5_blockslice_allocate(rci_t  nrows,  wi_t  width)
{
	
    m5d_t * block  = m1ri_calloc(nrows * width ,  sizeof(m5d_t * ) );
    return block;
}

m5d_t ** m5_rowslice_allocate(m5d_t * block,  wi_t width, rci_t nrows)
{
	int i;
    m5d_t **rows = m1ri_malloc( nrows * width * sizeof(m5d_t *));
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
 	c->block = m5_blockslice_allocate(c->nrows,   c->width);
    c->row = m5_rowslice_allocate(c->block , c->width, c->nrows);
    r = 0 ;
     
    for ( i = 0; i <  l;  i = i + slicesize)
    {       
    	for( f = 0; f <colroundeddown ; f++)
        {
        	c->row[(r * f) + f] = m5d_init_window(a,i , (f * slicesize), slicesize, slicesize);
        }
        
        if(extracols > 0)
        {
        	c->row[(r * f) + f] = m5d_init_window(a,i , (f * slicesize), slicesize, extracols);
		}
        r++;
        
   	}
    
    if(extrarows >0 )
    {
		for( f = 0; f <colroundeddown ; f++)
        {
           c->row[(r * f) + f] = m5d_init_window(a, i , (f * slicesize), extrarows, slicesize);
        }

    	if(extracols > 0)
    	{
           c->row[(r * f) + f] = m5d_init_window(a, i , (f * slicesize), extrarows, extracols);
     	}
    }
}

static inline vfd *  m5d_transpose_vfd(vfd  **a, vfd  **b  )
{
    int i, x;
    vfd temp;
    for(i = 0; i <64; i ++)
    {
      for(x = 0; x < 64; x ++)
       {
       
        temp.units =  (a[x][0].units & (rightbit << i ) );
        temp.sign =  (a[x][0].sign & (rightbit <<  i) );
        temp.middle =  (a[x][0].middle & (rightbit << i) );    
        b[i][0].units = (temp.units) ?  b[i][0].units | (rightbit << x) : b[i][0].units ;
        b[i][0].sign = (temp.sign) ? b[i][0].sign | (rightbit << x) : b[i][0].sign  ;
        b[i][0].middle = (temp.middle) ? b[i][0].middle | (rightbit <<  x) : b[i][0].middle;    
        }

    }
    
    return *b;
}

m5d_t *  m5d_transpose_sliced(m5d_t * a)
{
    int x, y;
    m5d_t * c;
    c = m5d_create(a->ncols, a->nrows);
    m5_slice * b, *d;
    d = malloc(sizeof(m5_slice));
    b = malloc(sizeof(m5_slice));
    m5d_slices(b, a, 1);
    m5d_slices(d, c, 1);
    for (x = 0; x < b->nrows; x++) {
        for (y = 0; y < b->ncols; y ++) {
			m5d_transpose_vfd(b->row[x][y].rows, d->row[y][x].rows);  
        }
    }
    return c;
}


/*
m5d_init_window without checks
	
*/
static inline m5d_t   * m5d_init_window_unshackled(const m5d_t *c,const rci_t strow, const rci_t stvfd, const rci_t sizerows, const  rci_t sizecols)
{
	
	m5d_t * submatrix = m1ri_malloc(sizeof(m5d_t));
    int i, f;
	f = strow * M1RI_RADIX;
    submatrix->nrows =   M1RI_RADIX * sizerows;
    submatrix->ncols =  M1RI_RADIX * sizecols;
    submatrix->flags = iswindowed;
    submatrix->width =  sizecols;
    submatrix->rows = m1ri_calloc(M1RI_RADIX * sizerows ,  sizecols * sizeof(vfd *));
    submatrix->lblock = ( (sizerows +  strow)  ==  c->width)? c->lblock:  0;
    submatrix->fcol   = 0;
    submatrix->svfd = stvfd;
    
    
    
    for(  i =   f; i < (f + (M1RI_RADIX * sizerows)) ; i++)
    {
        submatrix->rows[i - f] = c->rows[i] + stvfd;   
    }
    return submatrix;
}


m5_slice * m5d_quarter(const m5d_t * a)
{
	 m5_slice * c = m1ri_malloc(sizeof(m5_slice));
	 
	 
	
     c->row = m1ri_calloc( 4 , sizeof(m5d_t **));

     
     
     c->row[0] = m5d_init_window_unshackled(a,  0, 0 , a->nrows/128, a->ncols/128);
	 c->row[1] = m5d_init_window_unshackled(a, 0, a->ncols/128 , a->nrows/128, a->ncols/128);   
     c->row[2] = m5d_init_window_unshackled(a, a->nrows/128, 0 , a->nrows/128, a->ncols/128);
	 c->row[3] = m5d_init_window_unshackled(a, a->nrows/128,a->ncols/128,  a->nrows/128, a->ncols/128);
	 
	 
	 

    
    return c;
}

 
 
 /**
   Orders the elements computed in vfd_elem for the final part (addition)
 
 */
static inline void vfd_elem_order(vfd * a)
 {
   vec temp1, temp2;
   temp1 = a->units & ~(a->sign);
   temp2 = a->middle & ~(a->sign);
   a->middle = a->middle & a->sign;
   a->sign = temp1 | temp2 | a->sign;
   a->middle = a->middle | (temp1 & ~(temp2));
    a->units  = (temp1 & temp2) | (a->units & ~(temp1));
   
 }
 
 
void m5d_add_i(m5d_t * x, m5d_t *y) 
{
	 int i, j;
        for( i = 0; i < x->nrows; i++)
        {
            for(j = 0; j < (x->width ); j++)
            {
            	m5d_inc(&x->rows[i][j], &y->rows[i][j]);
        	}   
        }
	

} 
 
 
 
 void m5d_sub_i(m5d_t * x, m5d_t *y) 
{
	 int i, j;
        for( i = 0; i < x->nrows; i++)
        {
            for(j = 0; j < (x->width ); j++)
            {
            	m5d_dec(&x->rows[i][j], &y->rows[i][j]);
        	}   
        }
	

}
 
 
/*

	\brief incremental subtraction, but where the subtrahend is changed
	
*/


 

static inline void sub_m5d_r(vfd  const *r,vfd   *x)
{
	vfd y;
	y.units = x->units;
	y.sign = x->sign;
    
      vec c, d, e, f, g, h,  i, j, k, l, m, o, p, q;
    c = x->units ^ r->units;
    d = x->middle ^ r->middle;
    e = c ^ x->sign;
    f = e | x->units;
    g = f ^ r->sign;
    h = g | d;
    
	i = r->sign;
    x->sign =  h | c;
    j = x->sign ^ c;
    k = j & i;
    l = k | x->middle;
    m = l ^ g;
	q = r->middle;
    x->middle = m ^ r->middle;
    o = m | d;
    p = o ^ c;
    x->units = p ^ q;
         
}



void m5d_sub_r(m5d_t   *x , m5d_t   const *r)
{
  
	int n , i;
	for(i = 0; i < x->nrows; i++)
    {
    	for(n = 0; n < x->width; n++)
        {
        
       
		  sub_m5d_r(r->rows[i] + n,  x->rows[i] + n );
        }
    }
   



} 
 
static inline void vfd_elem(vfd * c, vfd const * a, vfd const * b)
{
  vfd one, two ,three, four;
  vfd tempone, temptwo;
  one.units = one.sign = one.middle = a->units & b->units;
  
  two.middle =  a->middle & b->middle;
  two.sign = a->middle & b->sign;
  two.units = a->middle & b->units;
  
  three.middle =  a->sign & b->middle;
  three.sign = a->sign & b->sign;
  three.units = a->sign & b->units;
  
  four.sign = four.middle = a->units & b->middle;
  
  four.units = a->units & b->sign;
  vfd_elem_order(&one);
  vfd_elem_order(&two);
  vfd_elem_order(&three);
  vfd_elem_order(&four);
    
    
  add_vfd(&tempone, &one, &two);
  add_vfd(&temptwo, &three, &four);
  add_vfd(c, &tempone, &temptwo);
  
  
  
  	  
}
m5d_t * m5d_hadamard(m5d_t * c, m5d_t const * a, m5d_t const * b )
{
    
    
     if (c == NULL)
	{
		c = m5d_create(a->nrows, b->ncols);

	} 
	else if( (c->nrows != a->nrows || c->ncols != b->ncols)) 
	{
		m1ri_die("m5d_hadamard: Provided return matrix has wrong dimensions.\n");	
	}
    if((a->nrows != b->nrows) || ( b->ncols != a->ncols))
    {
       
      m1ri_die("m5d_hadamard: Input Matrices must have same dimension.\n");
    }

    int i, j;

    for( i = 0; i < a->nrows; i++)
        {
            for(j = 0; j < (a->width ); j++)
            {  
               vfd_elem(c->rows[i] + j, a->rows[i]  + j, b->rows[i] + j);
            }  
          

    	}
        
        
    
    return c;

}
 
 
 
m5d_t *  m5d_copy(m5d_t * a, m5d_t const *b)
{
  if(a == NULL)
  {	
  	a = m5d_create( b->nrows, b->ncols);
  }
  
  if((a->ncols < b->ncols) || (a->nrows < b->nrows))
  {
  	m1ri_die("m5d_copy: Provided return matrix has wrong dimensions.\n");
  
  }
  
  for(int i = 0; i < b->nrows; i++)
  {
    for(int j = 0; j < b->width; j++)
    {
    
      a->rows[i][j] = b->rows[i][j];
    
    }
    
  
  
  }
  
     a->lblock = b->lblock; /*   first block pointed to in a window */
     a->fcol = b->fcol;  /* /column offset of first block */
     a->flags = b->flags;
  return a;
  


}
void  m5d_colswap(m5d_t *M, rci_t col_a, rci_t col_b)
{
    if((M->ncols >= (col_a ) && (M->ncols >= col_b)))
    {
        int i;
        vec block_a, block_b, dif_a, dif_b, a_place, b_place; 
        vfd tempa, tempb;
         block_a = (col_a-1)/M1RI_RADIX;
         block_b = (col_b-1)/M1RI_RADIX;
         dif_a = col_a%M1RI_RADIX;
         dif_b = col_b%M1RI_RADIX;
         a_place =  rightbit <<   dif_a ;
         b_place =  rightbit <<  dif_b ;
        if(block_a == block_b)
        { 

              
          for( i = 0; i < M->nrows; i++)
          {
		     
		  
		           tempa.units  = (b_place  & M->rows[i][block_b].units) ? (a_place  ): 0;
		     tempb.units  = (a_place  & M->rows[i][block_a].units) ? (b_place  ): 0; 
		     
		     tempa.middle  = (b_place  & M->rows[i][block_b].middle) ? (a_place  ): 0;
		     tempb.middle  = (a_place  & M->rows[i][block_a].middle) ? (b_place  ): 0; 
		     
		     tempa.sign  = (b_place  & M->rows[i][block_b].sign) ? (a_place  ): 0;
		     tempb.sign  = (a_place  & M->rows[i][block_a].sign) ? (b_place  ): 0; 
		       M->rows[i][block_a].units  = (tempa.units)  ? M->rows[i][block_a].units  | tempa.units :   M->rows[i][block_a].units  & ~a_place; 
		       M->rows[i][block_b].units  = (tempb.units)  ? M->rows[i][block_a].units  | tempb.units :   M->rows[i][block_b].units  & ~b_place;  
		     
		       M->rows[i][block_a].sign  = (tempa.sign)  ? (M->rows[i][block_a].sign  | tempa.sign) :   M->rows[i][block_a].sign  & ~a_place; 
		       M->rows[i][block_b].sign  = (tempb.sign)  ? (M->rows[i][block_a].sign  | tempb.sign) :   M->rows[i][block_b].sign  & ~b_place; 
		      
		       M->rows[i][block_a].middle  = (tempa.middle)  ? (M->rows[i][block_a].middle  | tempa.middle) :   M->rows[i][block_a].middle  & ~a_place; 
		       M->rows[i][block_b].middle  = (tempb.middle)  ? (M->rows[i][block_a].middle  | tempb.middle) :   M->rows[i][block_b].middle  & ~b_place; 
		     
		       
		       

		       
          }
    
        }
        
      
        
        
    }
    
}

void m5d_colswap_capped_row(m5d_t *M, rci_t col_a, rci_t col_b, rci_t start_row)
{
      if((M->ncols >= (col_a ) && (M->ncols >= col_b)))
    {
        int i;
        vec block_a, block_b, dif_a, dif_b, a_place, b_place; 
        vfd tempa, tempb;
         block_a = (col_a-1)/M1RI_RADIX;
         block_b = (col_b-1)/M1RI_RADIX;
         dif_a = col_a%M1RI_RADIX;
         dif_b = col_b%M1RI_RADIX;
         a_place =  rightbit <<  dif_a ;
         b_place =  rightbit <<  dif_b ;
        if(block_a == block_b)
        { 

          for( i = start_row; i < M->nrows; i++)
          {
		     
		  
		           tempa.units  = (b_place  & M->rows[i][block_b].units) ? (a_place  ): 0;
		     tempb.units  = (a_place  & M->rows[i][block_a].units) ? (b_place  ): 0; 
		     
		     tempa.middle  = (b_place  & M->rows[i][block_b].middle) ? (a_place  ): 0;
		     tempb.middle  = (a_place  & M->rows[i][block_a].middle) ? (b_place  ): 0; 
		     
		     tempa.sign  = (b_place  & M->rows[i][block_b].sign) ? (a_place  ): 0;
		     tempb.sign  = (a_place  & M->rows[i][block_a].sign) ? (b_place  ): 0; 
		       M->rows[i][block_a].units  = (tempa.units)  ? M->rows[i][block_a].units  | tempa.units :   M->rows[i][block_a].units  & ~a_place; 
		       M->rows[i][block_b].units  = (tempb.units)  ? M->rows[i][block_a].units  | tempb.units :   M->rows[i][block_b].units  & ~b_place;  
		     
		       M->rows[i][block_a].sign  = (tempa.sign)  ? (M->rows[i][block_a].sign  | tempa.sign) :   M->rows[i][block_a].sign  & ~a_place; 
		       M->rows[i][block_b].sign  = (tempb.sign)  ? (M->rows[i][block_a].sign  | tempb.sign) :   M->rows[i][block_b].sign  & ~b_place; 
		      
		       M->rows[i][block_a].middle  = (tempa.middle)  ? (M->rows[i][block_a].middle  | tempa.middle) :   M->rows[i][block_a].middle  & ~a_place; 
		       M->rows[i][block_b].middle  = (tempb.middle)  ? (M->rows[i][block_a].middle  | tempb.middle) :   M->rows[i][block_b].middle  & ~b_place; 
		     
		       
		       

		       
          }
    
        }
        
      
        
        
    } 

}


void  m5d_transpose(m5d_t   * a)
{

   
  	int x, y;
    m5d_t * c;
    c = m5d_create( a->ncols, a->nrows);
    m5_slice * b, *d;
    d = malloc(sizeof(m5_slice));
    b = malloc(sizeof(m5_slice));
    m5d_slices(b, a, 1);
    m5d_slices(d, c, 1);
    for (x = 0; x < b->nrows; x++) 
    {
        for (y = 0; y < b->ncols; y ++) 
        {
         	 m5d_transpose_vfd(b->row[x][y].rows, d->row[y][x].rows);
            
        }
    }


   
}


int m5d_is_zero(const m5d_t *A)
{
	if(A == NULL)
	{
	
		m1ri_die("m5d_is_zero: A cannot be null!\n");

	
	}
	int i, j;

	for(i = 0; i < (A->width ); i++)
   	{
        for(j = 0; j < (A->width ); j++)
        {
            if((A->rows[i][j].units != 0) || (A->rows[i][j].sign != 0) || (A->rows[i][j].middle != 0));
            return 0;
            
        }   
    
    }
  
  
  
  return 1;

}

inline void m5d_mul_zero(m5d_t * a)
{
	
	int i, j;
	for(i = 0; i < a->nrows; i++)
	{
		for(j = 0; j < a->width; j++)
		{
			a->rows[i][j].units = 0;
			a->rows[i][j].middle = 0;
			a->rows[i][j].sign  = 0;
		  
		
		}
	
	}

}


inline void m5d_mul_two(m5d_t * a, const m5d_t * b)
{
	int i, j;
	for(i = 0; i < a->nrows; i++)
	{
		for(j = 0; j < a->width; j++)
		{
			
					
			a->rows[i][j].sign = 	b->rows[i][j].sign;
			a->rows[i][j].middle = 	b->rows[i][j].middle;		
			a->rows[i][j].units = 	b->rows[i][j].units;
			m5d_mul2(a->rows[i] + j);
		
		}
	
	}


}


inline void m5d_mul_three(m5d_t * a, const m5d_t * b)
{
	int i, j;
	for(i = 0; i < a->nrows; i++)
	{
		for(j = 0; j < a->width; j++)
		{
					
			a->rows[i][j].sign = 	b->rows[i][j].sign;
			a->rows[i][j].middle = 	b->rows[i][j].middle;		
			a->rows[i][j].units = 	b->rows[i][j].units;
			m5d_mul3(a->rows[i] + j);
		}
	
	}


}


inline void m5d_mul_fourth(m5d_t * a, const m5d_t * b)
{
	int i, j;
	for(i = 0; i < a->nrows; i++)
	{
		for(j = 0; j < a->width; j++)
		{
			
					
			a->rows[i][j].sign = 	b->rows[i][j].sign;
			a->rows[i][j].middle = 	b->rows[i][j].middle;		
			a->rows[i][j].units = 	b->rows[i][j].units;
			m5d_mul4(a->rows[i] + j);
		}
	
	}


}

m5d_t *m5d_mul_scalar(m5d_t *C, const long a, const m5d_t *B)
{

	if(C == NULL)
    { 
      C = m5d_create( B->ncols, B->nrows);
    }
    
    else if(C->nrows != B->nrows || C->ncols !=  B->ncols)
    {
    	m1ri_die("m5d_mul_scalar: C has wrong dimensions!\n");
    
    }
    
    
	long m = a%5;
	switch(m)
	{
		case 0: m5d_mul_zero(C);
		break ;
		case 1: m5d_copy(C, B);
		break;
  		case 2: m5d_mul_two(C, B);
		break ;
  		case 3: m5d_mul_three(C, B);
		break ;
  		case 4: m5d_mul_fourth(C, B);
  		break;
  	}	
  
  return C;
}

m5d_t * m5d_create_rand(m5d_t * a, rci_t n)
{
	 
	 a = m5d_create( n, n);
	 m5d_rand(a);
	 return a;
	

}


void m5d_add_row(m5d_t *A, rci_t ar, const m5d_t *B, rci_t br, rci_t start_col)
{
  

}






