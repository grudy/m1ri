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
    bits = (spill <= 0) ? M->rows[x][block].middle << -spill : (M->rows[x][block + 1].sign << (M1RI_RADIX - spill)) | (M->rows[x][block].middle >> spill);
    return bits;
      
}


vec m5d_rs_bits(m5d_t *M, rci_t  x, rci_t  y, int  n) {
    
	wi_t  block = (y  ) / M1RI_RADIX;
    int  spill = (y  % M1RI_RADIX) + n - M1RI_RADIX;
    vec bits;
    bits = (spill <= 0) ? M->rows[x][block].sign << -spill : (M->rows[x][block + 1].sign << (M1RI_RADIX - spill)) | (M->rows[x][block].sign >> spill);

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
    bits = (spill <= 0) ? M->rows[x][block].units << -spill : (M->rows[x][block + 1].units<< (M1RI_RADIX - spill)) | (M->rows[x][block].units>> spill);
  
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


//unfinished
void *  m5d_write_elem( m5d_t * M,rci_t x, rci_t y, vec s, vec u , vec m)
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



vfd  * m5d_block_allocate(vfd * block, rci_t  nrows,  wi_t  width)
{
    
    block  = m1ri_malloc(nrows * width * sizeof(vfd) );
    return block;
     
}

/** 
 
 */

vfd ** m5d_row_alloc(vfd * block, vfd ** rows, wi_t width, rci_t nrows)
{
 
    int i;
    rows = m1ri_malloc( nrows * width * sizeof(vfd *));
    for ( i = 0; i <  nrows;  i++ )
    {
        rows[i]  = (block + (i * width));
        
        
    };
    
    return rows;
}

/** 
 
 */

m5d_t m5d_create( m5d_t * a, rci_t nrows, rci_t ncols)
{
    
    
    a->ncols = ncols;
    a->nrows = nrows;
    a->width = M1RI_DN(ncols, M1RI_RADIX);
    a->block = m5d_block_allocate(a->block,  a->nrows,    a->width);
    a->rows  = m5d_row_alloc(a->block, a->rows, a->width, a->nrows);
    a->flags = 0;
    a->fcol = 0;
    return *a;
    
}

/** 
 
 */

/** 
 windows in m1ri_word rows * m1ri_word column incriments
 stvfd = the vfd or width offset from the base matrix
 strow = row offset in increments of 64
 sizecol  = cols * M1RI_RADIX
 sizerow  = rows * M1RI_RADIX
 */

m5d_t   m5d_window(m5d_t *c, rci_t strow, rci_t stvfd, rci_t sizerows, rci_t sizecols)
{
    int i;
    m5d_t  submatrix;
     /** c->width should not be compared twice */
    if((strow + sizerows) > c->width)
    {    
        return submatrix;
    }
    
    if((stvfd + sizecols) > c->width)
    {   
        return submatrix;
    }
    
    submatrix.nrows =   M1RI_RADIX * sizerows;
    submatrix.ncols =  M1RI_RADIX * sizecols;
    submatrix.flags = iswindowed;
    submatrix.width =  sizecols;
    submatrix.block = &c->block[(stvfd * stvfd)];
    submatrix.rows = m1ri_calloc(M1RI_RADIX * sizerows ,  sizecols * sizeof(vfd *));
    submatrix.lblock = ( (sizerows +  strow)  ==  c->width)? c->lblock:  0;
    submatrix.fcol   = 0;
    submatrix.svfd = stvfd;

    for(  i =   strow; i < (strow + (M1RI_RADIX * sizerows)) ; i++)
    {
        submatrix.rows[i - strow] = c->rows[i];
    }
    
    return submatrix;
    
}

void   m5d_window_create(m5d_t *c, m5d_t * submatrix, rci_t strow, rci_t stvfd, rci_t sizerows, rci_t sizecols)
{
     /** c->width should not be compared twice */
    
    if((strow + sizerows) > c->width)
    {   
        return;
    }
    
    if((stvfd + sizecols) > c->width)
    {
        return;    
    }
    int f = strow * M1RI_RADIX;
    int i;
    submatrix->nrows =   M1RI_RADIX * sizerows;
    submatrix->ncols =  M1RI_RADIX * sizecols;
    submatrix->flags = iswindowed;
    submatrix->width =  sizecols;
    submatrix->block =      m1ri_calloc(sizecols * sizerows, sizeof(m5d_t));
    submatrix->block = &c->block[(strow * stvfd)];
    submatrix->rows = m1ri_calloc(M1RI_RADIX * sizerows ,  sizecols * sizeof(vfd *));
    submatrix->lblock = ( (sizerows +  strow)  ==  c->width)? c->lblock:  0;
    submatrix->fcol   = 0;
    submatrix->svfd = stvfd;
   
    for(  i =   f; i < (f + (M1RI_RADIX * sizerows)) ; i++)
    {    
        submatrix->rows[i - f] = c->rows[i] + stvfd;
    }
    
}

vfd * m5d_rand(m5d_t * a)
{
    int i;
    for( i = 0; i < (a->nrows * a->width); i++)
    {
        a->block[i].sign = m1ri_rand();
        a->block[i].middle = m1ri_rand() &  a->block[i].sign;
        a->block[i].units = m1ri_rand() &  a->block[i].sign;;
    }
    return a->block;
}


/** 
 Make an Identity Matrix
 a = Identity matrix
 n = matrix size (row length and column width)
*/

m5d_t  m5d_identity_set(m5d_t * a)

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

/** 
 Creates an Identity Matrix over GF(5)
 */
m5d_t   m5d_identity(m5d_t  *a, rci_t n)
{
    *a = m5d_create(a, n, n);
    *a = m5d_identity_set(a);
    
    return *a;   
}

/** 
  Releases a m5d_t into the wilderness.
 */

void m5d_free( m5d_t *  tofree)
{ 


    m1ri_free(tofree->rows);
    m1ri_free(tofree->block);   
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
    /** */r->units = q;
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
    /** */// = q;
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

/** *******************************************
 matrix r = (direct sum matrix r + matrix x)
 ********************************************/
void iadd_vfd(vfd *r,vfd *x)
{
    vec c, d, e, f, g, h, i, j, k, l, m, n ,o, p, q;
    c = x->units ^ r->units;
    d = x->middle ^ r->middle;
    e = x->sign ^ r->sign;
    f = d & c;
    g = f | x->middle;
    h = f ^ r->sign;
    i = h | e;
	r->sign = i;  /** */
    j = i ^ x->units;
    k = j ^ r->middle;
    l = k | c;
    m = l ^ e;
    n = m ^ g;
	r->middle = n; /** */
    o = m | d;
    p = o ^ c;
    q = p^n;
    r->sign = q; /** */
   
}

 //matrix r = (matrix r - matrix x)
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
vfd m5d_mul2(vfd a)
{
    vec temp = a.sign;
    a.sign = a.sign ^ temp;
    a.units = a.middle & temp;
    a.middle = temp;
    
    return a;
}
/** 
	Scalar Multiplication
*/
vfd m5d_mul3(vfd a)
{
    vec temp = a.middle ^ a.sign;
    a.sign = a.middle;
    a.middle = a.units | temp;
    a.units = temp;
    
    return a;
}
/** 
	Scalar Multiplication
*/
vfd m5d_mul4(vfd a)
{
    vec x = a.units ^ a.sign;
    a.units = x & a.sign;
    a.middle = x;
    
    return a; 
}
 //matrix r = (matrix r - matrix x)
void m5d_sub_i(vfd *r,vfd *x) 
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

void m5d_add_r(m5d_t * c, m5d_t  *a, m5d_t  *b)
{
	         


    
    if((a->nrows == b->nrows) && ( b->ncols == a->ncols))
    {
        int i, j;
         m5d_create(c, b->nrows, b->ncols);
        for( i = 0; i < a->nrows; i++)
        {
            for(j = 0; j < (a->width ); j++)
            {   
                //add_vfd(c->rows[i + j] , a->rows[i + j], b->rows[i + j]);
               add_vfd(&c->rows[i][j], &a->rows[i][j], &b->rows[i][j]);
               
               
            }
        }
        
    }
    
    
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

void m5d_copypadding(m5d_t  * r, m5d_t  const * x)
{
		int  s;
        for(int i = 0; i < x->nrows; i++)
        {
			for( s = 0; s < x->width; s++)
        	{
            	r->rows[i][s] = x->rows[i][s];
            }   
        }
}

void m5d_putpadding(m5d_t  * r, m5d_t  const * x)
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

void m5d_sub_64(m5d_t * c ,m5d_t  * a , m5d_t * b)
{
    /** todo: Test this functions */
        for(int i = 0; i < a->nrows; i++)
        {
            vfd_sub(c->rows[i], a->rows[i], b->rows[i]);
        }  
}



void m5d_sub(m5d_t * c ,m5d_t  * a , m5d_t * b)
{
    if((a->nrows == b->nrows) && ( b->ncols == a->ncols))
    {
        int i, j;
        m5d_create(c, b->nrows, b->ncols);
        for( i = 0; i < a->nrows; i++)
        {
            for(j = 0; j < (a->width ); j++)
            {   
                vfd_sub(&c->rows[i][j], &a->rows[i][j], &b->rows[i][j]);
            }
        }
        
    }
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
                m5d_sub_i(&a->rows[i][j], &b->rows[i][j]);
            }
        }
        
    }


}


void *  m5d_combine3(vfd *table, vfd *input )
{
    vfd t, a, b, c;
    t.sign = t.middle =  t.units = 0;
    a = input[0];
    b = input[1];
    c = input[2];
    
    table[0] = t;
    table[1] = a;
    table[2] = b;
    table[4] = c;
    
    add_vfd(&t, &a, &b);
    table[3] = t;
    iadd_vfd(&t, &c);
    table[7] = t;
    m5d_sub_i(&t, &a);
    table[6] = t;
    add_vfd((table + 5), &a , &b);
    
    return 0;
    
}


void m5d_combine4(vfd *table, vfd *input )
{
    vfd t, a, b, c , d;
    t.sign = t.middle = t.units = 0;
    a = input[0];
    b = input[1];
    c = input[2];
    d = input[3];
    
    table[0] = t;
    table[1] = a;
    table[2] = b;
    table[4] = c;
    table[8] = d;
    
    add_vfd(&t, &c, &d);
    
    table[12] = t;
    
    add_vfd(&t,&b,&c);
    table[6] = t;
    iadd_vfd(&t,&d);
    table[14] = t;
    m5d_sub_i(&t,&c);
    table[10] = t;
    
    add_vfd(&t,&b,&c);
    table[3] = t;
    iadd_vfd(&t, &d);
    
    table[11] = t;
    iadd_vfd(&t, &c);
    table[15] = t;
    m5d_sub_i(&t, &d);
    table[7] = t;
    m5d_sub_i(&t, &b);
    table[5] = t;
    iadd_vfd(&t, &d);
    table[13] = t;
    m5d_sub_i(&t, &c);
    table[9] = t;
  
}


void m5d_combine5(vfd *table, vfd *input )
{
	int i;
    vfd e, *t4;
    
    m5d_combine4(table, input);
    e = input[4];
    t4 = table + 16;
    table[16] = e;
    
    for ( i = 1; i < 16 ; i ++ ) {
        add_vfd(t4 + i, table + i, &e);
    }

}

void m5d_combine6(vfd *table, vfd *input )

{
    vfd f, *t5;
    int i;
    m5d_combine5(table, input);
    f = input[5];
    t5 = (table + 32);
    table [32] = f;
    
    for (i = 1; i < 32; i++)
        add_vfd((t5 + i), (table + i), &f);
    
}

void m5d_combine7(vfd *table, vfd *input )
{
    vfd g, *t6;
    int i;
   
    m5d_combine6(table, input);
    g = input[6];
    t6 = (table+64);
    table[64] = g;
    
    for (i = 1; i < 64; i = i +1) {
        add_vfd((t6 + i), (table + i), &g );
    }
    
}


void m5d_combine8(vfd *table, vfd *input)
{
    vfd h, *t7;
    int i;
    
    m5d_combine7(table, input);
    h = input[7];
    t7 = (table+128);
    table[128] = h;
    
    for (i = 1; i < 128; i++)
        add_vfd((t7 + i), (table+i), &h);
}
 
 
 
/** 
GF(5), base case  still needs to be tuned for 
optimization and checked for correctness 
*/  
void m5d_mul_64(vfd **R, vfd **A, vfd **B)
{
    int i;
    vfd t1, t2, t3,  r1, r2, r3 , a;
    vec v1, v2, v3;
    
    vfd  tables6[9][64];
    vfd tables5[2][32];
    
    for (i = 0; i < 9; i ++)
    {
        m5d_combine6(&tables6[i][0], &(B [6*i][0]));
    }
      
    for (i = 0; i < 2; i ++)
    {
        m5d_combine5(&tables5[i][0], &(B[54 + (5 * i)][0]));
    }
  
    
    for (i = 0; i < 64; i ++  )//i from 0 <= i < 64
    {
        a = A[i][0];
        
  		v3 = a.sign;
        v2 = a.sign & a.middle;
    	v1 = v2 & a.units;
        r1 = tables6[0][v1&63];
        v1 >>= 6;
        r2 = tables6[0][v2&63];
        v2 >>= 6;
        t1 = tables6[1][v1&63]; iadd_vfd(&r1, &t1);v1 >>= 6;
        t2 = tables6[1][v2&63]; iadd_vfd(&r2, &t2); v2 >>= 6;
        t3 = tables6[1][v3&63]; iadd_vfd(&r3, &t3); v3 >>= 6;
        t1 = tables6[2][v1&63]; iadd_vfd(&r1, &t1); v1 >>= 6;
        t2 = tables6[2][v2&63]; iadd_vfd(&r2, &t2); v2 >>= 6;
        t3 = tables6[1][v3&63]; iadd_vfd(&r3, &t3); v3 >>= 6;
        t1 = tables6[3][v1&63]; iadd_vfd(&r1, &t1); v1 >>= 6;
        t2 = tables6[3][v2&63]; iadd_vfd(&r2, &t2); v2 >>= 6;
        t3 = tables6[1][v3&63]; iadd_vfd(&r3, &t3); v3 >>= 6;
        t1 = tables6[4][v1&63]; iadd_vfd(&r1, &t1); v1 >>= 6;
        t2 = tables6[4][v2&63]; iadd_vfd(&r2, &t2); v2 >>= 6;
        t3 = tables6[1][v3&63]; iadd_vfd(&r3, &t3); v3 >>= 6;
        t1 = tables6[5][v1&63]; iadd_vfd(&r1, &t1); v1 >>= 6;
        t2 = tables6[5][v2&63]; iadd_vfd(&r2, &t2); v2 >>= 6;
        t3 = tables6[1][v3&63]; iadd_vfd(&r3, &t3); v3 >>= 6;
        t1 = tables6[6][v1&63]; iadd_vfd(&r1, &t1); v1 >>= 6;
        t2 = tables6[6][v2&63]; iadd_vfd(&r2, &t2); v2 >>= 6;
        t3 = tables6[1][v3&63]; iadd_vfd(&r3, &t3); v3 >>= 6;
        t1 = tables6[7][v1&63]; iadd_vfd(&r1, &t1); v1 >>= 6;
        t2 = tables6[7][v2&63]; iadd_vfd(&r2, &t2); v2 >>= 6;
        t3 = tables6[1][v3&63]; iadd_vfd(&r3, &t3); v3 >>= 6;
        t1 = tables6[8][v1&63]; iadd_vfd(&r1, &t1); v1 >>= 6;
        t2 = tables6[8][v2&63]; iadd_vfd(&r2, &t2); v2 >>= 6;
        t3 = tables6[1][v3&63]; iadd_vfd(&r3, &t3); v3 >>= 6;
        t1 = tables5[0][v1&31]; iadd_vfd(&r1, &t1); v1 >>= 5;
        t2 = tables5[0][v2&31]; iadd_vfd(&r2, &t2); v2 >>= 5;
        t1 = tables5[1][v1&31]; iadd_vfd(&r1, &t1);
        t2 = tables5[1][v2&31]; iadd_vfd(&r2, &t2);
        t3 = tables5[1][v3&31]; iadd_vfd(&r3, &t3);
        
        iadd_vfd(&r1, &r2);
        m5d_add2_i(&r1, &r3);
        R[i][0] = r1;
    }   
}

//32 * 64,2048 bit, 256 byte matrix(slice) multiplication
void m5d_mul_32(vfd *R, vfd *A, vfd *B)
{
    long i;
    vfd t1, t2, t3, r1, r2, r3, a;
    long v1, v2, v3;
    
    vfd tables5[4][32];
    vfd tables4[3][16];
    for (i = 1; i < 4; i ++)
        
        m5d_combine5(tables5[i], B + 0 + 5*i);
    for (i = 0; i < 3; i++)
        m5d_combine4(tables4[i], B + 20 + 4*i);
    
    for (i = 0;i < 32; i++)
    {
        a = A[i];
 		v3 = a.sign;
        v2 = a.sign & a.middle;
    	v1 = v2 & a.units;
        t1 = tables5[0][v1&31]; v1 >>= 5;
        t2 = tables5[0][v2&31]; v2 >>= 5;
        t3 = tables5[0][v3&31]; v3 >>= 5;
        t1 = tables5[1][v1&31]; iadd_vfd(&r1, &t1); v1 >>= 5;
        t2 = tables5[1][v2&31]; iadd_vfd(&r2, &t2); v2 >>= 5;
        t3 = tables5[0][v3&31]; iadd_vfd(&r3, &t3); v3 >>= 5;
        t1 = tables5[2][v1&31]; iadd_vfd(&r1, &t1); v1 >>= 5;
        t2 = tables5[2][v2&31]; iadd_vfd(&r2, &t2); v2 >>= 5;
        t3 = tables5[0][v3&31]; iadd_vfd(&r3, &t3); v3 >>= 5;
        t1 = tables5[3][v1&31]; iadd_vfd(&r1, &t1); v1 >>= 5;
        t2 = tables5[3][v2&31]; iadd_vfd(&r2, &t2); v2 >>= 5;
        t3 = tables5[0][v3&31]; iadd_vfd(&r3, &t3); v3 >>= 5;
        t1 = tables4[0][v1&15]; iadd_vfd(&r1, &t1); v1 >>= 4;
        t2 = tables4[0][v2&15]; iadd_vfd(&r2, &t2); v2 >>= 4;
        t3 = tables5[0][v3&31]; iadd_vfd(&r3, &t3); v3 >>= 5;
        t1 = tables4[1][v1&15]; iadd_vfd(&r1, &t1); v1 >>= 4;
        t2 = tables4[1][v2&15]; iadd_vfd(&r2, &t2); v2 >>= 4;
        t3 = tables5[0][v3&31]; iadd_vfd(&r3, &t3); v3 >>= 5;
        t1 = tables4[2][v1&15]; iadd_vfd(&r1, &t1);
        t2 = tables4[2][v2&15]; iadd_vfd(&r2, &t2);
        t3 = tables4[2][v3&15]; iadd_vfd(&r3, &t3);
        // m5d_sub_i(&r1, &r2);
        iadd_vfd(&r1, &r2);
        m5d_add2_i(&r1, &r3);
        
        R[i] = r1;
    }
    
}

//16 * 64,1024 bit, 128 byte matrix(slice) multiplication
void m5d_mul_16(vfd *R, vfd *A, vfd *B)
{
    long i;
    vfd t1, t2, t3, r1, r2, r3, a;
    long v1, v2, v3;
    
    vfd tables4[4][16];
    for (i = 0; i < 4; i++)
        m5d_combine4(tables4[i], B + (4*i));
    for (i = 0;  i < 16; i++)
    {
        a = A[i];
  		v3 = a.sign;
        v2 = a.sign & a.middle;
    	v1 = v2 & a.units;
        r1 = tables4[0][v1&15]; v1 >>= 4;
        r2 = tables4[0][v2&15]; v2 >>= 4;
        r3 = tables4[0][v3&15]; v3 >>= 4;
        t1 = tables4[1][v1&15]; iadd_vfd(&r1, &t1); v1 >>= 4;
        t2 = tables4[1][v2&15]; iadd_vfd(&r2, &t2); v2 >>= 4;
        t3 = tables4[1][v3&15]; iadd_vfd(&r3, &t3); v3 >>= 4;
        t1 = tables4[2][v1&15]; iadd_vfd(&r1, &t1); v1 >>= 4;
        t2 = tables4[2][v2&15]; iadd_vfd(&r2, &t2); v2 >>= 4;
        t3 = tables4[1][v3&15]; iadd_vfd(&r3, &t3); v3 >>= 4;
        t1 = tables4[3][v1&15]; iadd_vfd(&r1, &t1);
        t2 = tables4[3][v2&15]; iadd_vfd(&r2, &t2);
        t3 = tables4[3][v3&15]; iadd_vfd(&r3, &t3);
    
    	iadd_vfd(&r1, &r2);
        m5d_add2_i(&r1, &r3);
            
        R[i] = r1;
    }
}

//8 * 64,512 bit, m1ri_word byte matrix(slice) multiplication
void m5d_mul_8(vfd *R, vfd *A, vfd *B)
{

    int i;
    vfd t1, t2, t3 ,  r1, r2, r3,a;
    vec v1, v2, v3;
    
    vfd tables4[2][16];
    for (i = 0; i < 2; i++)
        m5d_combine4(tables4[i], B + (4*i));
    for (i = 0; i < 8; i++)
    {
    	a = A[i];
     	v3 = a.sign;
        v2 = a.sign & a.middle;
    	v1 = v2 & a.units;
    	r1 = tables4[0][v1&15]; v1 >>= 4;
    	r2 = tables4[0][v2&15]; v2 >>= 4;
    	r3 = tables4[0][v3&15]; v3 >>= 4;
   	 	t1 = tables4[1][v1&15]; iadd_vfd(&r1, &t1);
    	t2 = tables4[1][v2&15]; iadd_vfd(&r2, &t2);
    	t3 = tables4[1][v3&15]; iadd_vfd(&r3, &t3);
    	// m5d_sub_i((&r1, &r2);
    	
    	iadd_vfd(&r1, &r2);
		m5d_add2_i(&r1, &r3);
        
        R[i] = r1;
    }
}

//4 * 64,256 bit, 32 byte matrix(slice) multiplication
void m5d_mul_4(vfd *R, vfd *A, vfd *B)
{
    int i;
    vfd r1, r2, r3,  a;
    vec v1, v2, v3;
    
    vfd table4[16];
    for (i = 0; i < 1; i++)
        m5d_combine4(table4, B + (4*i));
    for(i = 0; i < 4; i++)
    {
        a = A[i];
        v3 = a.sign;
        v2 = a.sign & a.middle;
    	v1 = v2 & a.units;
		v2  = a.sign & a.middle;
		v3 =  v2 & a.units; 
        r1 = table4[v1&15];
        r2 = table4[v2&15];
        r3 = table4[v3&15];
     	
     	iadd_vfd(&r1, &r2);
        m5d_add2_i(&r1, &r3);
        
        R[i] = r1;
    }
    
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
m5d_t * m5d_hadamard(m5d_t const * a, m5d_t const * b )
{
    m5d_t  * c = malloc(sizeof(m5d_t));
    if((a->nrows == b->nrows) && ( b->ncols == a->ncols))
    {
       
        int i, j;
        m5d_create(c, a->nrows , b->ncols);
        if(a->ncols < 256)
        { 
          for( i = 0; i < a->nrows; i++)
          {
            for(j = 0; j < (a->width ); j++)
            {  
               vfd_elem(c->rows[i] + j, a->rows[i]  + j, b->rows[i] + j);
            }  
          }

        }
     }   
        
    
    return c;

}
 
 void m5d_copy(m5d_t * a, m5d_t const *b)
{
  m5d_create(a, b->ncols, b->nrows);
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

