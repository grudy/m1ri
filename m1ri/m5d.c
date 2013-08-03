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
 
 m5d.c
 */

#include "m5d.h"


/*
 Read n bits from a sign portion of an element
 x = rows
 y = columns
 M = Matrix read from
 */

vec m5d_rm_bits(m5d_t *M, rci_t  x, rci_t  y, int  n) {
    

    
    wi_t  block = (y  ) / 64;
    
    int  spill = (y  % 64) + n - 64;
    
    vec bits;
    
    bits = (spill <= 0) ? M->rows[x][block].middle << -spill : (M->rows[x][block + 1].sign << (64 - spill)) | (M->rows[x][block].middle >> spill);
    
    
    return bits;
    
    
}


vec m5d_rs_bits(m5d_t *M, rci_t  x, rci_t  y, int  n) {
    
    
    
    wi_t  block = (y  ) / 64;
    
    int  spill = (y  % 64) + n - 64;
    
    vec bits;
    
    bits = (spill <= 0) ? M->rows[x][block].sign << -spill : (M->rows[x][block + 1].sign << (64 - spill)) | (M->rows[x][block].sign >> spill);
    
    
    return bits;
    
    
}
/*
 Read n bits from units
 x = rows
 y = columns
 M = Matrix read from
 */

vec m5d_ru_bits(m5d_t *M, rci_t  x, rci_t  y, int  n) {
    
    
    wi_t  block = (y  ) / 64;
    
    int  spill = (y  % 64) + n - 64;
    
    vec bits;
    
    bits = (spill <= 0) ? M->rows[x][block].units << -spill : (M->rows[x][block + 1].units<< (64 - spill)) | (M->rows[x][block].units>> spill);
  
    return bits;
    
    
}





/*
 Read n elements
 x = rows
 y = columns
 M = Matrix read from
 */

vfd m5d_read_elems(m5d_t *M, rci_t  x, rci_t  y, int  n) {
    
    
    
    
    wi_t  block = (y  ) / 64;
    
    int  spill = (y  % 64) + n - 64;
    
    vfd elem;
    
    elem.units = (spill <= 0) ? M->rows[x][block].units << -spill : ((M->rows[x][(block + 1)].units<< (64 - spill)) | (M->rows[x][block].units >> spill));
    
    elem.sign = (spill <= 0) ?  (M->rows[x][block].sign << -spill) : (M->rows[x][block + 1].sign << (64 - spill)) | (M->rows[x][block].sign>> spill);
    
    elem.middle = (spill <= 0) ?  (M->rows[x][block].middle << -spill) : (M->rows[x][block + 1].middle << (64 - spill)) | (M->rows[x][block].middle>> spill);
    
    elem.middle = (elem.middle >> (64 - n));
    
    elem.units = (elem.units >> (64 - n));
    
    elem.sign = (elem.sign >> (64 - n));
    //a->frows
    
    
    return elem;
    
    
}





/*
 Swap rows in a matrix;
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
    
    
    {
        
        
    }
    return 0;
}


/*
 
 */


//unfinished
void *  m5d_write_elem( m5d_t * M,rci_t x, rci_t y, vec s, vec u , vec m)
{
    
    
    
    wi_t  block = (y  ) / 64;
    
    int   spill = (y  % 64) - 63;
    
    
    
    s = ~(s == 0);
    u = ~(u == 0);
    
    
    M->rows[x][block].units  = (u == 0) ? (~(rightbit << -spill) &  (M->rows[x][block].units))  : ((u << (64 - spill)) | (M->rows[x][block].units));
    
    M->rows[x][block].sign  = (s == 0) ? (~(rightbit << -spill) &  (M->rows[x][block].units))  : ((u << (64 - spill)) | (M->rows[x][block].units));
    M->rows[x][block].middle  = (m == 0) ? (~(rightbit << -spill) &  (M->rows[x][block].units))  : ((u << (64 - spill)) | (M->rows[x][block].units));
    
    return 0;
    
    
}


/*
 
 */



vfd  * m5d_block_allocate(vfd * block, rci_t  nrows,  wi_t  width)
{
    
    
    block  = m1ri_malloc(nrows * width * sizeof(vfd) );
    
    
    
    return block;
    
    
    
}

/*
 
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

/*
 
 */

m5d_t m5d_create( m5d_t * a, rci_t nrows, rci_t ncols)
{
    
    
    a->ncols = ncols;
    a->nrows = nrows;
    a->width =  DN(ncols, m1ri_word);
    a->block = m5d_block_allocate(a->block,  a->nrows,    a->width);
    a->rows  = m5d_row_alloc(a->block, a->rows, a->width, a->nrows);
    a->flags = notwindowed;
    a->fcol = 0;
    return *a;
    
}

/*
 
 */

/*
 windows in m1ri_word rows * m1ri_word column incriments
 stvfd = the vfd or width offset from the base matrix
 strow = row offset in increments of 64
 sizecol  = cols * 64
 sizerow  = rows * 64
 */

m5d_t   m5d_window(m5d_t *c, rci_t strow, rci_t stvfd, rci_t sizerows, rci_t sizecols)
{
    
    
    m5d_t  submatrix;
    
    if((strow + sizerows) > c->width)
    {
        
        
        return submatrix;
    }
    
    
    
    
    if((stvfd + sizecols) > c->width)
    {
        
        
        return submatrix;
    }
    
    
    
    submatrix.nrows =   m1ri_word * sizecols;
    submatrix.ncols =  m1ri_word * sizecols;
    submatrix.flags = iswindowed;
    submatrix.width =  sizecols;
    submatrix.block = &c->block[(stvfd * stvfd)];
    submatrix.rows = m1ri_calloc(m1ri_word * sizerows ,  sizecols * sizeof(vfd *));
    submatrix.lblock = ( (sizerows +  strow)  ==  c->width)? c->lblock:  0;
    submatrix.fcol   = 0;
    submatrix.svfd = stvfd;
    
    
    
    int i;
    
    for(  i =   strow; i < (strow + (m1ri_word * sizerows)) ; i++)
    {
        
        submatrix.rows[i - strow] = c->rows[i];
        
    }
    
    return submatrix;
    
}

void   m5d_window_create(m5d_t *c, m5d_t * submatrix, rci_t strow, rci_t stvfd, rci_t sizerows, rci_t sizecols)
{
    
    
    submatrix->nrows =   m1ri_word * sizecols;
    submatrix->ncols =  m1ri_word * sizecols;
    submatrix->flags = iswindowed;
    submatrix->width = 1 * sizecols;
    submatrix->block = &c->block[(stvfd * stvfd)];
    submatrix->rows = m1ri_calloc(m1ri_word * sizerows ,  sizecols * sizeof(vfd *));
    submatrix->lblock = ( (sizerows +  strow)  ==  c->width)? c->lblock:  0;
    submatrix->fcol   = 0;
    submatrix->svfd = stvfd;
    
    
    int i;
    
    for(  i =   strow; i < strow + (m1ri_word * sizerows) ; i++)
    {
        
        submatrix->rows[i - strow] = c->rows[i];
        
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


/*
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
            l =  ((i - 1) * 64);
            j = i - 1;
            for ( k = 0 ; k < 64; k++)
            {
                
                a->rows[l][j].units = lbit[k];
                l++;
                
            }
            
            
            
            
            
            
        }
        
        
        
        
        if((a->ncols%64) != 0)
        {
            l = a->ncols %64;
            k = ((a->width -1) * 64);
            l = m1ri_word - l;
            for(i = 0; i < (64 - l); i++)
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

/*
 
 */


m5d_t   m5d_identity(m5d_t  *a, rci_t n)
{
    *a = m5d_create(a, n, n);
    *a = m5d_identity_set(a);
    
    return *a;
    
    
}


/*
 
 */




/*
 
 Releases a m5d_t into the wilderness.
 */



void m5d_free( m5d_t *  tofree)
{
    
    
    m1ri_free(tofree->rows);
    m1ri_free(tofree->block);
    
}


void addgf5(vfd * r, vfd * a, vfd * b)

{
    vec c, d, e, f, g, h, i, j, k, l, m, n ,o, p, q;
    c = b->units ^ a->units;
    d = b->middle ^ a->middle;
    e = b->sign ^ a->sign;
    f = d & c;
    g = f | b->middle;
    h = f ^ a->sign;
    i = h | e;
/**/r->sign = i;
    j = i ^ b->units;
    k = j ^ a->middle;
    l = k | c;
    m = l ^ e;
    n = m ^ g;
/**/r->middle = n;
    o = m | d;
    p = o ^ c;
    q = p^n;
    /**/r->sign = q;

    
   /* def add(a,b):
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
/**/r->units = 
    j = h | e;
    k = h ^ a->sign;
    l = k | b->middle;
    m = l & j;
    n = m ^ c;
    o = n | c;
    p = o ^ f;
    q = p ^ e;

    




}


void m5d_sub( vfd *r, vfd *a, vfd *b)

{
    vec c, d, e, f, g, h, i, j, k, l, m, n ,o, p, q;
    c = b->units ^ a->units;
    d = b->middle ^ a->middle;
    e = c ^ b->sign;
    f = e ^ b->middle;
    g = f ^ a->sign;
    h = g | d;
    i = h | c;
    j = i ^ c;
    k = j & a->sign;
    l = k | b->middle;
    m = l ^ g;
    n = m ^ a->middle;
    o = m | d;
    p = o ^ c;
    q = p ^ a->middle;
    q = r->units;
    
    /*
    def sub(a,b):
    c = b[0] ^ a[0]
    d = b[1] ^ a[1]
    e = c ^ b[2]
    f = e | b[0]
    g = f ^ a[2]
    h = g | d
    i = h | c
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





/********************************************
 matrix r = (direct sum matrix r + matrix x)
 ********************************************/
void iaddgf5(vfd *r,vfd *x)
{
  
   
    
    
    
    
}

/*
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
/*
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
/*
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
void isubgf5(vfd *r,vfd *x) 
{
    
    
    
}


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

void add_m1riw_gff5(vfd *R, vfd *A, vfd *B)
{
    int i;
    for (i = 0; i < (sizeof(vec)); i++ )
    {
       addgf5(  &R[i], &A[i], &B[i]);
    }
    
}










