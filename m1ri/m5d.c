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
void *  m5d_write_elem( m5d_t * M,rci_t x, rci_t y, vec s, vec u )
{
    
    
    
    wi_t  block = (y  ) / 64;
    
    int   spill = (y  % 64) - 63;
    
    
    
    s = ~(s == 0);
    u = ~(u == 0);
    
    
    M->rows[x][block].units  = (u == 0) ? (~(rightbit << -spill) &  (M->rows[x][block].units))  : ((u << (64 - spill)) | (M->rows[x][block].units));
    
    M->rows[x][block].sign  = (s == 0) ? (~(rightbit << -spill) &  (M->rows[x][block].units))  : ((u << (64 - spill)) | (M->rows[x][block].units));
    M->rows[x][block].middle  = (s == 0) ? (~(rightbit << -spill) &  (M->rows[x][block].units))  : ((u << (64 - spill)) | (M->rows[x][block].units));
    
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
    
    
    
    
    rows = m1ri_malloc( nrows * width * sizeof(vfd *));
    
    
    for (int i = 0; i <  nrows;  i++ )
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
    a->width =  RU64(ncols);
    a->block = m5d_block_allocate(a->block,  a->nrows,    a->width);
    a->rows  = m5d_row_alloc(a->block, a->rows, a->width, a->nrows);
    a->flags = notwindowed;
    a->fcol = 0;
    return *a;
    
}

/*
 
 */

/*
 windows in 64 rows * 64 column incriments
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
    
    for(int i = 0; i < (a->nrows * a->width); i++)
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
            l = 64 - l;
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


void addgf5(vfd * r, vfd * x, vfd * y)

{
    S A;
    A.s0 =  x->sign ^ y->sign;
    A.s1 =  (x->middle ^ y->middle) ^ A.s0;
    A.s2 = (x->units ^ y->units)  ^ A.s1;
    A.s3 = (x->units & y->units);
   
    *r = fold5(A.s3, A.s2, A.s1, A.s0);
    
    
    
}

vfd m5d_add_r(vfd  *x, vfd *y)
{
   S A;
    A.s0 =  x->sign ^ y->sign;
    A.s1 =  (x->middle ^ y->middle) ^ A.s0;
    A.s2 = (x->units ^ y->units)  ^ A.s1;
    A.s3 = (x->units & y->units);
    
    vfd r = fold5(A.s3, A.s2, A.s1, A.s0);
    return r;
    
}

void m5d_sub( vfd *r, vfd *x, vfd *y)    

{
    
    
    
}

vfd m5d_sub_r(vfd x, vfd y)

{
    vfd r;
    return r;
}




/********************************************
 matrix r = (direct sum matrix r + matrix x)
 ********************************************/
void iaddgf5(vfd *r,vfd *x)
{
    
    S A;
    A.s0 =  x->sign ^ r->sign;
    A.s1 =  (x->middle ^ r->middle) ^ A.s0;
    A.s2 = (x->units ^ r->units)  ^ A.s1;
    A.s3 = (x->units & r->units);
    
     *r = fold5(A.s3, A.s2, A.s1, A.s0);
    
    
    
    
    
    
}

/*
 
 def v5mul2(a,b,c):
 return c&b, c, c^a
 
 def v5mul3(a,b,c):
 x = c ^ b
 return x, x|a, b
 
 def v5mul4(a,b,c):
 x = c ^ a
 y = y & c
 return y, x, y^b
 
 
 0 = 000
 1 = 001
 2 = 011
 3 = 101
 4 = 111
*/
vfd m5d_mul2(vfd a)
{
    vec temp = a.sign;
    a.sign = a.sign ^ temp;
    a.units = a.middle & temp;
    a.middle = temp;
    
    
    return a;

}

vfd m5d_mul3(vfd a)
{
    vec temp = a.middle ^ a.sign;
    a.sign = a.middle;
    a.middle = a.units | temp;
    a.units = temp;
    return a;
}

vfd m5d_mul4(vfd a)
{
    vec x = a.units ^ a.sign;
    a.units = x & a.sign;
    a.middle = x;
    
    return a;
    
}


void isubgf5(vfd *r,vfd *x)  //matrix r = (matrix r - matrix x)
{
    
    
    
}



void  m5d_mul( vfd *r, vfd *x, vfd *y)             //multiply matrix x by y assinging the output to r
{
    
}




//return the value of the matrix multiplied


m5d_t m5d_hadamard(m5d_t * a, m5d_t * b)
{
    m5d_t h_5;
    return h_5;
}



vfd m5d_mul_i(vfd x, vfd y)
{
    vfd r;
    r.units = x.units & y.units;
    r.sign  = (y.sign ^ x.sign) & (r.units);
    
    return r;
    
}







void sub_64gf5(vfd *R, vfd *A, vfd *B)
{
    int i;
    for (i = 0; i < (sizeof(vec)); i++ )
    {
        R[i] = m5d_sub_r(A[i], B[i]);
    }
}


void add_64gf5(vfd *R, vfd *A, vfd *B)
{
    int i;
    for (i = 0; i < (sizeof(vec)); i++ )
    {
        R[i] = m5d_add_r(&A[i], &B[i]);
    }
    
}










