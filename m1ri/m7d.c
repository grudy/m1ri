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
/** 
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





/** 
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





/** 
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


/** 
 
 */

void m7d_set_ui(m7d_t * a, rci_t value )

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
 \
}




m7d_t *   m7d_identity( rci_t n)
{
	m7d_t  *a = m7d_create(n, n);
    a = m7d_create( n, n);
    m7d_set_ui(a, 1);
    return a;
}



/* unfinished */
void *  m7d_write_elem( m7d_t * M,rci_t x, rci_t y, vec s, vec m,  vec u )
{
    wi_t  block = (y  ) / M1RI_RADIX;
    int   spill =  (y  % M1RI_RADIX) ;
    M->rows[x][block].units  = (u == 0) ? (~(leftbit >> spill) &  (M->rows[x][block].units))  : ((leftbit >> spill) | (M->rows[x][block].units));
    M->rows[x][block].middle  = (m == 0) ? (~(leftbit  >> spill) &  (M->rows[x][block].middle))  : ((leftbit  >> spill) | (M->rows[x][block].middle));
    M->rows[x][block].sign  = (s == 0) ? (~(leftbit  >> spill) &  (M->rows[x][block].sign))  : ((leftbit  >> spill) | (M->rows[x][block].sign));
    return 0;

}


/** 
 
 */



vtri  * m7d_block_allocate(vtri * block, rci_t  nrows,  wi_t  width)
{
    
    block  = m1ri_malloc(nrows * width * sizeof(vtri) );
    return block;
}

/** 
 
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

/** 
 
 */

m7d_t * m7d_create( rci_t nrows, rci_t ncols)
{

	m7d_t * a = m1ri_malloc(sizeof(m7d_t)); 
    a->ncols = ncols;
    a->nrows = nrows;
    a->width = M1RI_DN(ncols, M1RI_RADIX);
    a->block = m7d_block_allocate(a->block,  a->nrows,    a->width);
    a->rows  = m7d_row_alloc(a->block, a->rows, a->width, a->nrows);
    a->flags = notwindowed;
    
    
    return a;
    
}

/* Creates a window for m7d_t matrices   */

m7d_t  * m7d_init_window(const m7d_t *c,const rci_t strow, const rci_t svtri,const rci_t sizerows, const rci_t sizecols)
{
    
    
    m7d_t * submatrix = m1ri_malloc(sizeof(m7d_t));
    /** c->width should not be compared twice */
    if((strow + sizerows) > c->width)
    {    
        return  0;
    }
    
    if((svtri + sizecols) > c->ncols)
    {
        return  0;
    }
    int i, f;
	f = strow * M1RI_RADIX;
    submatrix->nrows =   M1RI_RADIX * sizerows;
    submatrix->ncols =  M1RI_RADIX * sizecols;
    submatrix->flags = iswindowed;
    submatrix->width =  sizecols;
    submatrix->rows = m1ri_calloc(M1RI_RADIX * sizerows ,  sizecols * sizeof(vtri *));
    submatrix->lblock = ( (sizerows +  strow)  ==  c->width)? c->lblock:  0;
    submatrix->fcol   = 0;
    submatrix->svtri = svtri;
    
    
    
    for(  i =   f; i < (f + (M1RI_RADIX * sizerows)) ; i++)
    {
        submatrix->rows[i - f] = c->rows[i] + svtri;   
    }
    return submatrix;
}


 
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
/** Negates input of a vtri*/
void vtri_negate(vtri * a)
{
	a->units = !a->units;
	a->middle = !a->middle;
	a->sign = !a->sign;
}
/** Fills a matrix over GF(7) with random Variables*/
m7d_t  m7d_rand(m7d_t * a)
{
    int i;
    
    for( i = 0; i < (a->nrows * a->width); i++)
    {
        
        a->block[i].sign = m1ri_rand();
        a->block[i].middle = m1ri_rand() ;
        a->block[i].units = m1ri_rand();
        vec temp = (~(a->block[i].sign) & ~(a->block[i].middle) & ~(a->block[i].units));
   		a->block[i].sign   = a->block[i].sign | temp;
   		a->block[i].middle = a->block[i].middle | temp;
   		a->block[i].units  =  a->block[i].units | temp;
   		
    }
    return *a;
}



/** 
 
 Releases a m7d_t into the wilderness.
 */



void m7d_free( m7d_t *  tofree)
{
    m1ri_free(tofree->rows);
    m1ri_free(tofree->block);
    
}


void add_vtri(vtri * r, vtri * x, vtri * y)

{
    
    
    /*
    s3s2s1s0 ← add(a2a1a0, b2b1b0)
    r2r1r0 ← add(s2s1s0, s3)

    */
    vec s;
    vec t;
    
    r->units = x->units ^ y->units;

    s = (x->units & y->units);
    r->middle = s^ x->middle ^ y->middle;
    t = ((s) & (x->middle | y->middle)) | (x->middle & y->middle);
    r->sign = x->sign ^ y->sign ^ t;
    /* to here I know is right */
    
    
    s = ((t) & (x->sign | y->sign)) | (x->sign & y->sign);
    
    t = s & r->units;
    r->units = s ^ r->units;
    s= t & r->middle;
    r->middle = r->middle ^ t;
    r->sign = r->sign | s;

    
    
}

void m7d_vtri_sub(vtri * r ,vtri * x, vtri * y)
{
   vtri temp;
   temp.units = ~y->units;
   temp.middle = ~y->middle;
   temp.sign = ~y->sign;
   add_vtri(r, x, &temp);
   
    

}


void m7d_sub_i(vtri  *r, vtri *y)
{
	/** 
	Subtraction function
	*/

}


inline void iadd_vtri(vtri  *x, vtri *y)
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
    /** Cleanup*/
      

}


void reduce_vtri( vtri * a)
{
    vtri b = *a ;
    a->units  = b.units ^ (b.units  | b.sign | b.middle) ;/*   ) */
    a->middle  = b.middle ^ (b.units  | b.sign | b.middle) ;
    a->sign  = b.units ^ (b.units  | b.sign | b.middle) ;
}





vtri vtri_mul_2(vtri a)
{
    vec temp;
    temp = a.units;
    a.units = a.middle;
    a.middle = a.sign;
    a.sign = temp;
    return a;
    
}
vtri vtri_mul_3(vtri a)
{
  
    vec z = a.units | a.middle | a.sign;
    vec temp = a.units;
    a.units =  a.sign ^ z;
    a.sign =   a.middle  ^ z;
    a.middle = temp ^ z;
   
    return a;
    
}
vtri vtri_mul_4(vtri a)
{
    vec temp = a.units;
    a.units = a.sign;
    a.sign = a.middle;
    a.middle = temp;
    
    return a;
    
    
}
vtri vtri_mul_5(vtri a)
{
    vec z = a.units| a.middle | a.sign;
    vec temp = a.units;
    a.units = a.middle ^ z;
    a.middle = a.sign ^ z;
    a.sign = z ^ temp;
    return a;
    
}
vtri vtri_mul_6(vtri a)
{
    vec z = a.units | a.middle | a.sign;
	a.units =  a.units ^ z;
    a.sign =   a.sign  ^ z;
    a.middle = a.middle ^ z;
    
    return a;
     
}

/** Summing and multiplying result by two, for method of four russian*/
void m7d_add_2r(vtri *x, vtri * y)
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


	/** Optimize later*/

    s = x->units;
    x->units = x->middle;
    x->middle = x->sign;
    x->sign = s;

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
    
    /* Optimize the multiplication later */
    r.units = x->units;
    x->units = x->sign;
    x->sign = x->middle;
    x->middle = r.units;
}





/** ************
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
m7d_t * m7d_sub(m7d_t * r, const   m7d_t  *x, const m7d_t  *y)
{
	
	int n , i;
  	if((x->nrows == y->nrows) && ( x->ncols == y->ncols))
  	{
  	
  	  r = m7d_create( x->nrows , y->ncols);
  	  for(i = 0; i < x->nrows; i++)
    	{
        
        	for(n = 0; n < x->width; n++)
        	{
		    	m7d_vtri_sub(&r->rows[i][n], &x->rows[i][n], &y->rows[i][n]);
        	}
    	}
    }
    return r;
}

int m7d_equal(m7d_t const *a, m7d_t const *b)
{
	u_int64_t temp = (a->ncols%64 == 0)? 0:  ((leftbit >> ((a->ncols%64 ) - 1)) -1) ;
	temp = ~temp;
    if ((a->nrows != b->nrows)    || ( a->ncols != b->ncols)  )
    {
        return 0;
    }
   
  
    
    for(int i = 0; i < a->nrows; i++)
    {
        
        for(int j = 0; j < (b->width -1); j++)
        {
            if((a->rows[i][j].sign != b->rows[i][j].sign) || (a->rows[i][j].units != b->rows[i][j].units ) || (a->rows[i][j].middle != b->rows[i][j].middle ))
            {
                return 0;
            }
            
        }
    }
	
	for(int i = 0; i < a->nrows; i++)
    {
        
        for(int j =  b->width - 1; j < b->width; j++)
        {
            if(((a->rows[i][j].sign & temp )!= (b->rows[i][j].sign & temp)) ||
             ((a->rows[i][j].units & temp) != (b->rows[i][j].units & temp)) ||
              ((a->rows[i][j].middle & temp) != (b->rows[i][j].middle & temp)))
            {
        
                return 0;
            
          
			}            
        }
    }
	
     
       
    return 1;
}




void m7d_copy_cutoff(m7d_t  * r, m7d_t  const * x)
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
       /*  R[i][0] = add_m7dr(A[i][0], B[i][0]); */
    }

}

m7d_t * m7d_add(m7d_t * c,const   m7d_t *a, const m7d_t *b)
{
	if (c == NULL)
	{
		c = m7d_create(a->nrows, b->ncols);

	} 
	else if( (c->nrows != a->nrows || c->ncols != b->ncols)) 
	{
		m1ri_die("m7d_add: Provided return matrix has wrong dimensions.\n");
    	
	
	}	

    if((a->nrows != b->nrows) || ( b->ncols != a->ncols))
    {   
      m1ri_die("m7d_add: Input Matrices must have same dimension.\n");
    }

        int i, j;
     	c = m7d_create(b->nrows, b->ncols);
        for( i = 0; i < a->nrows; i++)
        {
            for(j = 0; j < (a->width ); j++)
            {
                
                /*  */
                add_vtri(&c->rows[i][j], &a->rows[i][j], &b->rows[i][j]);    
            }
        }
        
    
    return c;
      
}
void *  m7d_combine3(vtri *table, vtri *input )
{
    vtri t, a, b, c;
    t.sign = t.middle = t.units = 0;
    a = input[0];
    b = input[1];
    c = input[2];
    
    table[0] = t;
    table[1] = a;
    table[2] = b;
    table[4] = c;
    
    add_vtri(&t, &a, &b);
    table[3] = t;
    iadd_vtri(&t, &c);
    table[7] = t;
    m7d_sub_i(&t, &a);
    table[6] = t;
    
    add_vtri((table + 5), &a , &b);

    return 0;
    
}

void m7d_combine4(vtri *table, vtri *input )
{
    vtri t, a, b, c , d;
    t.sign = t.middle =  t.units = 0;
    a = input[0];
    b = input[1];
    c = input[2];
    d = input[3];
    
    table[0] = t;
    table[1] = a;
    table[2] = b;
    table[4] = c;
    table[8] = d;
    
    add_vtri(&t, &c, &d);
    
    table[12] = t;
    
    add_vtri(&t,&b,&c);
    table[6] = t;
    iadd_vtri(&t,&d);
    table[14] = t;
    m7d_sub_i(&t,&c);
    table[10] = t;
    
    add_vtri(&t,&b,&c);
    table[3] = t;
    iadd_vtri(&t, &d);
     
    table[11] = t;
    iadd_vtri(&t, &c);
    table[15] = t;
    m7d_sub_i(&t, &d);
    table[7] = t;
    m7d_sub_i(&t, &b);
    table[5] = t;
    iadd_vtri(&t, &d);
    table[13] = t;
    m7d_sub_i(&t, &c);
    table[9] = t;
    
}

void m7d_combine5(vtri *table, vtri *input )
{
	int i;
    vtri e, *t4;
    m7d_combine4(table, input);
    e = input[4];
    t4 = table + 16;
    table[16] = e;
    
    for ( i = 1; i < 16 ; i ++ ) {
        add_vtri(t4 + i, table + i, &e);
    }
     
}


void m7d_combine6(vtri *table, vtri *input )
{
    vtri f, *t5;
    int i;
    m7d_combine5(table, input);
    f = input[5];
    t5 = (table + 32);
    table [32] = f;
    
    for (i = 1; i < 32; i++)
        add_vtri((t5 + i), (table + i), &f);
    
}

void m7d_combine7(vtri *table, vtri *input )

{
    
    vtri g, *t6;
    int i; 
    m7d_combine6(table, input);
    g = input[6];
    t6 = (table+64);
    table[64] = g;
    
    for (i = 1; i < 64; i = i +1)
	{
        add_vtri((t6 + i), (table + i), &g );
	}
}


void m7d_combine8(vtri *table, vtri *input)
{
    vtri h, *t7;
    int i;
    m7d_combine7(table, input);
    h = input[7];
    t7 = (table+128);
    table[128] = h;
    for (i = 1; i < 128; i++)
    add_vtri((t7 + i), (table+i), &h);
}
   



/*  Passed here is Multiplication base case for GF(7)*/

void m7d_mul_64(vtri **R, vtri **A, vtri **B)
{
    int i;
    vtri t1, t2, t3,  r1, r2, r3,  a;
    vec v1, v2, v3;
    vtri  tables6[9][64];
    vtri tables5[2][32];
     
    for (i = 0; i < 9; i ++)
    {
        m7d_combine6(&tables6[i][0], &(B [6*i][0]));
    }
     
    for (i = 0; i < 2; i ++)
    {
        m7d_combine5(&tables5[i][0], &(B[54 + (5 * i)][0]));
    }
   
    for (i = 0; i < 64; i ++  )/* i from 0 <= i < 64 */
    {
        a = A[i][0];
 		v2 = a.middle;
        v3 = a.sign;
        v1 = a.units;
        r1 = tables6[0][v1&63];
        v1 >>= 6;
        r2 = tables6[0][v2&63];
        v2 >>= 6;
        t1 = tables6[1][v1&63]; iadd_vtri(&r1, &t1);v1  >>= 6;
        t2 = tables6[1][v2&63]; iadd_vtri(&r2, &t2); v2 >>= 6;
        t3 = tables6[1][v3&63]; iadd_vtri(&r3, &t3); v3 >>= 6;
        t1 = tables6[2][v1&63]; iadd_vtri(&r1, &t1); v1 >>= 6;
        t2 = tables6[2][v2&63]; iadd_vtri(&r2, &t2); v2 >>= 6;
        t2 = tables6[1][v2&63]; iadd_vtri(&r2, &t2); v2 >>= 6;
        t3 = tables6[1][v3&63]; iadd_vtri(&r3, &t3); v3 >>= 6;
        t1 = tables6[3][v1&63]; iadd_vtri(&r1, &t1); v1 >>= 6;
        t2 = tables6[3][v2&63]; iadd_vtri(&r2, &t2); v2 >>= 6;
        t2 = tables6[1][v2&63]; iadd_vtri(&r2, &t2); v2 >>= 6;
        t3 = tables6[1][v3&63]; iadd_vtri(&r3, &t3); v3 >>= 6;
        t1 = tables6[4][v1&63]; iadd_vtri(&r1, &t1); v1 >>= 6;
        t2 = tables6[4][v2&63]; iadd_vtri(&r2, &t2); v2 >>= 6;
        t2 = tables6[1][v2&63]; iadd_vtri(&r2, &t2); v2 >>= 6;
        t3 = tables6[1][v3&63]; iadd_vtri(&r3, &t3); v3 >>= 6;
        t1 = tables6[5][v1&63]; iadd_vtri(&r1, &t1); v1 >>= 6;
        t2 = tables6[5][v2&63]; iadd_vtri(&r2, &t2); v2 >>= 6;
        t2 = tables6[1][v2&63]; iadd_vtri(&r2, &t2); v2 >>= 6;
        t3 = tables6[1][v3&63]; iadd_vtri(&r3, &t3); v3 >>= 6;
        t1 = tables6[6][v1&63]; iadd_vtri(&r1, &t1); v1 >>= 6;
        t2 = tables6[6][v2&63]; iadd_vtri(&r2, &t2); v2 >>= 6;
        t2 = tables6[1][v2&63]; iadd_vtri(&r2, &t2); v2 >>= 6;
        t3 = tables6[1][v3&63]; iadd_vtri(&r3, &t3); v3 >>= 6;
        t1 = tables6[7][v1&63]; iadd_vtri(&r1, &t1); v1 >>= 6;
        t2 = tables6[7][v2&63]; iadd_vtri(&r2, &t2); v2 >>= 6;
        t2 = tables6[1][v2&63]; iadd_vtri(&r2, &t2); v2 >>= 6;
        t3 = tables6[1][v3&63]; iadd_vtri(&r3, &t3); v3 >>= 6;
        t1 = tables6[8][v1&63]; iadd_vtri(&r1, &t1); v1 >>= 6;
        t2 = tables6[8][v2&63]; iadd_vtri(&r2, &t2); v2 >>= 6;
        t2 = tables6[1][v2&63]; iadd_vtri(&r2, &t2); v2 >>= 6;
        t3 = tables6[1][v3&63]; iadd_vtri(&r3, &t3); v3 >>= 6;
        t1 = tables5[0][v1&31]; iadd_vtri(&r1, &t1); v1 >>= 5;
        t2 = tables5[0][v2&31]; iadd_vtri(&r2, &t2); v2 >>= 5;
        t2 = tables6[1][v2&63]; iadd_vtri(&r2, &t2); v2 >>= 6;
        t3 = tables6[1][v3&63]; iadd_vtri(&r3, &t3); v3 >>= 6;
        t1 = tables5[1][v1&31]; iadd_vtri(&r1, &t1);
        t2 = tables5[1][v2&31]; iadd_vtri(&r2, &t2);
        t3 = tables6[1][v3&63]; iadd_vtri(&r3, &t3); 
       
        m7d_add_2r(&r1, &r2);
		m7d_add_4r(&r1 ,&r3);
        
        R[i][0] = r1;
       /*  */ 
    }
    
}

/* 32 * 64,2048 bit, 256 byte matrix(slice) multiplication */
void m7d_mul_32(vtri *R, vtri *A, vtri *B)
{
    long i;
    vtri t1, t2, t3,  r1, r2, r3,  a;
    long v1, v2, v3;
    
    vtri tables5[4][32];
    vtri tables4[3][16];
    for (i = 1; i < 4; i ++)
        
        m7d_combine5(tables5[i], B + 0 + 5*i);
    for (i = 0; i < 3; i++)
        m7d_combine4(tables4[i], B + 20 + 4*i);
    
    for (i = 0;i < 32; i++)
    {
        
        a = A[i];
        v3 = a.sign;
        v2 = a.middle;
        v1 = a.units ^ v2;
        t1 = tables5[0][v1&31]; v1 >>= 5;
        t2 = tables5[0][v2&31]; v2 >>= 5;
        t3 = tables5[0][v2&31]; v3 >>= 5;
        t1 = tables5[1][v1&31]; iadd_vtri(&r1, &t1); v1 >>= 5;
        t2 = tables5[1][v2&31]; iadd_vtri(&r2, &t2); v2 >>= 5;
        t3 = tables5[1][v3&31]; iadd_vtri(&r3, &t3); v2 >>= 5;
        t1 = tables5[2][v1&31]; iadd_vtri(&r1, &t1); v1 >>= 5;
        t2 = tables5[2][v2&31]; iadd_vtri(&r2, &t2); v2 >>= 5;
        t3 = tables5[2][v3&31]; iadd_vtri(&r3, &t3); v2 >>= 5;
        t1 = tables5[3][v1&31]; iadd_vtri(&r1, &t1); v1 >>= 5;
        t2 = tables5[3][v2&31]; iadd_vtri(&r2, &t2); v2 >>= 5;
        t3 = tables5[3][v3&31]; iadd_vtri(&r3, &t3); v2 >>= 5;
        t1 = tables4[0][v1&15]; iadd_vtri(&r1, &t1); v1 >>= 4;
        t2 = tables4[0][v2&15]; iadd_vtri(&r2, &t2); v2 >>= 4;
        t3 = tables4[0][v3&15]; iadd_vtri(&r3, &t3); v2 >>= 4;
        t1 = tables4[1][v1&15]; iadd_vtri(&r1, &t1); v1 >>= 4;
        t2 = tables4[1][v2&15]; iadd_vtri(&r2, &t2); v2 >>= 4;
        t3 = tables4[1][v3&15]; iadd_vtri(&r3, &t3); v2 >>= 4;
        t1 = tables4[2][v1&15]; iadd_vtri(&r1, &t1);
        t2 = tables4[2][v2&15]; iadd_vtri(&r2, &t2);
        t3 = tables4[3][v3&15]; iadd_vtri(&r3, &t3);
    
    	m7d_add_2r(&r1, &r2);
		m7d_add_4r(&r1 ,&r3);
        R[i] = r1;
    }
    
}

/**
	16 * 64,1024 bit, 128 byte matrix(slice) multiplication

*/
void m7d_mul_16(vtri *R, vtri *A, vtri *B)
{
    long i;
    vtri t1, t2, t3,  r1, r2, r3,  a;
    long v1, v2, v3;
    
    vtri tables4[4][16];
    for (i = 0; i < 4; i++)
        m7d_combine4(tables4[i], B + (4*i));
    for (i = 0;  i < 16; i++)
    {
        a = A[i];
        v2 = a.middle;
        v3 = a.sign;
        v1 = a.units;
        r1 = tables4[0][v1&15]; v1 >>= 4;
        r2 = tables4[0][v2&15]; v2 >>= 4;
        r3 = tables4[0][v3&15]; v2 >>= 4;
        t1 = tables4[1][v1&15]; iadd_vtri(&r1, &t1); v1 >>= 4;
        t2 = tables4[1][v2&15]; iadd_vtri(&r2, &t2); v2 >>= 4;
        t3 = tables4[1][v3&15]; iadd_vtri(&r3, &t3); v3 >>= 4;
        t1 = tables4[2][v1&15]; iadd_vtri(&r1, &t1); v1 >>= 4;
        t2 = tables4[2][v2&15]; iadd_vtri(&r2, &t2); v2 >>= 4;
        t3 = tables4[1][v3&15]; iadd_vtri(&r3, &t3); v3 >>= 4;
        t1 = tables4[3][v1&15]; iadd_vtri(&r1, &t1);
        t2 = tables4[3][v2&15]; iadd_vtri(&r2, &t2);
        t3 = tables4[1][v3&15]; iadd_vtri(&r3, &t3);
    	
     	m7d_add_2r(&r1, &r2);
		m7d_add_4r(&r1 ,&r3);
        R[i] = r1;
    }
}

/**
8 * 64,512 bit, m1ri_word byte matrix(slice) multiplication
*/
void m7d_mul_8(vtri *R, vtri *A, vtri *B)

{
    int i;
    vtri t1, t2, t3,  r1, r2, r3,  a;
    vec v1, v2, v3;
    
    vtri tables4[2][16];
    for (i = 0; i < 2; i++)
        m7d_combine4(tables4[i], B + (4*i));
    for (i = 0; i < 8; i++)
    {
        a = A[i];
    v3 = a.sign;
    v2 = a.middle;
    v1 = a.units;
    r1 = tables4[0][v1&15]; v1 >>= 4;
    r2 = tables4[0][v2&15]; v2 >>= 4;
    r3 = tables4[0][v3&15]; v3 >>= 4;
    t1 = tables4[1][v1&15]; iadd_vtri(&r1, &t1);
    t2 = tables4[1][v2&15]; iadd_vtri(&r2, &t2);
    t3 = tables4[1][v3&15]; iadd_vtri(&r3, &t3);
    
    m7d_add_2r(&r1, &r2);
	m7d_add_4r(&r1 ,&r3);
    R[i] = r1;
    }
}


/**
 4 * 64,256 bit, 32 byte matrix(slice) multiplication
*/
void m7d_mul_4(vtri *R, vtri *A, vtri *B)
{
    int i;
    vtri r1, r2, r3 ,  a;
    vec v1, v2, v3;
    
    vtri table4[16];
    for (i = 0; i < 1; i++)
        m7d_combine4(table4, B + (4*i));
    for(i = 0; i < 4; i++)
    {
        a = A[i];
        v2 = a.middle;
        v3 = a.sign;
        v1 = a.units;
        r1 = table4[v1&15];
        r2 = table4[v2&15];
        r3 = table4[v3&15];
    
    	m7d_add_2r(&r1, &r2);
		m7d_add_4r(&r1 ,&r3);
        
        R[i] = r1;
    }
    
}

/**

*/
m7d_t  * m7_blockslice_allocate(rci_t  nrows,  wi_t  width)
{
    m7d_t * block  = m1ri_calloc(nrows * width ,  sizeof(m7d_t  ) );
    return block;
}

m7d_t ** m7_rowslice_allocate(m7d_t * block,  wi_t width, rci_t nrows)
{
	int i;
    m7d_t ** rows = m1ri_malloc( nrows * width * sizeof(m7d_t **));
    for ( i = 0; i <  nrows;  i++ )
    {
        rows[i]  = &block[i * width];
    };
    return rows;
}



static inline  vtri *  m7d_transpose_vtri(vtri  **a, vtri  **b  )
{
    int i, x;
    vtri temp;
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

void  m7d_slices(m7_slice *  c, m7d_t * a, wi_t slicesize)
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
 	c->block = m7_blockslice_allocate(c->nrows,   c->width);
    c->row = m7_rowslice_allocate(c->block , c->width, c->nrows);
    r = 0 ;
     
    for ( i = 0; i <  l;  i = i + slicesize)
    {       
    	for( f = 0; f <colroundeddown ; f++)
        {
        	c->row[(r * f) + f] = m7d_init_window(a,i , (f * slicesize), slicesize, slicesize);
        }
        
        if(extracols > 0)
        {
        	c->row[(r * f) + f] = m7d_init_window(a,i , (f * slicesize), slicesize, extracols);
		}
        r++;
        
   	}
    
    if(extrarows >0 )
    {
		for( f = 0; f <colroundeddown ; f++)
        {
           c->row[(r * f) + f] = m7d_init_window(a, i , (f * slicesize), extrarows, slicesize);
        }

    	if(extracols > 0)
    	{
           c->row[(r * f) + f] = m7d_init_window(a, i , (f * slicesize), extrarows, extracols);
     	}
    }
   
}

m7d_t * m7d_transpose_sliced(m7d_t * a)
{
    int x, y;
    m7d_t *  c;
    c = m7d_create( a->ncols, a->nrows);
    m7_slice * b, *d;
    d = malloc(sizeof(m7_slice));
    b = malloc(sizeof(m7_slice));
    m7d_slices(b, a, 1);
    m7d_slices(d, c, 1);
    for (x = 0; x < b->nrows; x++)
     {
        for (y = 0; y < b->ncols; y ++)
         {
         
         	m7d_transpose_vtri(b->row[x][y].rows, d->row[y][x].rows);  
         	 
        }
    }
    
    return c;
}




m7_slice *  m7d_quarter(const  m7d_t * a)
{
	m7_slice * c = m1ri_malloc(sizeof(m7_slice));
	c->block = m7_blockslice_allocate( 2,   2);
    c->row = m7_rowslice_allocate(c->block,   2, 2);
    c->row[0] = m7d_init_window(a,  0, 0 , a->nrows/128, a->ncols/128);
	c->row[1] = m7d_init_window(a, 0, a->ncols/128 , a->nrows/128, a->ncols/128);   
    c->row[2] = m7d_init_window(a, a->nrows/128, 0 , a->nrows/128, a->ncols/128);
	c->row[3] = m7d_init_window(a,  a->nrows/128,a->ncols/128,  a->nrows/128, a->ncols/128);
    return c;
}




/*

Logic Friday results 
c->units = a->units a->mid b->units' b->mid' b->sign + a->units a->mid' a->sign' b->units b->mid  + a->units' a->mid a->sign' b->units b->sign + a->units' a->mid' a->sign b->mid  + a->units a->sign' b->units b->mid'  + a->units' a->mid b->units' b->sign + a->sign b->units' b->mid b->sign' + a->units b->units b->mid' b->sign';
c->mid = a->units a->mid' a->sign' b->mid b->sign + a->mid a->sign b->units b->mid' b->sign' + a->units a->mid' a->sign' b->units b->mid  + a->units' a->mid a->sign' b->units b->sign + a->mid' a->sign b->units' b->sign + a->units' a->sign b->mid' b->sign + a->mid a->sign' b->units b->sign' + a->units b->units' b->mid b->sign';
c->sign = a->mid' a->sign b->units  + a->units a->mid' a->sign' b->mid b->sign + a->mid a->sign b->units b->mid' b->sign' + a->units a->mid b->units' b->mid' b->sign + a->mid a->sign' b->units' b->mid  + a->units a->mid' b->mid' b->sign + a->units' a->mid b->mid b->sign';

*/

/**
  Inline elementwise multiplication of a single set of 64 bits
*/
static inline void vtri_elem(vtri * c, vtri const * a, vtri const * b)
{

  vtri one, two, three, four;
  one.units = a->units & b->units;
  two.units = a->sign & b->middle;
  three.units = a->middle & b->sign;
  one.middle = a->middle & b->units;
  two.middle = a->units & b->middle;
  three.middle =  a->sign & b->sign;
  one.sign = a->units & b->sign;
  two.sign = a->sign & b->units;
  three.sign = a->middle & b->middle;        
  /* three.middle =  */
  
vec temp = (~(one.sign) & ~(one.middle) & ~(one.units));
   		one.sign   = one.sign | temp;
   		one.middle = one.middle | temp;
   		one.units  =  one.units | temp;
   		
   		
temp = (~(two.sign) & ~(two.middle) & ~(two.units));
   		two.sign   = two.sign | temp;
   		two.middle = two.middle | temp;
   		two.units  =  two.units | temp;
   		
   		
 temp = (~(three.sign) & ~(three.middle) & ~(three.units));
   		three.sign   = three.sign | temp;
   		three.middle = three.middle | temp;
   		three.units  =  three.units | temp;
   		
  add_vtri(&four, &one, &two); 		
  add_vtri(c, &four, &three); 		  		
   		  		


}




m7d_t * m7d_hadamard(m7d_t * c, m7d_t const * a, m7d_t const * b )
{


    if (c == NULL)
	{
		c = m7d_create(a->nrows, b->ncols);

	} 
	else if( (c->nrows != a->nrows || c->ncols != b->ncols)) 
	{
		m1ri_die("m7d_hadamard: Provided return matrix has wrong dimensions.\n");	
	}
    if((a->nrows != b->nrows) || ( b->ncols != a->ncols))
    {
       
      m1ri_die("m7d_hadamard: Input Matrices must have same dimension.\n");
    }

    int i, j;
    if(a->ncols < 256)
    { 
    	for( i = 0; i < a->nrows; i++)
    	{
    		for(j = 0; j < (a->width ); j++)
    		{	  
        		vtri_elem(c->rows[i] + j, a->rows[i]  + j, b->rows[i] + j);

    		}  
    	}

    }
    
 
        
    
    return c;

}


void m7d_copy(m7d_t * a, m7d_t const *b)
{
  a = m7d_create(b->ncols, b->nrows);
  for(int i = 0; i < a->nrows; i++)
  {
    for(int j = 0; j < b->ncols; j++)
    {
    
      a->rows[i][j] = b->rows[i][j];
    
    }
    
     a->lblock = b->lblock; /*   first block pointed to in a window */
     a->fcol = b->fcol;  /* /column offset of first block */
     a->flags = b->flags;
  
  }
  
  


}


void  m7d_colswap(m7d_t *M, rci_t col_a, rci_t col_b)
{
    if((M->ncols >= (col_a ) && (M->ncols >= col_b)))
    {
        int i;
        vec block_a, block_b, dif_a, dif_b, a_place, b_place; 
        vtri tempa, tempb;
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


void m7d_colswap_capped_row(m7d_t *M, rci_t col_a, rci_t col_b, rci_t start_row)
{
   if((M->ncols >= (col_a ) && (M->ncols >= col_b)))
    {
        int i;
        vec block_a, block_b, dif_a, dif_b, a_place, b_place; 
        vtri tempa, tempb;
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



void  m7d_transpose(m7d_t   * a)
{

   
  int x, y;
     m7d_t * c;
    c = m7d_create( a->ncols, a->nrows);
    m7_slice * b, *d;
    d = malloc(sizeof(m7_slice));
    b = malloc(sizeof(m7_slice));
    m7d_slices(b, a, 1);
    m7d_slices(d, c, 1);
    for (x = 0; x < b->nrows; x++) {
        for (y = 0; y < b->ncols; y ++) {
         m7d_transpose_vtri(b->row[x][y].rows, d->row[y][x].rows);
            
        }
    }


   
}


int m7d_is_zero(const m7d_t *A)
{

	if(A == NULL)
	{
	
		m1ri_die("m7d_is_zero: A cannot be null!\n");

	
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




inline void m7d_mul_zero(m7d_t * a)
{
	
	int i, j;
	for(i = 0; i < a->nrows; i++)
	{
		for(j = 0; j < a->ncols; j++)
		{
			a->rows[i][j].units = 0;
			a->rows[i][j].middle = 0;
			a->rows[i][j].sign  = 0;
		  
		
		}
	
	}

}



inline void m7d_mul_two(m7d_t *a, const m7d_t * b)
{
	
	int i, j;
	for(i = 0; i < a->nrows; i++)
	{
		for(j = 0; j < a->ncols; j++)
		{
		
		
			a->rows[i][j] = vtri_mul_2(b->rows[i][j]);
		
		  
		
		}
	
	}


}


inline void m7d_mul_three(m7d_t *a, const m7d_t * b)
{
	int i, j;
	for(i = 0; i < a->nrows; i++)
	{
		for(j = 0; j < a->ncols; j++)
		{
		
			a->rows[i][j] = vtri_mul_3(b->rows[i][j]);		  
		
		}
	
	}


}


inline void m7d_mul_fourth(m7d_t *a, const m7d_t * b)
{
	int i, j;
	for(i = 0; i < a->nrows; i++)
	{
		for(j = 0; j < a->ncols; j++)
		{

		  	a->rows[i][j] = vtri_mul_4(b->rows[i][j]);
			
		  
		
		}
	
	}


}


inline void m7d_mul_five(m7d_t *a, const m7d_t * b)
{
	int i, j;
	for(i = 0; i < a->nrows; i++)
	{
		for(j = 0; j < a->ncols; j++)
		{
			a->rows[i][j] = vtri_mul_5(b->rows[i][j]);
		
		  
		
		}
	
	}


}


inline void m7d_mul_six(m7d_t *a, const m7d_t * b)
{
	int i, j;
	for(i = 0; i < a->nrows; i++)
	{
		for(j = 0; j < a->ncols; j++)
		{
					
		  	a->rows[i][j] = vtri_mul_6(b->rows[i][j]);
		}
	
	}


}

m7d_t *m7d_mul_scalar(m7d_t *C, const long a, const m7d_t *B)
{

	if(C == NULL)
    { 
      C = m7d_create( B->ncols, B->nrows);
    }
    
    else if(C->nrows != B->nrows || C->ncols !=  B->ncols)
    {
    	m1ri_die("m7d_mul_scalar: C has wrong dimensions!\n");
    
    }
    
	long m = a%7;
	switch(m)
	{
		case 0: m7d_mul_zero(C);
		break ;
  		case 2: m7d_mul_two(C, B);
		break ;
  		case 3: m7d_mul_three(C, B);
		break ;
  		case 4: m7d_mul_fourth(C, B);
		break ;
  		case 5: m7d_mul_five(C, B);
		break ;
  		case 6: m7d_mul_six(C, B);
		break ;
  	}	
  
  return C;
}


void m7d_add_row(m7d_t *A, rci_t ar, const m7d_t *B, rci_t br, rci_t start_col)
{
	rci_t s_width = startcol/M1RI_RADIX;
	if(startcol%64 != 0)
	{
		vbg temp = A->row[ar][s_width].units + 
	
	
	}
		

}

