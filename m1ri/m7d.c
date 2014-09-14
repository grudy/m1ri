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
        if ((a->ncols%M1RI_RADIX) == 0)
        {
            
            l = (a->width - 1) * M1RI_RADIX;
            for(i  = 0; i < M1RI_RADIX; i++)
                
            {
                a->rows[l][a->width -1].units = (rightbit)<<i;
                l++;
                
            }
            
        }
        
        
    }
 \
}




m7d_t *   m7d_identity(m7d_t * a,  rci_t n)
{

    a = m7d_create( n, n);
    m7d_set_ui(a, 1);
    return a;
}



/* unfinished */
void *  m7d_write_elem( m7d_t * M,rci_t x, rci_t y, vec s, vec m,  vec u )
{
    wi_t  block = (y  ) / M1RI_RADIX;
    int   spill =  (y  % M1RI_RADIX) ;
    M->rows[x][block].units  = (u == 0) ? (~(rightbit  << spill) &  (M->rows[x][block].units))  : ((rightbit << spill) | (M->rows[x][block].units));
    M->rows[x][block].middle  = (m == 0) ? (~(rightbit  << spill) &  (M->rows[x][block].middle))  : ((rightbit << spill) | (M->rows[x][block].middle));
    M->rows[x][block].sign  = (s == 0) ? (~(rightbit  << spill) &  (M->rows[x][block].sign))  : ((rightbit << spill) | (M->rows[x][block].sign));
    return 0;

}


/** 
 
 */



vtri  * m7d_block_allocate(vtri * block, rci_t  nrows,  wi_t  width)
{
    
    block  = m1ri_calloc((nrows * width) ,  sizeof(vtri) );
    return block;
}

/** 
 
 */




vtri ** m7d_row_alloc(vtri * block, vtri ** rows, wi_t width, rci_t nrows)
{
	int i;
    rows = m1ri_calloc( nrows ,  sizeof(vtri *));
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
    a->fcol = 0;
    
    
    return a;
    
}

/* Creates a window for m7d_t matrices   */

m7d_t  * m7d_init_window(const m7d_t *c,const rci_t strow, const rci_t svtri,const rci_t sizerows, const rci_t sizecols)
{
    
    
    m7d_t * submatrix = m1ri_malloc(sizeof(m7d_t));
    /** c->width should not be compared twice */
    if((strow + sizerows) > c->nrows)
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


static inline m7d_t  * m7d_init_window_unshackled(const m7d_t *c,const rci_t strow, const rci_t svtri,const rci_t sizerows, const rci_t sizecols)
{
    
    
    m7d_t * submatrix = m1ri_malloc(sizeof(m7d_t));
  
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
    submatrix->flags = iswindowed;
    
    
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
void  m7d_rand(m7d_t * a)
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

}



/** 
 
 Releases a m7d_t into the wilderness.
 */

void m7d_free( m7d_t *  a)
{ 		
    m1ri_free(a->rows);
    if(a->flags == notwindowed)
    {
	  m1ri_free(a->block); 	
	}
    
    m1ri_free(a);
    a == NULL;
    
}



void m7d_vtri_sub(vtri * r ,vtri * x, vtri * y)
{
   vtri temp;
   temp.units = ~y->units;
   temp.middle = ~y->middle;
   temp.sign = ~y->sign;
   add_vtri(r, x, &temp);
   
    

}






void reduce_vtri( vtri * a)
{
    vtri b = *a ;
    a->units  = b.units ^ (b.units  | b.sign | b.middle) ;/*   ) */
    a->middle  = b.middle ^ (b.units  | b.sign | b.middle) ;
    a->sign  = b.units ^ (b.units  | b.sign | b.middle) ;
}



void vtri_mul_2(vtri * a)
{
    vec temp = a->units;
    a->units = a->sign;
    a->sign = a->middle;
    a->middle = temp;
    
    
    
    
}


void vtri_mul_3(vtri * a)
{
    vec z = a->units| a->middle | a->sign;
    vec temp = a->units;
    a->units = a->middle ^ z;
    a->middle = a->sign ^ z;
    a->sign = z ^ temp;
    
}


void vtri_mul_4(vtri * a)
{
    vec temp;
    temp = a->units;
    a->units = a->middle;
    a->middle = a->sign;
    a->sign = temp;
    
    
}



void vtri_mul_5(vtri * a)
{
  
    vec z = a->units | a->middle | a->sign;
    vec temp = a->units;
    a->units =  a->sign ^ z;
    a->sign =   a->middle  ^ z;
    a->middle = temp ^ z;
   
    
    
}

void vtri_mul_6(vtri * a)
{
    vec z = a->units | a->middle | a->sign;
	a->units =  a->units ^ z;
    a->sign =   a->sign  ^ z;
    a->middle = a->middle ^ z;
    
     
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
m7d_t * m7d_sub(m7d_t * r,    m7d_t  const *x,  m7d_t const *y)
{
	
	if (r == NULL)
	{
		r= m7d_create(x->nrows, y->ncols);

	} 
	else if( (r->nrows != x->nrows || r->ncols != y->ncols)) 
	{
		m1ri_die("m7d_sub: Provided return matrix has wrong dimensions.\n");
    	
	
	}	
	
	int n , i;
  	if((x->nrows == y->nrows) && ( x->ncols == y->ncols))
  	{
  	
  	  
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

/**
  Checks if an m7d_t is equal to another.
*/
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
            if((a->rows[i][j].sign != b->rows[i][j].sign) || (a->rows[i][j].units != b->rows[i][j].units) ||(a->rows[i][j].middle != b->rows[i][j].middle) )
            {
                printf("row [%d][%d] not equal \n", i, j);
                return 0;
            }
             printf("row [%d][%d]  equal \n", i, j);
            
        }
    }
    return 1;
}


m7d_t *  m7d_copy_cutoff(m7d_t  * r, m7d_t  const * x)
{
		int i, s;
        for( i = 0; i < r->nrows; i++)
        {
        	 for( s = 0; s < r->width; s++)
        	 {
            r->rows[i][s] = x->rows[i][s];
            }
            
        }
	return r;
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
     	
        for( i = 0; i < a->nrows; i++)
        {
            for(j = 0; j < (a->width ); j++)
            {
                
        
                add_vtri(&c->rows[i][j], &a->rows[i][j], &b->rows[i][j]);    
            }
        }
        
    
    return c;
      
}



static inline void m7d_combine4(vtri *table, vtri  ** const input) 
{
    vtri  t,   a,  b,  c,  d;
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

    add_vtri(&t,&c,&d);
    table[12] = t;
    
    add_vtri(&t,&b,&c);
    table[6] = t;
    m7d_inc(&t,&d);
    table[14] = t;
    m7d_dec(&t,&c);
    table[10] = t;

    add_vtri(&t,&a,&b);
    table[3] = t;
    m7d_inc(&t,&d);
    
    table[11] = t;
    m7d_inc(&t,&c);
    table[15] = t;
    m7d_dec(&t,&d);
    
    table[7] = t;
    m7d_dec(&t,&b);
    table[5] = t;
    m7d_inc(&t,&d);
    table[13] = t;
	m7d_dec(&t,&c);
    table[9] = t;
}

static inline void m7d_combine5(vtri *table, vtri  **  const input) 
{
    vtri *e, *t4;
    int i;

    m7d_combine4(table, input);
    e = input[4];
    t4 = table+16;
    table[16] = *e;
    
    for(i=1;i<16;i++)
        add_vtri(t4 + i, table+i, e);
}
    
static inline void m7d_combine6(vtri *table, vtri  ** const input) 
{
    vtri *e, *t4;
    vtri * f, *t5;
    int i;
    
    m7d_combine4(table, input);
    e = input[4];
    t4 = table+16;
    table[16] = *e;

    f = input[5];
    t5 = table+32;
    table[32] = *f;
    
    for(i=1;i<16;i++)
        add_vtri(t4 + i, table+i, e);

    for(i=1;i<32;i++)
        add_vtri(t5 + i, table+i, f);

}



void m7d_mul_64(vtri **R, vtri  ** const A, vtri   ** const B)
{
   
    vtri t, r, rv, * a;
    vec v;
    int i;

    vtri tables6[4][64];
    vtri tables5[8][32];

    for(i=0;i<4;i++)
    {
        m7d_combine6(tables6[i], B + 6*i);
    }
    for(i=0;i<8;i++)
    {
        m7d_combine5(tables5[i], B + 24 + (5*i));
	}
    for(i=0;i<64;i++) 
    {
    /*
    	000 = 0    
		100  = 1
		010 = 2
		110 = 3	
		001 = 4
		101 = 5
		011 = 6
    
    */
    
    
    	/*   first part is ones */
        a = A[i];
        
		
		
		//one
		v = (a->units);

		
        
        
        r = tables6[0][v&63];                 v >>= 6;
        t = tables6[1][v&63]; m7d_inc(&r, &t); v >>= 6;
        t = tables6[2][v&63]; m7d_inc(&r, &t); v >>= 6;
        t = tables6[3][v&63]; m7d_inc(&r, &t); v >>= 6;
        
        t = tables5[0][v&31]; m7d_inc(&r, &t); v >>= 5;
        t = tables5[1][v&31]; m7d_inc(&r, &t); v >>= 5;
        t = tables5[2][v&31]; m7d_inc(&r, &t); v >>= 5;
        t = tables5[3][v&31]; m7d_inc(&r, &t); v >>= 5;
        t = tables5[4][v&31]; m7d_inc(&r, &t); v >>= 5;
        t = tables5[5][v&31]; m7d_inc(&r, &t); v >>= 5;
        t = tables5[6][v&31]; m7d_inc(&r, &t); v >>= 5;
        t = tables5[7][v&31]; m7d_inc(&r, &t);

      
      
      	//two
      	v =  a->middle;
		
		
		
        rv = tables6[0][v&63];                 v >>= 6;
        t = tables6[1][v&63]; m7d_inc(&rv, &t); v >>= 6;
        t = tables6[2][v&63]; m7d_inc(&rv, &t); v >>= 6;
        t = tables6[3][v&63]; m7d_inc(&rv, &t); v >>= 6;
        
        t = tables5[0][v&31]; m7d_inc(&rv, &t); v >>= 5;
        t = tables5[1][v&31]; m7d_inc(&rv, &t); v >>= 5;
        t = tables5[2][v&31]; m7d_inc(&rv, &t); v >>= 5;
        t = tables5[3][v&31]; m7d_inc(&rv, &t); v >>= 5;
        t = tables5[4][v&31]; m7d_inc(&rv, &t); v >>= 5;
        t = tables5[5][v&31]; m7d_inc(&rv, &t); v >>= 5;
        t = tables5[6][v&31]; m7d_inc(&rv, &t); v >>= 5;
        t = tables5[7][v&31]; m7d_inc(&rv, &t);

		

		vtri_mul_2(&rv);
		m7d_inc(&r, &rv);
		//three
		v = a->sign;


		
		rv = tables6[0][v&63];                 v >>= 6;
        t = tables6[1][v&63]; m7d_inc(&rv, &t); v >>= 6;
        t = tables6[2][v&63]; m7d_inc(&rv, &t); v >>= 6;
        t = tables6[3][v&63]; m7d_inc(&rv, &t); v >>= 6;
        
        t = tables5[0][v&31]; m7d_inc(&rv, &t); v >>= 5;
        t = tables5[1][v&31]; m7d_inc(&rv, &t); v >>= 5;
        t = tables5[2][v&31]; m7d_inc(&rv, &t); v >>= 5;
        t = tables5[3][v&31]; m7d_inc(&rv, &t); v >>= 5;
        t = tables5[4][v&31]; m7d_inc(&rv, &t); v >>= 5;
        t = tables5[5][v&31]; m7d_inc(&rv, &t); v >>= 5;
        t = tables5[6][v&31]; m7d_inc(&rv, &t); v >>= 5;
        t = tables5[7][v&31]; m7d_inc(&rv, &t);

    	
		
		vtri_mul_4(&rv);
		

    	m7d_inc(&r, &rv);
    	
    	
		R[i][0] = r;
        
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


m7_slice *  m7d_quarter(const  m7d_t * a)
{

 	m7_slice * c = m1ri_malloc(sizeof(m7_slice));
	c->row = m1ri_calloc( 4 , sizeof(m7d_t **));
	c->row[0] = m7d_init_window_unshackled(a,  0, 0 , a->nrows/128, a->ncols/128);
	c->row[1] = m7d_init_window_unshackled(a, 0, a->ncols/128 , a->nrows/128, a->ncols/128);   
    c->row[2] = m7d_init_window_unshackled(a, a->nrows/128, 0 , a->nrows/128, a->ncols/128);
	c->row[3] = m7d_init_window_unshackled(a, a->nrows/128,a->ncols/128,  a->nrows/128, a->ncols/128);
	 
	 
	 
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
  
    	for( i = 0; i < a->nrows; i++)
    	{
    		for(j = 0; j < (a->width ); j++)
    		{	  
        		vtri_elem(c->rows[i] + j, a->rows[i]  + j, b->rows[i] + j);

    		}  
    	}

    
    
 
        
    
    return c;

}


m7d_t *  m7d_copy(m7d_t * a, m7d_t const *b)
{
   if(a == NULL)
  {	
  	a = m7d_create( b->nrows, b->ncols);
  }
  
  if((a->ncols < b->ncols) || (a->nrows < b->nrows))
  {
  	m1ri_die("m7d_copy: Provided return matrix has wrong dimensions.\n");
  
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
         a_place =  rightbit <<   dif_a ;
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
		for(j = 0; j < a->width; j++)
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
		for(j = 0; j < a->width; j++)
		{
			
			a->rows[i][j].sign = 	b->rows[i][j].sign;
			a->rows[i][j].middle = 	b->rows[i][j].middle;		
			a->rows[i][j].units = 	b->rows[i][j].units;	
			vtri_mul_2(a->rows[i] + j);
		
		  
		
		}
	
	}


}


inline void m7d_mul_three(m7d_t *a, const m7d_t * b)
{
	int i, j;
	for(i = 0; i < a->nrows; i++)
	{
		for(j = 0; j < a->width; j++)
		{
			a->rows[i][j].sign = 	b->rows[i][j].sign;
			a->rows[i][j].middle = 	b->rows[i][j].middle;		
			a->rows[i][j].units = 	b->rows[i][j].units;	
			vtri_mul_3(a->rows[i] + j);	  
		
		}
	
	}


}


inline void m7d_mul_fourth(m7d_t *a, const m7d_t * b)
{
	int i, j;
	for(i = 0; i < a->nrows; i++)
	{
		for(j = 0; j < a->width; j++)
		{
			a->rows[i][j].sign = 	b->rows[i][j].sign;
			a->rows[i][j].middle = 	b->rows[i][j].middle;		
			a->rows[i][j].units = 	b->rows[i][j].units;	
		  	 vtri_mul_4(a->rows[i] + j);
			
		  
		
		}
	
	}


}


inline void m7d_mul_five(m7d_t *a, const m7d_t * b)
{
	int i, j;
	for(i = 0; i < a->nrows; i++)
	{
		for(j = 0; j < a->width; j++)
		{
		
		
			a->rows[i][j].sign = 	b->rows[i][j].sign;
			a->rows[i][j].middle = 	b->rows[i][j].middle;		
			a->rows[i][j].units = 	b->rows[i][j].units;	
			vtri_mul_5(a->rows[i] + j);
		
		  
		
		}
	
	}


}


inline void m7d_mul_six(m7d_t *a, const m7d_t * b)
{
	int i, j;
	for(i = 0; i < a->nrows; i++)
	{
		for(j = 0; j < a->width; j++)
		{
			a->rows[i][j].sign = 	b->rows[i][j].sign;
			a->rows[i][j].middle = 	b->rows[i][j].middle;		
			a->rows[i][j].units = 	b->rows[i][j].units;	
		  	 vtri_mul_6(a->rows[i] + j);
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

void m7d_add_i(m7d_t * x, m7d_t *y) 
{
	 int i, j;
        for( i = 0; i < x->nrows; i++)
        {
            for(j = 0; j < (x->width ); j++)
            {
            	m7d_inc(&x->rows[i][j], &y->rows[i][j]);
        	}   
        }
	

}

void m7d_sub_i(m7d_t * x, m7d_t *y) 
{
	 int i, j;
        for( i = 0; i < x->nrows; i++)
        {
            for(j = 0; j < (x->width ); j++)
            {
            	m7d_dec(&x->rows[i][j], &y->rows[i][j]);
        	}   
        }
	

}




static inline void sub_m7d_r(vtri  const *r,vtri   *y)
{
	vtri a;
	vtri b;
   	a.units = ~y->units;
   	a.middle = ~y->middle;
   	a.sign = ~y->sign;
   	
   	b.units = r->units;
   	b.middle = r->middle;
   	b.sign = r->sign;
   	
   	
   	m7d_inc( &b, &a);
   	y->units = b.units;
   	y->middle = b.middle;
   	y->sign  = b.sign;
   

         
}



void m7d_sub_r(m7d_t   *x , m7d_t   const *r)
{
  
	int n , i;
	for(i = 0; i < x->nrows; i++)
    {
    	for(n = 0; n < x->width; n++)
        {
        
       
		  sub_m7d_r(r->rows[i] + n,  x->rows[i] + n );
        }
    }
   



} 
void m7d_add_row(m7d_t *A, rci_t ar, const m7d_t *B, rci_t br, rci_t start_col)
{
	rci_t s_width = start_col/M1RI_RADIX;
	if(start_col%64 != 0)
	{
		add_vtri(A->rows[ar] + s_width, A->rows[ar]+ s_width, B->rows[ar] + s_width);
	
	
	}
		

}

