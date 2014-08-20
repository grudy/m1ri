
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
 
 m1ri_strassen.c
 */
 
 

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "m1ri_strassen.h"
#include <math.h>
#include <stdlib.h>

#if __M1RI_HAVE_OPENMP
#include <omp.h>
#endif

/**
	Multiplies m3d_t matrices in squares.  
	Not to be used directly, but called by m3d_strassen
*/


/**
	Recursive Matrix Multiplication over GF(3), on a square matrix.
*/
static inline void m3d_mul_naive_square(m3d_t *c, const m3d_t *a, const m3d_t *b)
{
	
	m3d_t  * x1, * x2; 
  	m3_slice *  a_slice, *  b_slice, *  c_slice;

   	a_slice = m3d_quarter( a);
   /* m3d_print(a_slice->row[0][0]); */
    b_slice = m3d_quarter( b);
    c_slice = m3d_quarter(c);
   
    if((c_slice->row[0]->ncols) > M1RI_RADIX)
    {
       
    	x1 = m3d_create( c_slice->row[0]->nrows, c_slice->row[0]->ncols);		    
     	x2 = m3d_create( c_slice->row[0]->nrows, c_slice->row[0]->ncols);	
     	    	   
    	m3d_mul_naive_square(x1, a_slice->row[0], b_slice->row[0]);
    	m3d_mul_naive_square(x2, a_slice->row[1], b_slice->row[2]);
	    
	    c_slice->row[0] = m3d_add(c_slice->row[0],  x1, x2) ; 
		m3d_mul_naive_square(x1, a_slice->row[0], b_slice->row[1]);
       	m3d_mul_naive_square(x2, a_slice->row[1], b_slice->row[3]);
    	c_slice->row[1] = m3d_add(c_slice->row[1],  x1, x2) ;
	    m3d_mul_naive_square(x1, a_slice->row[2], b_slice->row[0]);
      	m3d_mul_naive_square(x2, a_slice->row[3], b_slice->row[2]);
       	c_slice->row[2] = m3d_add(c_slice->row[2], x1, x2); 
	    m3d_mul_naive_square(x1, a_slice->row[2], b_slice->row[1]);
		m3d_mul_naive_square(x2, a_slice->row[3], b_slice->row[3]); 
 	    c_slice->row[3] = m3d_add(c_slice->row[3], x1, x2) ;
 	 	m3d_free(x1);
		m3d_free(x2);
		
	
    }
   
    else if((c->ncols ) == (M1RI_RADIX  << 1))
    {
    	x1 = m3d_create(M1RI_RADIX,M1RI_RADIX);		    
		x2 = m3d_create( M1RI_RADIX, M1RI_RADIX);
   	 	m3d_mul_64(x1->rows[0], a_slice->row[0]->rows[0], b_slice->row[0]->rows[0]);
    	m3d_mul_64(x2->rows[0], a_slice->row[1]->rows[0], b_slice->row[2]->rows[0]);
  	    m3d_add_64(c_slice->row[0]->rows[0], x1->rows[0], x2->rows[0]) ;
        m3d_mul_64(x1->rows[0],  a_slice->row[0]->rows[0], b_slice->row[1]->rows[0] );
        m3d_mul_64(x2->rows[0], a_slice->row[1]->rows[0], b_slice->row[3]->rows[0]);
		m3d_add_64(c_slice->row[1]->rows[0], x1->rows[0], x2->rows[0]) ; 
       	m3d_mul_64(x1->rows[0], a_slice->row[2]->rows[0], b_slice->row[0]->rows[0]);
       	m3d_mul_64(x2->rows[0], a_slice->row[3]->rows[0], b_slice->row[2]->rows[0]);
     	m3d_add_64(c_slice->row[2]->rows[0], x1->rows[0], x2->rows[0]); 
        m3d_mul_64(x1->rows[0], a_slice->row[2]->rows[0], b_slice->row[1]->rows[0]);
        m3d_mul_64(x2->rows[0], a_slice->row[3]->rows[0], b_slice->row[3]->rows[0]); 
		m3d_add_64(c_slice->row[3]->rows[0], x1->rows[0], x2->rows[0]); 
		m3d_free(x1);
		m3d_free(x2);
	
    }

}


static inline void m3d_qrt_mul(m3d_t * c,const m3d_t *   a,const m3d_t *   b )
{
  	m3d_t * x1;
    m3d_t * x2;
	 	
    
	m3_slice *  a_slice, *  b_slice, *  c_slice;     	
 
    a_slice = m3d_quarter( a);
   	b_slice = m3d_quarter(b);
   	c_slice = m3d_quarter( c);
	
	x1 = m3d_create( c->ncols/2, c->nrows/2);
    x2 = m3d_create(c->ncols/2, c->nrows/2);
	if((c->ncols) > (M1RI_RADIX << 1))
    {
     
       
        x1 = m3d_sub(x1, a_slice->row[0], a_slice->row[2]);      // 1 
        
    	
        x2 = m3d_sub(x2, b_slice->row[1],b_slice->row[1]);      // 2 
        
        
        
    	m3d_qrt_mul(c_slice->row[2], x1, x2);      // 3 
        x1 = m3d_add(x1, a_slice->row[2],a_slice->row[1]);      // 4 
        x2 = m3d_sub(x2, b_slice->row[1],b_slice->row[0]);      // 5 
        
    	
        m3d_qrt_mul(c_slice->row[3], x1, x2);        // 6 
        x1 = m3d_sub(x1, x1,a_slice->row[0]);    // 7 
        x2 = m3d_sub(x2, b_slice->row[1],x2);      // 8 
        m3d_qrt_mul(c_slice->row[1],x1,x2);     // 9 
        x1 = m3d_sub(x1, a_slice->row[1],x1);        // 10 
        m3d_qrt_mul(c_slice->row[0],x1,b_slice->row[1]);       // 11 
        m3d_qrt_mul( x1 , a_slice->row[1], b_slice->row[1]);      // 12 
        c_slice->row[1] = m3d_add(c_slice->row[1], x1 , c_slice->row[1]);       // 13 
        c_slice->row[2] = m3d_add(c_slice->row[2], c_slice->row[1] , c_slice->row[2]);       // 14 
        c_slice->row[1] = m3d_add(c_slice->row[1], c_slice->row[1] , c_slice->row[3]);       // 15 
        c_slice->row[3] = m3d_add(c_slice->row[3], c_slice->row[2] , c_slice->row[3]);        // 16 
        c_slice->row[3] = m3d_add(c_slice->row[3], c_slice->row[2] , c_slice->row[3]);      // 17 
        x2 = m3d_sub(x2, x2, b_slice->row[2]);                // 18 
        m3d_qrt_mul(c_slice->row[2], a_slice->row[1], x2);                // 19 
        c_slice->row[2] = m3d_sub(c_slice->row[2], c_slice->row[2], c_slice->row[0]);      // 20 
		
        x1 = m3d_sub(x1 ,  a_slice->row[0], a_slice->row[2]);      // 1 
        x2 = m3d_sub(x2, b_slice->row[1],b_slice->row[1]);      // 2 
        m3d_qrt_mul(c_slice->row[2], x1, x2);      // 3 
        x1 = m3d_add(x1, a_slice->row[2],a_slice->row[1]);      // 4 
        x2 = m3d_sub(x2, b_slice->row[1],b_slice->row[0]);      // 5 
        m3d_qrt_mul(c_slice->row[3], x1, x2);        // 6 
        x1 = m3d_sub(x1, x1,a_slice->row[0]);    // 7 
        x2 = m3d_sub(x2, b_slice->row[1],x2);      // 8 
        m3d_qrt_mul(c_slice->row[1],x1,x2);     // 9 
        x1 = m3d_sub(x1, a_slice->row[1],x1);        // 10 
        m3d_qrt_mul(c_slice->row[0],x1,b_slice->row[1]);       // 11 
        m3d_qrt_mul( x1 , a_slice->row[1], b_slice->row[1]);      // 12 
        c_slice->row[1] = m3d_add(c_slice->row[1], x1 , c_slice->row[1]);       // 13 
        c_slice->row[2] = m3d_add(c_slice->row[2], c_slice->row[1] , c_slice->row[2]);       // 14 
        c_slice->row[1] = m3d_add(c_slice->row[1], c_slice->row[1] , c_slice->row[3]);       // 15 
        c_slice->row[3] = m3d_add(c_slice->row[3], c_slice->row[2] , c_slice->row[3]);        // 16 
        c_slice->row[3] = m3d_add(c_slice->row[3], c_slice->row[2] , c_slice->row[3]);      // 17 
        x2 = m3d_sub(x2, x2, b_slice->row[2]);                // 18 
        m3d_qrt_mul(c_slice->row[2], a_slice->row[1], x2);                // 19 
        c_slice->row[2] = m3d_sub( c_slice->row[2], c_slice->row[2], c_slice->row[0]);      // 20 

        m3d_qrt_mul(c_slice->row[0], a_slice->row[1], b_slice->row[2]);
        c_slice->row[0] = m3d_add(c_slice->row[0], x1,c_slice->row[0] );
    	
       	
       	
    }
    
    else if((c->ncols ) == (M1RI_RADIX  << 1))
    {
		m3d_sub_64(x1->rows, a_slice->row[0]->rows, a_slice->row[2]->rows); //   1 * /
        m3d_sub_64(x2->rows,b_slice->row[3]->rows,b_slice->row[1]->rows) ; //   2  
        m3d_mul_64(c_slice->row[2]->rows[0], x1->rows[0], x2->rows[0]); //   3         
        m3d_add_64(x1->rows[0],a_slice->row[2]->rows[0],a_slice->row[3]->rows[0]) ; //   4 
        
        m3d_sub_64(x2->rows,b_slice->row[1]->rows,b_slice->row[0]->rows) ; //   5 

     	m3d_mul_64(c_slice->row[3]->rows[0], x1->rows[0], x2->rows[0]); //     6 
        m3d_sub_64(x1->rows,x1->rows,a_slice->row[0]->rows) ; // 7 
        m3d_sub_64(x2->rows,b_slice->row[3]->rows,x2->rows); //   8 
        m3d_mul_64(c_slice->row[1]->rows[0],x1->rows[0],x2->rows[0]); //  9 
        m3d_sub_64(x1->rows,a_slice->row[1]->rows,x1->rows); //     10 

        m3d_mul_64(c_slice->row[0]->rows[0],x1->rows[0],b_slice->row[3]->rows[0]); //    11 
       
		
  	    m3d_mul_64(x1->rows[0], a_slice->row[3]->rows[0], b_slice->row[3]->rows[0]); //   12 
  	    
  	    
        m3d_add_64(c_slice->row[1]->rows[0],x1->rows[0] , c_slice->row[1]->rows[0]) ; //    13 
        m3d_add_64(c_slice->row[2]->rows[0],c_slice->row[1]->rows[0] , c_slice->row[2]->rows[0]) ; //    14 
        m3d_add_64(c_slice->row[1]->rows[0],c_slice->row[1]->rows[0] , c_slice->row[3]->rows[0]) ; //    15 
        m3d_add_64(c_slice->row[3]->rows[0],c_slice->row[2]->rows[0] , c_slice->row[3]->rows[0]) ; //     16 
        m3d_add_64(c_slice->row[3]->rows[0],c_slice->row[2]->rows[0] , c_slice->row[3]->rows[0]) ; //   17 
    	m3d_sub_64(x2->rows, x2->rows, b_slice->row[2]->rows) ; //             18 
        m3d_mul_64(c_slice->row[2]->rows[0], a_slice->row[3]->rows[0], x2->rows[0]); //             19 
        m3d_sub_64(c_slice->row[2]->rows, c_slice->row[2]->rows,c_slice->row[0]->rows); //   20 
        m3d_mul_64(c_slice->row[0]->rows[0], a_slice->row[1]->rows[0],b_slice->row[2]->rows[0]); //
        m3d_add_64(c_slice->row[0]->rows[0], x1->rows[0],c_slice->row[0]->rows[0]) ; // 
 			
  
		
    	
    }
    
    else if(c->ncols  == M1RI_RADIX )
    {
    	 m3d_mul_64(c->rows[0], a->rows[0], b->rows[0]);
    
    } 
    
    
    m1ri_free(a_slice);
    m1ri_free(b_slice);
    m1ri_free(c_slice);
        m3d_free(x1);
    	m3d_free(x2);
}

/**
  Strassen  algorithm on an m3d_t	
*/
m3d_t *  m3d_strassen(m3d_t *c,const m3d_t  *a,const m3d_t   *b)
{
	if (c == NULL)
	{
		c = m3d_create(a->nrows, b->ncols);

	} 

	else 
	{
	if (c->nrows != a->nrows || c->ncols != b->ncols) 
	{
		m1ri_die("m3d_strassen: Provided return matrix has wrong dimensions.\n");
    	}
	
	}
	
    if(a->ncols == b->nrows)
    {
 
      	/*  These hold the padded matrix sizes */
    	
	u_int32_t  arcr, acbr, bccc;
	arcr = a->nrows;
	acbr = a->ncols;
	bccc = b->ncols;

    
	arcr =  powerof2(arcr);
	acbr =  powerof2(acbr);
	bccc =  powerof2(bccc);
	
	c = m3d_create( a->nrows, b->ncols); 
	int lasta, lastb, lastboth;
	lasta = 64 - a->nrows%64;
	lastb = 64 -  b->ncols%64;  
	lastboth = 64 - a->nrows; 

	if((arcr != a->nrows) || (acbr != a->ncols) || (bccc) != (b->ncols))
	{
		m3d_t * padded_a,   * padded_b , * padded_c;
	
	
		padded_a = m3d_create(arcr, acbr);
		padded_b = m3d_create(acbr, bccc);
		padded_c = m3d_create(arcr, bccc);
		padded_a = m3d_copy(padded_a, a);
		/*padded_b = m3d_copy(padded_b, b);
		padded_c = m3d_copy(padded_c, c);
	
		m3d_qrt_mul(padded_c, padded_a, padded_b); 
		c  = m3d_copy_cutoff(c, padded_c);
		m3d_free(padded_a);
		m3d_free(padded_b);
		m3d_free(padded_c);
		*/
		
	}
	
	
	else
	{
	
		m3d_qrt_mul(c, a, b); 
	}
	

          
    }
    
    return c;
    
}

/**
	This handles the arithmetic of m5d_strassen
*/

void m5d_qrt_mul(m5d_t * c,const m5d_t *   a, const m5d_t *   b )
{
  	m5d_t * x1;
    m5d_t * x2;
    	
    
	m5_slice  * a_slice, * b_slice,*  c_slice;
    	
    	
    	
    	x1 = m5d_create( c->nrows/2, c->ncols/2);
    	x2 = m5d_create( c->nrows/2, c->ncols/2);
    	a_slice = m5d_quarter( a);
   	b_slice = m5d_quarter(b);
   	c_slice = m5d_quarter( c);
   	

	if((c->ncols) > (M1RI_RADIX << 1))
    {
     
       {

        x1 = m5d_sub(x1, a_slice->row[0], a_slice->row[2]);  /* 1 */
        x2 = m5d_sub(x2, b_slice->row[1],b_slice->row[1]);  /* 2 */
        m5d_qrt_mul(c_slice->row[2], x1, x2);  /* 3 */
        x1 = m5d_add(x1, a_slice->row[2],a_slice->row[1]);  /* 4 */
        x2 = m5d_sub(x2, b_slice->row[1],b_slice->row[0]);  /* 5 */
        m5d_qrt_mul(c_slice->row[3], x1, x2);    /* 6 */
        x1 = m5d_sub(x1, x1,a_slice->row[0]);/* 7 */
        x2 = m5d_sub(x2, b_slice->row[1],x2);  /* 8 */
        m5d_qrt_mul(c_slice->row[1],x1,x2); /* 9 */
        x1 = m5d_sub(x1, a_slice->row[1],x1);    /* 10 */
        m5d_qrt_mul(c_slice->row[0],x1,b_slice->row[1]);   /* 11 */
        m5d_qrt_mul( x1 , a_slice->row[1], b_slice->row[1]);  /* 12 */
        c_slice->row[1] = m5d_add(c_slice->row[1], x1 , c_slice->row[1]);   /* 13 */
        c_slice->row[2] = m5d_add(c_slice->row[2], c_slice->row[1] , c_slice->row[2]);   /* 14 */
        c_slice->row[1] = m5d_add(c_slice->row[1], c_slice->row[1] , c_slice->row[3]);   /* 15 */
        c_slice->row[3] = m5d_add(c_slice->row[3] , c_slice->row[2] , c_slice->row[3]);    /* 16 */
        c_slice->row[3] = m5d_add(c_slice->row[3] , c_slice->row[2] , c_slice->row[3]);  /* 17 */
        x2 = m5d_sub(x2, x2, b_slice->row[2]);            /* 18 */
        m5d_qrt_mul(c_slice->row[2], a_slice->row[1], x2);            /* 19 */
        c_slice->row[2] = m5d_sub(c_slice->row[2] , c_slice->row[2], c_slice->row[0]);  /* 20 */

        x1 = m5d_sub(x1,  a_slice->row[0], a_slice->row[2]);  /* 1 */
        x2 = m5d_sub(x2, b_slice->row[1],b_slice->row[1]);  /* 2 */
        m5d_qrt_mul(c_slice->row[2], x1, x2);  /* 3 */
        x1 = m5d_add(x1 ,a_slice->row[2],a_slice->row[1]);  /* 4 */
        x2 = m5d_sub(x2, b_slice->row[1],b_slice->row[0]);  /* 5 */
        m5d_qrt_mul(c_slice->row[3], x1, x2);    /* 6 */
        x1 = m5d_sub(x1, x1,a_slice->row[0]);/* 7 */
        x2 = m5d_sub(x2, b_slice->row[1],x2);  /* 8 */
        m5d_qrt_mul(c_slice->row[1],x1,x2); /* 9 */
        x1 = m5d_sub(x1, a_slice->row[1],x1);    /* 10 */
        m5d_qrt_mul(c_slice->row[0],x1,b_slice->row[1]);   /* 11 */
        m5d_qrt_mul( x1 , a_slice->row[1], b_slice->row[1]);  /* 12 */
        c_slice->row[1] = m5d_add(c_slice->row[1] , x1 , c_slice->row[1]);   /* 13 */
        c_slice->row[2] = m5d_add(c_slice->row[2] , c_slice->row[1] , c_slice->row[2]);   /* 14 */
        c_slice->row[1] = m5d_add(c_slice->row[1] , c_slice->row[1] , c_slice->row[3]);   /* 15 */
        c_slice->row[3] = m5d_add(c_slice->row[3], c_slice->row[2] , c_slice->row[3]);    /* 16 */
        c_slice->row[3] = m5d_add(c_slice->row[3], c_slice->row[2] , c_slice->row[3]);  /* 17 */
        x2 = m5d_sub(x2, x2, b_slice->row[2]);            /* 18 */
        m5d_qrt_mul(c_slice->row[2], a_slice->row[1], x2);            /* 19 */
        c_slice->row[2] = m5d_sub(c_slice->row[2],  c_slice->row[2], c_slice->row[0]);  /* 20 */
        m5d_qrt_mul(c_slice->row[0], a_slice->row[1], b_slice->row[2]);
        c_slice->row[0] = m5d_add(c_slice->row[0], x1,c_slice->row[0] );
        m5d_free(x1);
    	m5d_free(x2);
       }	
    }
    
    else if((c->ncols ) == (M1RI_RADIX  << 1))
    {
	
		m5d_sub_64(x1->rows, a_slice->row[0]->rows, a_slice->row[2]->rows);  /* 1 */
        m5d_sub_64(x2->rows,b_slice->row[3]->rows,b_slice->row[1]->rows) ;  /* 2 */
        m5d_mul_64(c_slice->row[2]->rows, x1->rows, x2->rows);  /* 3 */
        m5d_add_64(x1->rows,a_slice->row[2]->rows,a_slice->row[3]->rows) ;  /* 4 */
        m5d_sub_64(x2->rows,b_slice->row[1]->rows,b_slice->row[0]->rows) ;  /* 5 */

     	m5d_mul_64(c_slice->row[3]->rows, x1->rows, x2->rows);    /* 6 */
        m5d_sub_64(x1->rows,x1->rows,a_slice->row[0]->rows) ;/* 7 */
        m5d_sub_64(x2->rows,b_slice->row[3]->rows,x2->rows);  /* 8 */
        m5d_mul_64(c_slice->row[1]->rows,x1->rows,x2->rows); /* 9 */
        m5d_sub_64(x1->rows,a_slice->row[1]->rows,x1->rows);    /* 10 */

        m5d_mul_64(c_slice->row[0]->rows,x1->rows,b_slice->row[3]->rows);   /* 11 */
        m5d_mul_64(x1->rows, a_slice->row[3]->rows, b_slice->row[3]->rows);  /* 12 */
        m5d_add_64(c_slice->row[1]->rows,x1->rows , c_slice->row[1]->rows) ;   /* 13 */
        m5d_add_64(c_slice->row[2]->rows,c_slice->row[1]->rows , c_slice->row[2]->rows) ;   /* 14 */
        m5d_add_64(c_slice->row[1]->rows,c_slice->row[1]->rows , c_slice->row[3]->rows) ;   /* 15 */
        m5d_add_64(c_slice->row[3]->rows,c_slice->row[2]->rows , c_slice->row[3]->rows) ;    /* 16 */
        m5d_add_64(c_slice->row[3]->rows,c_slice->row[2]->rows , c_slice->row[3]->rows) ;  /* 17 */

    	m5d_sub_64(x2->rows, x2->rows, b_slice->row[2]->rows) ;            /* 18 */
        m5d_mul_64(c_slice->row[2]->rows, a_slice->row[3]->rows, x2->rows);            /* 19 */
        m5d_sub_64(c_slice->row[2]->rows, c_slice->row[2]->rows,c_slice->row[0]->rows);  /* 20 */
        m5d_mul_64(c_slice->row[0]->rows, a_slice->row[1]->rows,b_slice->row[2]->rows);
        m5d_add_64(c_slice->row[0]->rows, x1->rows,c_slice->row[0]->rows) ; 
 
	   


        m5d_free(x1);
    	m5d_free(x2);
    }
    
  
}


/**
	 m5d_strassen
*/

m5d_t *  m5d_strassen(m5d_t *c, m5d_t const  *a, const  m5d_t   *b)
{
    
	if (c == NULL)
	{
	c = m5d_create(a->nrows, b->ncols);

	} 
	
	else 
	{
	if (c->nrows != a->nrows || c->ncols != b->ncols) 
	{
		m1ri_die("m5d_strassen: Provided return matrix has wrong dimensions.\n");
    	}
	
	}
	
    if(a->ncols == b->nrows)
    {
 
      	/*  These hold the padded matrix sizes */
	u_int32_t  arcr, acbr, bccc;
	arcr = a->nrows;
	acbr = a->ncols;
	bccc = b->ncols;

    
	arcr =  powerof2(arcr);
	acbr =  powerof2(acbr);
	bccc =  powerof2(bccc);
	
	c = m5d_create( a->nrows, b->ncols); 
	int lasta, lastb, lastboth;
	lasta = 64 - a->nrows%64;
	lastb = 64 -  b->ncols%64;  
	lastboth = 64 - a->nrows; 

	if((arcr != a->nrows) || (acbr != a->ncols) || (bccc) != (b->ncols))
	{
	m5d_t * padded_a,   * padded_b , * padded_c;
	padded_a = m5d_create(arcr, acbr);
	padded_b = m5d_create(acbr, bccc);
	padded_a = m5d_create(arcr, bccc);
	m5d_copy(padded_a, a);
	m5d_copy(padded_b, b);
	m5d_copy(padded_c, c);
	
	
	m5d_qrt_mul(padded_c, padded_a, padded_b); 

	m5d_free(padded_a);
	m5d_free(padded_b);
	m5d_free(padded_c);
	}
	
	
	else
	{
	m5d_qrt_mul(c, a, b); 
	}
	

          
    }
    return c;
    
    
}




/**
	This handles the arithmetic of m7d_strassen 
	
*/


void m7d_qrt_mul(m7d_t * c,const m7d_t *   a,const m7d_t *   b )
{
  	m7d_t * x1;
    	m7d_t * x2;
    	
    
m7_slice  * a_slice, * b_slice,*  c_slice;
    	
    	
    	
    	x1 = m7d_create( c->nrows/2, c->ncols/2);
    	x2 = m7d_create( c->nrows/2, c->ncols/2);
    	a_slice = m7d_quarter( a);
   	b_slice = m7d_quarter(b);
   	c_slice = m7d_quarter( c);
   	

	if((c->ncols) > (M1RI_RADIX << 1))
    {
     
       {

        x1 = m7d_sub(x1, a_slice->row[0], a_slice->row[2]);  /* 1 */
        x2 = m7d_sub(x2, b_slice->row[1],b_slice->row[1]);  /* 2 */
        m7d_qrt_mul(c_slice->row[2], x1, x2);  /* 3 */
        x1 = m7d_add(x1, a_slice->row[2],a_slice->row[1]);  /* 4 */
        x2 = m7d_sub(x2, b_slice->row[1],b_slice->row[0]);  /* 5 */
        m7d_qrt_mul(c_slice->row[3], x1, x2);    /* 6 */
        x1 = m7d_sub(x1, x1,a_slice->row[0]);/* 7 */
        x2 = m7d_sub(x2,  b_slice->row[1],x2);  /* 8 */
        m7d_qrt_mul(c_slice->row[1],x1,x2); /* 9 */
        x1 = m7d_sub(x1, a_slice->row[1],x1);    /* 10 */
        m7d_qrt_mul(c_slice->row[0],x1,b_slice->row[1]);   /* 11 */
        m7d_qrt_mul( x1 , a_slice->row[1], b_slice->row[1]);  /* 12 */
        c_slice->row[1] = m7d_add(c_slice->row[1], x1 , c_slice->row[1]);   /* 13 */
        c_slice->row[2] = m7d_add(c_slice->row[2], c_slice->row[1] , c_slice->row[2]);   /* 14 */
        c_slice->row[1] = m7d_add(c_slice->row[1], c_slice->row[1] , c_slice->row[3]);   /* 15 */
        c_slice->row[3] = m7d_add(c_slice->row[3], c_slice->row[2] , c_slice->row[3]);    /* 16 */
        c_slice->row[3] = m7d_add(c_slice->row[3], c_slice->row[2] , c_slice->row[3]);  /* 17 */
        x2 = m7d_sub(x2, x2, b_slice->row[2]);            /* 18 */
        m7d_qrt_mul(c_slice->row[2], a_slice->row[1], x2);            /* 19 */
        c_slice->row[2] = m7d_sub( c_slice->row[2], c_slice->row[2], c_slice->row[0]);  /* 20 */

        x1 = m7d_sub( x1, a_slice->row[0], a_slice->row[2]);  /* 1 */
        x2 = m7d_sub(x2, b_slice->row[1],b_slice->row[1]);  /* 2 */
        m7d_qrt_mul(c_slice->row[2], x1, x2);  /* 3 */
        x1 = m7d_add(x1, a_slice->row[2],a_slice->row[1]);  /* 4 */
        x2 = m7d_sub(x2, b_slice->row[1],b_slice->row[0]);  /* 5 */
        m7d_qrt_mul(c_slice->row[3], x1, x2);    /* 6 */
        x1 = m7d_sub(x1, x1,a_slice->row[0]);/* 7 */
        x2 = m7d_sub(x2, b_slice->row[1],x2);  /* 8 */
        m7d_qrt_mul(c_slice->row[1],x1,x2); /* 9 */
        x1 = m7d_sub(x1, a_slice->row[1],x1);    /* 10 */
        m7d_qrt_mul(c_slice->row[0],x1,b_slice->row[1]);   /* 11 */
        m7d_qrt_mul( x1 , a_slice->row[1], b_slice->row[1]);  /* 12 */
        c_slice->row[1] = m7d_add(c_slice->row[1] , x1 , c_slice->row[1]);   /* 13 */
        c_slice->row[2] = m7d_add( c_slice->row[2], c_slice->row[1] , c_slice->row[2]);   /* 14 */
        c_slice->row[1] = m7d_add(c_slice->row[1] , c_slice->row[1] , c_slice->row[3]);   /* 15 */
        c_slice->row[3] = m7d_add(c_slice->row[3], c_slice->row[2] , c_slice->row[3]);    /* 16 */
        c_slice->row[3] = m7d_add(c_slice->row[3] ,c_slice->row[2] , c_slice->row[3]);  /* 17 */
        x2 = m7d_sub(x2, x2, b_slice->row[2]);            /* 18 */
        m7d_qrt_mul(c_slice->row[2], a_slice->row[1], x2);            /* 19 */
        c_slice->row[2] = m7d_sub(c_slice->row[2] ,  c_slice->row[2], c_slice->row[0]);  /* 20 */

        m7d_qrt_mul(c_slice->row[0], a_slice->row[1], b_slice->row[2]);
        c_slice->row[0] = m7d_add(c_slice->row[0], x1,c_slice->row[0] );
        m7d_free(x1);
    	m7d_free(x2);
       }	
    }
    
    else if((c->ncols ) == (M1RI_RADIX  << 1))
    {

	
		m7d_sub_64(x1->rows, a_slice->row[0]->rows, a_slice->row[2]->rows);  /* 1 */
        m7d_sub_64(x2->rows,b_slice->row[3]->rows,b_slice->row[1]->rows) ;  /* 2 */
        m7d_mul_64(c_slice->row[2]->rows, x1->rows, x2->rows);  /* 3 */
        m7d_add_64(x1->rows,a_slice->row[2]->rows,a_slice->row[3]->rows) ;  /* 4 */
        m7d_sub_64(x2->rows,b_slice->row[1]->rows,b_slice->row[0]->rows) ;  /* 5 */

     	m7d_mul_64(c_slice->row[3]->rows, x1->rows, x2->rows);    /* 6 */
        m7d_sub_64(x1->rows,x1->rows,a_slice->row[0]->rows) ;/* 7 */
        m7d_sub_64(x2->rows,b_slice->row[3]->rows,x2->rows);  /* 8 */
        m7d_mul_64(c_slice->row[1]->rows,x1->rows,x2->rows); /* 9 */
        m7d_sub_64(x1->rows,a_slice->row[1]->rows,x1->rows);    /* 10 */

        m7d_mul_64(c_slice->row[0]->rows,x1->rows,b_slice->row[3]->rows);   /* 11 */
        m7d_mul_64(x1->rows, a_slice->row[3]->rows, b_slice->row[3]->rows);  /* 12 */
        m7d_add_64(c_slice->row[1]->rows,x1->rows , c_slice->row[1]->rows) ;   /* 13 */
        m7d_add_64(c_slice->row[2]->rows,c_slice->row[1]->rows , c_slice->row[2]->rows) ;   /* 14 */
        m7d_add_64(c_slice->row[1]->rows,c_slice->row[1]->rows , c_slice->row[3]->rows) ;   /* 15 */
        m7d_add_64(c_slice->row[3]->rows,c_slice->row[2]->rows , c_slice->row[3]->rows) ;    /* 16 */
        m7d_add_64(c_slice->row[3]->rows,c_slice->row[2]->rows , c_slice->row[3]->rows) ;  /* 17 */

    	m7d_sub_64(x2->rows, x2->rows, b_slice->row[2]->rows) ;            /* 18 */
        m7d_mul_64(c_slice->row[2]->rows, a_slice->row[3]->rows, x2->rows);            /* 19 */
        m7d_sub_64(c_slice->row[2]->rows, c_slice->row[2]->rows,c_slice->row[0]->rows);  /* 20 */
        m7d_mul_64(c_slice->row[0]->rows, a_slice->row[1]->rows,b_slice->row[2]->rows);
        m7d_add_64(c_slice->row[0]->rows, x1->rows,c_slice->row[0]->rows) ; 
 


        m7d_free(x1);
    	m7d_free(x2);
    }
    
  
}

/**
	Strassen  algorithm on an m7d_t
*/

/**
  Strassen  algorithm on an m7d_t	
*/
m7d_t *  m7d_strassen(m7d_t *c,const m7d_t  *a,const m7d_t   *b)
{
	if (c == NULL)
	{
	c = m7d_create(a->nrows, b->ncols);

	} 
	
	else 
	{
	if (c->nrows != a->nrows || c->ncols != b->ncols) 
	{
		m1ri_die("m7d_strassen: Provided return matrix has wrong dimensions.\n");
    	}
	
	}
	
    if(a->ncols == b->nrows)
    {
 
      	/*  These hold the padded matrix sizes */
	u_int32_t  arcr, acbr, bccc;
	arcr = a->nrows;
	acbr = a->ncols;
	bccc = b->ncols;

    
	arcr =  powerof2(arcr);
	acbr =  powerof2(acbr);
	bccc =  powerof2(bccc);
	
	c = m7d_create( a->nrows, b->ncols); 
	int lasta, lastb, lastboth;
	lasta = 64 - a->nrows%64;
	lastb = 64 -  b->ncols%64;  
	lastboth = 64 - a->nrows; 

	if((arcr != a->nrows) || (acbr != a->ncols) || (bccc) != (b->ncols))
	{
	m7d_t * padded_a,   * padded_b , * padded_c;
	padded_a = m7d_create(arcr, acbr);
	padded_b = m7d_create(acbr, bccc);
	padded_a = m7d_create(arcr, bccc);
	m7d_copy(padded_a, a);
	m7d_copy(padded_b, b);
	m7d_copy(padded_c, c);
	
	
	m7d_qrt_mul(padded_c, padded_a, padded_b); 

	m7d_free(padded_a);
	m7d_free(padded_b);
	m7d_free(padded_c);
	}
	
	
	else
	{
	m7d_qrt_mul(c, a, b); 
	}
	

          
    }
    return c;
    
}




/**
	Classic O(N)^3 algorithm for Matrix Multiplication over GF(3) this function
	handles padding and calls  m3d_mul_naive_square for the multiplication portion.
*/
m3d_t * m3d_classic_mul(m3d_t *c, const m3d_t  *a, const m3d_t  *b)
{
	if (c == NULL)
	{
	c = m3d_create(a->nrows, b->ncols);

	} 
	
	else 
	{
	if (c->nrows != a->nrows || c->ncols != b->ncols) 
	{
		m1ri_die("m3d_mul_naive: Provided return matrix has wrong dimensions.\n");
    	}
	
	}
	
	
	c = m3d_create(  a->nrows, b->ncols); 
	/** * arcr, acbr, bccc hold the padded matrix sizes*/
	
	
	u_int32_t  arcr, acbr, bccc, g;
	arcr = a->nrows;
	acbr = a->ncols;
	bccc  = b->ncols;
	arcr =  powerof2(arcr);
	acbr =  powerof2(acbr);
	bccc =  powerof2(bccc);
	g = (M1RI_RADIX  << 1);
	while (arcr <  g)
	{
		arcr = arcr << 1;
	
	}
	while (acbr < g )
	{
	acbr = 	acbr << 1;
	
	}
	
	while (bccc <  g)
	{
		bccc = bccc << 1;
	
	}
	

	m3d_t * padded_a  = m1ri_malloc(sizeof(m3d_t));
	m3d_t  * padded_b  = m1ri_malloc(sizeof(m3d_t));
	m3d_t * padded_c = m1ri_malloc(sizeof(m3d_t));;
	
	if((arcr != a->nrows) || (acbr != a->ncols) || (bccc != b->ncols))
	{
		padded_a = m3d_create( arcr, acbr);
		padded_b = m3d_create( acbr, bccc);
		padded_c = m3d_create( arcr, bccc);
		m3d_copy(padded_a, a);
		m3d_copy(padded_b, b);
		m3d_mul_naive_square(padded_c, padded_a, padded_b); 
		m3d_copy_cutoff(c, padded_c);
		m3d_free(padded_a);
		m3d_free(padded_b);
		m3d_free(padded_c);
	}
		
	else
	{
		m3d_mul_naive_square(c, a, b); 
	}
	
	return c;
	
}

/*

	Classic O(N)^3 algorithm for Matrix Multiplication over GF(5)
*/
m5d_t * m5d_classic_mul(m5d_t *c, const m5d_t  *a, const m5d_t  *b)
{
	if (c == NULL)
	{
		c = m5d_create(a->nrows, b->ncols);

	} 
	
	else 
	{
	if (c->nrows != a->nrows || c->ncols != b->ncols) 
	{
		m1ri_die("m5d_mul_naive: Provided return matrix has wrong dimensions.\n");
    	}
	
	}
	
	
	c = m5d_create(  a->nrows, b->ncols); 
	/** * arcr, acbr, bccc hold the padded matrix sizes*/
	u_int32_t  arcr, acbr, bccc, g;
	arcr = a->nrows;
	acbr = a->ncols;
	bccc  = b->ncols;
	arcr =  powerof2(arcr);
	acbr =  powerof2(acbr);
	bccc =  powerof2(bccc);
	g = (M1RI_RADIX  << 1);
	while (arcr <  g)
	{
		arcr = arcr << 1;
	
	}
	while (acbr < g )
	{
	acbr = 	acbr << 1;
	
	}
	
	while (bccc <  g)
	{
		bccc = bccc << 1;
	
	}
	

	m5d_t * padded_a  = m1ri_malloc(sizeof(m5d_t));
	m5d_t  * padded_b  = m1ri_malloc(sizeof(m5d_t));
	m5d_t * padded_c = m1ri_malloc(sizeof(m5d_t));;

	if((arcr != a->nrows) || (acbr != a->ncols) || (bccc != b->ncols))
	{
		padded_a = m5d_create( arcr, acbr);
		padded_b = m5d_create( acbr, bccc);
		padded_c = m5d_create( arcr, bccc);
		m5d_copy(padded_a, a);
		m5d_copy(padded_b, b);
		m5d_mul_naive_square(padded_c, padded_a, padded_b); 
		m5d_copy_cutoff(c, padded_c);
		m5d_free(padded_a);
		m5d_free(padded_b);
		m5d_free(padded_c);
	}
		
	else
	{
		//m5d_mul_naive_square(c, a, b); 
	}
	
	return c;
	
	
}



m7d_t * m7d_classic_mul(m7d_t *c, const m7d_t  *a, const m7d_t  *b)
{
	if (c == NULL)
	{
	c = m7d_create(a->nrows, b->ncols);

	} 
	
	else 
	{		
		if (c->nrows != a->nrows || c->ncols != b->ncols) 
		{
			m1ri_die("m7d_mul_naive: Provided return matrix has wrong dimensions.\n");
	    	
		
		}
		
		
		c = m7d_create(  a->nrows, b->ncols); 
		/** * arcr, acbr, bccc hold the padded matrix sizes*/
		
		
		u_int32_t  arcr, acbr, bccc, g;
		arcr = a->nrows;
		acbr = a->ncols;
		bccc  = b->ncols;
		arcr =  powerof2(arcr);
		acbr =  powerof2(acbr);
		bccc =  powerof2(bccc);
		g = (M1RI_RADIX  << 1);
		while (arcr <  g)
		{
			arcr = arcr << 1;
		
		}
		while (acbr < g )
		{
		acbr = 	acbr << 1;
		
		}
		
		while (bccc <  g)
		{
			bccc = bccc << 1;
		
		}
		
	
		m7d_t * padded_a  = m1ri_malloc(sizeof(m7d_t));
		m7d_t  * padded_b  = m1ri_malloc(sizeof(m7d_t));
		m7d_t * padded_c = m1ri_malloc(sizeof(m7d_t));;
		
		if((arcr != a->nrows) || (acbr != a->ncols) || (bccc != b->ncols))
		{
			padded_a = m7d_create( arcr, acbr);
			padded_b = m7d_create( acbr, bccc);
			padded_c = m7d_create( arcr, bccc);
			m7d_copy(padded_a, a);
			m7d_copy(padded_b, b);
			m7d_mul_naive_square(padded_c, padded_a, padded_b); 
			m7d_copy_cutoff(c, padded_c);
			m7d_free(padded_a);
			m7d_free(padded_b);
			m7d_free(padded_c);
		}
			
		else
		{
			m7d_mul_naive_square(c, a, b); 
		}
	}
	
	return c;
	
}








