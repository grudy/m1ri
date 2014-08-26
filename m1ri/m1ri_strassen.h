
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
 
 m1ri_strassen.h
 */


#ifndef M1RIGF3_STRASSEN_H
#define M1RIGF3_STRASSEN_H
#include <m1ri/m3d.h>
#include <m1ri/m5d.h>
#include <m1ri/m7d.h>
#include <m1ri/m1ri_io.h>


/**
	\brief Strassen Matrix Multiplication over GF(3), on a matrix.
	\param c product matrix, may be null
	\param b multiplicand matrix, must not be null
	\param a multiplier matrix, must not be null 	
	\	
	\	c = a * b
	
*/

m3d_t *  m3d_strassen(m3d_t *, const m3d_t *,const m3d_t*);


/**
	\brief Strassen Matrix Multiplication over GF(3), on a matrix.
	\param c product matrix, may be null
	\param b multiplicand matrix, must not be null
	\param a multiplier matrix, must not be null 	
	\	
	\	c = a * b
*/

m5d_t *  m5d_strassen(m5d_t * ,const m5d_t *, const m5d_t *);


/**
	\brief Strassen Matrix Multiplication over GF(3), on a matrix.
	\param c product matrix, may be null
	\param b multiplicand matrix, must not be null
	\param a multiplier matrix, must not be null 	
	\	
	\	c = a * b
	
*/
m7d_t *  m7d_strassen(m7d_t * ,const m7d_t *, const m7d_t *);



/**
	\brief Classical Matrix Multiplication over GF(3), on a matrix.
	\param c product matrix, may be null
	\param b multiplicand matrix, must not be null
	\param a multiplier matrix, must not be null 	
	\	
	\	c = a * b
	
*/

m3d_t *  m3d_classic_mul(m3d_t *,const   m3d_t  * , const m3d_t  *);

/**
	\brief Classical Matrix Multiplication over GF(5), on a matrix.
	\param c product matrix, may be null
	\param b multiplicand matrix, must not be null
	\param a multiplier matrix, must not be null 	
	\	
	\	c = a * b
	
*/

m5d_t *  m5d_classic_mul(m5d_t *, const m5d_t  * , const m5d_t  *);


/**
	\brief Classical Matrix Multiplication over GF(7), on a matrix.
	\param c product matrix, may be null
	\param b multiplicand matrix, must not be null
	\param a multiplier matrix, must not be null 	
	\	
	\	c = a * b
	
*/

m7d_t *  m7d_classic_mul(m7d_t *, const  m7d_t  * , const m7d_t  *);

/**
	Recursive Matrix Multiplication over GF(5), on a square matrix.
*/
static inline void m5d_mul_naive_square(m5d_t *c, const m5d_t *a, const m5d_t *b)
{
	
	m5d_t  * x1, * x2; 
  	m5_slice *  a_slice, *  b_slice, *  c_slice;
	a_slice = m5d_quarter( a);
    b_slice = m5d_quarter( b);
    c_slice = m5d_quarter(c);
    
    if((c_slice->row[0]->ncols) > M1RI_RADIX)
    {
       
    	x1 = m5d_create( c_slice->row[0]->nrows, c_slice->row[0]->ncols);		    
     	x2 = m5d_create( c_slice->row[0]->nrows, c_slice->row[0]->ncols);	    	   
		
		m5d_mul_naive_square(x1, a_slice->row[0], b_slice->row[0]);
    	
    	m5d_mul_naive_square(x2, a_slice->row[1], b_slice->row[2]);
	    c_slice->row[0] = m5d_add(c_slice->row[0],  x1, x2) ; 
		m5d_mul_naive_square(x1, a_slice->row[0], b_slice->row[1]);
       	m5d_mul_naive_square(x2, a_slice->row[1], b_slice->row[3]);
    	
    	c_slice->row[1] = m5d_add(c_slice->row[1],  x1, x2) ;
	    m5d_mul_naive_square(x1, a_slice->row[2], b_slice->row[0]);
      	m5d_mul_naive_square(x2, a_slice->row[3], b_slice->row[2]);
       	c_slice->row[2] = m5d_add(c_slice->row[2], x1, x2); 
	    m5d_mul_naive_square(x1, a_slice->row[2], b_slice->row[1]);
		m5d_mul_naive_square(x2, a_slice->row[3], b_slice->row[3]); 
 	    c_slice->row[3] = m5d_add(c_slice->row[3], x1, x2) ;
 	 	m5d_free(x1);
		m5d_free(x2);
		
    }
   
    else if((c->ncols ) == (M1RI_RADIX  << 1))
    {
       
    	x1 = m5d_create(M1RI_RADIX,M1RI_RADIX);		    
		x2 = m5d_create( M1RI_RADIX, M1RI_RADIX);
   	 	
   	 	m5d_mul_64(x1->rows, a_slice->row[0]->rows, b_slice->row[0]->rows);
    	m5d_mul_64(x2->rows, a_slice->row[1]->rows, b_slice->row[2]->rows);
  	    m5d_add_64(c_slice->row[0]->rows, x1->rows, x2->rows) ;
        m5d_mul_64(x1->rows,  a_slice->row[0]->rows, b_slice->row[1]->rows );
        m5d_mul_64(x2->rows, a_slice->row[1]->rows, b_slice->row[3]->rows);
		m5d_add_64(c_slice->row[1]->rows, x1->rows, x2->rows) ; 
       	m5d_mul_64(x1->rows, a_slice->row[2]->rows, b_slice->row[0]->rows);
       	m5d_mul_64(x2->rows, a_slice->row[3]->rows, b_slice->row[2]->rows);
     	m5d_add_64(c_slice->row[2]->rows, x1->rows, x2->rows); 
        m5d_mul_64(x1->rows, a_slice->row[2]->rows, b_slice->row[1]->rows);
        m5d_mul_64(x2->rows, a_slice->row[3]->rows, b_slice->row[3]->rows); 
		m5d_add_64(c_slice->row[3]->rows, x1->rows, x2->rows); 
		
		m5d_free(x1);
		m5d_free(x2);
		
		
    }
 	
 	else if(c->ncols  == M1RI_RADIX )
    {
    	  m5d_mul_64(c->rows, a->rows, b->rows);
    
    } 
 	  m5d_quarter_free(a_slice);
      m5d_quarter_free(b_slice);
      m5d_quarter_free(c_slice);  
}


/**
	\brief Recursive Matrix Multiplication over GF(7), on a square matrix.
*/

static inline void m7d_mul_naive_square(m7d_t *c, const m7d_t *a, const m7d_t *b)
{
	
	m7d_t  * x1, * x2; 
  	m7_slice *  a_slice, *  b_slice, *  c_slice;

   	a_slice = m7d_quarter( a);
   /* m7d_print(a_slice->row[0][0]); */
    b_slice = m7d_quarter( b);
    c_slice = m7d_quarter(c);
   
    if((c_slice->row[0]->ncols) > M1RI_RADIX)
    {
       
    	x1 = m7d_create( c_slice->row[0]->nrows, c_slice->row[0]->ncols);		    
     	x2 = m7d_create( c_slice->row[0]->nrows, c_slice->row[0]->ncols);	
     	    	   
    	m7d_mul_naive_square(x1, a_slice->row[0], b_slice->row[0]);
    	m7d_mul_naive_square(x2, a_slice->row[1], b_slice->row[2]);
	    c_slice->row[0] = m7d_add(c_slice->row[0],  x1, x2) ; 
		m7d_mul_naive_square(x1, a_slice->row[0], b_slice->row[1]);
       	m7d_mul_naive_square(x2, a_slice->row[1], b_slice->row[3]);
    	c_slice->row[1] = m7d_add(c_slice->row[1],  x1, x2) ;
	    m7d_mul_naive_square(x1, a_slice->row[2], b_slice->row[0]);
      	m7d_mul_naive_square(x2, a_slice->row[3], b_slice->row[2]);
       	c_slice->row[2] = m7d_add(c_slice->row[2], x1, x2); 
	    m7d_mul_naive_square(x1, a_slice->row[2], b_slice->row[1]);
		m7d_mul_naive_square(x2, a_slice->row[3], b_slice->row[3]); 
 	    c_slice->row[3] = m7d_add(c_slice->row[3], x1, x2) ;
 	 	m7d_free(x1);
		m7d_free(x2);
	
    }
   
    else if((c->ncols ) == (M1RI_RADIX  << 1))
    {
    	x1 = m7d_create(M1RI_RADIX,M1RI_RADIX);		    
		x2 = m7d_create( M1RI_RADIX, M1RI_RADIX);
   	 	m7d_mul_64(x1->rows, a_slice->row[0]->rows, b_slice->row[0]->rows);
    	m7d_mul_64(x2->rows, a_slice->row[1]->rows, b_slice->row[2]->rows);
  	    m7d_add_64(c_slice->row[0]->rows, x1->rows, x2->rows) ;
        m7d_mul_64(x1->rows,  a_slice->row[0]->rows, b_slice->row[1]->rows );
        m7d_mul_64(x2->rows, a_slice->row[1]->rows, b_slice->row[3]->rows);
		m7d_add_64(c_slice->row[1]->rows, x1->rows, x2->rows) ; 
       	m7d_mul_64(x1->rows, a_slice->row[2]->rows, b_slice->row[0]->rows);
       	m7d_mul_64(x2->rows, a_slice->row[3]->rows, b_slice->row[2]->rows);
     	m7d_add_64(c_slice->row[2]->rows, x1->rows, x2->rows); 
        m7d_mul_64(x1->rows, a_slice->row[2]->rows, b_slice->row[1]->rows);
        m7d_mul_64(x2->rows, a_slice->row[3]->rows, b_slice->row[3]->rows); 
		m7d_add_64(c_slice->row[3]->rows, x1->rows, x2->rows); 
		m7d_free(x1);
		m7d_free(x2);
		
    }

}



#endif

