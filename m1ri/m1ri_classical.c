
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
 
 m1ri_classical.c
 */

#include <stdio.h>
#include <m1ri/m1ri_classical.h>
#include <m1ri/m1ri_cubes.h>
#include <m1ri/m7d.h>
void m3d_mul_naive_square(m3d_t *c, m3d_t *a, m3d_t *b)
{
	
	m3d_t  x1, x2; 
  	m3_slice *  a_slice, *  b_slice, *  c_slice;
	a_slice =  m1ri_malloc(sizeof(m3_slice));
	b_slice =   m1ri_malloc(sizeof(m3_slice));
	c_slice = m1ri_malloc(sizeof(m3_slice));
   	m3d_quarter(a_slice, a);
   //m3d_print(a_slice->row[0][0]);
    m3d_quarter(b_slice, b);
    m3d_quarter(c_slice, c);
   
    if((c_slice->row[0][0].ncols) > M1RI_RADIX)
    {
     
    	m3d_create(&x1, c_slice->row[0][0].nrows, c_slice->row[0][0].ncols);			    
     	m3d_create(&x2, c_slice->row[0][0].nrows, c_slice->row[0][0].ncols);			    	   
    	m3d_mul_naive_square(&x1, &a_slice->row[0][0], &b_slice->row[0][0]);
       	m3d_mul_naive_square(&x2, &a_slice->row[0][1], &b_slice->row[1][0]);
	    m3d_add_r(&c_slice->row[0][0], &x1, &x2) ; 
		m3d_mul_naive_square(&x1, &a_slice->row[0][0], &b_slice->row[0][1]);
       	m3d_mul_naive_square(&x2, &a_slice->row[0][1], &b_slice->row[1][1]);
    	m3d_add_r(&c_slice->row[0][1], &x1, &x2) ;
	    m3d_mul_naive_square(&x1, &a_slice->row[1][0], &b_slice->row[0][0]);
      	m3d_mul_naive_square(&x2, &a_slice->row[1][1], &b_slice->row[1][0]);
       	m3d_add_r(&c_slice->row[1][0], &x1, &x2); 
	    m3d_mul_naive_square(&x1, &a_slice->row[1][0], &b_slice->row[0][1]);
		m3d_mul_naive_square(&x2, &a_slice->row[1][1], &b_slice->row[1][1]); 
 	    m3d_add_r(&c_slice->row[1][1], &x1, &x2) ;
	
    }
   
    else if((c_slice->row[0][0].ncols ) <= M1RI_RADIX)
    {
    	m3d_create(&x1, M1RI_RADIX,M1RI_RADIX);			    
		m3d_create(&x2, M1RI_RADIX, M1RI_RADIX);
   	 	m3d_mul_64(x1.rows, a_slice->row[0][0].rows, b_slice->row[0][0].rows);
    	m3d_mul_64(x2.rows, a_slice->row[0][1].rows, b_slice->row[1][0].rows);
  	    m3d_add_64(c_slice->row[0][0].rows, x1.rows, x2.rows) ;
        m3d_mul_64(x1.rows,  a_slice->row[0][0].rows, b_slice->row[0][1].rows );
        m3d_mul_64(x2.rows, a_slice->row[0][1].rows, b_slice->row[1][1].rows);
		m3d_add_64(c_slice->row[0][1].rows, x1.rows, x2.rows) ; 
       	m3d_mul_64(x1.rows, a_slice->row[1][0].rows, b_slice->row[0][0].rows);
       	m3d_mul_64(x2.rows, a_slice->row[1][1].rows, b_slice->row[1][0].rows);
     	m3d_add_64(c_slice->row[1][0].rows, x1.rows, x2.rows); 
        m3d_mul_64(x1.rows, a_slice->row[1][0].rows, b_slice->row[0][1].rows);
        m3d_mul_64(x2.rows, a_slice->row[1][1].rows, b_slice->row[1][1].rows); 
		m3d_add_64(c_slice->row[1][1].rows, x1.rows, x2.rows); 
	
    }

}

void m3d_classic_mul(m3d_t *c, m3d_t  *a, m3d_t  *b)
{
	
	if(a->ncols == b->nrows);
	{
		
		m3d_create( c, a->nrows, b->ncols); 
		/* arcr, acbr, bccc hold the padded matrix sizes*/
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
			m3d_create(padded_a, arcr, acbr);
			m3d_create(padded_b, acbr, bccc);
			m3d_create(padded_c, arcr, bccc);
			m3d_copypadding(padded_a, a);
			m3d_copypadding(padded_b, b);
			m3d_mul_naive_square(padded_c, padded_a, padded_b); 
			m3d_putpadding(c, padded_c);
			m3d_free(padded_a);
			m3d_free(padded_b);
			m3d_free(padded_c);
		}
				
		else
		{
			m3d_mul_naive_square(c, a, b); 
		}
		
	}
	
}


//m5d and m7d multiplication are commented out because many  of the functions are 
//unwritten, untested, or unfinished

void m5d_mul_naive_square(m5d_t *c, m5d_t *a, m5d_t *b)
{

  	/*m5_slice *  a_slice, *  b_slice, *  c_slice;
	a_slice =  m1ri_malloc(sizeof(m5_slice));
	b_slice =   m1ri_malloc(sizeof(m5_slice));
	c_slice = m1ri_malloc(sizeof(m5_slice));
   	m5d_quarter(a_slice, a);
    m5d_quarter(b_slice, b);
    m5d_quarter(c_slice, c);
    m5d_t  x1; 
    m5d_t  x2; 

    if((c_slice->row[0][0].ncols) > M1RI_RADIX)
    {
     
		m5d_create(&x1, c_slice->row[0][0].nrows, c_slice->row[0][0].ncols);			    
     	m5d_create(&x2, c_slice->row[0][0].nrows, c_slice->row[0][0].ncols);			    
		m5d_mul_naive_square(&x1, &a_slice->row[0][0], &b_slice->row[0][0]);
       	m5d_mul_naive_square(&x2, &a_slice->row[0][1], &b_slice->row[1][0]);
        m5d_add_r(&c_slice->row[0][0], &x1, &x2) ; 
		m5d_mul_naive_square(&x1, &a_slice->row[0][0], &b_slice->row[0][1]);
       	m5d_mul_naive_square(&x2, &a_slice->row[0][1], &b_slice->row[1][1]);
       	m5d_add_r(&c_slice->row[0][1], &x1, &x2) ;
		m5d_mul_naive_square(&x1, &a_slice->row[1][0], &b_slice->row[0][0]);
       	m5d_mul_naive_square(&x2, &a_slice->row[1][1], &b_slice->row[1][0]);
       	m5d_add_r(&c_slice->row[1][0], &x1, &x2); 
		m5d_mul_naive_square(&x1, &a_slice->row[1][0], &b_slice->row[0][1]);
		m5d_mul_naive_square(&x2, &a_slice->row[1][1], &b_slice->row[1][1]); 
		m5d_add_r(&c_slice->row[1][1], &x1, &x2) ;
	
    }
   
    else if((c_slice->row[0][0].ncols ) <= M1RI_RADIX)
    {
      	m5d_create(&x1, M1RI_RADIX,M1RI_RADIX);			    
      	m5d_create(&x2, M1RI_RADIX, M1RI_RADIX);
   	  	m5d_mul_64(x1.rows, a_slice->row[0][0].rows, b_slice->row[0][0].rows);
		m5d_mul_64(x2.rows, a_slice->row[0][1].rows, b_slice->row[1][0].rows);
		m5d_add_r(&c_slice->row[0][0], &x1, &x2) ;
        m5d_mul_64(x1.rows,  a_slice->row[0][0].rows, b_slice->row[0][1].rows );
        m5d_mul_64(x2.rows, a_slice->row[0][1].rows, b_slice->row[1][1].rows);
		m5d_add_r(&c_slice->row[0][1], &x1, &x2) ;
       	m5d_mul_64(x1.rows, a_slice->row[1][0].rows, b_slice->row[0][0].rows);
      	m5d_mul_64(x2.rows, a_slice->row[1][1].rows, b_slice->row[1][0].rows);
		m5d_add_r(&c_slice->row[1][0], &x1, &x2) ;
		m5d_mul_64(x1.rows, a_slice->row[1][0].rows, b_slice->row[0][1].rows);
       	m5d_mul_64(x2.rows, a_slice->row[1][1].rows, b_slice->row[1][1].rows); 
		m5d_add_r(&c_slice->row[1][1],&x1, &x2) ;
		
    }
    */
}

void m5d_classic_mul(m5d_t *c, m5d_t  *a, m5d_t  *b)
{
	
	if(a->ncols == b->nrows);
	{
		m5d_create( c, a->nrows, b->ncols); 
		/* arcr, acbr, bccc hold the padded matrix sizes*/
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

			m5d_create(padded_a, arcr, acbr);
			m5d_create(padded_b, acbr, bccc);
			m5d_create(padded_c, arcr, bccc);
			m5d_copypadding(padded_a, a);
			m5d_copypadding(padded_b, b);
			m5d_mul_naive_square(padded_c, padded_a, padded_b); 
			m5d_putpadding(c, padded_c);
			m5d_free(padded_a);
			m5d_free(padded_b);
			m5d_free(padded_c);
		}
		
		else
		{
			m5d_mul_naive_square(c, a, b); 
		}	
		 	
	}
	   
}

void m7d_mul_naive_square(m7d_t *c, m7d_t *a, m7d_t *b)
{
	m7d_t  x1, x2;
 	m7_slice *  a_slice, *  b_slice, *  c_slice;
	a_slice =  m1ri_malloc(sizeof(m7_slice));
	b_slice =   m1ri_malloc(sizeof(m7_slice));
	c_slice = m1ri_malloc(sizeof(m7_slice));
   	m7d_quarter(a_slice, a);
    m7d_quarter(b_slice, b);
    m7d_quarter(c_slice, c);

    if((c_slice->row[0][0].ncols) > M1RI_RADIX)
    {
     
    	m7d_create(&x1, c_slice->row[0][0].nrows, c_slice->row[0][0].ncols);			    
     	m7d_create(&x2, c_slice->row[0][0].nrows, c_slice->row[0][0].ncols);			    
     	m7d_mul_naive_square(&x1, &a_slice->row[0][0], &b_slice->row[0][0]);
     	m7d_mul_naive_square(&x2, &a_slice->row[0][1], &b_slice->row[1][0]);
    	m7d_add_r(&c_slice->row[0][0], &x1, &x2) ; 
		m7d_mul_naive_square(&x1, &a_slice->row[0][0], &b_slice->row[0][1]);
      	m7d_mul_naive_square(&x2, &a_slice->row[0][1], &b_slice->row[1][1]);
		m7d_add_r(&c_slice->row[0][1], &x1, &x2) ;
		m7d_mul_naive_square(&x1, &a_slice->row[1][0], &b_slice->row[0][0]);
       	m7d_mul_naive_square(&x2, &a_slice->row[1][1], &b_slice->row[1][0]);
       	m7d_add_r(&c_slice->row[1][0], &x1, &x2); 
        m7d_mul_naive_square(&x1, &a_slice->row[1][0], &b_slice->row[0][1]);
        m7d_mul_naive_square(&x2, &a_slice->row[1][1], &b_slice->row[1][1]); 
	    m7d_add_r(&c_slice->row[1][1], &x1, &x2) ;
	
    }
   
    else if((c_slice->row[0][0].ncols ) <= M1RI_RADIX)
    {
        m7d_create(&x1, M1RI_RADIX,M1RI_RADIX);			    
      	m7d_create(&x2, M1RI_RADIX, M1RI_RADIX);
     	m7d_mul_64(x1.rows, a_slice->row[0][0].rows, b_slice->row[0][0].rows);
		m7d_mul_64(x2.rows, a_slice->row[0][1].rows, b_slice->row[1][0].rows);
		m7d_add_64(c_slice->row[0][0].rows, x1.rows, x2.rows) ;
		m7d_mul_64(x1.rows,  a_slice->row[0][0].rows, b_slice->row[0][1].rows );
        m7d_mul_64(x2.rows, a_slice->row[0][1].rows, b_slice->row[1][1].rows);
	    m7d_add_64(c_slice->row[0][1].rows, x1.rows, x2.rows) ; 
		m7d_mul_64(x1.rows, a_slice->row[1][0].rows, b_slice->row[0][0].rows);
       	m7d_mul_64(x2.rows, a_slice->row[1][1].rows, b_slice->row[1][0].rows);
      	m7d_add_64(c_slice->row[1][0].rows, x1.rows, x2.rows); 
		m7d_mul_64(x1.rows, a_slice->row[1][0].rows, b_slice->row[0][1].rows);
		m7d_mul_64(x2.rows, a_slice->row[1][1].rows, b_slice->row[1][1].rows); 
		m7d_add_64(c_slice->row[1][1].rows, x1.rows, x2.rows); 
		
    }
     
}


void m7d_classic_mul(m7d_t *c, m7d_t  *a, m7d_t  *b)
{

	if(a->ncols == b->nrows);
	{
		m7d_create( c, a->nrows, b->ncols); 
		/* arcr, acbr, bccc hold the padded matrix sizes*/
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
		
		//This creates a padded matrix for recursive multiplication
		if((arcr != a->nrows) || (acbr != a->ncols) || (bccc != b->ncols))
		{
		
			m7d_create(padded_a, arcr, acbr);
			m7d_create(padded_b, acbr, bccc);
			m7d_create(padded_c, arcr, bccc);
			m7d_copypadding(padded_a, a);
			m7d_copypadding(padded_b, b);
			m7d_mul_naive_square(padded_c, padded_a, padded_b); 
			m7d_putpadding(c, padded_c);
			m7d_free(padded_a);
			m7d_free(padded_b);
			m7d_free(padded_c);
			
		}
		else
		{
			m7d_mul_naive_square(c, a, b); 
		}

	}
   
}

