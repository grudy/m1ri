
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
#include "m1ri_classical.h"




void m3d_mul_naive_square(m3d_t *c, m3d_t *a, m3d_t *b)
{
  m3_slice  a_slice,  b_slice,  c_slice;
   m3d_create(c, a->nrows, b->ncols);
   m3d_quarter(&a_slice, a);
    m3d_quarter(&b_slice, b);
    m3d_quarter(&c_slice, c);
   
			    
    if((c_slice.row[0][0].ncols) > cutoff)
    {
       
      m3d_t * x1 = m1ri_malloc(sizeof(m3d_t)); 
     m3d_t * x2 = m1ri_malloc(sizeof(m3d_t)); 
     m3d_create(x1, c_slice.row[0][0].nrows, c_slice.row[0][0].ncols);			    
     m3d_create(x2, c_slice.row[0][0].nrows, c_slice.row[0][0].ncols);			    
			    
       m3d_mul_naive_square(x1, &a_slice.row[0][0], &b_slice.row[0][0]);
       m3d_mul_naive_square(x2, &a_slice.row[0][1], &a_slice.row[1][0]);
      
       m3d_add_r(&c_slice.row[0][0], x1, x2) ; 

       m3d_mul_naive_square(x1, &a_slice.row[0][0], &b_slice.row[0][1]);
       m3d_mul_naive_square(x2, &a_slice.row[0][1], &b_slice.row[1][1]);

       m3d_add_r(&c_slice.row[0][1], x1, x2) ;

       m3d_mul_naive_square(x1, &a_slice.row[1][0], &b_slice.row[0][0]);
       m3d_mul_naive_square(x2, &a_slice.row[1][1], &b_slice.row[1][0]);
       m3d_add_r(&c_slice.row[1][0], x1, x2); 

       m3d_mul_naive_square(x1, &a_slice.row[1][0], &b_slice.row[0][1]);
       m3d_mul_naive_square(x2, &a_slice.row[1][1], &b_slice.row[1][1]); 
 
      m3d_add_r(&c_slice.row[1][1], x1, x2) ;
	
        m3d_free(x1);
        m3d_free(x2);  
    
    }
    
    else if((c_slice.row[0][0].ncols ) <= cutoff)
    {
        
       mul_64_m3d(c_slice.row[0][0].rows, b_slice.row[0][0].rows, a_slice.row[0][0].rows);
        mul_64_m3d(c_slice.row[0][1].rows, b_slice.row[0][1].rows, a_slice.row[0][1].rows);
        mul_64_m3d(c_slice.row[1][0].rows, b_slice.row[1][0].rows, a_slice.row[1][0].rows);
         mul_64_m3d(c_slice.row[1][1].rows, b_slice.row[1][1].rows, a_slice.row[1][1].rows);
        
     
    }
     
}


void m3d_classic_mul(m3d_t *c, m3d_t *a, m3d_t *b)
{
	if(a->ncols == b->nrows);
	{
		m3d_create( c, a->nrows, b->ncols); 
		int lrows, lcols;
		lrows = a->nrows%64;
		lcols = b->ncols%64;  
		m3d_t * c_main_partition = malloc(sizeof(m3d_t));
		m3d_t * a_main_partition = malloc(sizeof(m3d_t));
		m3d_t * b_main_partition = malloc(sizeof(m3d_t));
		m3d_window_create(c, c_main_partition, 0, 0, c->nrows -1, c->ncols -1); 
		m3d_window_create(a, a_main_partition, 0, 0, a->nrows -1, a->ncols -1   ); 
		m3d_window_create(b, b_main_partition,  0, 0,b->nrows -1 , b->ncols -1); 
		m3d_mul_naive_square(c_main_partition, b_main_partition, a_main_partition); 
		
		
		
		
	
	
	}
	
	


}




