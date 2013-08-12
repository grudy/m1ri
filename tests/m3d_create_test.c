/* M1RI
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
 
 m3d_arith_test.c
 */ 





#include <m1ri/m1ri.h>
#include "time.h"
int main(int argc, const char * argv[])
{
 	 m3d_t *a, *b, *c;
    a = malloc(sizeof(m3d_t));
    b = malloc(sizeof(m3d_t));
    c =  malloc(sizeof(m3d_t));
 	 m3d_create(a, 128, 128);
  m3d_create(b, 128, 128);
   m3d_create(c, 128, 128);
 // m3d_rand(a);
    //m3d_rand(b);
    	m3_slice * a_slice , * b_slice, * c_slice;
    	a_slice = m1ri_malloc(sizeof(m3_slice));
    	b_slice = m1ri_malloc(sizeof(m3_slice));
    	c_slice = m1ri_malloc(sizeof(m3_slice));
	//m3d_add_r(c, a ,b);
	// m3d_create(a, 64, 64);
  // m3d_create(b, 64, 64);
    /// m3d_create(c, 64, 64);
	//a->rows[0][0].units = leftbit;
 	//m3d_print(a);
	//m3d_print(b);
	//m3d_slices(a_slice, a, 1);
	//m3d_identity_set(&a_slice->row[0][0]);
	m3d_identity_set(a);
 
  m3d_print(a);
  		//m3d
  
  
  ///	m3d_print(c);
	//m3d_classic_mul(c, a, b);
  //	m3d_mul_naive_square(c, a, b);
	
	//m3d_print(c);
    
    return 0;
    
    
}

