/** * M1RI
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
 
 m1ri_hadamard.c
 */


#include <m1ri/m1ri.h>

#include "time.h"
int main(int argc, const char * argv[])
{
    
    
    
     //A*(B+C) == A*B + A*C
     
     //m3d
    m3d_t * a_m3d,  * b_m3d,  * c_m3d, * t1_m3d,* t2_m3d, *t3_m3d, * r1_m3d, *  r2_m3d;
    a_m3d = malloc(sizeof(m3d_t));
    b_m3d = malloc(sizeof(m3d_t));
    c_m3d = malloc(sizeof(m3d_t));
    t1_m3d = malloc(sizeof(m3d_t));
	t2_m3d = malloc(sizeof(m3d_t));
	t3_m3d = malloc(sizeof(m3d_t));
	r1_m3d = malloc(sizeof(m3d_t));
	r2_m3d = malloc(sizeof(m3d_t));
	m3d_create(a_m3d, 4, 4);
	m3d_create(b_m3d, 4, 4);
	m3d_create(c_m3d, 4, 4);

	
	m3d_rand(a_m3d);
	m3d_rand(b_m3d);
	m3d_rand(c_m3d);
	 /*
	     a_m3d * (b_m3d  + c_m3d)   == (a_m3d * b_m3d) + (a_m3d *  c_m3d);
	 */
	 
	 m3d_add_r(t1_m3d, b_m3d, c_m3d);
	 r1_m3d = m3d_hadamard(a_m3d, t1_m3d);
	 printf("\n\t\tM3D \n********************************************\n");
	 printf("\n a_m3d\n");
	 m3d_print(a_m3d);
	 printf("\n b_m3d\n");
	 m3d_print(b_m3d);
	 printf("\n c_m3d\n");
	 m3d_print(c_m3d);
	 printf("\n t1_m3d\n");
	 m3d_print(t1_m3d);
	 printf("\n r1_m3d\n");
	 m3d_print(r1_m3d);
	 t2_m3d = m3d_hadamard(a_m3d , b_m3d);
	 t3_m3d = m3d_hadamard(a_m3d ,  c_m3d);
	 m3d_add_r(r2_m3d, t2_m3d, t3_m3d);
	 if(!(m3d_equal(r1_m3d, r2_m3d)))
	 {
	 	printf("Hadamard m3d test failed");
	    return 1;
	 
	 }   
	 printf("Hadamard m3d test passed"); 
	m3d_free(a_m3d);
    m3d_free(b_m3d);
    m3d_free(c_m3d);
    m3d_free(t1_m3d);
    m3d_free(t2_m3d);
    m3d_free(t3_m3d);
	m3d_free(r1_m3d);
    m3d_free(r2_m3d);
    
       
    
     //A*(B+C) == A*B + A*C
     
     //m5d
     
     
    m5d_t * a_m5d,  * b_m5d,  * c_m5d, * t1_m5d,* t2_m5d, *t3_m5d, * r1_m5d, *  r2_m5d;
    a_m5d = malloc(sizeof(m5d_t));
    b_m5d = malloc(sizeof(m5d_t));
    c_m5d = malloc(sizeof(m5d_t));
    t1_m5d = malloc(sizeof(m5d_t));
	t2_m5d = malloc(sizeof(m5d_t));
	t3_m5d = malloc(sizeof(m5d_t));
	r1_m5d = malloc(sizeof(m5d_t));
	r2_m5d = malloc(sizeof(m5d_t));
	m5d_create(a_m5d, 4, 4);
	m5d_create(b_m5d, 4, 4);
	m5d_create(c_m5d, 4, 4);

	
	m5d_rand(a_m5d);
	m5d_rand(b_m5d);
	m5d_rand(c_m5d);
//	 /*
//	     a_m5d * (b_m5d  + c_m5d)   == (a_m5d * b_m5d) + (a_m5d *  c_m5d);
//	 * /
	 
	 
	 
	 m5d_add_r(t1_m5d, b_m5d, c_m5d);
	 r1_m5d = m5d_hadamard(a_m5d, t1_m5d);
	 
	 printf("\n\t\tm5d \n********************************************\n");
	 printf("\n a_m5d\n");
	 m5d_print(a_m5d);
	 printf("\n b_m5d\n");
	 m5d_print(b_m5d);
	 
	
	 t2_m5d = m5d_hadamard(a_m5d , b_m5d);
	  printf("\n t2_m5d\n");
	 m5d_print(t2_m5d);
	 
	 
	  printf("\n c_m5d\n");
	 m5d_print(c_m5d);
	 printf("\n t1_m5d\n");
	 m5d_print(t1_m5d);
	 printf("\n r1_m5d\n");
	 m5d_print(r1_m5d);
	 
	 t3_m5d = m5d_hadamard(a_m5d ,  c_m5d);
	 m5d_add_r(r2_m5d, t2_m5d, t3_m5d);
	  
	 printf("\n r2_m5d\n");
	  m5d_print(r2_m5d);
	  
	 if(!(m5d_equal(r1_m5d, r2_m5d)))
	 {
	 	printf("Hadamard m5d test failed");
	    return 1;
	 
	 }   
	 
	 
	 
	 
	 printf("Hadamard m5d test passed"); 
	m5d_free(a_m5d);
    m5d_free(b_m5d);
    m5d_free(c_m5d);
    
    m5d_free(t1_m5d);
    m5d_free(t2_m5d);
    m5d_free(t3_m5d);
	m5d_free(r1_m5d);
    m5d_free(r2_m5d);
      
      //A*(B+C) == A*B + A*C
     
     //m7d
    m7d_t * a_m7d,  * b_m7d,  * c_m7d, * t1_m7d,* t2_m7d, *t3_m7d, * r1_m7d, *  r2_m7d;
    a_m7d = malloc(sizeof(m7d_t));
    b_m7d = malloc(sizeof(m7d_t));
    c_m7d = malloc(sizeof(m7d_t));
    t1_m7d = malloc(sizeof(m7d_t));
	t2_m7d = malloc(sizeof(m7d_t));
	t3_m7d = malloc(sizeof(m7d_t));
	r1_m7d = malloc(sizeof(m7d_t));
	r2_m7d = malloc(sizeof(m7d_t));
	m7d_create(a_m7d, 4, 4);
	m7d_create(b_m7d, 4, 4);
	m7d_create(c_m7d, 4, 4);

	
	m7d_rand(a_m7d);
	m7d_rand(b_m7d);
	m7d_rand(c_m7d);
	 /*
	     a_m7d * (b_m7d  + c_m7d)   == (a_m7d * b_m7d) + (a_m7d *  c_m7d);
	 */
	 
	 m7d_add_r(t1_m7d, b_m7d, c_m7d);
	 r1_m7d = m7d_hadamard(a_m7d, t1_m7d);
	
	 m7d_print(r2_m7d);
	 t2_m7d = m7d_hadamard(a_m7d , b_m7d);
	 t3_m7d = m7d_hadamard(a_m7d ,  c_m7d);
	 m7d_add_r(r2_m7d, t2_m7d, t3_m7d);
	  printf("\n\t\tm7d \n********************************************\n");
	 printf("\n a_m7d\n");
	 m7d_print(a_m7d);
	 printf("\n b_m7d\n");
	 m7d_print(b_m7d);
	 printf("\n c_m7d\n");
	 m7d_print(c_m7d);
	 printf("\n t1_m7d\n");
	 m7d_print(t1_m7d);
	 printf("\n r1_m7d\n");
	 m7d_print(r1_m7d);
	 printf("\n r2_m7d\n");
	 m7d_print(r2_m7d);
	 if(!(m7d_equal(r1_m7d, r2_m7d)))
	 {
	 	printf("Hadamard m7d test failed");
	    return 1;
	 
	 }   
	 printf("Hadamard m7d test passed"); 
	m7d_free(a_m7d);
    m7d_free(b_m7d);
    m7d_free(c_m7d);
    m7d_free(t1_m7d);
    m7d_free(t2_m7d);
    m7d_free(t3_m7d);
	m7d_free(r1_m7d);
    m7d_free(r2_m7d);
    
    
	return 0;
}


