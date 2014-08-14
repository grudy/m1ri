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
 
 test_hadamard.c
 */


#include <m1ri/m1ri.h>

#include "time.h"


void m3d_hadamard_test(int m, int n)
{
	    
     /* A*(B+C) == A*B + A*C */
     
     /* m3d */
    m3d_t * a_m3d,  * b_m3d,  * c_m3d, * t1_m3d,* t2_m3d, *t3_m3d, * r1_m3d, *  r2_m3d;
 
	a_m3d = m3d_create(m, n);
	b_m3d = m3d_create(m, n);
	c_m3d = m3d_create(m, n);

	/*
	m3d_rand(a_m3d);
	m3d_rand(b_m3d);
	m3d_rand(c_m3d);
	 /*
	     a_m3d * (b_m3d  + c_m3d)   == (a_m3d * b_m3d) + (a_m3d *  c_m3d);
	 */
	 
	 t1_m3d = m3d_add(t1_m3d, b_m3d, c_m3d);
	 r1_m3d = m3d_hadamard(r1_m3d, a_m3d, t1_m3d);
	 printf("\n\t\tM3D \n******************************\n");
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
	 t2_m3d = m3d_hadamard(t2_m3d, a_m3d , b_m3d);
	 t3_m3d = m3d_hadamard(t3_m3d, a_m3d ,  c_m3d);
	 r2_m3d = m3d_add(r2_m3d, t2_m3d, t3_m3d);
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
    


}

void m5d_hadamard_test(int m, int n)
{
	    
     /* A*(B+C) == A*B + A*C */
     
     /* m5d */
    m5d_t * a_m5d,  * b_m5d,  * c_m5d, * t1_m5d,* t2_m5d, *t3_m5d, * r1_m5d, *  r2_m5d;
 
	a_m5d = m5d_create( m, n);
	b_m5d = m5d_create(m, n);
	c_m5d = m5d_create(m, n);

	/*
	m5d_rand(a_m5d);
	m5d_rand(b_m5d);
	m5d_rand(c_m5d);
	 /*
	     a_m5d * (b_m5d  + c_m5d)   == (a_m5d * b_m5d) + (a_m5d *  c_m5d);
	 */
	 
	 t1_m5d = m5d_add(t1_m5d, b_m5d, c_m5d);
	 r1_m5d = m5d_hadamard(r1_m5d, a_m5d, t1_m5d);
	 printf("\n\t\tm5D \n******************************\n");
	 printf("\n a_m5d\n");
	 m5d_print(a_m5d);
	 printf("\n b_m5d\n");
	 m5d_print(b_m5d);
	 printf("\n c_m5d\n");
	 m5d_print(c_m5d);
	 printf("\n t1_m5d\n");
	 m5d_print(t1_m5d);
	 printf("\n r1_m5d\n");
	 m5d_print(r1_m5d);
	 t2_m5d = m5d_hadamard(t2_m5d, a_m5d , b_m5d);
	 t3_m5d = m5d_hadamard(t3_m5d, a_m5d ,  c_m5d);
	 r2_m5d = m5d_add(r2_m5d, t2_m5d, t3_m5d);
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
    


}



void m7d_hadamard_test(int m, int n)
{
	    
     /* A*(B+C) == A*B + A*C */
     
     /* m7d */
    m7d_t * a_m7d,  * b_m7d,  * c_m7d, * t1_m7d,* t2_m7d, *t3_m7d, * r1_m7d, *  r2_m7d;
 
	a_m7d = m7d_create(m, n);
	b_m7d = m7d_create(m, n);
	c_m7d = m7d_create(m, n);

	/*
	m7d_rand(a_m7d);
	m7d_rand(b_m7d);
	m7d_rand(c_m7d);
	 /*
	     a_m7d * (b_m7d  + c_m7d)   == (a_m7d * b_m7d) + (a_m7d *  c_m7d);
	 */
	 
	 t1_m7d = m7d_add(t1_m7d, b_m7d, c_m7d);
	 r1_m7d = m7d_hadamard(r1_m7d, a_m7d, t1_m7d);
	 printf("\n\t\tm7D \n******************************\n");
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
	 t2_m7d = m7d_hadamard(t2_m7d, a_m7d , b_m7d);
	 t3_m7d = m7d_hadamard(t3_m7d, a_m7d ,  c_m7d);
	 r2_m7d = m7d_add(r2_m7d, t2_m7d, t3_m7d);
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
    


}
int main(int argc, const char * argv[])
{
    
    
	
       
    
    
    
    
	return 0;
}


