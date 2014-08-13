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
 
 test_arithmatic.c
 */


#include <m1ri/m1ri.h>

#include "time.h"


 void m3d_test_addition(int m,int  n)
 {

    m3d_t * a, * b, *c, *d;
    a = m3d_create(m, n);
    b = m3d_create( m, n);
    m3d_rand(a);
    m3d_rand(b);
    
    c = m3d_add(c, a, b);
    
    d = m3d_sub(d, c, b);
   
	m3d_print(d);

	if(!(m3d_equal(d, a)))
    {
     	{
          m1ri_die("Error in m3d addition and subtraction of size %d by %d,\n", m, n);
    
         }
    
    }
   
  m3d_free(a);
  m3d_free(b);
  m3d_free(c);
  m3d_free(d);
  
 
 }


 void m5d_test_addition(int m,int  n)
 {

	  m5d_t * e, * f, *g, * h;
  
  

 
    e = m5d_create( m, n);
    f = m5d_create( m, n);

    m5d_rand(e);
    m5d_rand(f);
    
     g = m5d_add(g, e, f);
     h = m5d_sub(h,   g, f);



	m5d_print(f);
	
	if(!(m5d_equal(h, e)))
    {
     	{
          m1ri_die("Error in m5d addition and subtraction of size %d by %d,\n", m, n);
    
         }
    
    }
    
    
  m5d_free(e);
  m5d_free(f);
  m5d_free(g);
  m5d_free(h);

 
 }


void m7d_test_addition(int m,int  n)
 {

	  m7d_t * e, * f, *g, * h;
  
  

 	
    e = m7d_create( m, n);
    f = m7d_create( m, n);

    m7d_rand(e);
    m7d_rand(f);
    
     g = m7d_add(g, e, f);
     h = m7d_sub(h,   g, f);


	m7d_print(f);
	
	if(!(m7d_equal(h, e)))
    {
     	{
          m1ri_die("Error in m7d addition and subtraction of size %d by %d,\n", m, n);
    
         }
    
    }
    
    
  m7d_free(e);
  m7d_free(f);
  m7d_free(g);
  m7d_free(h);

 
 }


int main(int argc, const char * argv[])
{
    
    m3d_test_addition(4, 4);
    m3d_test_addition(64, 64);
    m3d_test_addition(64, 4);
    m3d_test_addition(4, 800);
    m3d_test_addition(14, 294);
    m3d_test_addition(342, 64);
    m3d_test_addition(64, 44);
    m3d_test_addition(142, 181);
    printf("\nm3d addition and subtraction test passed\n");

        
    m5d_test_addition(4, 4);
    m5d_test_addition(64, 64);
    m5d_test_addition(64, 4);
    m5d_test_addition(4, 800);
    m5d_test_addition(14, 294);
    m5d_test_addition(342, 64);
    m5d_test_addition(64, 44);
    m5d_test_addition(142, 181);
    printf("\nm5d addition and subtraction test passed\n");

    
            
    m7d_test_addition(4, 4);
    m7d_test_addition(64, 64);
    m7d_test_addition(64, 4);
    m7d_test_addition(4, 800);
	m7d_test_addition(14, 294);
    m7d_test_addition(342, 64);
    m7d_test_addition(64, 44);
    m7d_test_addition(142, 181);
    printf("\nm7d addition and subtraction test passed\n");

  printf("\n \n ****************************************\n \t\t\t\t\tm5d\n");
  
  
  
    printf("\n\n****************************************\n \t\t\t\t\tm7d\n");
     m7d_t * i, * j, *k,  * m;
  
    i = m7d_create( 4, 4);
   
    j = m7d_create( 4, 4);

    m7d_rand(i);
    m7d_rand(j);
   /*  m7d_specs(i); */

     k = m7d_add(k,  i, j);
     m = m7d_sub(m,  k,  j);


    printf("Matrix i\n");
    m7d_print(i);
    printf("Matrix j\n");
    m7d_print(j);
    printf("Matrix k\n");
    m7d_print(k);
	printf("Matrix m\n");
    m7d_print(m); 
	
	
	if(!(m7d_equal(i, m)))
    {
     	{
     	    m7d_free(i);
  			m7d_free(j);
  			m7d_free(k);
  			m7d_free(m);
  	
     	
     	
          printf("Error in m7d addition and subtraction\n");
          return 1; 
    
         }
    
    }
    /*  */
	printf("m7d addition and subtraction test passed");
  
  m7d_free(i);
  m7d_free(j);
  m7d_free(k);
  m7d_free(m);
  
  
  return 0;
  

}


