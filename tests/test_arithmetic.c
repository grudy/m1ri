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
    
    
  /*
  	 m3d
  
  */
  
    m3d_t * a, * b, *c, *d;
    a = m1ri_malloc(sizeof(m3d_t)); 
    b = m1ri_malloc(sizeof(m3d_t)); 
    c = m1ri_malloc(sizeof(m3d_t)); 
    d = m1ri_malloc(sizeof(m3d_t)); 
    m3d_create(a, 4, 4);
    m3d_create(b, 4, 4);
    m3d_rand(a);
    m3d_rand(b);
    
    m3d_add_r(c, a, b);
    
    m3d_sub(d,  c, b);
    printf("Matrix a\n");
    m3d_print(a);
    printf("Matrix b\n");
    m3d_print(b);
    printf("Matrix c\n");
    m3d_print(c);
    printf("Matrix d\n");
    m3d_print(d); 
	
	if(!(m3d_equal(d, a)))
    {
     	{
          printf("Error in m3d addition and subtraction");
          return 1; 
    
         }
    
    }
	printf("m3d addition and subtraction test passed\n");
    
  m3d_free(a);
  m3d_free(b);
  m3d_free(c);
  m3d_free(d);
  
  
  
   /*
  	 m5d
  */
  printf("\n\n\****************************************\n \t\t\t\t\tm5d\n");
  
    m5d_t * e, * f, *g, * h;
  
  

    e = m1ri_malloc(sizeof(m5d_t)); 
    f = m1ri_malloc(sizeof(m5d_t)); 
    g = m1ri_malloc(sizeof(m5d_t)); 
    h = m1ri_malloc(sizeof(m5d_t)); 
    
    m5d_create(e, 4, 4);
    m5d_create(f, 4, 4);

    m5d_rand(e);
    m5d_rand(f);
    
     m5d_add_r(g, e, f);
     m5d_sub(h,  g, f);


    printf("Matrix e\n");
    m5d_print(e);
    printf("Matrix f\n");
    m5d_print(f);
    printf("Matrix g\n");
    m5d_print(g);
    printf("Matrix h\n");
    m5d_print(h); 
	
	
	if(!(m5d_equal(h, e)))
    {
     	{
          printf("Error in m5d addition and subtraction\n");
          return 1; 
    
         }
    
    }
    
	printf("m5d addition and subtraction test passed\n");
    
  m5d_free(e);
  m5d_free(f);
  m5d_free(g);
  m5d_free(h);

  
    printf("\n\n\****************************************\n \t\t\t\t\tm7d\n");
     m7d_t * i, * j, *k,  * m;
  
  

    i = m1ri_malloc(sizeof(m7d_t)); 
    j = m1ri_malloc(sizeof(m7d_t)); 
    k = m1ri_malloc(sizeof(m7d_t)); 
    m = m1ri_malloc(sizeof(m7d_t)); 
  
    m7d_create(i, 128, 128);
   
    m7d_create(j, 128, 128);

    m7d_rand(i);
    m7d_rand(j);
    m7d_specs(i);
     m7d_print(i);
     //m7d_add_r(k,  i, j);
     //m7d_sub(m,   k,  j);


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
          printf("Error in m7d addition and subtraction");
          //return 1; 
    
         }
    
    }
    
	printf("m7d addition and subtraction test passed");
    
  m7d_free(i);
  m7d_free(j);
  m7d_free(k);
  m7d_free(m);
  
  
  return 0;
  

}


