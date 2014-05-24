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
	printf("m3d addition and subtraction test passed");
	
  
  
  return 0;
  

}


