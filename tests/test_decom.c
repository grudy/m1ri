/** * M1RI
TOMAS J. 
 AND ROBERT W. BRADSHAW "BITSLICING AND THE METHOD OF FOUR
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

test_decom
 */ 

#include <m1ri/m1ri.h>
#include "time.h"




void m3d_transpose_test(int y, int z)
{
	m3d_t * a, *b, *c;
	a = m3d_create(y, z);
	b = m3d_create(z, y);
	c = m3d_create(y, z);
	 
	m3d_rand(a);
	b = m3d_transpose(b, a);
	c = m3d_transpose(c, b);
	
	/*
		testing if 
		(a * b)  * c  == a * (b * c)
		 d 	     * c  == a *   e  
				f     ==   g
 		
		where * is matrix multiplication 
	
	*/
	

	/*
	if((y <= 512) && (z <= 512))
	{
		
		printf("\n matrix a \n ");
		m3d_print(a);
		printf("\n matrix b \n");
		m3d_print(b);
		printf("\n matrix c \n");
		m3d_print(c);
		
	}
	
	*/
	if(!(m3d_equal(a, c)))
	{
	
	
		printf("\nm3d_transpose on two %d by %d matrix matrices not associative ", y, z  );
		m1ri_die("");
	
	}
	
	

	
    printf("----------------------------------------------------------------------");
    
 	m3d_print(a);
 	m3d_print(c);
	
    
    m3d_free(a);
    m3d_free(b);
    m3d_free(c);
    
 

}

void m5d_transpose_test(int y, int z)
{
	m5d_t * a, *b, *c;
	a = m5d_create(y, z);
	b = m5d_create(z, y);
	c = m5d_create(y, z);
	
	 
	m5d_rand(a);
	b = m5d_transpose(b, a);
	c = m5d_transpose(c, b);
	
	/*
		testing if 
		(a * b)  * c  == a * (b * c)
		 d 	     * c  == a *   e  
				f     ==   g
 		
		where * is matrix multiplication 
	
	*/
	

	/*
	if((y <= 256) && (z <= 256))
	{
		
		printf("\n matrix a \n ");
		m5d_print(a);
		printf("\n matrix b \n");
		m5d_print(b);
		printf("\n matrix c \n");
		m5d_print(c);
		
	}
	*/
	
	if(!(m5d_equal(a, c)))
	{
	
	
		printf("\nm5d_transpose on two %d by %d matrix matrices not associative ", y, z  );
		m1ri_die("");
	
	}
	
	

	
    printf("----------------------------------------------------------------------");
    

	
    
    m5d_free(a);
    m5d_free(b);
    m5d_free(c);
 
}


void m7d_transpose_test(int y, int z)
{
	m7d_t * a, *b, *c;
	a = m7d_create(y, z);
	b = m7d_create(z, y);
	c = m7d_create(y, z);
	 
	m7d_rand(a);
	b = m7d_transpose(b, a);
	c = m7d_transpose(c, b);
	
	/*
		testing if 
		(a * b)  * c  == a * (b * c)
		 d 	     * c  == a *   e  
				f     ==   g
 		
		where * is matrix multiplication 
	
	*/
	

	/*
	if((y <= 256) && (z <= 256))
	{
		
		printf("\n matrix a \n ");
		m7d_print(a);
		printf("\n matrix b \n");
		m7d_print(b);
		printf("\n matrix c \n");
		m7d_print(c);
		
	}
	
	
	if(!(m7d_equal(a, c)))
	{
	
	
		printf("\nm7d_transpose on two %d by %d matrix matrices not associative ", y, z  );
		m1ri_die("");
	
	}
	
	

	
    printf("----------------------------------------------------------------------");
    

	
    
    m7d_free(a);
    m7d_free(b);
    m7d_free(c);
    */
 
}


 int main(int argc, const char * argv[])
{
 	
	m3d_transpose_test(64, 64);
	
	m3d_transpose_test(128, 128);


	m3d_transpose_test(40, 120);

	m3d_transpose_test(512, 512);
	
	
	m5d_transpose_test(64, 64);
	
	m5d_transpose_test(128, 128);


	m5d_transpose_test(40, 120);
	
	m5d_transpose_test(512, 512);
	
	m5d_transpose_test(40, 120);
	
	m7d_transpose_test(64, 64);
	
	m7d_transpose_test(128, 128);



	m7d_transpose_test(59, 59);

	m7d_transpose_test(40, 120);

	m7d_transpose_test(512, 512);


	return 0;
    
    
}


  
  

