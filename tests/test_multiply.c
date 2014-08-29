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
 
 */ 

#include <m1ri/m1ri.h>




#include "time.h"




void m3d_mul_associative_test(int y, int z)
{
	m3d_t * a, *b, *c, *d, * e, *f, *g;
	a = m3d_create(y, z);
	b = m3d_create(y, z);
	c  = m3d_create(y, z);
	
	m3d_rand(a);
	m3d_rand(b);
	m3d_rand(c);

	/*
		testing if 
		(a * b)  * c  == a * (b * c)
		 d 	     * c  == a *   e  
				f     ==   g
 		
		where * is matrix multiplication 
	
	*/
	

	d = m3d_strassen(d, a, b);
	f = m3d_strassen(f, d, c);
	
	
	e = m3d_strassen(e, b, c);
	g = m3d_strassen(g, a, e);
	
	
	
	
	
	if((y <= 256) && (z <= 256))
	{
		printf("\n matrix f \n");
		m3d_print(f);
		printf("\n matrix g \n");
		m3d_print(g);
	
	}
	
	
	if(!(m3d_equal(f, g)))
	{
	
		printf("\nm3d_strassen on two %d by %d matrix matrices not associative ", y, z  );
		m1ri_die("");
	
	}
	
	
	
	
    printf("----------------------------------------------------------------------");
    
    
    
    m3d_free(a);
    m3d_free(b);
    m3d_free(c);
    m3d_free(d);
    m3d_free(e);
    m3d_free(f);
    m3d_free(g);
    
    

}


void m3d_strassen_test(int y, int z)
{
	m3d_t * a, *b, *c;
	a = m3d_create(y, z);
    b = m3d_create(y, z);
    
    m3d_rand(a);
    m3d_rand(b);
    
   // printf("\nOutput of first matrix a of size %d, by %d\n", y, z);
   //	m3d_print(a);
    
    m3d_print(a);
    m3d_print(b);
    
	clock_t begin, end;
	double time_spent;
	begin = clock();
	c = m3d_strassen(c, a, b);
	//printf("\n %d by %d matrix \n", c->nrows, c->ncols);
	m3d_print(c);
    time_t after;
    time(&after);
	end = clock();
	time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    
    printf("----------------------------------------------------------------------");
    printf("\nm3d_strassen on two %d by %d matrix matrices.", y, z  );
    printf(" \n------------------->Runs in%9f seconds. \n", time_spent);
    printf("----------------------------------------------------------------------");
    m3d_free(a);
    m3d_free(b);
    m3d_free(c);
    
    

}


int main(int argc, const char * argv[])
{

	
 	//m3d_strassen_test(64, 64);
   	//m3d_strassen_test(512, 512);


   m3d_mul_associative_test(
   4, 4);
 
   
    
   return 0 ; 
}    