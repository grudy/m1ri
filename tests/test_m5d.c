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
 
 test_m5d.c
 */


#include <m1ri/m1ri.h>
#include <time.h>
int main(int argc, const char * argv[])
{
	
    int isequal;
    m5d_t * a, * b,  *d, *e, *f, *g, *h, *i, *j, *k;
    a = m5d_create(128, 128);
    b = m5d_create( 128, 128);
    m5d_rand(a);
    m5d_rand(b);
    isequal = m5d_equal(a, b);
    //m5d_print(a);
    //m5d_print(b);
    
    
     if(isequal)
    {
        m1ri_die("Equaltest: failed ");
        
    }
    
        printf("Equaltest: passed ");
        
    
    
   

    
    m5d_t *  test_m5d_output  = m5d_create( 3   ,3);
    
    m5d_rand(test_m5d_output);
    
    m5d_print(test_m5d_output);
    
    m5d_write_elem(test_m5d_output, 1, 1, 1, 1, 1);
    
    
    d = m5d_identity(64);
    e = m5d_identity(64);
	f = m5d_identity(64);
	

    isequal = m5d_equal(d, e);
    
    if(isequal)
    {
        printf("Equaltest: passed ");
        
    }
    
    if(!isequal)
    {
        printf("Equaltest: failed ");
        return 1;
        
    }
    
    
    m5d_free(f);
    
    
    
   f = m5d_create( 256, 256);
   g = m5d_create( 256, 256);
   h = m5d_create( 256, 256);
   i = m5d_create( 256, 256);
   j = m5d_create( 256, 256);
   k = m5d_create( 256, 256);
    
    
    m5d_rand(f);
    m5d_rand(g);
    m5d_rand(i);
    m5d_rand(j);
    
    
    k = m5d_classic_mul(k, i, j);
	m5d_print(k);
    
	m5d_t * y0, *y1,* y2,* y3, *y4,* y5,* y6;
	y0 = m5d_create(512, 512);
	y1 = m5d_create(512, 512);
	y2 = m5d_create(512, 512);
    m5d_rand(y0);
  	m5d_rand(y1);
  	m5d_rand(y2);
  	
	y3 = m5d_strassen(y3, y0, y1);
   	//m5d_print(y3);
   	/*
   	y4 = m5d_strassen(y4, y1, y2);
  	y5 = m5d_strassen(y5, y3, y2);
  	y6 = m5d_strassen(y6, y0, y4);
    m5d_print(y6);

	
	
	m5d_free(y0);
	m5d_free(y1);
	m5d_free(y2);
	m5d_free(y3);
	m5d_free(y4);
	m5d_free(y5);
	m5d_free(y6);   
	


	k = m5d_classic_mul(k, i, j);
	m5d_print(k); 

	m5d_classic_mul(k, i, j);
	/* m5d_print(k); * /

	m5d_free(h); 

/* 	m5d_free(h); * /

/* 	m5d_free(h); * /

	m5d_free(i);
	m5d_free(j);
	m5d_free(k);
	*/
	
    return 0;
}
