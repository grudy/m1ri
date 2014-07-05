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

 test_m3d.c
 */ 

#include <m1ri/m1ri.h>
#include "time.h"
 int main(int argc, const char * argv[])
{
  /*
   [a][b]
   [c][d]
   
  */
  m3d_t * o, *ab, *cd, *abcd, *v; 
  m3_slice * s;  
 
  o   = m3d_create( 128, 128);
  m3d_rand(o);
		//  
  
   s = m3d_quarter( o);
  
  
  m3d_specs(s->row[0]);
  m3d_print(s->row[0]);
  m3d_print(s->row[1]);
  m3d_print(s->row[2]);
  m3d_print(s->row[3]);

  printf("\n\n s done \n");
  ab = m3d_concat(s->row[0], s->row[1]);
  cd = m3d_concat( s->row[2], s->row[3]);
  
  m3d_specs(s->row[0]);
  abcd = m3d_stack(abcd,  ab, cd);
  m3d_specs(o);
  m3d_specs(abcd);
  
  printf("\nab\n");
  m3d_print(ab);
  printf("\ncd\n");
  m3d_print(cd);
  printf("\n abcd \n");
  m3d_print(abcd);
  printf("\no\n");
  m3d_print(o);
  if(!m3d_equal(o, abcd))
  {
     printf("\n o and abcd not equal \n");
  
  }
  
  v = m3d_create( 4, 4);
  
  m3d_rand(v);
  m3d_print(v);

  m3d_colswap(v, 1, 2);
  m3d_print(v);


  /*
    Multiplication (classic) Test
    (wo*w1)w2 = wo*(w1*w2);
    w3 = wo*w1;
    w4 = w1*w2;
    w5 = w3*w2;
    w6 = w0*w4;
  */ 
  
  m3d_t * w0, *w1,* w2,* w3, *w4,* w5,* w6;

  

  w0 = m3d_create(512, 512);
  w1 = m3d_create(512, 512);
  w2 = m3d_create(512, 512);
  
  m3d_rand(w0);
  m3d_rand(w1);
  m3d_rand(w2);
  
  w3 = m3d_classic_mul(w3, w0, w1);
  w4 = m3d_classic_mul(w4, w1, w2);
  w5 = m3d_classic_mul(w5, w3, w2);
  w6 = m3d_classic_mul(w6, w0, w4);
   	

  /* 	
  if(!(m3d_equal(w5, w6)))
	 {
	 	printf("Classic Multiplication m3d test failed");
	//return 1;
	 
	 }   
  
  
  printf("Classic Multiplication m3d test passed"); 
  */
  
	m3d_t * y0, *y1,* y2,* y3, *y4,* y5,* y6;
	y0 = m3d_create(512, 512);
	y1 = m3d_create(512, 512);
	y2 = m3d_create(512, 512);
  
	m3d_rand(y0);
	m3d_rand(y1);
	m3d_rand(y2);
  
	y3 = m3d_strassen(y3, y0, y1);
  
	y4 = m3d_strassen(y4, y1, y2);
	y5 = m3d_strassen(y5, y3, y2);
	y6 = m3d_strassen(y6, y0, y4);
 
  
	m3d_free(w0);
	m3d_free(w1);
	m3d_free(w2);
	m3d_free(w3);
	m3d_free(w4);
	m3d_free(w5);
	m3d_free(w6);
	
	m3d_free(y0);
	m3d_free(y1);
	m3d_free(y2);
	m3d_free(y3);
	m3d_free(y4);
	m3d_free(y5);
	m3d_free(y6);
	
	m3d_free(v);	
	m3d_free(o);
	m3d_free(ab);
	m3d_free(cd);
	m3d_free(abcd);
	m1ri_free(s);	
  
  
    return 0;
}

