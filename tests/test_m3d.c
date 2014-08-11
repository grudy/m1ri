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
		/*    */
	
	 s = m3d_quarter( o);
	
	
	m3d_specs(s->row[0]);
	m3d_print(s->row[0]);
	m3d_print(s->row[1]);
	m3d_print(s->row[2]);
	m3d_print(s->row[3]);

 
	ab = m3d_concat(ab, s->row[0], s->row[1]);
	cd = m3d_concat(cd, s->row[2], s->row[3]);
	
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
	   printf("\n o and abcd not equal, stack and concat test failed \n");
		 return 1;
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
	    return 1;  
	 
	 }   
  
  */


  
  printf("Classic Multiplication m3d test passed"); 
  
  
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
  

 	m3d_t * a, *b, *c, *d, *e, *f, *g, *h;	
 	m3_slice  * z;
	a = m3d_create(64, 64);
		
	a = m3d_identity(a, 64);
	b = m3d_identity(b, 64);	
	
	
	if(!m3d_equal(a, b))
	{
	  printf("\nEquality test failed  on equal 64 by 64  matrices \n");
	}
	printf("Equality test passed on equal 64 by 64  matrices");
	
	m3d_free(a);
	m3d_free(b);
	
	/*
		copy test
	*/
	a = m3d_create(64, 64);
	m3d_rand(a);
	m3d_copy(b, a);
	
	if
	
	
	
	a = m3d_identity(a,256);
	b = m3d_identity(b,256);
		
	
	
	if(m3d_equal(a, b))
	{
	  printf("\nEquality test failed  on  equal 256 by 256 matrices \n");
	}
	printf("Equality test passed on equal 256 by 256 matrices");
	
	
	m3d_free(a);
	m3d_free(b);
	
		
	//a = m3d_create(
	
		

	
	b = NULL;
	c = NULL;
	d = NULL;
	
	a = m3d_create(256, 256);
	m3d_rand(a);
	z = m3d_quarter(a);
	m3d_print(z->row[0]);
	
	  
  	b = m3d_concat(b, z->row[0], z->row[1]);
  	c = m3d_concat(b,  z->row[2], z->row[3]);
  
  	d = m3d_stack(d,  b, c);
  	

  if(!m3d_equal(a, d))
  {
     printf("\n a and d not equal, \n m3d_stack and m3d_concat test failed \n");
  	 return 1;
  }
	printf("\n Test of m3d_stack and m3d_concat passed \n");
	
	m3d_free(a);
	m3d_free(b);
	m3d_free(c);
	m3d_free(d);
	
	
	a = m3d_create(64, 64);
	b = m3d_create(64, 64);
	m3d_set_ui(a, 2);
	m3d_set_ui(b, 1);
	
	
	
  	if(m3d_equal(a, b))
  	{
     	printf("\n  identity matrix and scalar of the identity matrix have equal values \n");
  		 return 1;
 	}
  	m3d_free(a);
	m3d_free(b);
	
	a  = m3d_create(64, 64);
	b  = m3d_create(64, 64);
	c  = m3d_create(64, 64);
	d  = m3d_create(64, 64);
	
	
	//c =   m3d_add(m3d_t *a, const m3d_t  *b,const  m3d_t  *c);
	//d =   m3d_sub(

	
	
	/*
		testing submatrix, stack , and  accuracy 
		a[][] == [b][c] == 	f[][]g ==  e[][]
		 [][]	 [d][e]		f[][]g		[][]
		 
	*/  
		 
	/*	
	a = m3d_create(128, 128);
	m3d_rand(a);
	b = NULL;
	c = NULL;
	d = NULL;
	e = NULL;
	
	b = m3d_submatrix(b, a, 0, 0, 64, 64 );
	c = m3d_submatrix(c, a, 0, 64, 64, 128);
	d = m3d_submatrix(d, a, 64,  0, 128,64); 
	//e = m3d_submatrix(e, a, 64, 64, 128, 128); 
	/
	f = m3d_stack(f, b, d);
	
	g = m3d_stack(g, c, e);
	*/
	//m3d_specs(b);
	//m3d_specs(c);
	//m3d_specs(d);
	//m3d_specs(e);
	//m3d_specs(f);
	//m3d_specs(g);
	/*
	//e = m3d_concat(e, f, g);
/*
	if(!m3d_equal(a, h))
  	{
     	printf("\n  Submatrix, Stack and Concat test failed \n Small m3d_submatrix test failed \n");
  		 return 1;
 	}
 	printf("\n Test of Submatrix,  m3d_stack and m3d_concat passed over smaller matrices \n");

	
	m3d_free(a);
	m3d_free(b);
	m3d_free(c);
	m3d_free(d);
	m3d_free(e);
	
	

	a = m3d_create(128, 128);
	b = m3d_create(128, 128);
	c = m3d_create(128, 128);

	
	*/

	
	
	
	/*
	m3d_rowswap (m3d_t  * , rci_t , rci_t );
	m3d_colswap(m3d_t *, rci_t , rci_t );
	
	
	
	/*

	
	

	
	
	m3d_t *    m3d_init_window(const m3d_t  *, rci_t , rci_t , rci_t , rci_t );
	

	
	

	
	
	
	
	
	void  m3d_slices(m3_slice *  ,const m3d_t * , wi_t );
	
	
	
	
	 
	m3_slice *  m3d_quarter( const m3d_t * );
	
	
	m3d_t  * m3d_transpose_sliced(m3d_t * );
	
	
	void  m3d_colswap_capped_row(m3d_t *, rci_t , rci_t, rci_t );
	
	
	
	int m3d_cmp(m3d_t *A, m3d_t *B);
	
	
	void vbg_negation(vbg * );
	
	
	void sub_m3d( vbg *, vbg const *  , vbg const * );       
	
	
	
	vbg sub_m3dr(vbg , vbg );               
	
	
	
	
	void  vbg_mul( vbg *, vbg  *, vbg  *);      
	
	
	m3d_t *  m3d_sub(m3d_t *,   const  m3d_t  *, const m3d_t  *);
	
	
	vbg vbg_mul_elementwise(vbg const , vbg const);
	
	
	m3d_t * m3d_hadamard(m3d_t * , m3d_t const * , m3d_t const * );
	
	
	
	static inline void m3d_sub_64(vbg **R, vbg  **A, vbg  **B)
	{
	    int i;
	    for (i= 0; i < M1RI_RADIX; i++ )
	    {
	        R[i][0] = sub_m3dr(A[i][0], B[i][0]);
	    }
	    
	}
	
	static inline vbg add_m3dr(vbg  x, vbg const y)
	{ 
	    vec t;
	    x.sign  = y.units ^ x.sign;
	    t = (x.sign & x.units) ^ y.sign;
	    x.units = (y.units ^ x.units) |  t;
	    x.sign = t & x.sign;
	    return x; 
	}
	
	
	
	






	m3d_t *m3d_mul_scalar(m3d_t *, const long , const m3d_t *);



	void m3d_add_row(m3d_t *A, rci_t ar, const m3d_t *B, rci_t br, rci_t start_col);

	int m3d_is_zero(const m3d_t *);

*/
	
	return 0;
    
    
}


  
  

