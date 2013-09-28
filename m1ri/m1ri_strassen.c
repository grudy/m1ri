
/** 
 Matrix Represenations and basic operations
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
 
 m1ri_strassen.c
 */

#include "m1ri_strassen.h"
#include <math.h>
#include <stdlib.h>


/**
	Multiplies m3d_t matrices in squares.  
	Not to be used directly, but called by m3d_strassen
*/
void m3d_qrt_mul(m3d_t * c,m3d_t *a, m3d_t * b )
{
    m3d_t * x1;
    m3d_t * x2;
    x1 = m1ri_malloc(sizeof(m3d_t));
    x2 = m1ri_malloc(sizeof(m3d_t));
    m3_slice  a_slice, b_slice, c_slice;
    m3d_create(x1, c->nrows/2, c->ncols/2);
    m3d_create(x2, c->nrows/2, c->ncols/2);
    m3d_quarter(&a_slice, a);
    m3d_quarter(&b_slice, b);
    m3d_quarter(&c_slice, c);
    
    if((c_slice.row[0][0].ncols) > M1RI_RADIX)
    {
       
        m3d_sub(x1, &a_slice.row[0][0], &a_slice.row[1][0]);  //1
        m3d_sub(x2,&b_slice.row[1][1],&b_slice.row[0][1]);  //2
        m3d_qrt_mul(&c_slice.row[1][0], x1, x2);  //3
        m3d_add_r(x1,&a_slice.row[1][0],&a_slice.row[1][1]);  //4
        m3d_sub(x2,&b_slice.row[0][1],&b_slice.row[0][0]);  //5
        m3d_qrt_mul(&c_slice.row[1][1], x1, x2);    //6
        m3d_sub(x1,x1,&a_slice.row[0][0]);//7
        m3d_sub(x2,&b_slice.row[1][1],x2);  //8
        m3d_qrt_mul(&c_slice.row[0][1],x1,x2); //9
        m3d_sub(x1,&a_slice.row[0][1],x1);    //10
        m3d_qrt_mul(&c_slice.row[0][0],x1,&b_slice.row[1][1]);   //11
        m3d_qrt_mul(x1, &a_slice.row[1][1], &b_slice.row[1][1]);  //12
        m3d_add_r(&c_slice.row[0][1],x1 , &c_slice.row[0][1]);   //13
        m3d_add_r(&c_slice.row[1][0],&c_slice.row[0][1] , &c_slice.row[1][0]);   //14
        m3d_add_r(&c_slice.row[0][1],&c_slice.row[0][1] , &c_slice.row[1][1]);   //15
        m3d_add_r(&c_slice.row[1][1],&c_slice.row[1][0] , &c_slice.row[1][1]);    //16
        m3d_add_r(&c_slice.row[1][1],&c_slice.row[1][0] , &c_slice.row[1][1]);  //17
        m3d_sub(x2, x2, &b_slice.row[1][0]);            //18
        m3d_qrt_mul(&c_slice.row[1][0], &a_slice.row[1][1], x2);            //19
        m3d_sub(&c_slice.row[1][0], &c_slice.row[1][0], &c_slice.row[0][0]);  //20
        m3d_qrt_mul(&c_slice.row[0][0], &a_slice.row[0][1], &b_slice.row[1][0]);
        m3d_add_r(&c_slice.row[0][0], x1,&c_slice.row[0][0] );
        
    }
    
    else if((c_slice.row[0][0].ncols ) == M1RI_RADIX)
    {
		m3d_sub_64(x1->rows, a_slice.row[0][0].rows, a_slice.row[1][0].rows);  //1
        m3d_sub_64(x2->rows,b_slice.row[1][1].rows,b_slice.row[0][1].rows) ;  //2
        m3d_mul_64(c_slice.row[1][0].rows, x1->rows, x2->rows);  //3
        m3d_add_64(x1->rows,a_slice.row[1][0].rows,a_slice.row[1][1].rows) ;  //4
        m3d_sub_64(x2->rows,b_slice.row[0][1].rows,b_slice.row[0][0].rows) ;  //5
        m3d_mul_64(c_slice.row[1][1].rows, x1->rows, x2->rows);    //6
        m3d_sub_64(x1->rows,x1->rows,a_slice.row[0][0].rows) ;//7
        m3d_sub_64(x2->rows,b_slice.row[1][1].rows,x2->rows);  //8
        m3d_mul_64(c_slice.row[0][1].rows,x1->rows,x2->rows); //9
        m3d_sub_64(x1->rows,a_slice.row[0][1].rows,x1->rows);    //10
        m3d_mul_64(c_slice.row[0][0].rows,x1->rows,b_slice.row[1][1].rows);   //11
        m3d_mul_64(x1->rows, a_slice.row[1][1].rows, b_slice.row[1][1].rows);  //12
        m3d_add_64(c_slice.row[0][1].rows,x1->rows , c_slice.row[0][1].rows) ;   //13
        m3d_add_64(c_slice.row[1][0].rows,c_slice.row[0][1].rows , c_slice.row[1][0].rows) ;   //14
        m3d_add_64(c_slice.row[0][1].rows,c_slice.row[0][1].rows , c_slice.row[1][1].rows) ;   //15
        m3d_add_64(c_slice.row[1][1].rows,c_slice.row[1][0].rows , c_slice.row[1][1].rows) ;    //16
        m3d_add_64(c_slice.row[1][1].rows,c_slice.row[1][0].rows , c_slice.row[1][1].rows) ;  //17
        m3d_sub_64(x2->rows, x2->rows, b_slice.row[1][0].rows) ;            //18
        m3d_mul_64(c_slice.row[1][0].rows, a_slice.row[1][1].rows, x2->rows);            //19
        m3d_sub_64(c_slice.row[1][0].rows, c_slice.row[1][0].rows,c_slice.row[0][0].rows);  //20
        m3d_mul_64(c_slice.row[0][0].rows, a_slice.row[0][1].rows,b_slice.row[1][0].rows);
        m3d_add_64(c_slice.row[0][0].rows, x1->rows,c_slice.row[0][0].rows) ; 
    }
    
    m1ri_free(x1);
    m1ri_free(x2);
}

/**
  Strassen  algorithm on an m3d_t	
*/
void  m3d_strassen(m3d_t *c, m3d_t  *a, m3d_t   *b)
{
    if(a->ncols == b->nrows)
    {
 
      	// These hold the padded matrix sizes
		u_int32_t  arcr, acbr, bccc;
		arcr = a->nrows;
		acbr = a->ncols;
		bccc = b->ncols;
		
		arcr =  powerof2(arcr);
		acbr =  powerof2(acbr);
		bccc =  powerof2(bccc);
		
		m3d_create( c, a->nrows, b->ncols); 
		int lasta, lastb, lastboth;
		lasta = 64 - a->nrows%64;
		lastb = 64 -  b->ncols%64;  
		lastboth = 64 - a->nrows; 

		if((arcr != a->nrows) || (acbr != a->ncols) || (bccc) != (bccc))
		{
		m3d_t * padded_a  = m1ri_malloc(sizeof(m3d_t));
		m3d_t  * padded_b  = m1ri_malloc(sizeof(m3d_t));
		m3d_t * padded_c = m1ri_malloc(sizeof(m3d_t));;
		m3d_create(padded_a, arcr, acbr);
		m3d_create(padded_b, acbr, bccc);
		m3d_create(padded_c, arcr, bccc);
		m3d_copypadding(padded_a, a);
		m3d_copypadding(padded_b, b);
		m3d_copypadding(padded_c, c);
		
		
		m3d_qrt_mul(padded_c, padded_a, padded_b); 

		m3d_free(padded_a);
		m3d_free(padded_b);
		m3d_free(padded_c);
		}
		
	
		else
		{
		m3d_qrt_mul(c, a, b); 
		}
		

          
    }
    
}

/**
	This handles the arithmetic of m5d_strassen
*/

void m5d_qrt_mul(m5d_t * c,m5d_t *a, m5d_t * b )
{
    m5d_t * x1, *x2;
    
    x1 = m1ri_malloc(sizeof(m5d_t));
    x2 = m1ri_malloc(sizeof(m5d_t));
    m5_slice  a_slice, b_slice, c_slice;
    m5d_create(x1, c->nrows/2, c->ncols/2);
    m5d_create(x2, c->nrows/2, c->ncols/2);
    m5d_quarter(&a_slice, a);
    m5d_quarter(&b_slice, b);
    m5d_quarter(&c_slice, c);
    if((c_slice.row[0][0].ncols) > M1RI_RADIX)
    {
        m5d_sub_r(x1, &a_slice.row[0][0], &a_slice.row[1][0]);  //1
        m5d_sub_r(x2,&b_slice.row[1][1],&b_slice.row[0][1]);  //2
        m5d_qrt_mul(&c_slice.row[1][0], x1, x2);  //3
        m5d_add_r(x1,&a_slice.row[1][0],&a_slice.row[1][1]);  //4
        m5d_sub_r(x2,&b_slice.row[0][1],&b_slice.row[0][0]);  //5
        m5d_qrt_mul(&c_slice.row[1][1], x1, x2);    //6
        m5d_sub_r(x1,x1,&a_slice.row[0][0]);//7
        m5d_sub_r(x2,&b_slice.row[1][1],x2);  //8
    	m5d_qrt_mul(&c_slice.row[0][1],x1,x2); //9
		m5d_sub_d(x1,&a_slice.row[0][1]);    //10  //todo is this the correct order????
        m5d_qrt_mul(&c_slice.row[0][0],x1,&b_slice.row[1][1]);   //11
        m5d_qrt_mul(x1, &a_slice.row[1][1], &b_slice.row[1][1]);  //12
        m5d_add_r(&c_slice.row[0][1],x1 , &c_slice.row[0][1]);   //13
        m5d_add_r(&c_slice.row[1][0],&c_slice.row[0][1] , &c_slice.row[1][0]);   //14
        m5d_add_r(&c_slice.row[0][1],&c_slice.row[0][1] , &c_slice.row[1][1]);   //15
        m5d_add_r(&c_slice.row[1][1],&c_slice.row[1][0] , &c_slice.row[1][1]);    //16
        m5d_add_r(&c_slice.row[1][1],&c_slice.row[1][0] , &c_slice.row[1][1]);  //15
        m5d_sub_r(x2, x2, &b_slice.row[1][0]);            //18
        m5d_qrt_mul(&c_slice.row[1][0], &a_slice.row[1][1], x2);            //19
        m5d_sub_r(&c_slice.row[1][0], &c_slice.row[1][0], &c_slice.row[0][0]);  //20
        m5d_qrt_mul(&c_slice.row[0][0], &a_slice.row[0][1], &b_slice.row[1][0]);
        m5d_add_r(&c_slice.row[0][0], x1,&c_slice.row[0][0] );
        
    }
  
    else if((c_slice.row[0][0].ncols ) == M1RI_RADIX)
    {	
    
		m5d_sub_64(x1, &a_slice.row[0][0], &a_slice.row[1][0]);  //1
   		m5d_sub_64(x2,&b_slice.row[1][1],&b_slice.row[0][1]) ;  //2
     	m5d_mul_64(c_slice.row[1][0].rows, x1->rows, x2->rows);  //3
        m5d_add_64(x1->rows,a_slice.row[1][0].rows,a_slice.row[1][1].rows) ;  //4
        m5d_sub_64(x2,&b_slice.row[0][1],&b_slice.row[0][0]) ;  //5
        m5d_mul_64(c_slice.row[1][1].rows, x1->rows, x2->rows);    //6
        m5d_sub_64(x1,x1,&a_slice.row[0][0]) ;//5
        m5d_sub_64(x2,&b_slice.row[1][1],x2);  //8
        m5d_mul_64(c_slice.row[0][1].rows,x1->rows,x2->rows); //9
        m5d_sub_64(x1,&a_slice.row[0][1],x1);    //10
        m5d_mul_64(c_slice.row[0][0].rows,x1->rows,b_slice.row[1][1].rows);   //11
        m5d_mul_64(x1->rows, a_slice.row[1][1].rows, b_slice.row[1][1].rows);  //12
        m5d_add_64(c_slice.row[0][1].rows,x1->rows , c_slice.row[0][1].rows) ;   //15
        m5d_add_64(c_slice.row[1][0].rows,c_slice.row[0][1].rows , c_slice.row[1][0].rows) ;   //14
        m5d_add_64(c_slice.row[0][1].rows,c_slice.row[0][1].rows , c_slice.row[1][1].rows) ;   //15
        m5d_add_64(c_slice.row[1][1].rows,c_slice.row[1][0].rows , c_slice.row[1][1].rows) ;    //16
        m5d_add_64(c_slice.row[1][1].rows,c_slice.row[1][0].rows , c_slice.row[1][1].rows) ;  //15
        m5d_sub_64(x2, x2, &b_slice.row[1][0]) ;            //18
        m5d_mul_64(c_slice.row[1][0].rows, a_slice.row[1][1].rows, x2->rows);            //19
        m5d_sub_64(&c_slice.row[1][0], &c_slice.row[1][0],&c_slice.row[0][0]);  //20
        m5d_mul_64(c_slice.row[0][0].rows, a_slice.row[0][1].rows,b_slice.row[1][0].rows);
        m5d_add_64(c_slice.row[0][0].rows, x1->rows,c_slice.row[0][0].rows) ;
    }
    

    m1ri_free(x1);
    m1ri_free(x2);
}

/**
	This handles the arithmetic of m5d_strassen
*/

void  m5d_strassen(m5d_t *c, m5d_t  *a, m5d_t   *b)
{
    if(a->ncols == b->nrows)
    {
        
      	// These hold the padded matrix sizes
		u_int32_t  arcr, acbr, bccc;
		arcr = a->nrows;
		acbr = a->ncols;
		bccc = b->ncols;
	

		arcr =  powerof2(arcr);
		acbr =  powerof2(acbr);
		bccc =  powerof2(bccc);
		
		m5d_create( c, a->nrows, b->ncols); 
		int lasta, lastb, lastboth;
		lasta = 64 - a->nrows%64;
		lastb = 64 -  b->ncols%64;  
		lastboth = 64 - a->nrows; 
	
		
		if((arcr != a->nrows) || (acbr != a->ncols) || (bccc) != (bccc))
		{
		
		m5d_t * padded_a  = m1ri_malloc(sizeof(m5d_t));
		m5d_t  * padded_b  = m1ri_malloc(sizeof(m5d_t));
		m5d_t * padded_c = m1ri_malloc(sizeof(m5d_t));;
		m5d_create(padded_a, arcr, acbr);
		m5d_create(padded_b, acbr, bccc);
		m5d_create(padded_c, arcr, bccc);
		m5d_copypadding(padded_a, a);
		m5d_copypadding(padded_b, b);
		m5d_copypadding(padded_c, c);
		
		
		m5d_qrt_mul(padded_c, padded_a, padded_b); 

		m5d_free(padded_a);
		m5d_free(padded_b);
		m5d_free(padded_c);
		}

		else
		{
		m5d_qrt_mul(c, a, b); 
		}
		

          
    }
    
}


/**
	This handles the arithmetic of m7d_strassen 
	
*/

void m7d_qrt_mul(m7d_t * c,m7d_t *a, m7d_t * b )
{
    m7d_t * x1, *x2;
    
    x1 = m1ri_malloc(sizeof(m7d_t));
    x2 = m1ri_malloc(sizeof(m7d_t));
    m7_slice  a_slice, b_slice, c_slice;
    m7d_create(x1, c->nrows, c->ncols);
    m7d_create(x2, c->nrows, c->ncols);
    m7d_quarter(&a_slice, a);
    m7d_quarter(&b_slice, b);
    m7d_quarter(&c_slice, c);
    
    if((c_slice.row[0][0].ncols) > M1RI_RADIX)
    {
       
        m7d_sub(x1, &a_slice.row[0][0], &a_slice.row[1][0]);  //1
        m7d_sub(x2,&b_slice.row[1][1],&b_slice.row[0][1]);  //2
        m7d_qrt_mul(&c_slice.row[1][0], x1, x2);  //5
        m7d_add_r(x1,&a_slice.row[1][0],&a_slice.row[1][1]);  //4
        m7d_sub(x2,&b_slice.row[0][1],&b_slice.row[0][0]);  //5
        m7d_qrt_mul(&c_slice.row[1][1], x1, x2);    //6
        m7d_sub(x1,x1,&a_slice.row[0][0]);//7
        m7d_sub(x2,&b_slice.row[1][1],x2);  //8
        m7d_qrt_mul(&c_slice.row[0][1],x1,x2); //9
        m7d_sub(x1,&a_slice.row[0][1],x1);    //10
        m7d_qrt_mul(&c_slice.row[0][0],x1,&b_slice.row[1][1]);   //11
        m7d_qrt_mul(x1, &a_slice.row[1][1], &b_slice.row[1][1]);  //12
        m7d_add_r(&c_slice.row[0][1],x1 , &c_slice.row[0][1]);   //15
        m7d_add_r(&c_slice.row[1][0],&c_slice.row[0][1] , &c_slice.row[1][0]);   //14
        m7d_add_r(&c_slice.row[0][1],&c_slice.row[0][1] , &c_slice.row[1][1]);   //15
        m7d_add_r(&c_slice.row[1][1],&c_slice.row[1][0] , &c_slice.row[1][1]);    //16
        m7d_add_r(&c_slice.row[1][1],&c_slice.row[1][0] , &c_slice.row[1][1]);  //17
        m7d_sub(x2, x2, &b_slice.row[1][0]);            //18
        m7d_qrt_mul(&c_slice.row[1][0], &a_slice.row[1][1], x2);            //19
        m7d_sub(&c_slice.row[1][0], &c_slice.row[1][0], &c_slice.row[0][0]);  //20
        m7d_qrt_mul(&c_slice.row[0][0], &a_slice.row[0][1], &b_slice.row[1][0]);
        m7d_add_r(&c_slice.row[0][0], x1,&c_slice.row[0][0] );
        
    }
    
    else if((c_slice.row[0][0].ncols ) == M1RI_RADIX)
    {
		m7d_sub_64(x1->rows, a_slice.row[0][0].rows, a_slice.row[1][0].rows);  //1
        m7d_sub_64(x2->rows,b_slice.row[1][1].rows,b_slice.row[0][1].rows) ;  //2
        m7d_mul_64(c_slice.row[1][0].rows, x1->rows, x2->rows);  //3
        m7d_add_64(x1->rows,a_slice.row[1][0].rows,a_slice.row[1][1].rows) ;  //4
        m7d_sub_64(x2->rows,b_slice.row[0][1].rows,b_slice.row[0][0].rows) ;  //5
        m7d_mul_64(c_slice.row[1][1].rows, x1->rows, x2->rows);    //6
        m7d_sub_64(x1->rows,x1->rows,a_slice.row[0][0].rows) ;//7
        m7d_sub_64(x2->rows,b_slice.row[1][1].rows,x2->rows);  //8
        m7d_mul_64(c_slice.row[0][1].rows,x1->rows,x2->rows); //9
        m7d_sub_64(x1->rows,a_slice.row[0][1].rows,x1->rows);    //10
        m7d_mul_64(c_slice.row[0][0].rows,x1->rows,b_slice.row[1][1].rows);   //11
        m7d_mul_64(x1->rows, a_slice.row[1][1].rows, b_slice.row[1][1].rows);  //12
        m7d_add_64(c_slice.row[0][1].rows,x1->rows , c_slice.row[0][1].rows) ;   //15
        m7d_add_64(c_slice.row[1][0].rows,c_slice.row[0][1].rows , c_slice.row[1][0].rows) ;   //14
        m7d_add_64(c_slice.row[0][1].rows,c_slice.row[0][1].rows , c_slice.row[1][1].rows) ;   //15
        m7d_add_64(c_slice.row[1][1].rows,c_slice.row[1][0].rows , c_slice.row[1][1].rows) ;    //16
        m7d_add_64(c_slice.row[1][1].rows,c_slice.row[1][0].rows , c_slice.row[1][1].rows) ;  //17
        m7d_sub_64(x2->rows, x2->rows, b_slice.row[1][0].rows) ;            //18
        m7d_mul_64(c_slice.row[1][0].rows, a_slice.row[1][1].rows, x2->rows);            //19
        m7d_sub_64(c_slice.row[1][0].rows, c_slice.row[1][0].rows,c_slice.row[0][0].rows);  //20
        m7d_mul_64(c_slice.row[0][0].rows, a_slice.row[0][1].rows,b_slice.row[1][0].rows);
        m7d_add_64(c_slice.row[0][0].rows, x1->rows,c_slice.row[0][0].rows) ;
    }
    
    m1ri_free(x1);
    m1ri_free(x2);
}
/**
	Strassen  algorithm on an m7d_t
*/
void  m7d_strassen(m7d_t *c, m7d_t  *a, m7d_t   *b)
{
    if(a->ncols == b->nrows)
    {
        m7d_create( c, a->nrows, b->ncols); 
      	// These hold the padded matrix sizes
		u_int32_t  arcr, acbr, bccc;
		a->nrows = arcr;
		a->ncols = acbr;
		b->ncols = bccc;

		arcr =  powerof2(arcr);
		acbr =  powerof2(acbr);
		bccc =  powerof2(bccc);
		
		m7d_create( c, a->nrows, b->ncols); 
		int lasta, lastb, lastboth;
		lasta = 64 - a->nrows%64;
		lastb = 64 -  b->ncols%64;  
		lastboth = 64 - a->nrows; 
		
		

		if((arcr != a->nrows) || (acbr != a->ncols) || (bccc) != (bccc))
		{
		m7d_t * padded_a  = m1ri_malloc(sizeof(m7d_t));
		m7d_t  * padded_b  = m1ri_malloc(sizeof(m7d_t));
		m7d_t * padded_c = m1ri_malloc(sizeof(m7d_t));;
		m7d_create(padded_a, arcr, acbr);
		m7d_create(padded_b, acbr, bccc);
		m7d_create(padded_c, arcr, bccc);
		m7d_copypadding(padded_a, a);
		m7d_copypadding(padded_b, b);
		m7d_copypadding(padded_c, c);
		
		
		m7d_qrt_mul(padded_c, padded_a, padded_b); 

		m7d_free(padded_a);
		m7d_free(padded_b);
		m7d_free(padded_c);
		}
		
		
		else
		{
		m7d_qrt_mul(c, a, b); 
		}
		

          
    }
    
}

