
/*
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


void m3d_qrt_mul(m3d_t * c,m3d_t *a, m3d_t * b )
{
    m3d_t * x1, *x2;
    x1 = x2 = m1ri_malloc(sizeof(m3d_t));
    m3_slice  a_slice, b_slice, c_slice;
    m3d_create(x1, c->nrows, c->ncols);
    m3d_create(x2, c->nrows, c->ncols);
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
        
       
        

    
     	
    }
    
}




void  m3d_strassen(m3d_t *c, m3d_t  *a, m3d_t   *b)
{
    if(a->ncols == b->nrows)
    {
      		// These hold the padded matrix sizes
		u_int32_t  arcr, acbr, bccc;
		a->nrows = arcr;
		a->ncols = acbr;
		b->ncols = bccc;

		arcr =  powerof2(arcr);
		acbr =  powerof2(acbr);
		bccc =  powerof2(bccc);
		
		m3d_create( c, a->nrows, b->ncols); 
		int lasta, lastb, lastboth;
		lasta = 64 - a->nrows%64;
		lastb = 64 -  b->ncols%64;  
		lastboth = 64 - a->nrows; 
		
		
		m3d_t * c_main_partition = m1ri_malloc(sizeof(m3d_t));
		m3d_t * a_main_partition = m1ri_malloc(sizeof(m3d_t));
		m3d_t * b_main_partition = m1ri_malloc(sizeof(m3d_t));
		
		
		m3d_window_create(c, c_main_partition, 0, 0, c->nrows -1, c->ncols -1); 
		m3d_window_create(a, a_main_partition, 0, 0, a->nrows -1, a->ncols -1   ); 
		m3d_window_create(b, b_main_partition,  0, 0,b->nrows -1 , b->ncols -1); 
		
		
		
		if((arcr != a->nrows) || (acbr != a->ncols) || (bccc) != (bccc))
		{
		m3d_t * padded_a  = m1ri_malloc(sizeof(m3d_t));
		m3d_t  * padded_b  = m1ri_malloc(sizeof(m3d_t));
		m3d_t * padded_c = m1ri_malloc(sizeof(m3d_t));;
		m3d_create(padded_a, arcr, acbr);
		m3d_create(padded_b, acbr, bccc);
		m3d_create(padded_c, arcr, bccc);
		m3d_copypadding(padded_a, a_main_partition);
		m3d_copypadding(padded_b, b_main_partition);
		m3d_copypadding(padded_c, c_main_partition);
		
		
		m3d_qrt_mul(padded_c, padded_a, padded_b); 

		m3d_free(padded_a);
		m3d_free(padded_b);
		m3d_free(padded_c);
		}
		
		//m3d_create(padded_c
		
		else
		{
		m3d_qrt_mul(c_main_partition, a_main_partition, b_main_partition); 
		}
		

		m3d_free(a_main_partition);
		m3d_free(b_main_partition);
		m3d_free(c_main_partition);
        
           
    }
    
}




