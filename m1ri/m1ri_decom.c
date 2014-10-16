 
/** *
 
 TOMAS J. BOOTHBY AND ROBERT W. BRADSHAW "BITSLICING AND THE METHOD OF FOUR
 RUSSIANS OVER LARGER FINITE FIELDS"
  Rank-profile revealing Gaussian elimination and the CUP matrix decomposition
  Claude-Pierre Jeannerod, Clément Pernet, Arne Storjohann
 
 http://arxiv.org/abs/1112.5717
 
 
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
 
 
 Matrix decompositions


 m3d_decom.c
 
 */
#include <m1ri/m1ri_decom.h>
 
 

void   m3d_transpose_64(vbg  **b, vbg  **a  )
{

    int i, x;
    vbg temp;
    for(i = 0; i <64; i ++)
    {
        for(x = 0; x < 64; x ++)
        {
        	temp.units =  (a[x][0].units & (rightbit << i) );
        	temp.sign =  (a[x][0].sign & (rightbit << i ) );
       	   	b[i][0].units = (temp.units) ?  b[i][0].units | (rightbit << x) : b[i][0].units ;
       		b[i][0].sign = (temp.sign) ? b[i][0].sign | (rightbit << x) : b[i][0].sign  ;
        }
    }
 



}


static inline void m3d_transpose_il(m3d_t * c, m3d_t  const *  a )
{

    
	

	if((c->ncols) > (M1RI_RADIX << 1))
    {
     	m3_slice *  a_slice, *  c_slice;    	
 
    	a_slice = m3d_quarter( a);
    	c_slice = m3d_quarter(c);

    	m3d_transpose_il( c_slice->row[0], a_slice->row[0]); 
    	m3d_transpose_il( c_slice->row[1], a_slice->row[2]);  
    	m3d_transpose_il( c_slice->row[2], a_slice->row[1]);
    	m3d_transpose_il( c_slice->row[3], a_slice->row[3]);  
  

 
    	
    	
    	m3d_quarter_free(a_slice);
    	m3d_quarter_free(c_slice);
       	
       	
    }
    
    else if((c->ncols ) == (M1RI_RADIX  << 1))
    {
    	m3_slice *  a_slice, *  c_slice;    	
 
   		 a_slice = m3d_quarter( a);
   		 c_slice = m3d_quarter(c);

	  	m3d_transpose_64( c_slice->row[0]->rows, a_slice->row[0]->rows); 
    	m3d_transpose_64( c_slice->row[1]->rows, a_slice->row[2]->rows);  
    	m3d_transpose_64( c_slice->row[2]->rows, a_slice->row[1]->rows);
    	m3d_transpose_64( c_slice->row[3]->rows, a_slice->row[3]->rows);  
  	   	
   	 	m3d_quarter_free(a_slice);
    	m3d_quarter_free(c_slice);
    }
    
    else if(c->ncols  == M1RI_RADIX )
    {
    	 m3d_transpose_64(c->rows, a->rows);
    
    } 
    
    
   
     
}


 

/**
  /brief transposing  arbitrary algorithm on an m3d_t	
*/
m3d_t *  m3d_transpose(m3d_t *c, m3d_t  const *a)
{
	if (c == NULL)
	{
		c = m3d_create(a->ncols, a->nrows);

	} 

	else 
	{
		if (c->nrows != a->ncols || c->ncols != a->nrows) 
		{
			m1ri_die("m3d_transpose: Provided return matrix has wrong dimensions.\n");
    	}
	
	}
	
    
    
 
      	/*  These hold the padded matrix sizes */
    	
	u_int32_t  arcr, acbr, bccc;
	arcr = a->nrows;
	acbr = a->ncols;

    
	arcr =  powerof2(arcr);
	acbr =  powerof2(acbr);
	
	arcr = MAX(arcr, 64);
	acbr = MAX(acbr, 64);
	
	arcr = MAX(arcr, acbr);

	
	if((arcr != a->nrows) || (acbr != a->ncols)  )
	{
		m3d_t * padded_a,   * padded_c;
	
	
		padded_a = m3d_create(arcr, arcr);
		padded_c = m3d_create(arcr, arcr);
		padded_a = m3d_copy(padded_a, a);

	
	
		m3d_transpose_il(padded_c, padded_a); 
		c  = m3d_copy_cutoff(c, padded_c);
		m3d_free(padded_a);
		m3d_free(padded_c);
		
		
	}
	
	
	
	else
	{
	
		m3d_transpose_il(c, a); 
		
	}
	

          
    
    
    return c;
    
}



void   m5d_transpose_64(vfd  **b, vfd  **a  )
{

    int i, x;
    vfd temp;
    for(i = 0; i <64; i ++)
    {
        for(x = 0; x < 64; x ++)
        {
 			temp.units =  (a[x][0].units & (rightbit << i) );
        	temp.middle =  (a[x][0].middle & (rightbit << i ) );
        	temp.sign =  (a[x][0].sign & (rightbit << i ) );

       	   	b[i][0].units = (temp.units) ?  b[i][0].units | (rightbit << x) : b[i][0].units ;
       	   	b[i][0].middle = (temp.middle) ?  b[i][0].middle | (rightbit << x) : b[i][0].middle ;
			b[i][0].sign = (temp.sign) ? b[i][0].sign | (rightbit << x) : b[i][0].sign  ;
        }
    }
 



}


static inline void m5d_transpose_il(m5d_t * c, m5d_t  const *  a )
{


	if((c->ncols) > (M1RI_RADIX << 1))
    {
     	m5_slice *  a_slice, *  c_slice;    	
 
    	a_slice = m5d_quarter( a);
    	c_slice = m5d_quarter(c);

    	m5d_transpose_il( c_slice->row[0], a_slice->row[0]); 
    	m5d_transpose_il( c_slice->row[1], a_slice->row[2]);  
    	m5d_transpose_il( c_slice->row[2], a_slice->row[1]);
    	m5d_transpose_il( c_slice->row[3], a_slice->row[3]);  
  

 
    	
    	
    	m5d_quarter_free(a_slice);
    	m5d_quarter_free(c_slice);
       	
       	
    }
    
    else if((c->ncols ) == (M1RI_RADIX  << 1))
    {
    	m5_slice *  a_slice, *  c_slice;    	
 
   		 a_slice = m5d_quarter( a);
   		 c_slice = m5d_quarter(c);

	  	m5d_transpose_64( c_slice->row[0]->rows, a_slice->row[0]->rows); 
    	m5d_transpose_64( c_slice->row[1]->rows, a_slice->row[2]->rows);  
    	m5d_transpose_64( c_slice->row[2]->rows, a_slice->row[1]->rows);
    	m5d_transpose_64( c_slice->row[3]->rows, a_slice->row[3]->rows);  
  	   	
   	 	m5d_quarter_free(a_slice);
    	m5d_quarter_free(c_slice);
    }
    
    else if(c->ncols  == M1RI_RADIX )
    {
    	 m5d_transpose_64(c->rows, a->rows);
    
    } 
    
    
   
     
}


 

/**
  /brief transposing  arbitrary algorithm on an m5d_t	
*/
m5d_t *  m5d_transpose(m5d_t *c, m5d_t  const *a)
{
	if (c == NULL)
	{
		c = m5d_create(a->ncols, a->nrows);

	} 

	else 
	{
		if (c->nrows != a->ncols || c->ncols != a->nrows) 
		{
			m1ri_die("m5d_transpose: Provided return matrix has wrong dimensions.\n");
    	}
	
	}
	
    
    
 
      	/*  These hold the padded matrix sizes */
    	
	u_int32_t  arcr, acbr, bccc;
	arcr = a->nrows;
	acbr = a->ncols;

    
	arcr =  powerof2(arcr);
	acbr =  powerof2(acbr);
	
	arcr = MAX(arcr, 64);
	acbr = MAX(acbr, 64);
	
	arcr = MAX(arcr, acbr);

	
	if((arcr != a->nrows) || (acbr != a->ncols)  )
	{
		m5d_t * padded_a,   * padded_c;
	
	
		padded_a = m5d_create(arcr, arcr);
		padded_c = m5d_create(arcr, arcr);
		padded_a = m5d_copy(padded_a, a);

	
	
		m5d_transpose_il(padded_c, padded_a); 
		c  = m5d_copy_cutoff(c, padded_c);
		m5d_free(padded_a);
		m5d_free(padded_c);
		
		
	}
	
	else
	{
	
		m5d_transpose_il(c, a); 
		
	}
	
    
    return c;
    
}

void   m7d_transpose_64(vtri  **b, vtri  **a  )
{

    int i, x;
    vtri temp;
    for(i = 0; i <64; i ++)
    {
        for(x = 0; x < 64; x ++)
        {
        	temp.units =  (a[x][0].units & (rightbit << i) );
        	temp.middle =  (a[x][0].middle & (rightbit << i ) );
        	temp.sign =  (a[x][0].sign & (rightbit << i ) );

       	   	b[i][0].units = (temp.units) ?  b[i][0].units | (rightbit << x) : b[i][0].units ;
       	   	b[i][0].middle = (temp.middle) ?  b[i][0].middle | (rightbit << x) : b[i][0].middle ;
			b[i][0].sign = (temp.sign) ? b[i][0].sign | (rightbit << x) : b[i][0].sign  ;
        }
    }
 



}


static inline void m7d_transpose_il(m7d_t * c, m7d_t  const *  a )
{

    
	

	if((c->ncols) > (M1RI_RADIX << 1))
    {
     	m7_slice *  a_slice, *  c_slice;    	
 
    	a_slice = m7d_quarter( a);
    	c_slice = m7d_quarter(c);

    	m7d_transpose_il( c_slice->row[0], a_slice->row[0]); 
    	m7d_transpose_il( c_slice->row[1], a_slice->row[2]);  
    	m7d_transpose_il( c_slice->row[2], a_slice->row[1]);
    	m7d_transpose_il( c_slice->row[3], a_slice->row[3]);  
  

 
    	
    	
    	m7d_quarter_free(a_slice);
    	m7d_quarter_free(c_slice);
       	
       	
    }
    
    else if((c->ncols ) == (M1RI_RADIX  << 1))
    {
    	m7_slice *  a_slice, *  c_slice;    	
 
   		 a_slice = m7d_quarter( a);
   		 c_slice = m7d_quarter(c);

	  	m7d_transpose_64( c_slice->row[0]->rows, a_slice->row[0]->rows); 
    	m7d_transpose_64( c_slice->row[1]->rows, a_slice->row[2]->rows);  
    	m7d_transpose_64( c_slice->row[2]->rows, a_slice->row[1]->rows);
    	m7d_transpose_64( c_slice->row[3]->rows, a_slice->row[3]->rows);  
  	   	
   	 	m7d_quarter_free(a_slice);
    	m7d_quarter_free(c_slice);
    }
    
    else if(c->ncols  == M1RI_RADIX )
    {
    	 m7d_transpose_64(c->rows, a->rows);
    
    } 
    
    
   
     
}


 

/**
  /brief transposing  arbitrary algorithm on an m7d_t	
*/
m7d_t *  m7d_transpose(m7d_t *c, m7d_t  const *a)
{
	if (c == NULL)
	{
		c = m7d_create(a->ncols, a->nrows);

	} 

	else 
	{
		if (c->nrows != a->ncols || c->ncols != a->nrows) 
		{
			m1ri_die("m7d_transpose: Provided return matrix has wrong dimensions.\n");
    	}
	
	}
	
    
    
 
      	/*  These hold the padded matrix sizes */
    	
	u_int32_t  arcr, acbr, bccc;
	arcr = a->nrows;
	acbr = a->ncols;

    
	arcr =  powerof2(arcr);
	acbr =  powerof2(acbr);
	
	arcr = MAX(arcr, 64);
	acbr = MAX(acbr, 64);
	
	arcr = MAX(arcr, acbr);

	
	if((arcr != a->nrows) || (acbr != a->ncols)  )
	{
		m7d_t * padded_a,   * padded_c;
	
	
		padded_a = m7d_create(arcr, arcr);
		padded_c = m7d_create(arcr, arcr);
		padded_a = m7d_copy(padded_a, a);

	
	
		m7d_transpose_il(padded_c, padded_a); 
		c  = m7d_copy_cutoff(c, padded_c);

		
		
	}
	
	
	
	else
	{
	
		m7d_transpose_il(c, a); 
		
	}
	

          
    
    
    return c;
    
}


#include "m1ri.h"


void m3d_tsrm_ur(m3d_t * L, m3d_t B)
{  


}
void m5d_tsrm_ur(m5d_t * L, m5d_t B)
{  


}
void m7d_tsrm_ur(m7d_t * L, m7d_t * B)
{  


}
 
 
 
 
 
 /*
   Solves U X = B with X and B matrices and U upper left triangular
   X replaces B
 */

void m3d_tsrm_ul(m3d_t const *U, m3d_t *B  )
{  


}
void m5d_tsrm_ul(m5d_t const *U, m5d_t *B )
{  


}
void m7d_tsrm_ul(m7d_t const *U, m7d_t *B)
{  


}

 
 /*
   Solves U X = B with X and B matrices and U lower lefttriangular
   X replaces B
 */


void m3d_tsrm_ll(m3d_t const *U, m3d_t *B)
{  


}
void m5d_tsrm_ll(m5d_t const *U, m5d_t *B)
{  


}
void m7d_tsrm_ll(m7d_t const *U, m7d_t *B)
{  


}
 
  /*
   Solves U X = B with X and B matrices and U  lower right triangular
   X replaces B
 */
void m3d_tsrm_lr(m3d_t const *U, m3d_t *B )
{  


}
void m5d_lower_r_triangular(m5d_t const *U, m5d_t *B)
{  


} 
void m7d_tsrm_lr(m7d_t const *U, m7d_t *B)
{  
  

}
 
 
/*
   Invert triangular matrix a
 
*/
 
 void m3d_trtri_upper(m3d_t * a)
{  
  

}
 void m5d_trtri_upper(m5d_t * a)
{  


}
 void m7d_trtri_upper(m7d_t * a)
{  


}
 

