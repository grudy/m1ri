
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
void m3d_strassen_total16(m3_slice *c, m3_slice const *a, m3_slice const *b)
{
    m3d_t * x1 = m1ri_calloc(64, sizeof(m3d_t));
    m3d_t * x2 = m1ri_calloc(64, sizeof(m3d_t));
    m3d_create(x2, 64, 64);
    m3d_create(x1,64,  64);
    m3d_sub(x1, &a->row[0][0], &a->row[1][0]);
    m3d_sub(x2,&b->row[1][1],&b->row[0][1]);
    mul_64_m3d(c->row[1][0].rows ,x1->rows,x2->rows);
    
    
    
    
    add_64_m3d(x1->rows,a->row[1][0].rows,a->row[1][1].rows);
    m3d_sub(x2,&b->row[0][1],&b->row[0][0]);
    mul_64_m3d(c->row[1][1].rows,x1->rows,x2->rows);
    
    m3d_sub(x1,x1,&a->row[0][0]);
    m3d_sub(x2,&b->row[1][1],x2);
    mul_64_m3d(c->row[0][1].rows,x1->rows,x2->rows);
    
    
    
    
    m3d_sub(x1,&a->row[0][1],x1);
    mul_64_m3d(c->row[0][0].rows,x1->rows,b->row[1][1].rows);
    mul_64_m3d(x1->rows,a->row[0][0].rows,b->row[0][0].rows);
    
    add_64_m3d(c->row[0][1].rows,x1->rows,c->row[0][1].rows);
    add_64_m3d(c->row[1][0].rows,c->row[0][1].rows,c->row[1][0].rows);
    add_64_m3d(c->row[0][1].rows,c->row[0][1].rows,c->row[1][1].rows);
    add_64_m3d(c->row[1][1].rows,c->row[1][0].rows,c->row[1][1].rows);
    add_64_m3d(c->row[0][1].rows,c->row[0][1].rows,c->row[0][0].rows);
    
    m3d_sub(x2,x2,&b->row[1][0]);
    mul_64_m3d(c->row[0][0].rows,a->row[1][1].rows,x2->rows);
    
    m3d_sub(&c->row[1][0],&c->row[1][0],&c->row[0][0]);
    mul_64_m3d(c->row[0][0].rows,a->row[0][1].rows,b->row[1][0].rows);
    add_64_m3d(c->row[0][0].rows,x1->rows,c->row[0][0].rows);
    
    //
    m1ri_free(x1);
    m1ri_free(x2);
    
    
    
}


//Strassen winograd arithmatic on slices of size row
void m3d_mul_slicerow(m3_slice *  c, m3_slice  *  a, m3_slice *   b, rci_t * rownum )
{
    int colnum;
    m3d_t * x1 = m1ri_calloc(64, sizeof(m3d_t));
    m3d_t * x2 = m1ri_calloc(64, sizeof(m3d_t));
    m3d_create(x2, 64, 64);
    m3d_create(x1,64,  64);
    for( colnum = 0; colnum < c->nrows; colnum = colnum +2)
    {
        
        
        m3d_sub(x1, &a->row[*rownum][colnum], &a->row[*rownum + 1][colnum]);
        m3d_sub(x2,&b->row[*rownum + 1][colnum + 1],&b->row[*rownum][colnum + 1]);
        mul_64_m3d(c->row[*rownum + 1][colnum].rows ,x1->rows,x2->rows);
        
        
        
        
        add_64_m3d(x1->rows,a->row[*rownum + 1][colnum].rows,a->row[*rownum + 1][colnum + 1].rows);
        m3d_sub(x2,&b->row[*rownum][colnum + 1],&b->row[*rownum][colnum]);
        mul_64_m3d(c->row[*rownum + 1][colnum + 1].rows,x1->rows,x2->rows);
        
        m3d_sub(x1,x1,&a->row[*rownum][colnum]);
        m3d_sub(x2,&b->row[*rownum + 1][colnum + 1],x2);
        mul_64_m3d(c->row[*rownum][colnum + 1].rows,x1->rows,x2->rows);
        
        
        
        
        m3d_sub(x1,&a->row[*rownum][colnum + 1],x1);
        mul_64_m3d(c->row[*rownum][colnum].rows,x1->rows,b->row[*rownum + 1][colnum + 1].rows);
        mul_64_m3d(x1->rows,a->row[*rownum][colnum].rows,b->row[*rownum][colnum].rows);
        
        add_64_m3d(c->row[*rownum][colnum + 1].rows,x1->rows,c->row[*rownum][colnum + 1].rows);
        add_64_m3d(c->row[*rownum + 1][colnum].rows,c->row[*rownum][colnum + 1].rows,c->row[*rownum + 1][colnum].rows);
        add_64_m3d(c->row[*rownum][colnum + 1].rows,c->row[*rownum][colnum + 1].rows,c->row[*rownum + 1][colnum + 1].rows);
        add_64_m3d(c->row[*rownum + 1][colnum + 1].rows,c->row[*rownum + 1][colnum].rows,c->row[*rownum + 1][colnum + 1].rows);
        add_64_m3d(c->row[*rownum][colnum + 1].rows,c->row[*rownum][colnum + 1].rows,c->row[*rownum][colnum].rows);
        
        m3d_sub(x2,x2,&b->row[*rownum + 1][colnum]);
        mul_64_m3d(c->row[*rownum][colnum].rows,a->row[*rownum + 1][colnum + 1].rows,x2->rows);
        
        m3d_sub(&c->row[*rownum + 1][colnum],&c->row[*rownum + 1][colnum],&c->row[*rownum][colnum]);
        mul_64_m3d(c->row[*rownum][colnum].rows,a->row[*rownum][colnum + 1].rows,b->row[*rownum + 1][colnum].rows);
        add_64_m3d(c->row[*rownum][colnum].rows,x1->rows,c->row[*rownum][colnum].rows);
        
    }
    
}


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
    
    if((c_slice.row[0][0].ncols) > cutoff)
    {
       
        
        m3d_qrt_mul(&c_slice.row[0][0], &b_slice.row[0][0], &a_slice.row[0][0]);
        m3d_qrt_mul(&c_slice.row[0][1], &b_slice.row[0][1], &a_slice.row[0][1]);
        m3d_qrt_mul(&c_slice.row[1][0], &b_slice.row[1][0], &a_slice.row[1][0]);
        m3d_qrt_mul(&c_slice.row[1][1], &b_slice.row[1][1], &a_slice.row[1][1]);
        
        
        
    
    }
    
    else if((c_slice.row[0][0].ncols ) <= cutoff)
    {
        
        m3d_strassen_total16(&c_slice, &a_slice, &c_slice);
        /*
        m3d_sub(x1, &a_slice.row[0][0], &a_slice.row[1][0]);
        m3d_sub(x2,&b_slice.row[1][1],&b_slice.row[0][1]);
        m3d_qrt_mul(&c_slice.row[1][0], x1, x2);
        m3d_add_r(x1,&a_slice.row[1][0],&a_slice.row[1][1]);
        m3d_sub(x2,&b_slice.row[0][1],&b_slice.row[0][0]);  //5
        m3d_qrt_mul(&c_slice.row[1][0], x1, x2);    //6
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
        m3d_sub(&c_slice.row[0][0], &c_slice.row[1][0], &c_slice.row[0][0]);  //20
        m3d_qrt_mul(&c_slice.row[0][0], &a_slice.row[0][1], &b_slice.row[1][0]);
        m3d_add_r(&c_slice.row[0][0], x1,&c_slice.row[0][0] );
        */
     
    }
    
    
    
    
}



void  m3d_strassen_window_directly(m3d_t *c, m3d_t  *a, m3d_t   *b)
{
    if(a->ncols == b->nrows)
    {
        int i;
        int rowsextra, colsextra;
        m3_slice  * a_slices ,  *  b_slices  ,  * c_slices ;
        a_slices = m1ri_malloc(sizeof(m3_slice));
        b_slices = m1ri_malloc(sizeof(m3_slice));
        c_slices = m1ri_malloc(sizeof(m3_slice));
        m3d_create(c, a->nrows   , b->ncols);
        m3d_slices(a_slices, a, 1);
        m3d_slices(b_slices, b, 1);
        m3d_slices(c_slices, c, 1);
        rowsextra = c->nrows%2;
        colsextra = c->width%2;
        
        for( i = 0;i < c_slices->nrows -rowsextra; i = i + 2)
        {
            m3d_mul_slicerow( c_slices, b_slices, a_slices, &i);
         
            
            
            
        }
        
     
        if(rowsextra)
        {
            
            
            
        }
        
        
    }
    
    
    
    
}


void  m3d_strassen(m3d_t *c, m3d_t  *a, m3d_t   *b)
{
    if(a->ncols == b->nrows)
    {
       // int maxdiv, divfac, ln_cut, slicecount;
 
        m3d_create(c, a->nrows   , b->ncols);
       /* if((c->nrows  >= 128)&& (c->width  >= 2))
        {
            int ecol, erow;
            erow = a->nrows%
            
            while()
            m3d_t * a_fwindow, * b_fwindow, * c_fwindow;
            m3d_window_create(a, a_fwindow, 0, 0, a->nrows - , <#rci_t#>)
            
        
        
        }*/
    
        
        
        
        
        
    }
    
    
    
}




