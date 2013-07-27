
/*
 Matrix Represenations and basic operations
 TOMAS J. BOOTHBY AND ROBERT W. BRADSHAW "BITSLICING AND THE METHOD OF FOUR
 RUSSIANS OVER L/Users/grudy/Documents/c folders/m1rigf3/m1rigf3/m3d_tests.cARGER FINITE FIELDS"
 
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
#import <math.h>
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

void  m3d_strassen(m3d_t *c, m3d_t  *a, m3d_t   *b)
{
    if(a->ncols == b->nrows)
    {
        m3d_create(c, a->nrows   , b->ncols);
        
        
        if(((a->ncols%m1ri_word) && (b->nrows%m1ri_word)) == 0 )
        {
            
            int t, x,recursions,  slrow, slwidth, rowdiv, coldiv, themax ;
            slrow =   DN(a->nrows, cutoff);
            slwidth = b->width;
            
            t = 64;
            rowdiv = 1;
            
            while (t < slrow) {
                t = t * 2;
                rowdiv++;
                
            }
            
            t = 64;
            coldiv = 1;
            while (t < slwidth) {
                
                
                t = t * 2;
                rowdiv++;
                
            }
            themax = MAX(rowdiv, coldiv);
            m3_slice a_sliced[themax], b_sliced[themax], c_sliced[themax];
            
            if(themax > 1)
            {
                recursions = 1;
                for(x = themax; x > 1 ; x =  DN(x, 2))
                {
                    
                    m3d_slices(&a_sliced[recursions], &a_sliced[recursions -1].row[0][1], x);
                    m3d_slices(&b_sliced[recursions], &b_sliced[recursions -1].row[0][1], x);
                    m3d_slices(&c_sliced[recursions], &c_sliced[recursions -1].row[0][1], x);
                    
                    m3d_slices(&a_sliced[recursions], &a_sliced[recursions -1].row[1][1], x);
                    m3d_slices(&b_sliced[recursions], &b_sliced[recursions -1].row[1][1], x);
                    m3d_slices(&c_sliced[recursions], &c_sliced[recursions -1].row[1][1], x);
                    
                    m3d_slices(&a_sliced[recursions], &a_sliced[recursions -1].row[0][1], x);
                    m3d_slices(&b_sliced[recursions], &b_sliced[recursions -1].row[0][1], x);
                    m3d_slices(&c_sliced[recursions], &c_sliced[recursions -1].row[0][1], x);
                    
                    m3d_slices(&a_sliced[recursions], &a_sliced[recursions -1].row[0][1], x);
                    m3d_slices(&b_sliced[recursions], &b_sliced[recursions -1].row[0][1], x);
                    m3d_slices(&c_sliced[recursions], &c_sliced[recursions -1].row[0][1], x);
                    recursions ++;
                }
                
                if(themax == 1)
                {
                    mul_64_m3d(c->rows, a->rows, b->rows);
                    
                }
                
                
            }
            
            
            
        }
        
        
        
    }
    
    
}

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

void  m3d_strassen_window_directly(m3d_t *c, m3d_t  *a, m3d_t   *b)
{
    if(a->ncols == b->nrows)
    {
        int i;
        bool rowsextra, colsextra;
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





