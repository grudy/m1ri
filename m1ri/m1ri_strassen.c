
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

void  m3d_strassen_(m3d_t *c, m3d_t *a, m3d_t*b)
{
   if(a->ncols == b->nrows)
    {
        m3d_create(c, a->nrows   , b->ncols);
        
        
     if(((a->ncols%m1ri_word) && (b->nrows%m1ri_word)) == 0 )
        {
          /*  m3_slice * a_slices = malloc(sizeof(m3_slice));
            m3_slice * b_slices = malloc(sizeof(m3_slice));
            int strop, brush, cutoff;
            
             strop =   a->width;
             brush =   DN(a->nrows, 64);
             cutoff = 2;  //When to stop slicing
            
            
            m3d_slices(a_slices , a, 2);
            m3d_slices(b_slices, b, 2);
            
            */
            
            //a->nrows;
            
            
            //m3_slice a = m3d_slices(&b, <#m3d_t *#>, <#wi_t#>)
            
            
            
            ///even width square 64 * 64
            if(a->ncols == m1ri_word)
            {
                             
                
            }
            
        
                
                
             ///even width smt
          
            }
        
        
        }
        
     if((a->ncols%m1ri_word) && (b->nrows%m1ri_word))
     {
     
     
     
     
     
     }
      
    
        
        
        
        
        
     else
     {
        //can't be multiplied
    
    
    }
    
    
    
    
    /*
    m3d_qrt cs =  m3d_qtrwindows(C);
    m3d_qrt bs =  m3d_qtrwindows(B);
    m3d_qrt as =  m3d_qtrwindows(A);
    
   */

}



/*
 Naive multiplication
*/





/*
void  m3d_mul_naive(m3d_t *c, m3d_t *a, m3d_t*b)
{

    if (a->ncols == b->nrows)
    {
        int i, x,y,  lastvecs;
        m3d_create(c, a->nrows, b->ncols);
        lastvecs = c->ncols%m1ri_word;
        for(i = 0; i < c->nrows; i++)
        {
            
            for(x = 0;x  < a->width; x++)
            {
                    
                for(y = 0; y < 64; y++)
                {
                    
                
                   c->rows[i][x] =  a->rows[i][  b->rows[x + y][i];
                    
                    
                    
                    
                    
                    
                }
            
            
            
            
            
            }
            
            
            
            
        
        
        
        }
    
    
    
    
    
    
    
    }
        
        





}

*/








