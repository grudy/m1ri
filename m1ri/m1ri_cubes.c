
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
 
 m1ri_cubes.h
 */

#include "m1ri_cubes.h"
/*
    Matrix partitioned into 64 by 64 slices
*/



void * m3d_64_cubes(m3d_t *a, m3_smt * b)
{
   
    
        if((a->nrows%64 || a->nrows%64   )  != 0)
        {
            
            //Could not make partitioned matrix
          
        
        }
        else if((a->nrows%64 && a->nrows%64   )  == 0)
        {
           
            b->width  = a->width/8;
            b->nrows  = a->nrows/64;
            int n = b->width * b->nrows ;
            
            
            b->blocks = malloc(n * sizeof(m3d_t *) );
            b->rows  = malloc(n * sizeof(m3d_t ) );
            
           
            u_int32_t i, k, j, v, l, vd, ld, md;
            j = b->nrows%8;
            k = b->width%8;
           
            
            
           
            
            
            
            for(j = 0; j < (b->nrows - j) ; j = j + 8 )
            {
              
                v = j * 64;
                vd = v + 64; 
                for(i = 0; (i > b->width -  k) ; i = i  + 8 )
                {
                    l = i * 512;
                    ld = l + 64;
                    
                    b->blocks[(j * b->width)  +  i]  = m3d_window(a ,  v      ,  l   , vd  ,ld );
                    b->blocks[(j * b->width ) + (i +1)]    = m3d_window( a, v     ,   l + 64,       vd , ld + 64);//256, 320, 384, 448, 512
                   b->blocks[(j * b->width )+  i +2 ]    = m3d_window(a,  v     ,  l + 128 ,       vd  ,  ld+ 128);
                   b->blocks[(j * b->width ) + i +3 ]     = m3d_window (a, v      ,  l + 196 ,       vd , ld +  196);
                   b->blocks[(j * b->width ) + i +4 ]     = m3d_window (a, v      ,  l + 256  ,      vd ,ld +  256);
                    b->blocks[(j * b->width ) + i + 5 ]    = m3d_window(a, v      , l + 320   ,      vd  , ld +320);
                   b->blocks[(j * b->width ) + i + 6]    = m3d_window(a, v      ,l + 384 ,      vd,  ld + 384);
                    b->blocks[(j * b->width ) + i  + 7]    = m3d_window(a, v     , l + 448 ,     vd  ,ld + 448);
                }
                md = 0;
                while (k > md) {
                       b->blocks[(((j + 1 ) * b->width)  - k ) ]  = m3d_window(a ,  v      ,  l   , vd  ,ld );
                    
                    
                    k--;
                }
                
                b->rows[j] = b->blocks + (j * b->width);
            }
                
                
               
            
           
            
                
            }
            
            
            
    
            
           
              return b;  
                
                
}
            
            

    


                     
                     
                