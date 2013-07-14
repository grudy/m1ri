
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



void * m3d_64_cubes(m3_smt * b, m3d_t  *a )
{
   
    
    if((a->nrows%m1ri_word || a->ncols%m1ri_word   )  != 0)
    {
        
        b->width  = a->width;
        b->nrows  = DN(a->width, 64);
        int n = b->width * b->nrows ;
        b->blocks = malloc(n * sizeof(m3d_t ) );
        b->rows  = malloc(n * sizeof(m3d_t *) );
        u_int32_t  k, j;
        j = b->nrows%8;
        k = b->width%8;
        /*
         The last column and row in the structure
         */
        signed short rlast, clast;
        rlast = 63 - a->nrows%m1ri_word;
        clast = 63 - a->ncols%m1ri_word;
            
    for(k  = 0; k < (b->nrows - 1) ; k ++ )
        {
            
            for(j = 0; j < (b->width - 1) ; j++  )
            {
                    
                    
            b->blocks[(k*b->width) + j ] = m3d_window(a, (k * m1ri_word), (j * m1ri_word), ((k * m1ri_word) + 63), ((j * m1ri_word) + 63));
                    
                    
            }
                b->blocks[(k*b->width) + j ] = m3d_window(a, (k * m1ri_word), (j * m1ri_word), ((k * m1ri_word) + 63), ((j * m1ri_word) + clast));
                b->rows[j] = b->blocks + (k * b->width);
        }
        
        while(j < b->width )
        {
                    
        b->blocks[(k*b->width) + j ] = m3d_window(a, (k * m1ri_word), (j * m1ri_word), ((k * m1ri_word) + rlast), ((j * m1ri_word) + 63));
            j++;
        }
        b->blocks[(k*b->width) + j ] = m3d_window(a, (k * m1ri_word), (j * m1ri_word), ((k * m1ri_word) + rlast), ((j * m1ri_word) + clast));
        b->rows[j] = b->blocks + (k * b->width);
                
        
    
            
        
        }
    else if((a->nrows%m1ri_word && a->ncols%m1ri_word   )  == 0)
        {
           
            b->width  = a->width;
            b->nrows  = a->nrows/64;
            int n = b->width * b->nrows ;
            
            
            b->blocks = malloc(n * sizeof(m3d_t ) );
            b->rows  = malloc(n * sizeof(m3d_t *) );
            
           
            u_int32_t  k, j;
            j = b->nrows%8;
            k = b->width%8;
            
        for(k  = 0; k < b->nrows; k ++ )
        {
        
        for(j = 0; j < (b->width) ; j ++ )
        {
                
            b->blocks[(k*b->width) + j ] = m3d_window(a, (k * m1ri_word), (j * m1ri_word), ((k * m1ri_word) + 63), ((j * m1ri_word) + 63));
            
        }
                   b->rows[j] = b->blocks + (k * b->width);
            
        }
            
    }
            return b;  
                
                
}
            



m3d_t m3d_transpose(m3d_t   * a)
{
    m3_smt *b = malloc(sizeof(m3_smt));
    m3d_64_cubes(b, a);
    int i, j, k;
    
    for(i = 0; i < b->nrows; i++)
    {
        
        for(j = 0; j < b->width ; j = j + 1)
        {
            
            for(k = 0; k < 64; k =  k + 8)
            {
                
                // b.rows[i][j].units =   b.rows[i][j].units | ((a->rows[j][i].units  & (leftbit >> (k + 3))));
                b->rows[i][j].rows[0][0].units =    b->rows[i][j].rows[0][0].units  | (a->rows[j][i].units  & (leftbit >> (k )));
                b->rows[i][j].rows[0][0].units  =    b->rows[i][j].rows[0][0].units  | (a->rows[j][i].units & (leftbit >> (k +1))) ;
                b->rows[i][j].rows[0][0].units  =    b->rows[i][j].rows[0][0].units  | ((a->rows[j][i].units & (leftbit >> (k + 2)))) ;
                b->rows[i][j].rows[0][0].units  =    b->rows[i][j].rows[0][0].units  | ((a->rows[j][i].units  & (leftbit >> (k + 3))));
                b->rows[i][j].rows[0][0].units  =    b->rows[i][j].rows[ 0][0].units  | (a->rows[j][i].units & (leftbit >> (k + 4))) ;
                b->rows[i][j].rows[0][0].units  =    b->rows[i][j].rows[0][0].units  | ((a->rows[j][i].units & (leftbit >> (k + 5)))) ;
                b->rows[i][j].rows[0][0].units  =    b->rows[i][j].rows[0][0].units  | ((a->rows[j][i].units & (leftbit >> (k + 6)))) ;
                b->rows[i][j].rows[0][0].units  =    b->rows[i][j].rows[0][0].units  | ((a->rows[j][i].units & (leftbit >> (k + 7)))) ;
                b->rows[i][j].rows[0][0].sign =   b->rows[i][j].rows[0][0].sign  | (a->rows[j][i].sign & (leftbit >> k)) ;
                b->rows[i][j].rows[0][0].sign =   b->rows[i][j].rows[0][0].sign  | ((a->rows[j][i].sign & (leftbit >> (k + 1))));
                b->rows[i][j].rows[0][0].sign =   b->rows[i][j].rows[0][0].sign  | (a->rows[j][i].sign  & (leftbit >> (k + 2)));
                b->rows[i][j].rows[0][0].sign =   b->rows[i][j].rows[0][0].sign  | ((a->rows[j][i].sign  & (leftbit >> (k + 3))));
                b->rows[i][j].rows[0][0].sign =   b->rows[i][j].rows[0][0].sign  | (a->rows[j][i].sign & (leftbit >> (k + 4))) ;
                b->rows[i][j].rows[0][0].sign =   b->rows[i][j].rows[0][0].sign  | ((a->rows[j][i].sign & (leftbit >> (k + 5)))) ;
                b->rows[i][j].rows[0][0].sign =   b->rows[i][j].rows[0][0].sign  | ((a->rows[j][i].sign & (leftbit >> (k + 6)))) ;
                b->rows[i][j].rows[0][0].sign =   b->rows[i][j].rows[0][0].sign  | ((a->rows[j][i].sign & (leftbit >> (k + 7)))) ;
                
                
            
            
                
            }
            
            
            
            
        }
        
    }
    
    
    
    
    
    return  *a;
    
    
}





                