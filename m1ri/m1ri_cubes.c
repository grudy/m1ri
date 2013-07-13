
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
            
void * m3d_64_cubes_take2(m3_slice * b, m3d_t  *a )
{
    
    
    if((a->nrows%m1ri_word || a->ncols%m1ri_word   )  != 0)
    {
        
        b->width  = a->width;
        b->nrows  = DN(a->width, 64);
        int n = b->width * b->nrows ;
        b->blocks = malloc(a->nrows * a->width * sizeof(vec));
        b->rows  = malloc(n * sizeof(a->nrows * a->width * sizeof(vec)));
        u_int32_t  k, j;
        j = b->nrows%8;
        k = b->width%8;
        /*
         The last column and row in the structure
         */
        signed short rlast, clast;
        rlast = 63 - a->nrows%m1ri_word;
        clast = 63 - a->ncols%m1ri_word;
        
        for(k  = 0; k < (b->nrows - 1) ; k++ )
        {
            
            for(j = 0; j < (a->ncols - (a->ncols%64) ) ; j = j + 64 )
            {
                
            
                b->blocks[0] = a->rows[ m1ri_word + 0][ m1ri_word + 0];
                b->blocks[0] = a->rows[ m1ri_word + 1][ m1ri_word + 1];
                b->blocks[0] = a->rows[ m1ri_word + 2][ m1ri_word + 2];
                b->blocks[0] = a->rows[ m1ri_word + 3][ m1ri_word + 3];
                b->blocks[0] = a->rows[ m1ri_word + 4][ m1ri_word + 4];
                b->blocks[0] = a->rows[ m1ri_word + 5][ m1ri_word + 5];
                b->blocks[0] = a->rows[ m1ri_word + 6][ m1ri_word + 6];
                b->blocks[0] = a->rows[ m1ri_word + 7][ m1ri_word + 7];
                b->blocks[0] = a->rows[ m1ri_word + 8][ m1ri_word + 8];
                b->blocks[0] = a->rows[ m1ri_word + 9][ m1ri_word + 9];
                b->blocks[0] = a->rows[ m1ri_word + 10][ m1ri_word + 10];
                b->blocks[0] = a->rows[ m1ri_word + 11][ m1ri_word + 11];
                b->blocks[0] = a->rows[ m1ri_word + 12][ m1ri_word + 12];
                b->blocks[0] = a->rows[ m1ri_word + 13][ m1ri_word + 13];
                b->blocks[0] = a->rows[ m1ri_word + 14][ m1ri_word + 14];
                b->blocks[0] = a->rows[ m1ri_word + 15][ m1ri_word + 15];
                b->blocks[0] = a->rows[ m1ri_word + 16][ m1ri_word + 16];
                b->blocks[0] = a->rows[ m1ri_word + 17][ m1ri_word + 17];
                b->blocks[0] = a->rows[ m1ri_word + 18][ m1ri_word + 18];
                b->blocks[0] = a->rows[ m1ri_word + 19][ m1ri_word + 19];
                b->blocks[0] = a->rows[ m1ri_word + 20][ m1ri_word + 20];
                b->blocks[0] = a->rows[ m1ri_word + 21][ m1ri_word + 21];
                b->blocks[0] = a->rows[ m1ri_word + 22][ m1ri_word + 22];
                b->blocks[0] = a->rows[ m1ri_word + 23][ m1ri_word + 23];
                b->blocks[0] = a->rows[ m1ri_word + 24][ m1ri_word + 24];
                b->blocks[0] = a->rows[ m1ri_word + 25][ m1ri_word + 25];
                b->blocks[0] = a->rows[ m1ri_word + 26][ m1ri_word + 26];
                b->blocks[0] = a->rows[ m1ri_word + 27][ m1ri_word + 27];
                b->blocks[0] = a->rows[ m1ri_word + 28][ m1ri_word + 28];
                b->blocks[0] = a->rows[ m1ri_word + 29][ m1ri_word + 29];
                b->blocks[0] = a->rows[ m1ri_word + 30][ m1ri_word + 30];
                b->blocks[0] = a->rows[ m1ri_word + 31][ m1ri_word + 31];
                b->blocks[0] = a->rows[ m1ri_word + 32][ m1ri_word + 32];
                b->blocks[0] = a->rows[ m1ri_word + 33][ m1ri_word + 33];
                b->blocks[0] = a->rows[ m1ri_word + 34][ m1ri_word + 34];
                b->blocks[0] = a->rows[ m1ri_word + 35][ m1ri_word + 35];
                b->blocks[0] = a->rows[ m1ri_word + 36][ m1ri_word + 36];
                b->blocks[0] = a->rows[ m1ri_word + 37][ m1ri_word + 37];
                b->blocks[0] = a->rows[ m1ri_word + 38][ m1ri_word + 38];
                b->blocks[0] = a->rows[ m1ri_word + 39][ m1ri_word + 39];
                b->blocks[0] = a->rows[ m1ri_word + 40][ m1ri_word + 40];
                b->blocks[0] = a->rows[ m1ri_word + 41][ m1ri_word + 41];
                b->blocks[0] = a->rows[ m1ri_word + 42][ m1ri_word + 42];
                b->blocks[0] = a->rows[ m1ri_word + 43][ m1ri_word + 43];
                b->blocks[0] = a->rows[ m1ri_word + 44][ m1ri_word + 44];
                b->blocks[0] = a->rows[ m1ri_word + 45][ m1ri_word + 45];
                b->blocks[0] = a->rows[ m1ri_word + 46][ m1ri_word + 46];
                b->blocks[0] = a->rows[ m1ri_word + 47][ m1ri_word + 47];
                b->blocks[0] = a->rows[ m1ri_word + 48][ m1ri_word + 48];
                b->blocks[0] = a->rows[ m1ri_word + 49][ m1ri_word + 49];
                b->blocks[0] = a->rows[ m1ri_word + 50][ m1ri_word + 50];
                b->blocks[0] = a->rows[ m1ri_word + 51][ m1ri_word + 51];
                b->blocks[0] = a->rows[ m1ri_word + 52][ m1ri_word + 52]; 
                b->blocks[0] = a->rows[ m1ri_word + 53][ m1ri_word + 53]; 
                b->blocks[0] = a->rows[ m1ri_word + 54][ m1ri_word + 54]; 
                b->blocks[0] = a->rows[ m1ri_word + 55][ m1ri_word + 55]; 
                b->blocks[0] = a->rows[ m1ri_word + 56][ m1ri_word + 56]; 
                b->blocks[0] = a->rows[ m1ri_word + 57][ m1ri_word + 57]; 
                b->blocks[0] = a->rows[ m1ri_word + 58][ m1ri_word + 58]; 
                b->blocks[0] = a->rows[ m1ri_word + 59][ m1ri_word + 59]; 
                b->blocks[0] = a->rows[ m1ri_word + 60][ m1ri_word + 60]; 
                b->blocks[0] = a->rows[ m1ri_word + 61][ m1ri_word + 61]; 
                b->blocks[0] = a->rows[ m1ri_word + 62][ m1ri_word + 62]; 
                b->blocks[0] = a->rows[ m1ri_word + 63][ m1ri_word + 63];
                
               // b->blocks[(k*b->width) + j ] = m3d_window(a, (k * m1ri_word), (j * m1ri_word), ((k * m1ri_word) + 63), ((j * m1ri_word) + 63));
                
                
            }
            
             for(j = 0; j < (a->ncols - (a->ncols%64) ) ; j++  )
             {
             
             
             
             }
           
        }
        
        while(j < b->width )
        {
            
           // b->blocks[(k*b->width) + j ] = m3d_window(a, (k * m1ri_word), (j * m1ri_word), ((k * m1ri_word) + rlast), ((j * m1ri_word) + 63));
            j++;
        }
       // b->blocks[(k*b->width) + j ] = m3d_window(a, (k * m1ri_word), (j * m1ri_word), ((k * m1ri_word) + rlast), ((j * m1ri_word) + clast));
       // b->rows[j] = b->blocks + (k * b->width);
        
        
        
        
        
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
                
             //   b->blocks[(k*b->width) + j ] = m3d_window(a, (k * m1ri_word), (j * m1ri_word), ((k * m1ri_word) + 63), ((j * m1ri_word) + 63));
                
            }
          //  b->rows[j] = b->blocks + (k * b->width);
            
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





                