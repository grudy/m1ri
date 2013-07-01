/*
//  m1ri_classical.h
//  m1riproject
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

//  Copyright (c) 2013 William Alumbaugh. All rights reserved.
*/

#ifndef M1RIPROJECT_CLASSICAL_H
#define M1RIPROJECT_CLASSICAL_H
#include "m1ri_3dt.h"
#include "m1riwrappers.h"
#include "m1rielarith.h"




/*
 m3d_t m3d_mul(m3d_t *a,m3d_t *b )
 
 {
 m3d_t c;
 
     m3d_t btranspose = m3d_transpose(b);
     if (a->ncols == b->nrows)
     {
 
            int i, j;
            vec temp;
 
            for (i = 0; i < a->nrows; i++)
            {
 
            for(j = 0; j < a->width; i++)
            {
                
                for(i = 0; i < 63; i++)
               // c[i][j] =
 
 
            }
                
 
            }
 
 
 
 
 
 
 }
 
 
 return c;
 
 }

 


vbg mul_128_inner(vbg a, vbg r1, vbg tables6[9][64], vbg tables5[9][32]){
    vec v1, v2;
    vbg t1, t2, r2;
    v2 = a.sign;
    v1 = a.units ^ v2;
    t1 = tables6[0][v1&63]; iaddgf3(&r1, &t1); v1 >>= 6;
    r2 = tables6[0][v2&63];                    v2 >>= 6;
    t1 = tables6[1][v1&63]; iaddgf3(&r1, &t1); v1 >>= 6;
    t2 = tables6[1][v2&63]; iaddgf3(&r2, &t2); v2 >>= 6;
    t1 = tables6[2][v1&63]; iaddgf3(&r1, &t1); v1 >>= 6;
    t2 = tables6[2][v2&63]; iaddgf3(&r2, &t2); v2 >>= 6;
    t1 = tables6[3][v1&63]; iaddgf3(&r1, &t1); v1 >>= 6;
    t2 = tables6[3][v2&63]; iaddgf3(&r2, &t2); v2 >>= 6;
    t1 = tables6[4][v1&63]; iaddgf3(&r1, &t1); v1 >>= 6;
    t2 = tables6[4][v2&63]; iaddgf3(&r2, &t2); v2 >>= 6;
    t1 = tables6[5][v1&63]; iaddgf3(&r1, &t1); v1 >>= 6;
    t2 = tables6[5][v2&63]; iaddgf3(&r2, &t2); v2 >>= 6;
    t1 = tables6[6][v1&63]; iaddgf3(&r1, &t1); v1 >>= 6;
    t2 = tables6[6][v2&63]; iaddgf3(&r2, &t2); v2 >>= 6;
    t1 = tables6[7][v1&63]; iaddgf3(&r1, &t1); v1 >>= 6;
    t2 = tables6[7][v2&63]; iaddgf3(&r2, &t2); v2 >>= 6;
    t1 = tables6[8][v1&63]; iaddgf3(&r1, &t1); v1 >>= 6;
    t2 = tables6[8][v2&63]; iaddgf3(&r2, &t2); v2 >>= 6;
    t1 = tables5[0][v1&31]; iaddgf3(&r1, &t1); v1 >>= 5;
    t2 = tables5[0][v2&31]; iaddgf3(&r2, &t2); v2 >>= 5;
    t1 = tables5[1][v1&31]; iaddgf3(&r1, &t1);
    t2 = tables5[1][v2&31]; iaddgf3(&r2, &t2);
    isubgf3(&r1, &r2);
    return r1;
}



inline void mul_128(vbg *RR, vbg *AA, vbg *BB)
{
    vec i, j;
   
    vbg *A, *B, *R;
    
    
    vbg tables6[9][64]	;
    vbg tables5[2][32]  ;
    
    for ( i  = 0 ; i < 255 ; i++)
    {
        RR[i].units = RR[i].sign = 0;
    }
    for ( j  = 0 ; i < 3 ; i++)
    {
        
        B = BB + j * 64;
        R = RR + (j/2)*128;
        A = AA + (j%2)*128;
        
        for ( i  = 0 ; i < 8 ; i++)
        {
            combine6(tables6[i], B + 0 + 6*i);
        }
        
        for ( i  = 0 ; i < 1 ; i++)
        {
            combine5(tables5[i], B + 54 + 5*i);
            
        }
        for ( i  = 0 ; i < 127 ; i++)
        {
            R[i] = mul_128_inner(A[i], R[i], tables6, tables5);
        }
    }
    
}

*/









#endif