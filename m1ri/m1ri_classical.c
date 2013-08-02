
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
 
 m1ri_classical.c
 */

#include <stdio.h>
#include "m1ri_classical.h"




m3d_t m3d_mul_naive(m3d_t *c, m3d_t *a, m3d_t *b)
{
    
    //vbg table[64];
    
   // m3d_t * tb = m1ri_malloc(sizeof(m3d_t));
   
    *c = m3d_create(c, a->nrows, b->ncols);
    if (a->ncols == b->nrows)
        // temp;
    {
        int j, k, u, i,  remainder;
        u = 0;
        remainder = c->nrows%8;
        
        for ( i = 0 ;  i  < c->nrows - remainder ; i++) {
            
            k = i * c->width;
            for (j = 0; j < c->width ; j++)
            {
                for(u = 0; u < 64;u ++)
                {
                    
                    
                    
                }
         
            }
            
        }
        
    }
    
    return *c;
    
}

m3d_t m3d_mul_transposed(m3d_t *c, m3d_t *a, m3d_t *b)
{
   
      vbg table[64];
     
    m3d_t * tb = m1ri_malloc(sizeof(m3d_t));
   // m3d_transpose(b, tb);
    *c = m3d_create(c, a->nrows, b->ncols);
    if (a->ncols == b->nrows)
       // temp;
    {
        int j, k, u, i, remainder;
        u = 0;
        remainder = c->nrows%8;
        
        for ( i = 0 ;  i  < c->nrows - remainder ; i = i + 8) {
            
            k = i * c->width;
            for (j = 0; j < c->width ; j++)
            {
                for(u = 0; u < 64;u = + 8 )
                {
                    
                    vbg_mul(&table[0] , &a->rows[i][j+0], &tb->rows[i][j+0]); vbg_mul(&table[1] , &a->rows[i][j+1], &tb->rows[i][j+1]);
                    vbg_mul(&table[2] , &a->rows[i][j+2], &tb->rows[i][j+2]); vbg_mul(&table[3] , &a->rows[i][j+3], &tb->rows[i][j+3]);
                    vbg_mul(&table[4] , &a->rows[i][j+4], &tb->rows[i][j+4]); vbg_mul(&table[5] , &a->rows[i][j+5], &tb->rows[i][j+5]);
                    vbg_mul(&table[6] , &a->rows[i][j+6], &tb->rows[i][j+6]); vbg_mul(&table[7] , &a->rows[i][j+7], &tb->rows[i][j+7]);
                    
                    vbg_mul(&table[8] , &a->rows[i][j+0], &tb->rows[i + 1][j+8]); vbg_mul(&table[9] , &a->rows[i][j+1], &tb->rows[i+1][j+1]);
                    vbg_mul(&table[10] , &a->rows[i][j+2], &tb->rows[i + 1][j+10]); vbg_mul(&table[11] , &a->rows[i][j+3], &tb->rows[i+1][j+3]);
                    vbg_mul(&table[12] , &a->rows[i][j+4], &tb->rows[i + 1][j+12]); vbg_mul(&table[13] , &a->rows[i][j+5], &tb->rows[i+1][j+5]);
                    vbg_mul(&table[14] , &a->rows[i][j+6], &tb->rows[i + 1][j+14]); vbg_mul(&table[15] , &a->rows[i][j+7], &tb->rows[i+1][j+7]);
                    
                    vbg_mul(&table[16] , &a->rows[i][j+0], &tb->rows[i + 2][j+16]); vbg_mul(&table[17] , &a->rows[i][j+1], &tb->rows[i+2][j+1]);
                    vbg_mul(&table[18] , &a->rows[i][j+2], &tb->rows[i + 2][j+18]); vbg_mul(&table[19] , &a->rows[i][j+3], &tb->rows[i+2][j+3]);
                    vbg_mul(&table[20] , &a->rows[i][j+4], &tb->rows[i + 2][j+20]); vbg_mul(&table[21] , &a->rows[i][j+5], &tb->rows[i+2][j+5]);
                    vbg_mul(&table[22] , &a->rows[i][j+6], &tb->rows[i + 2][j+22]); vbg_mul(&table[23] , &a->rows[i][j+7], &tb->rows[i+2][j+7]);
                    
                    vbg_mul(&table[24] , &a->rows[i][j+0], &tb->rows[i + 3][j+24]); vbg_mul(&table[25] , &a->rows[i][j+1], &tb->rows[i+3][j+1]);
                    vbg_mul(&table[26] , &a->rows[i][j+2], &tb->rows[i + 3][j+26]); vbg_mul(&table[27] , &a->rows[i][j+3], &tb->rows[i+3][j+3]);
                    vbg_mul(&table[28] , &a->rows[i][j+4], &tb->rows[i + 3][j+28]); vbg_mul(&table[29] , &a->rows[i][j+5], &tb->rows[i+3][j+5]);
                    vbg_mul(&table[30] , &a->rows[i][j+6], &tb->rows[i + 3][j+30]); vbg_mul(&table[31] , &a->rows[i][j+7], &tb->rows[i+3][j+7]);
                    
                    vbg_mul(&table[32] , &a->rows[i][j+0], &tb->rows[i + 4][j+32]); vbg_mul(&table[33] , &a->rows[i][j+1], &tb->rows[i+4][j+1]);
                    vbg_mul(&table[34] , &a->rows[i][j+2], &tb->rows[i + 4][j+34]); vbg_mul(&table[35] , &a->rows[i][j+3], &tb->rows[i+4][j+3]);
                    vbg_mul(&table[36] , &a->rows[i][j+4], &tb->rows[i + 4][j+36]); vbg_mul(&table[37] , &a->rows[i][j+5], &tb->rows[i+4][j+5]);
                    vbg_mul(&table[38] , &a->rows[i][j+6], &tb->rows[i + 4][j+38]); vbg_mul(&table[39] , &a->rows[i][j+7], &tb->rows[i+4][j+7]);
                    
                    vbg_mul(&table[40] , &a->rows[i][j+0], &tb->rows[i + 5][j+40]); vbg_mul(&table[41] , &a->rows[i][j+1], &tb->rows[i+5][j+1]);
                    vbg_mul(&table[42] , &a->rows[i][j+2], &tb->rows[i + 5][j+42]); vbg_mul(&table[43] , &a->rows[i][j+3], &tb->rows[i+5][j+3]);
                    vbg_mul(&table[44] , &a->rows[i][j+4], &tb->rows[i + 5][j+44]); vbg_mul(&table[45] , &a->rows[i][j+5], &tb->rows[i+5][j+5]);
                    vbg_mul(&table[46] , &a->rows[i][j+6], &tb->rows[i + 5][j+46]); vbg_mul(&table[47] , &a->rows[i][j+7], &tb->rows[i+5][j+7]);
                    
                    vbg_mul(&table[48] , &a->rows[i][j+0], &tb->rows[i + 6][j+48]); vbg_mul(&table[49] , &a->rows[i][j+1], &tb->rows[i+6][j+1]);
                    vbg_mul(&table[50] , &a->rows[i][j+2], &tb->rows[i + 6][j+50]); vbg_mul(&table[51] , &a->rows[i][j+3], &tb->rows[i+6][j+3]);
                    vbg_mul(&table[52] , &a->rows[i][j+4], &tb->rows[i + 6][j+52]); vbg_mul(&table[53] , &a->rows[i][j+5], &tb->rows[i+6][j+5]);
                    vbg_mul(&table[54] , &a->rows[i][j+6], &tb->rows[i + 6][j+54]); vbg_mul(&table[55] , &a->rows[i][j+7], &tb->rows[i+6][j+7]);
                    
                    vbg_mul(&table[56] , &a->rows[i][j+0], &tb->rows[i + 7][j+56]); vbg_mul(&table[57] , &a->rows[i][j+1], &tb->rows[i+7][j+1]);
                    vbg_mul(&table[58] , &a->rows[i][j+2], &tb->rows[i + 7][j+58]); vbg_mul(&table[59] , &a->rows[i][j+3], &tb->rows[i+7][j+3]);
                    vbg_mul(&table[60] , &a->rows[i][j+4], &tb->rows[i + 7][j+60]); vbg_mul(&table[61] , &a->rows[i][j+5], &tb->rows[i+7][j+5]);
                    vbg_mul(&table[62] , &a->rows[i][j+6], &tb->rows[i + 7][j+62]); vbg_mul(&table[63] , &a->rows[i][j+7], &tb->rows[i+7][j+7]);


              
                }
             
            }
            
        }
        
        
    }
    
    
    return  *c;
}




/*
 Tom Boothys method

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

