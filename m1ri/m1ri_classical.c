
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
   
      vbg table[64];
    m3d_t * tb = m1ri_malloc(sizeof(m3d_t));
    m3d_transpose(b, tb);
    *c = m3d_create(c, a->nrows, b->ncols);
    if (a->ncols == b->nrows)
       // temp;
    {
        int j, k, u;
        u = 0;
        
        
        for (int i = 0 ;  i  < c->nrows; i ++) {
            
            k = i * c->width;
            for (j = 0; j < c->width ; j++)
            {
                for(u = 0; u < 64;u++ )
                {
                    vbg_mul(&table[0] , &a->rows[i][j], &tb->rows[i][j]); vbg_mul(&table[1] , &a->rows[i][j+1], &tb->rows[i][j+1]);
                    vbg_mul(&table[2] , &a->rows[i][j+2], &tb->rows[i][j+2]); vbg_mul(&table[3] , &a->rows[i][j+3], &tb->rows[i][j+3]);
                    vbg_mul(&table[4] , &a->rows[i][j+4], &tb->rows[i][j+4]); vbg_mul(&table[5] , &a->rows[i][j+5], &tb->rows[i][j+5]);
                    vbg_mul(&table[6] , &a->rows[i][j+6], &tb->rows[i][j+6]); vbg_mul(&table[7] , &a->rows[i][j+7], &tb->rows[i][j+7]);
                    vbg_mul(&table[8] , &a->rows[i][j+8], &tb->rows[i][j+8]); vbg_mul(&table[9] , &a->rows[i][j+9], &tb->rows[i][j+9]);
                    vbg_mul(&table[10] , &a->rows[i][j+10], &tb->rows[i][j+10]); vbg_mul(&table[11] , &a->rows[i][j+11], &tb->rows[i][j+11]);
                    vbg_mul(&table[12] , &a->rows[i][j+12], &tb->rows[i][j+12]); vbg_mul(&table[13] , &a->rows[i][j+13], &tb->rows[i][j+13]);
                    vbg_mul(&table[14] , &a->rows[i][j+14], &tb->rows[i][j+14]); vbg_mul(&table[15] , &a->rows[i][j+15], &tb->rows[i][j+15]);
                    vbg_mul(&table[16] , &a->rows[i][j+16], &tb->rows[i][j+16]); vbg_mul(&table[17] , &a->rows[i][j+17], &tb->rows[i][j+17]);
                    vbg_mul(&table[18] , &a->rows[i][j+18], &tb->rows[i][j+18]); vbg_mul(&table[19] , &a->rows[i][j+19], &tb->rows[i][j+19]);
                    vbg_mul(&table[20] , &a->rows[i][j+20], &tb->rows[i][j+20]); vbg_mul(&table[21] , &a->rows[i][j+21], &tb->rows[i][j+21]);
                    vbg_mul(&table[22] , &a->rows[i][j+22], &tb->rows[i][j+22]); vbg_mul(&table[23] , &a->rows[i][j+23], &tb->rows[i][j+23]);
                    vbg_mul(&table[24] , &a->rows[i][j+24], &tb->rows[i][j+24]); vbg_mul(&table[25] , &a->rows[i][j+25], &tb->rows[i][j+25]);
                    vbg_mul(&table[26] , &a->rows[i][j+26], &tb->rows[i][j+26]); vbg_mul(&table[27] , &a->rows[i][j+27], &tb->rows[i][j+27]);
                    vbg_mul(&table[28] , &a->rows[i][j+28], &tb->rows[i][j+28]); vbg_mul(&table[29] , &a->rows[i][j+29], &tb->rows[i][j+29]);
                    vbg_mul(&table[30] , &a->rows[i][j+30], &tb->rows[i][j+30]); vbg_mul(&table[31] , &a->rows[i][j+31], &tb->rows[i][j+31]);
                    vbg_mul(&table[32] , &a->rows[i][j+32], &tb->rows[i][j+32]); vbg_mul(&table[33] , &a->rows[i][j+33], &tb->rows[i][j+33]);
                    vbg_mul(&table[34] , &a->rows[i][j+34], &tb->rows[i][j+34]); vbg_mul(&table[35] , &a->rows[i][j+35], &tb->rows[i][j+35]);
                    vbg_mul(&table[36] , &a->rows[i][j+36], &tb->rows[i][j+36]); vbg_mul(&table[37] , &a->rows[i][j+37], &tb->rows[i][j+37]);
                    vbg_mul(&table[38] , &a->rows[i][j+38], &tb->rows[i][j+38]); vbg_mul(&table[39] , &a->rows[i][j+39], &tb->rows[i][j+39]);
                    vbg_mul(&table[40] , &a->rows[i][j+40], &tb->rows[i][j+40]); vbg_mul(&table[41] , &a->rows[i][j+41], &tb->rows[i][j+41]);
                    vbg_mul(&table[42] , &a->rows[i][j+42], &tb->rows[i][j+42]); vbg_mul(&table[43] , &a->rows[i][j+43], &tb->rows[i][j+43]);
                    vbg_mul(&table[44] , &a->rows[i][j+44], &tb->rows[i][j+44]); vbg_mul(&table[45] , &a->rows[i][j+45], &tb->rows[i][j+45]);
                    vbg_mul(&table[46] , &a->rows[i][j+46], &tb->rows[i][j+46]); vbg_mul(&table[47] , &a->rows[i][j+47], &tb->rows[i][j+47]);
                    vbg_mul(&table[48] , &a->rows[i][j+48], &tb->rows[i][j+48]); vbg_mul(&table[49] , &a->rows[i][j+49], &tb->rows[i][j+49]);
                    vbg_mul(&table[50] , &a->rows[i][j+50], &tb->rows[i][j+50]); vbg_mul(&table[51] , &a->rows[i][j+51], &tb->rows[i][j+51]);
                    vbg_mul(&table[52] , &a->rows[i][j+52], &tb->rows[i][j+52]); vbg_mul(&table[53] , &a->rows[i][j+53], &tb->rows[i][j+53]);
                    vbg_mul(&table[54] , &a->rows[i][j+54], &tb->rows[i][j+54]); vbg_mul(&table[55] , &a->rows[i][j+55], &tb->rows[i][j+55]);
                    vbg_mul(&table[56] , &a->rows[i][j+56], &tb->rows[i][j+56]); vbg_mul(&table[57] , &a->rows[i][j+57], &tb->rows[i][j+57]);
                    vbg_mul(&table[58] , &a->rows[i][j+58], &tb->rows[i][j+58]); vbg_mul(&table[59] , &a->rows[i][j+59], &tb->rows[i][j+59]);
                    vbg_mul(&table[60] , &a->rows[i][j+60], &tb->rows[i][j+60]); vbg_mul(&table[61] , &a->rows[i][j+61], &tb->rows[i][j+61]);
                    vbg_mul(&table[62] , &a->rows[i][j+62], &tb->rows[i][j+62]); vbg_mul(&table[63] , &a->rows[i][j+63], &tb->rows[i][j+63]);
                    
                    
                    
                    
                //table[64]add_m3dr(&a);
                
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

