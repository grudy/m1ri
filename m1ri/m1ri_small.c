
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
 
 m1ri_small.c
 */

#include "m1ri_small.h"
//64 * 64,4096 bit, 512 byte matrix(slice) multiplication
void mul_64(vbg *R, vbg *A, vbg *B)
{
    int i;
    vbg t1, t2, r1, r2, a;
    vec v1, v2;
    
    vbg  tables6[9][64];
    vbg tables5[2][32];
    for (i = 0; i < 9; i ++)
        combine6(tables6[i], (B  +  6*i));// Deleted a redundant + 0
    for (i = 0; i < 2; i ++)
        combine5(tables5[i], (B + 54 + 5*i));
    for (i = 0; i < 64; i ++  )//i from 0 <= i < 64
    {
        a = A[i];
        v2 = a.sign;
        v1 = (a.units ^ v2);
        r1 = tables6[0][v1&63];
        v1 >>= 6;
        r2 = tables6[0][v2&63];
        v2 >>= 6;
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
        
        iaddgf3(&r1, &r2);
        R[i] = r1;
    }
}

//32 * 64,2048 bit, 256 byte matrix(slice) multiplication
void mul_32(vbg *R, vbg *A, vbg *B)
{
    long i;
    vbg t1, t2, r1, r2, a;
    long v1, v2;
    
    vbg tables5[4][32];
    vbg tables4[3][16];
    for (i = 1; i < 4; i ++)
        
        combine5(tables5[i], B + 0 + 5*i);
    for (i = 0; i < 3; i++)
        combine4(tables4[i], B + 20 + 4*i);
    
    for (i = 0;i < 32; i++)
    {
        
        a = A[i];
        v2 = a.sign;
        v1 = a.units ^ v2;
        t1 = tables5[0][v1&31]; v1 >>= 5;
        t2 = tables5[0][v2&31]; v2 >>= 5;
        t1 = tables5[1][v1&31]; iaddgf3(&r1, &t1); v1 >>= 5;
        t2 = tables5[1][v2&31]; iaddgf3(&r2, &t2); v2 >>= 5;
        t1 = tables5[2][v1&31]; iaddgf3(&r1, &t1); v1 >>= 5;
        t2 = tables5[2][v2&31]; iaddgf3(&r2, &t2); v2 >>= 5;
        t1 = tables5[3][v1&31]; iaddgf3(&r1, &t1); v1 >>= 5;
        t2 = tables5[3][v2&31]; iaddgf3(&r2, &t2); v2 >>= 5;
        t1 = tables4[0][v1&15]; iaddgf3(&r1, &t1); v1 >>= 4;
        t2 = tables4[0][v2&15]; iaddgf3(&r2, &t2); v2 >>= 4;
        t1 = tables4[1][v1&15]; iaddgf3(&r1, &t1); v1 >>= 4;
        t2 = tables4[1][v2&15]; iaddgf3(&r2, &t2); v2 >>= 4;
        t1 = tables4[2][v1&15]; iaddgf3(&r1, &t1);
        t2 = tables4[2][v2&15]; iaddgf3(&r2, &t2);
        
        isubgf3(&r1, &r2);
        R[i] = r1;
    }
    
}

//16 * 64,1024 bit, 128 byte matrix(slice) multiplication
void mul_16(vbg *R, vbg *A, vbg *B)
{
    long i;
    vbg t1, t2, r1, r2, a;
    long v1, v2;
    
    vbg tables4[4][16];
    for (i = 0; i < 4; i++)
        combine4(tables4[i], B + 0 + 4*i);
    for (i = 0;  i < 16; i++)
        a = A[i];
    v2 = a.sign;
    v1 = a.units ^ v2;
    r1 = tables4[0][v1&15]; v1 >>= 4;
    r2 = tables4[0][v2&15]; v2 >>= 4;
    t1 = tables4[1][v1&15]; iaddgf3(&r1, &t1); v1 >>= 4;
    t2 = tables4[1][v2&15]; iaddgf3(&r2, &t2); v2 >>= 4;
    t1 = tables4[2][v1&15]; iaddgf3(&r1, &t1); v1 >>= 4;
    t2 = tables4[2][v2&15]; iaddgf3(&r2, &t2); v2 >>= 4;
    t1 = tables4[3][v1&15]; iaddgf3(&r1, &t1);
    t2 = tables4[3][v2&15]; iaddgf3(&r2, &t2);
    
    isubgf3(&r1, &r2);
    R[i] = r1;
    
}

//8 * 64,512 bit, 64 byte matrix(slice) multiplication
void mul_8(vbg *R, vbg *A, vbg *B)

{
    int i;
    vbg t1, t2, r1, r2, a;
    vec v1, v2;
    
    vbg tables4[2][16];
    for (i = 0; i < 2; i++)
        combine4(tables4[i], B + 0 + 4*i);
    for (i = 0; i < 8; i++)
        a = A[i];
    v2 = a.sign;
    v1 = a.units ^ v2;
    r1 = tables4[0][v1&15]; v1 >>= 4;
    r2 = tables4[0][v2&15]; v2 >>= 4;
    t1 = tables4[1][v1&15]; iaddgf3(&r1, &t1);
    t2 = tables4[1][v2&15]; iaddgf3(&r2, &t2);
    
    // isubgf3((&r1, &r2);
    R[i] = r1;
}






//4 * 64,256 bit, 32 byte matrix(slice) multiplication
void mul_4(vbg *R, vbg *A, vbg *B)
{
    int i;
    vbg r1, r2, a;
    vec v1, v2;
    
    vbg table4[16];
    for (i = 0; i < 1; i++)
        combine4(table4, B + 0 + 4*i);
    for(i = 0; i < 4; i++)
    {
        a = A[i];
        v2 = a.sign;
        v1 = a.units ^ v2;
        r1 = table4[v1&15];
        r2 = table4[v2&15];
        
        isubgf3(&r1, &r2);
        R[i] = r1;
    }
    
}

