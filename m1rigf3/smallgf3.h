// smallfg3.h
//// computions for small matrices over gf3
// TOMAS J. BOOTHBY AND ROBERT W. BRADSHAW "BITSLICING AND THE METHOD OF FOUR
// RUSSIANS OVER LARGER FINITE FIELDS"
//
/*Copyright 2013 William Andrew Alumbaugh <williamandrewalumbaugh@gmail.com>
 
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
 

 
 */
#ifndef m1rigf3_smallgf3_h
#define m1rigf3_smallgf3_h
#include"m1rigf3.h"


void mul_64(matrixgf3 *R, matrixgf3 *A, matrixgf3 *B)
{
    int i;
matrixgf3 t1, t2, r1, r2, a
long v1, v2

matrixgf3 tables6[9][64]
matrixgf3 tables5[2][32]
for i from 0 <= i < 9:
combine6(tables6[i], B + 0 + 6*i)
for i from 0 <= i < 2:
combine5(tables5[i], B + 54 + 5*i)
    for (i = 0; i <  )//i from 0 <= i < 64
        a = A[i]
        v2 = a.s
        v1 = a.u ^ v2
        r1 = tables6[0][v1&63];                    v1 >>= 6;
        r2 = tables6[0][v2&63];                    v2 >>= 6;
        t1 = tables6[1][v1&63]; v3_iadd(&r1, &t1); v1 >>= 6;
        t2 = tables6[1][v2&63]; v3_iadd(&r2, &t2); v2 >>= 6;
        t1 = tables6[2][v1&63]; v3_iadd(&r1, &t1); v1 >>= 6;
        t2 = tables6[2][v2&63]; v3_iadd(&r2, &t2); v2 >>= 6;
        t1 = tables6[3][v1&63]; v3_iadd(&r1, &t1); v1 >>= 6;
        t2 = tables6[3][v2&63]; v3_iadd(&r2, &t2); v2 >>= 6;
        t1 = tables6[4][v1&63]; v3_iadd(&r1, &t1); v1 >>= 6;
        t2 = tables6[4][v2&63]; v3_iadd(&r2, &t2); v2 >>= 6;
        t1 = tables6[5][v1&63]; v3_iadd(&r1, &t1); v1 >>= 6;
        t2 = tables6[5][v2&63]; v3_iadd(&r2, &t2); v2 >>= 6;
        t1 = tables6[6][v1&63]; v3_iadd(&r1, &t1); v1 >>= 6;
        t2 = tables6[6][v2&63]; v3_iadd(&r2, &t2); v2 >>= 6;
        t1 = tables6[7][v1&63]; v3_iadd(&r1, &t1); v1 >>= 6;
        t2 = tables6[7][v2&63]; v3_iadd(&r2, &t2); v2 >>= 6;
        t1 = tables6[8][v1&63]; v3_iadd(&r1, &t1); v1 >>= 6;
        t2 = tables6[8][v2&63]; v3_iadd(&r2, &t2); v2 >>= 6;
        t1 = tables5[0][v1&31]; v3_iadd(&r1, &t1); v1 >>= 5;
        t2 = tables5[0][v2&31]; v3_iadd(&r2, &t2); v2 >>= 5;
        t1 = tables5[1][v1&31]; v3_iadd(&r1, &t1);
        t2 = tables5[1][v2&31]; v3_iadd(&r2, &t2);

v3_isub(&r1, &r2)
R[i] = r1
}
inline  void mul_32(matrixgf3 *R, matrixgf3 *A, matrixgf3 *B)
{
    long i
    matrixgf3 t1, t2, r1, r2, a;
    long v1, v2;

    matrixgf3 tables5[4][32];
    matrixgf3 tables4[3][16];
for i from 0 <= i < 4:

combine5(tables5[i], B + 0 + 5*i)
for i from 0 <= i < 3:
combine4(tables4[i], B + 20 + 4*i)

for i from 0 <= i < 32
{
    
    a = A[i];
    v2 = a.s;
    v1 = a.u ^ v2;
    t1 = tables5[0][v1&31];                    v1 >>= 5;
    t2 = tables5[0][v2&31];                    v2 >>= 5;
    t1 = tables5[1][v1&31]; v3_iadd(&r1, &t1); v1 >>= 5;
    t2 = tables5[1][v2&31]; v3_iadd(&r2, &t2); v2 >>= 5;
    t1 = tables5[2][v1&31]; v3_iadd(&r1, &t1); v1 >>= 5;
    t2 = tables5[2][v2&31]; v3_iadd(&r2, &t2); v2 >>= 5;
    t1 = tables5[3][v1&31]; v3_iadd(&r1, &t1); v1 >>= 5;
    t2 = tables5[3][v2&31]; v3_iadd(&r2, &t2); v2 >>= 5;
    t1 = tables4[0][v1&15]; v3_iadd(&r1, &t1); v1 >>= 4;
    t2 = tables4[0][v2&15]; v3_iadd(&r2, &t2); v2 >>= 4;
    t1 = tables4[1][v1&15]; v3_iadd(&r1, &t1); v1 >>= 4;
    t2 = tables4[1][v2&15]; v3_iadd(&r2, &t2); v2 >>= 4;
    t1 = tables4[2][v1&15]; v3_iadd(&r1, &t1);
    t2 = tables4[2][v2&15]; v3_iadd(&r2, &t2);
    
    v3_isub(&r1, &r2);
    R[i] = r1;
}
inline  void mul_16(matrixgf3 *R, matrixgf3 *A, matrixgf3 *B)

{
    long i;
    matrixgf3 t1, t2, r1, r2, a
    long v1, v2
    
    matrixgf3 tables4[4][16]
    for i from 0 <= i < 4:
        combine4(tables4[i], B + 0 + 4*i)
        for i from 0 <= i < 16:
            a = A[i];
    v2 = a.s
    v1 = a.u ^ v2
    r1 = tables4[0][v1&15];                    v1 >>= 4;
    r2 = tables4[0][v2&15];                    v2 >>= 4;
    t1 = tables4[1][v1&15]; v3_iadd(&r1, &t1); v1 >>= 4;
    t2 = tables4[1][v2&15]; v3_iadd(&r2, &t2); v2 >>= 4;
    t1 = tables4[2][v1&15]; v3_iadd(&r1, &t1); v1 >>= 4;
    t2 = tables4[2][v2&15]; v3_iadd(&r2, &t2); v2 >>= 4;
    t1 = tables4[3][v1&15]; v3_iadd(&r1, &t1);
    t2 = tables4[3][v2&15]; v3_iadd(&r2, &t2);
    
    v3_isub(&r1, &r2)
    R[i] = r1
    
}

void mul_8(matrixgf3 *R, matrixgf3 *A, matrixgf3 *B)

{
    long i;
    matrixgf3 t1, t2, r1, r2, a;
    long v1, v2;
    
    matrixgf3 tables4[2][16]
    for i from 0 <= i < 2:
        combine4(tables4[i], B + 0 + 4*i)
        for i from 0 <= i < 8:
            a = A[i];
    v2 = a.s;
    v1 = a.u ^ v2;
    r1 = tables4[0][v1&15];                    v1 >>= 4;
    r2 = tables4[0][v2&15];                    v2 >>= 4;
    t1 = tables4[1][v1&15]; v3_iadd(&r1, &t1);
    t2 = tables4[1][v2&15]; v3_iadd(&r2, &t2);
    
    isubgf3((&r1, &r2);
            R[i] = r1;
            }
            
            void mul_4(matrixgf3 *R, matrixgf3 *A, matrixgf3 *B)
    {
        long i;
        matrixgf3 r1, r2, a;
        long v1, v2;
        
        matrixgf3 table4[16];
        for i from 0 <= i < 1:
            combine4(table4, B + 0 + 4*i)
            for i from 0 <= i < 4:
                a = A[i];
        v2 = a.s;
        v1 = a.u ^ v2;
        r1 = table4[v1&15];
        r2 = table4[v2&15];
        
        isubgf3(&r1, &r2);
        R[i] = r1s;
        
        
    }
            
            
#endif
