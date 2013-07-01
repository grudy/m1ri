//
// strassen.h
//
//// computions for small matrices over gf3
// TOMAS J. BOOTHBY AND ROBERT W. BRADSHAW "BITSLICING AND THE METHOD OF FOUR
// RUSSIANS OVER LARGER FINITE FIELDS"
//
/*Copyright 2013 William Andrew Alumbaugh <williamandrewalumbaugh@gmail.com>
<<<<<<< HEAD
<<<<<<< HEAD
=======
 
>>>>>>> 
=======


>>>>>>> ff5b2b342a1c60afd5b26258e5d1070644cc496c
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
<<<<<<< HEAD
 */



<<<<<<< HEAD
=======
*/
=======
>>>>>>> ff5b2b342a1c60afd5b26258e5d1070644cc496c

#ifndef m1rigf3_strassen_h
#define m1rigf3_strassen_h
#include "m1rigf3.h"
#include "m1rielarith.h"

void mul_128(vbg *C, vbg *A, vbg *B)
   
    {
 
    int i;
    
    vbg *A11, *C12, *A21, *A22;
    vbg *B11, *B12, *B21, *B22;
    vbg *C11, *A12, *C21, *C22;
 
    
    vbg X1[64],X2[64];
    
    for (i = 0; i < 64; i = i +1)
        X1[i].units = X1[i].sign = X2[i].units = X2[i].sign = 0;
 
    for (i = 0; i < 256; i = i +1)
    {
        C[i].units = C[i].sign = 0;
    }
 
    A11 = A;
    B11 = B;
    C11 = C;
    A21 = A+64;
    B21 = B+64;
    C21 = C+64;
    A12 = A+128;
    B12 = B+128;
    C12 = C+128;
    A22 = A+192;
    B22 = B+192;
    C22 = C+192;
 
    sub_64gf3(X1,A11,A21);
    sub_64gf3(X2,B22,B12);
    mul_64(C21,X1,X2);
    add_64gf3(X1,A21,A22);
    sub_64gf3(X2,B12,B11);
    mul_64(C22,X1,X2);
    sub_64gf3(X1,X1,A11);
    sub_64gf3(X2,B22,X2);
    mul_64(C12,X1,X2);
    sub_64gf3(X1,A12,X1);
    mul_64(C11,X1,B22);
    mul_64(X1,A11,B11);
    add_64gf3(C12,X1,C12);
    add_64gf3(C21,C12,C21);
    add_64gf3(C12,C12,C22);
    add_64gf3(C22,C21,C22);
    add_64gf3(C12,C12,C11);
    sub_64gf3(X2,X2,B21);
    mul_64(C11,A22,X2);
    sub_64gf3(C21,C21,C11);
    mul_64(C11,A12,B21);
    add_64gf3(C11,X1,C11);
 
    }
 


#endif
