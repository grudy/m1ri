//
//  strassen.h
//  m1rigf3
//
//  Created by grudy on 5/25/13.
//  Copyright (c) 2013 William Alumbaugh. All rights reserved.
//

#ifndef m1rigf3_strassen_h
#define m1rigf3_strassen_h


 mul_128(vec3 *C, vec3 *A, vec3 *B)
   
    {
 
    int i;
    
    matrixgf3 *A11, *C12, *A21, *A22;
    matrixgf3 *B11, *B12, *B21, *B22;
    matrixgf3 *C11, *A12, *C21, *C22;
 
    
    matrixgf3 X1[64],X2[64];
    
    for (i = 0; i < 64; i = i +1)
        X1[i].u = X1[i].s = X2[i].u = X2[i].s = 0;
 
    for (i = 0; i < 256; i = i +1)
    {
        C[i].u = C[i].s = 0;
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
 
    sub_64(X1,A11,A21);
    sub_64(X2,B22,B12);
    mul_64(C21,X1,X2);
    add_64(X1,A21,A22);
    sub_64(X2,B12,B11);
    mul_64(C22,X1,X2);
    sub_64(X1,X1,A11);
    sub_64(X2,B22,X2);
    mul_64(C12,X1,X2);
    sub_64(X1,A12,X1);
    mul_64(C11,X1,B22);
    mul_64(X1,A11,B11);
    add_64(C12,X1,C12);
    add_64(C21,C12,C21);
    add_64(C12,C12,C22);
    add_64(C22,C21,C22);
    add_64(C12,C12,C11);
    sub_64(X2,X2,B21);
    mul_64(C11,A22,X2);
    sub_64(C21,C21,C11);
    mul_64(C11,A12,B21);
    add_64(C11,X1,C11);
 
    }
 


#endif
