//
//  main.c
//  proofofconcept of  the addition logic used   in  the paper  published by
//  TOMAS J. BOOTHBY AND ROBERT W. BRADSHAW "BITSLICING AND THE METHOD OF FOUR
//    RUSSIANS OVER LARGER FINITE FIELDS"
//
// William Andrew Alumbaugh 


#include <stdio.h>
#include<stddef.h>
#define true  1 
#define false 0
#define fn(a, b, c, d)  ((a^b)&(c^d))   //for finding R[0] (the first half of the value representingthe sum of vectory and vectorx, vectorr)
#define st(a, b , c)    (a^b^c)   //performing the  (S=  x[0]  XOR  y[1] XOR [x1]) and (T = x[1] XOR Y[0] XOR Y[1]) operations of addition



typedef union vector{       //defines a 32-bit  bit vector
    struct {
        unsigned int v1  :1;
        unsigned int v2  :1;
        unsigned int v3  :1;
        unsigned int v4  :1;
        unsigned int v5  :1;
        unsigned int v6  :1;
        unsigned int v7  :1;
        unsigned int v8  :1;
        unsigned int v9  :1;
        unsigned int v10 :1;
        unsigned int v11 :1;
        unsigned int v12 :1;
        unsigned int v13 :1;
        unsigned int v14 :1;
        unsigned int v15 :1;
        unsigned int v16 :1;
        unsigned int v17 :1;
        unsigned int v18 :1;
        unsigned int v19 :1;
        unsigned int v20 :1;
        unsigned int v21 :1;
        unsigned int v22 :1;
        unsigned int v23 :1;
        unsigned int v24 :1;
        unsigned int v25 :1;
        unsigned int v26 :1;
        unsigned int v27 :1;
        unsigned int v28 :1;
        unsigned int v29 :1;
        unsigned int v30 :1;
        unsigned int v31 :1;
        unsigned int v32 :1;
    } bit;
    
} vec;


int print(vec a)
{
    printf("\n \n");
    printf("[%d %d][%d %d][%d %d][%d %d] \n" , a.bit.v1, a.bit.v2, a.bit.v3, a.bit.v4, a.bit.v5, a.bit.v6, a.bit.v7, a.bit.v8  );
    printf("[%d %d][%d %d][%d %d][%d %d] \n", a.bit.v9, a.bit.v10, a.bit.v11, a.bit.v12, a.bit.v13, a.bit.v14, a.bit.v15, a.bit.v16);
    printf("[%d %d][%d %d][%d %d][%d %d] \n", a.bit.v17, a.bit.v18, a.bit.v19, a.bit.v20, a.bit.v21, a.bit.v22, a.bit.v23, a.bit.v24);
    printf("[%d %d][%d %d][%d %d][%d %d] \n" ,a.bit.v25, a.bit.v26, a.bit.v27, a.bit.v28,  a.bit.v29, a.bit.v30, a.bit.v31, a.bit.v32);

    printf("\n \n");
    return 0;
   
}



int main(int argc, const char * argv[])
{

    vec x;  //vector x
    vec y;  //vector y 
    vec r;  //vector r
    
    // setting variabes for vector x
    x.bit.v1 = true;
    x.bit.v2 = true;
    x.bit.v3 = false;
    x.bit.v4 = true;
    x.bit.v5 = true;
    x.bit.v6 = false;
    x.bit.v7 = true;
    x.bit.v8 = true;
    x.bit.v9 = true;
    x.bit.v10 = true;
    x.bit.v11 = true;
    x.bit.v12 = true;
    x.bit.v13 = true;
    x.bit.v14 = true;
    x.bit.v15 = true;
    x.bit.v16 = true;
    x.bit.v17 = false;
    x.bit.v18 = false;
    x.bit.v19 = true;
    x.bit.v20 = true;
    x.bit.v21 = false;
    x.bit.v22 = true;
    x.bit.v23 = true;
    x.bit.v24 = false;
    x.bit.v25 = false;
    x.bit.v26 = true;
    x.bit.v27 = false;
    x.bit.v28 = true;
    x.bit.v29 = false;
    x.bit.v30 = false;
    x.bit.v31 = true,
    x.bit.v32 = false;
    
    
    // setting variabes for vector y
    y.bit.v1 = true;  
    y.bit.v2 = true;
    y.bit.v3 = false;
    y.bit.v4 = true;
    y.bit.v5 = true;
    y.bit.v6 = false;
    y.bit.v7 = true;
    y.bit.v8 = true;
    y.bit.v9 = true;
    y.bit.v10 = true;
    y.bit.v11 = true;
    y.bit.v12 = true;
    y.bit.v13 = true;
    y.bit.v14 = true;
    y.bit.v15 = true;
    y.bit.v16 = true;
    y.bit.v17 = false;
    y.bit.v18 = false;
    y.bit.v19 = true;
    y.bit.v20 = true;
    y.bit.v21 = false;
    y.bit.v22 = true;
    y.bit.v23 = true;
    y.bit.v24 = false;
    y.bit.v25 = false;
    y.bit.v26 = true;
    y.bit.v27 = false;
    y.bit.v28 = true;
    y.bit.v29 = false;
    y.bit.v30 = false;
    y.bit.v31 = true;
    y.bit.v32 = false;
    
    
    
    
    
    
    
    
    // setting variabes for vector r
    
    r.bit.v1 = false;
    r.bit.v2 = false;
    r.bit.v3 = false;
    r.bit.v4 = false;
    r.bit.v5 = false;
    r.bit.v6 = false;
    r.bit.v7 = false;
    r.bit.v8 = false;
    r.bit.v9 = false;
    r.bit.v10 = false;
    r.bit.v11 = false;
    r.bit.v12 = false;      
    r.bit.v13 = false;
    r.bit.v14 = false;
    r.bit.v15 = false;
    r.bit.v16 = false;
    r.bit.v17 = false;
    r.bit.v18 = false;
    r.bit.v19 = false;
    r.bit.v20 = false;
    r.bit.v21 = false;
    r.bit.v22 = false;
    r.bit.v23 = false;
    r.bit.v24 = false;
    r.bit.v25 = false;
    r.bit.v26 = false;
    r.bit.v27 = false;
    r.bit.v28 = false;
    r.bit.v29 = false;
    r.bit.v30 = false;
    r.bit.v31 = false;
    r.bit.v32 = false;
    
    
    
    
    r.bit.v1 = fn(x.bit.v1, y.bit.v2, x.bit.v2, y.bit.v1);  //for representing  
    r.bit.v3 = fn(x.bit.v3, y.bit.v4, x.bit.v4, y.bit.v3);///r0 ← (x0 ⊕y1)∧(x1 ⊕y0)
    r.bit.v5 = fn(x.bit.v5, y.bit.v6, x.bit.v6, y.bit.v5);
    r.bit.v7 = fn(x.bit.v7, y.bit.v8, x.bit.v8, y.bit.v7);
    r.bit.v9 = fn(x.bit.v9, y.bit.v10, x.bit.v10, y.bit.v9);
    r.bit.v11 = fn(x.bit.v11, y.bit.v12, x.bit.v12, y.bit.v11);
    r.bit.v13 = fn(x.bit.v13, y.bit.v14, x.bit.v14, y.bit.v13);
    r.bit.v15 = fn(x.bit.v15, y.bit.v16, x.bit.v16, y.bit.v15);
    r.bit.v17 = fn(x.bit.v17, y.bit.v18, x.bit.v18, y.bit.v17);
    r.bit.v19 = fn(x.bit.v19, y.bit.v20, x.bit.v20, y.bit.v19);
    r.bit.v21 = fn(x.bit.v21, y.bit.v22, x.bit.v22, y.bit.v21);
    r.bit.v23 = fn(x.bit.v23, y.bit.v24, x.bit.v24, y.bit.v23);
    r.bit.v25 = fn(x.bit.v25, y.bit.v26, x.bit.v26, y.bit.v25);
    r.bit.v27 = fn(x.bit.v27, y.bit.v28, x.bit.v28, y.bit.v27);
    r.bit.v29 = fn(x.bit.v29, y.bit.v30, x.bit.v30, y.bit.v29);
    r.bit.v31 = fn(x.bit.v31, y.bit.v32, x.bit.v32, y.bit.v31);
    
    
    
    r.bit.v2 = (st(x.bit.v1, y.bit.v2, x.bit.v2)      |   st(x.bit.v2, y.bit.v1, y.bit.v2));
    r.bit.v4 = (st(x.bit.v3, y.bit.v4, x.bit.v4)      |   st(x.bit.v4, y.bit.v3, y.bit.v4)); //for representing
    r.bit.v6 = (st(x.bit.v5, y.bit.v6, x.bit.v6)      |   st(x.bit.v6, y.bit.v5, y.bit.v6));// r1 ← s XOR t.
    r.bit.v8 = (st(x.bit.v7, y.bit.v8, x.bit.v8)      |   st(x.bit.v8, y.bit.v7, y.bit.v8));
    r.bit.v10 = (st(x.bit.v9, y.bit.v10, x.bit.v10)   |   st(x.bit.v10, y.bit.v9, y.bit.v10));
    r.bit.v12 = (st(x.bit.v11, y.bit.v12, x.bit.v12)  |   st(x.bit.v12, y.bit.v11, y.bit.v12));
    r.bit.v14 = (st(x.bit.v13, y.bit.v14, x.bit.v14)  |   st(x.bit.v14, y.bit.v13, y.bit.v14));
    r.bit.v16 = (st(x.bit.v15, y.bit.v16, x.bit.v16)  |   st(x.bit.v16, y.bit.v15, y.bit.v16));
    r.bit.v18 = (st(x.bit.v17, y.bit.v18, x.bit.v18)  |   st(x.bit.v18, y.bit.v17, y.bit.v18));
    r.bit.v20 = (st(x.bit.v19, y.bit.v20, x.bit.v20)  |   st(x.bit.v20, y.bit.v19, y.bit.v20));
    r.bit.v22 = (st(x.bit.v21, y.bit.v22, x.bit.v22)  |   st(x.bit.v22, y.bit.v21, y.bit.v22));
    r.bit.v24 = (st(x.bit.v23, y.bit.v24, x.bit.v24)  |   st(x.bit.v24, y.bit.v23, y.bit.v24));
    r.bit.v26 = (st(x.bit.v25, y.bit.v26, x.bit.v26)  |   st(x.bit.v26, y.bit.v25, y.bit.v26));
    r.bit.v28 = (st(x.bit.v27, y.bit.v28, x.bit.v28)  |   st(x.bit.v28, y.bit.v27, y.bit.v28));
    r.bit.v30 = (st(x.bit.v29, y.bit.v30, x.bit.v30)  |   st(x.bit.v30, y.bit.v29, y.bit.v30));
    r.bit.v32 = (st(x.bit.v31, y.bit.v32, x.bit.v32)      |   st(x.bit.v32, y.bit.v31, y.bit.v32));
    
    
    
    
    // testing if bit packing was done properly
    int xsize = sizeof(x);
    int ysize = sizeof(y);
    int rsize  = sizeof(r);
    printf("x is %d bytes wide, y is %d bytes wide, r is %d bytes wide \n \n" ,xsize, ysize, rsize);//prints the size of the bit vectors, should be 4 bytes
    
    
    
    printf("-------------------");  //print the 3 matrices
    printf("\n");
    printf("    Matrix x \n");
    printf("-------------------");
    print(x);
    printf("-------------------");
    printf("\n");
    printf("    Matrix y \n");
    printf("-------------------");
    print(y);
    printf("-------------------");
    printf("\n");
    printf("    Matrix r \n");
    printf("-------------------");
    print(r);

    
    
    return 0;
}

