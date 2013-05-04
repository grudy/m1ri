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
#define fn(a, b, c, d)  (a^b^c^d)   //for finding R[0] (the first half of the value representingthe sum of vectory and vectorx, vectorr)
#define st(a, b , c)    (a^b^c)   //performing the  (S=  x[0]  XOR  y[1] XOR [x1]) and (T = x[1] XOR Y[0] XOR Y[1]) operations of addition

union vectory{  //Bit vector of 32 bits
    
    unsigned int y1 :1;
    unsigned int y2 :1;
    unsigned int y3 :1;
    unsigned int y4 :1;
    unsigned int y5 :1;
    unsigned int y6 :1;
    unsigned int y7 :1;
    unsigned int y8 :1;
    unsigned int y9 :1;
    unsigned int y10 :1;
    unsigned int y11 :1;
    unsigned int y12 :1;
    unsigned int y13 :1;
    unsigned int y14 :1;
    unsigned int y15 :1;
    unsigned int y16 :1;
    unsigned int y17 :1;
    unsigned int y18 :1;
    unsigned int y19 :1;
    unsigned int y20 :1;
    unsigned int y21 :1;
    unsigned int y22 :1;
    unsigned int y23 :1;
    unsigned int y24 :1;
    unsigned int y25 :1;
    unsigned int y26 :1;
    unsigned int y27 :1;
    unsigned int y28 :1;
    unsigned int y29 :1;
    unsigned int y30 :1;
    unsigned int y31 :1;
    unsigned int y32 :1;
    
    
    
}yvec;   //tag of vectory


union vectorr{
    
    unsigned int r1  :1;
    unsigned int r2  :1;
    unsigned int r3  :1;
    unsigned int r4  :1;
    unsigned int r5  :1;
    unsigned int r6  :1;
    unsigned int r7  :1;
    unsigned int r8  :1;
    unsigned int r9  :1;
    unsigned int r10 :1;
    unsigned int r11 :1;
    unsigned int r12 :1;
    unsigned int r13 :1;
    unsigned int r14 :1;
    unsigned int r15 :1;
    unsigned int r16 :1;
    unsigned int r17 :1;
    unsigned int r18 :1;
    unsigned int r19 :1;
    unsigned int r20 :1;
    unsigned int r21 :1;
    unsigned int r22 :1;
    unsigned int r23 :1;
    unsigned int r24 :1;
    unsigned int r25 :1;
    unsigned int r26 :1;
    unsigned int r27 :1;
    unsigned int r28 :1;
    unsigned int r29 :1;
    unsigned int r30 :1;
    unsigned int r31 :1;
    unsigned int r32 :1;
    
    
   }rvec; //tag of vectorr
    
    
    union vectorx {
        
        unsigned int x1 :1;
        unsigned int x2  :1;
        unsigned int x3  :1;
        unsigned int x4  :1;
        unsigned int x5  :1;
        unsigned int x6  :1;
        unsigned int x7  :1;
        unsigned int x8  :1;
        unsigned int x9  :1;
        unsigned int x10 :1;
        unsigned int x11 :1;
        unsigned int x12 :1;
        unsigned int x13 :1;
        unsigned int x14 :1;
        unsigned int x15 :1;
        unsigned int x16 :1;
        unsigned int x17 :1;
        unsigned int x18 :1;
        unsigned int x19 :1;
        unsigned int x20 :1;
        unsigned int x21 :1;
        unsigned int x22 :1;
        unsigned int x23 :1;
        unsigned int x24 :1;
        unsigned int x25 :1;
        unsigned int x26 :1;
        unsigned int x27 :1;
        unsigned int x28 :1;
        unsigned int x29 :1;
        unsigned int x30 :1;
        unsigned int x31 :1;
        unsigned int x32 :1;
        
    }xvec;   //tag of vectorx 
    
    



int main(int argc, const char * argv[])
{
    
    //declaring  values for the individual bits 
    
    xvec.x1 = true;
    xvec.x2 = true;
    xvec.x3 = false;
    xvec.x4 = true;
    xvec.x5 = true;
    xvec.x6 = false;
    xvec.x7 = true;
    xvec.x8 = true;
    xvec.x9 = true;
    xvec.x10 = true;
    xvec.x11 = true;
    xvec.x12 = true;
    xvec.x13 = true;
    xvec.x14 = true;
    xvec.x15 = true;
    xvec.x16 = true;
    xvec.x17 = false;
    xvec.x18 = false;
    xvec.x19 = true;
    xvec.x20 = true;
    xvec.x21 = false;
    xvec.x22 = true;
    xvec.x23 = true;
    xvec.x24 = false;
    xvec.x25 = false;
    xvec.x26 = true;
    xvec.x27 = false;
    xvec.x28 = true;
    xvec.x29 = false;
    xvec.x30 = false;
    xvec.x31 = true,
    xvec.x32 = false;
    

    yvec.y1 = true;
    yvec.y2 = true;
    yvec.y3 = false;
    yvec.y4 = true;
    yvec.y5 = true;
    yvec.y6 = false;
    yvec.y7 = true;
    yvec.y8 = true;
    yvec.y9 = true;
    yvec.y10 = true;
    yvec.y11 = true;
    yvec.y12 = true;
    yvec.y13 = true;
    yvec.y14 = true;
    yvec.y15 = true;
    yvec.y16 = true;
    yvec.y17 = false;
    yvec.y18 = false;
    yvec.y19 = true;
    yvec.y20 = true;
    yvec.y21 = false;
    yvec.y22 = true;
    yvec.y23 = true;
    yvec.y24 = false;
    yvec.y25 = false;
    yvec.y26 = true;
    yvec.y27 = false;
    yvec.y28 = true;
    yvec.y29 = false;
    yvec.y30 = false;
    yvec.y31 = true;
    yvec.y32 = false;
    
    
    
    
   
    
    

    
    
    rvec.r1 = false;
    rvec.r2 = false;
    rvec.r3 = false;
    rvec.r4 = false;
    rvec.r5 = false;
    rvec.r6 = false;
    rvec.r7 = false;
    rvec.r8 = false;
    rvec.r9 = false;
    rvec.r10 = false;
    rvec.r11 = false;
    rvec.r12 = false;
    rvec.r13 = false;
    rvec.r14 = false;
    rvec.r15 = false;
    rvec.r16 = false;
    rvec.r17 = false;
    rvec.r18 = false;
    rvec.r19 = false;
    rvec.r20 = false;
    rvec.r21 = false;
    rvec.r22 = false;
    rvec.r23 = false;
    rvec.r24 = false;
    rvec.r25 = false;
    rvec.r26 = false;
    rvec.r27 = false;
    rvec.r28 = false;
    rvec.r29 = false;
    rvec.r30 = false;
    rvec.r31 = false;
    rvec.r32 = false;
    
    
    
    //Here the arithmatic starts
    // r[0]= (x[0], y[1], xvec.x2, yvec.y1);
    rvec.r1 = fn(xvec.x1, yvec.y2, xvec.x2, yvec.y1);    
    rvec.r3 = fn(xvec.x3, yvec.y4, xvec.x4, yvec.y3);
    rvec.r5 = fn(xvec.x5, yvec.y6, xvec.x6, yvec.y5);
    rvec.r7 = fn(xvec.x7, yvec.y8, xvec.x8, yvec.y7);
    rvec.r9 = fn(xvec.x9, yvec.y10, xvec.x10, yvec.y9);
    rvec.r11 = fn(xvec.x11, yvec.y12, xvec.x12, yvec.y11);
    rvec.r13 = fn(xvec.x13, yvec.y14, xvec.x14, yvec.y13);
    rvec.r15 = fn(xvec.x15, yvec.y16, xvec.x16, yvec.y15);
    rvec.r17 = fn(xvec.x17, yvec.y18, xvec.x18, yvec.y17);
    rvec.r19 = fn(xvec.x19, yvec.y20, xvec.x20, yvec.y19);
    rvec.r21 = fn(xvec.x21, yvec.y22, xvec.x22, yvec.y21);
    rvec.r23 = fn(xvec.x23, yvec.y24, xvec.x24, yvec.y23);
    rvec.r25 = fn(xvec.x25, yvec.y26, xvec.x26, yvec.y25);
    rvec.r27 = fn(xvec.x27, yvec.y28, xvec.x28, yvec.y27);
    rvec.r29 = fn(xvec.x29, yvec.y30, xvec.x30, yvec.y29);
    rvec.r31 = fn(xvec.x31, yvec.y32, xvec.x32, yvec.y31);
    
    
    
    
    
    //    r[1] = s OR t
    
    rvec.r2 = (st(xvec.x1, yvec.y2, xvec.x2)      ||   st(xvec.x2, yvec.y1, yvec.y2));
    rvec.r4 = (st(xvec.x3, yvec.y4, xvec.x4)      ||   st(xvec.x4, yvec.y3, yvec.y4)); //
    rvec.r6 = (st(xvec.x5, yvec.y6, xvec.x6)      ||   st(xvec.x6, yvec.y5, yvec.y6));
    rvec.r8 = (st(xvec.x7, yvec.y8, xvec.x8)      ||   st(xvec.x8, yvec.y7, yvec.y8));
    rvec.r10 = (st(xvec.x9, yvec.y10, xvec.x10)   ||   st(xvec.x10, yvec.y9, yvec.y10));
    rvec.r12 = (st(xvec.x11, yvec.y12, xvec.x12)  ||   st(xvec.x12, yvec.y11, yvec.y12));
    rvec.r14 = (st(xvec.x13, yvec.y14, xvec.x14)  ||   st(xvec.x14, yvec.y13, yvec.y14));
    rvec.r16 = (st(xvec.x15, yvec.y16, xvec.x16)  ||   st(xvec.x16, yvec.y15, yvec.y16));
    rvec.r18 = (st(xvec.x17, yvec.y18, xvec.x18)  ||   st(xvec.x18, yvec.y17, yvec.y18));
    rvec.r20 = (st(xvec.x19, yvec.y20, xvec.x20)  ||   st(xvec.x20, yvec.y19, yvec.y20));
    rvec.r22 = (st(xvec.x21, yvec.y22, xvec.x22)  ||   st(xvec.x22, yvec.y21, yvec.y22));
    rvec.r24 = (st(xvec.x23, yvec.y24, xvec.x24)  ||   st(xvec.x24, yvec.y23, yvec.y24));
    rvec.r26 = (st(xvec.x25, yvec.y26, xvec.x26)  ||   st(xvec.x26, yvec.y25, yvec.y26));
    rvec.r28 = (st(xvec.x27, yvec.y28, xvec.x28)  ||   st(xvec.x28, yvec.y27, yvec.y28));
    rvec.r30 = (st(xvec.x29, yvec.y30, xvec.x30)  ||   st(xvec.x30, yvec.y29, yvec.y30));
    rvec.r32 = (st(xvec.x31, yvec.y32, xvec.x32)  ||   st(xvec.x32, yvec.y31, yvec.y32));
    
    
    // testing if bit packing was done properly 
    int xsize = sizeof(xvec);
    int ysize = sizeof(yvec);
    int rsize  = sizeof(rvec);    
    printf("xvec is %d bytes wide, yvec is %d bytes wide, rvec is %d bytes wide" ,xsize, ysize, rsize);




}



