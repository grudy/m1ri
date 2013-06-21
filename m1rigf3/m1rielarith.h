// GF3 arithmatic
//  m1rielarith.h
//  m1riproject
//
//  Created by grudy on 6/11/13.
//  Copyright (c) 2013 William Alumbaugh. All rights reserved.
//

#ifndef m1riproject_m1rielarith_h
#define m1riproject_m1rielarith_h

#include "m1rigf3.h"

void addgf3(vbg * r, vbg * x, vbg * y)

{
    r->units = (x->units ^ y->sign) & (x->sign ^ y->units); // ///r0 ← (x0 ⊕y->1)∧(x1 ⊕y->0);
    r->sign = (ST(x->units, y->sign, x->sign ) | ST(x->units, y->units, y->sign)); //// r1 ← s XOR t.
    
    
    
    
}



vbg addgf3r(vbg  x, vbg y)
{
    vec t;
    x.sign  = y.units ^ x.sign;
    t = (x.sign & x.units) ^ y.sign;
    x.units = (y.units ^ x.units) |  t;
    x.sign = t & x.sign;
    return x;
    
    
    
}

void subgf3( vbg *r, vbg *x, vbg *y)               //multiply matrix x by by matrix y.   The product is matrix r.

{
    r->units = ((x->units^y->units) | (x->sign^y->sign));
    r->sign = (((x->units^y->units)^x->sign)&(y->units ^ x->sign));
    
    
    
}



vbg subgf3r(vbg x, vbg y)               //multiply matrix x by by matrix y.   The product is matrix r.

{
    vbg r;
    r.units = ((x.units^y.units) | (x.sign^y.sign));
    r.sign = (((x.units^y.units)^x.sign)&(y.units ^ x.sign));
    
    return r;
    
}








void  mplygf3( vbg *r, vbg *x, vbg *y)             //multiply matrix x by y assinging the output to r
{
    r->units = y->units ^ x->units ;
    r->sign = (y->sign ^ x->sign) & (r->units);
    
}






vbg mplygf3r(vbg x, vbg y)    //return the value of the matrix
{
    
    vbg r;
    r.units = y.units & y.units;
    r.sign  = (y.sign ^ x.sign) & (r.units);
    
    return r;
    
}



/********************************************
 matrix r = (direct sum matrix r + matrix x)
 ********************************************/
void iaddgf3(vbg *r,vbg *x)
{
    
    vec t;
    
    t = x->units ^ r->sign;
    r->sign = x->units ^ r->units;
    r->units = x->units ^ r->units;
    r->sign = r->sign & t;
    t = t ^ x->sign;
    r->units = t | r->units;
    
    
    
    
    
}

void isubgf3(vbg *r,vbg *x)  //matrix r = (matrix r - matrix x)
{
    vec t;
    
    r->units = x->units ^ r->units;
    t  = r->units | r->sign;
    t = t ^ x->sign;
    r->sign = r->units;
    r->sign = r->sign & t;
    r->units = t | r->units;
    
    
    
    
}


/* * * * * * * * * * * * * * * * * * * *
 Subtract a 1 Megabyte Matrix from another
 1 megabyte Matrix
 * * * * * * * * * * * * * * * * * * * * */

void isub_64gf3(vbg * R, vbg * A)

{
    int i;
    for (i = 0; i < (sizeof(vec)); i++ )
    {
        isubgf3((R+i), (A+i));
    }
}

void sub_64gf3(vbg *R, vbg *A, vbg *B)
{
    int i;
    for (i = 0; i < (sizeof(vec)); i++ )
    {
        R[i] = subgf3r(A[i], B[i]);
    }
}


void add_64gf3(vbg *R, vbg *A, vbg *B)
{
    int i;
    for (i = 0; i < (sizeof(vec)); i++ )
    {
        R[i] = addgf3r(A[i], B[i]);
    }
    
}

//Adds two 64x matrix blocks
void iadd_64gf3(vbg *R, vbg *A)
{
    int i;
    for (i = 0; i < (sizeof(vec)); i++ )
    {
        iaddgf3((R+i), (A+i));
    }
    
}







void isub_128(vbg *R, vbg *A)
{
    int i;
    for (i = 0; i < (4 * (sizeof(vec))); i++ )
    {
        isubgf3((R+i), (A+i));
    }
}



#endif
