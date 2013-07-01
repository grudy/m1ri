// GF3 arithmatic
//  m1rielarith.h
//  m1riproject

//  Copyright (c) 2013 William Alumbaugh. All rights reserved.
//

#ifndef M1RIPROJECT_M1RIELARITH_H
#define M1RIPROJECT_M1RIELARITH_H

#include "m1ri_3dt.h"

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



void  vbg_mul( vbg *r, vbg *x, vbg *y)             //multiply matrix x by y assinging the output to r
{
    r->units = y->units ^ x->units ;
    r->sign = (y->sign ^ x->sign) & (r->units);
    
}




//return the value of the matrix multiplied


vbg vbg_mul_i(vbg x, vbg y)
{
    
    vbg r;
    r.units = x.units & y.units;
    r.sign  = (y.sign ^ x.sign) & (r.units);
    
    return r;
    
}



/*
    Hadamard multiplication
*/
m3d_t m3d_hadamard(m3d_t *a, m3d_t *b)
{
    
    
    m3d_t c;
    if((a->nrows == b->nrows) && ( b->ncols == a->ncols))
    {
         c = *a;
        int i, j;
        for( i = 0; i < a->nrows; i++)
        {
            for(j = 0; j < a->width; i++)
            {
            
                c.rows[i][j] = vbg_mul_i(a->rows[i][j], b->rows[i][j]);
            }
                
                
        }

        }
    
    return c;

    
}


/* * * * * * * * * * * * * * * * * * * *
 Subtract a 1 kilobyte Matrix from another
 1 kilobyte Matrix
 * * * * * * * * * * * * * * * * * * * * */



void sub_64gf3(vbg *R, vbg *A, vbg *B)
{
    int i;
    for (i = 0; i < (sizeof(vec)); i++ )
    {
        R[i] = subgf3r(A[i], B[i]);
    }
}


/* * * * * * * * * * * * * * * * * * * * * *
 Add a 1 kilobyte Matrix from another
 1 kilobyte Matrix
 * * * * * * * * * * * * * * * * * * * * * * */

void add_64gf3(vbg *R, vbg *A, vbg *B)
{
    int i;
    for (i = 0; i < (sizeof(vec)); i++ )
    {
        R[i] = addgf3r(A[i], B[i]);
    }
    
    
    
    
}






#endif
