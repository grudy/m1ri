//// GF3 arithmatic
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
 
    m1rigf3.h
 */

#ifndef m1rigf3_m1rigf3_h
#define m1rigf3_m1rigf3_h
#define fn(a, b, c, d) (a^b)&(c^d) //for finding R[0]# (the first half of the value representingthe sum of vectory and vectorx, vectorr)
#define st(a, b , c) (a^b^c) //performing the (S= x[0] XOR y[1] XOR [x1]) and (T = x[1] XOR Y[0] XOR Y[1]) operations of addition
#include <stdlib.h>





typedef union vector{ //defines a 64-bit bit vector
    
    unsigned long long v;
    
    
    struct vectorbits{
        unsigned int v0 :1;
        unsigned int v1 :1;
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
        unsigned int v33 :1;
        unsigned int v34 :1;
        unsigned int v35 :1;
        unsigned int v36 :1;
        unsigned int v37 :1;
        unsigned int v38 :1;
        unsigned int v39 :1;
        unsigned int v40 :1;
        unsigned int v41 :1;
        unsigned int v42 :1;
        unsigned int v43 :1;
        unsigned int v44 :1;
        unsigned int v45 :1;
        unsigned int v46 :1;
        unsigned int v47 :1;
        unsigned int v48 :1;
        unsigned int v49 :1;
        unsigned int v50 :1;
        unsigned int v51 :1;
        unsigned int v52 :1;
        unsigned int v53 :1;
        unsigned int v54 :1;
        unsigned int v55 :1;
        unsigned int v56 :1;
        unsigned int v57 :1;
        unsigned int v58 :1;
        unsigned int v59 :1;
        unsigned int v60 :1;
        unsigned int v61 :1;
        unsigned int v62 :1;
        unsigned int v63 :1;
      
    } bit;
    
} vec;

typedef union{  //calls a union of 128 bits
    
    vec units;
    vec sign;
} vbg;






/*
vbg ** creatematrix(nrows, ncols, gfs)   //rows, columns, field size
{
    vbg ** matrix;
    matrix = malloc(nrows *  sizeof(vbg *));
    if(matrix == 0)
    {
        //error message to implement later 
        return 0;
    }

for(int i = 0; i < nrows; i++)
    {
        matrix[i] = malloc(ncols * sizeof(vbg));
        
    }

    return matrix;

}
*/  //first design of a matrix







void addgf3(vbg * r, vbg * x, vbg * y)

{
    r->units.v = (x->units.v ^ y->sign.v) & (x->sign.v ^ y->units.v); // ///r0 ← (x0 ⊕y->1)∧(x1 ⊕y->0);
    r->sign.v = (st(x->units.v, y->sign.v, x->sign.v ) | st(x->units.v, y->units.v, y->sign.v)); //// r1 ← s XOR t.
    
    
   
}

/*
 
 cdef inline vec3 v3_addc(vec3 a, vec3 b):
 cdef long t
 
 a.s = b.u ^ a.s
 t = a.s & a.u
 t = t ^ b.s
 a.u = b.u ^ a.u
 a.u = a.u | t
 a.s = t & a.s
 
 return a
 */

vbg addgf3r(vbg  x, vbg y)
{
    vec t;
    x.sign.v  = y.units.v ^ x.sign.v;
    t.v = (x.sign.v & x.units.v) ^ y.sign.v;
    x.units.v = (y.units.v ^ x.units.v) |  t.v;
    x.sign.v = t.v & x.sign.v;
    return x;
    
    
    
}

void subgf3( vbg *r, vbg *x, vbg *y)               //multiply matrix x by by matrix y.   The product is matrix r.

{
    r->units.v = ((x->units.v^y->units.v) | (x->sign.v^y->sign.v));
    r->sign.v = (((x->units.v^y->units.v)^x->sign.v)&(y->units.v ^ x->sign.v));
    
    

}



vbg subgf3r(vbg x, vbg y)               //multiply matrix x by by matrix y.   The product is matrix r.

{
    vbg r;
    r.units.v = ((x.units.v^y.units.v) | (x.sign.v^y.sign.v));
    r.sign.v = (((x.units.v^y.units.v)^x.sign.v)&(y.units.v ^ x.sign.v));
    
    return r;
    
}








void  mplygf3( vbg *r, vbg *x, vbg *y)             //multiply matrix x by y assinging the output to r
{
    r->units.v = y->units.v ^ x->units.v ;
    r->sign.v = (y->sign.v ^ x->sign.v) & (r->units.v);

}






vbg mplygf3r(vbg x, vbg y)    //return the value of the matrix 
{
    
    vbg r;
    r.units.v = y.units.v & y.units.v;
    r.sign.v  = (y.sign.v ^ x.sign.v) & (r.units.v);
    
    return r;
    
}


void iaddgf3(vbg *r,vbg *x)  ////matrix r = (direct sum matrix r + matrix x) 
{
    
    vec t;
    
    t.v = x->units.v ^ r->sign.v;
    r->sign.v = x->units.v ^ r->units.v;
    r->units.v = x->units.v ^ r->units.v;
    r->sign.v = r->sign.v & t.v;
    t.v = t.v ^ x->sign.v;
    r->units.v = t.v | r->units.v;
    
    
    
    
    
}


void isubgf3(vbg *r,vbg *x)  //matrix r = (matrix r - matrix x)
{
    vec t;
    
    r->units.v = x->units.v ^ r->units.v;
    t .v = r->units.v | r->sign.v;
    t.v = t.v ^ x->sign.v;
    r->sign.v = r->units.v;
    r->sign.v = r->sign.v & t.v;
    r->units.v = t.v | r->units.v;
    
    
    
    
}








#endif
