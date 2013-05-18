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






typedef union vector{ //defines a 64-bit bit vector
    
    unsigned long long v;
    
    
    struct vectorbits{
        unsigned long v1 :1;
        unsigned int v2 :1;
        unsigned int v3 :1;
        unsigned int v4 :1;
        unsigned int v5 :1;
        unsigned int v6 :1;
        unsigned int v7 :1;
        unsigned int v8 :1;
        unsigned int v9 :1;
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
        unsigned int v64 :1;
    } bit;
    
} vec;

typedef union{  //calls a union of 128 bits
    
    vec units;
    vec sign;
} matrixgf3;



matrixgf3 addgf3(matrixgf3 x, matrixgf3 y, matrixgf3 r)

{
    r.units.v = (x.units.v ^ y.sign.v) & (x.sign.v ^ y.units.v); // ///r0 ← (x0 ⊕y1)∧(x1 ⊕y0);
    r.sign.v = (st(x.units.v, y.sign.v, x.sign.v ) | st(x.units.v, y.units.v, y.sign.v)); //// r1 ← s XOR t.
    
    
    return r;
}


void subgf3(matrixgf3 *x, matrixgf3 *y, matrixgf3 *r)               //multiply matrix x by by matrix y.   The product is matrix r.

{
    r->units.v = ((x->units.v^y->units.v) | (x->sign.v^y->sign.v));
    r->sign.v = (((x->units.v^y->units.v)^x->sign.v)&(y->units.v ^ x->sign.v));
    
    

}
/*
 
 cdef inline vec3 v3_subc(vec3 a, vec3 b):
 cdef long t
 cdef vec3 r
 
 r.s = b.u ^ a.u
 r.u = b.s ^ a.s
 r.u = r.u | r.s
 r.s = r.s ^ b.s
 t = b.u ^ a.s
 r.s = t & r.s
 
 return r

 
*/

matrixgf3 subgf3r(matrixgf3 x, matrixgf3 y)               //multiply matrix x by by matrix y.   The product is matrix r.

{
    matrixgf3 r;
    r.units.v = ((x.units.v^y.units.v) | (x.sign.v^y.sign.v));
    r.sign.v = (((x.units.v^y.units.v)^x.sign.v)&(y.units.v ^ x.sign.v));
    
    return r;
    
}




/*
 r.s = b.s ^ a.s
 r.s = r.s & r.u
*/



void  mplygf3(matrixgf3 *x, matrixgf3 *y, matrixgf3 *r)             //multiply matrix x by y saving the output to r
{
    r->units.v = y->units.v ^ x->units.v ;
    r->sign.v = (y->sign.v ^ x->sign.v) & (r->units.v);

}













#endif
