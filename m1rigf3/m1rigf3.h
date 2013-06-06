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
#include <stdio.h>




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

typedef struct{  //calls a union of 128 bits
    
    vec units;
    vec sign;
} vbg;


int print(vec a, vec b)
{
    printf("\n \n");
    printf("[%d %d][%d %d][%d %d][%d %d][%d %d][%d %d][%d %d][%d %d] " , a.bit.v0, b.bit.v0, a.bit.v1, b.bit.v1, a.bit.v2 , b.bit.v2, a.bit.v3,b.bit.v3, a.bit.v4,b.bit.v4, a.bit.v5,b.bit.v5, a.bit.v6,b.bit.v6, a.bit.v7, b.bit.v7 );
    printf("[%d %d][%d %d][%d %d][%d %d][%d %d][%d %d][%d %d][%d %d] " ,  a.bit.v8, b.bit.v8,  a.bit.v9, b.bit.v9, a.bit.v10 , b.bit.v10, a.bit.v11,b.bit.v11, a.bit.v12,b.bit.v12, a.bit.v13,b.bit.v13, a.bit.v14,b.bit.v14, a.bit.v15, b.bit.v15);
    printf("[%d %d][%d %d][%d %d][%d %d][%d %d][%d %d][%d %d][%d %d] ", a.bit.v16, b.bit.v16, a.bit.v17, b.bit.v17, a.bit.v18 , b.bit.v18, a.bit.v19,b.bit.v19, a.bit.v20,b.bit.v20, a.bit.v21,b.bit.v21, a.bit.v22,b.bit.v22, a.bit.v23, b.bit.v23);
    printf("[%d %d][%d %d][%d %d][%d %d][%d %d][%d %d][%d %d][%d %d] ", a.bit.v24, b.bit.v24 , a.bit.v25, b.bit.v25, a.bit.v26 , b.bit.v26, a.bit.v27,b.bit.v27, a.bit.v28,b.bit.v28, a.bit.v29,b.bit.v29, a.bit.v30,b.bit.v30, a.bit.v31, b.bit.v31);
    printf("[%d %d][%d %d][%d %d][%d %d][%d %d][%d %d][%d %d][%d %d] ", a.bit.v32, b.bit.v32 , a.bit.v33, b.bit.v33, a.bit.v34 , b.bit.v34, a.bit.v35,b.bit.v35, a.bit.v36,b.bit.v36, a.bit.v37,b.bit.v37, a.bit.v38,b.bit.v38, a.bit.v39, b.bit.v39);
    printf("[%d %d][%d %d][%d %d][%d %d][%d %d][%d %d][%d %d][%d %d] " , a.bit.v40, b.bit.v40 , a.bit.v41, b.bit.v41, a.bit.v42 , b.bit.v42, a.bit.v43,b.bit.v43, a.bit.v44,b.bit.v44, a.bit.v45,b.bit.v45, a.bit.v46,b.bit.v46, a.bit.v47, b.bit.v47 );
    printf("[%d %d][%d %d][%d %d][%d %d][%d %d][%d %d][%d %d][%d %d] ", a.bit.v48, b.bit.v48, a.bit.v49, b.bit.v49, a.bit.v50 , b.bit.v50, a.bit.v51,b.bit.v51, a.bit.v52,b.bit.v52, a.bit.v53,b.bit.v53, a.bit.v54,b.bit.v54, a.bit.v55, b.bit.v55);
    printf("[%d %d][%d %d][%d %d][%d %d][%d %d][%d %d][%d %d][%d %d] ", a.bit.v56, b.bit.v56 , a.bit.v57, b.bit.v57, a.bit.v58 , b.bit.v58, a.bit.v59,b.bit.v59, a.bit.v60, b.bit.v60, a.bit.v61, b.bit.v61, a.bit.v62, b.bit.v62, a.bit.v63, b.bit.v63  );
    
    
    printf("\n \n");
    return 0;
    
}



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
