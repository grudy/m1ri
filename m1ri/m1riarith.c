
/*
 Matrix Represenations and basic operations
 TOMAS J. BOOTHBY AND ROBERT W. BRADSHAW "BITSLICING AND THE METHOD OF FOUR
 RUSSIANS OVER LARGER FINITE FIELDS"
 
 Copyright 2013 William Andrew Alumbaugh <williamandrewalumbaugh@gmail.com>
 
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
 
 m1ri_arith.c
 */

#include "m1riarith.h"
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

void subgf3( vbg *r, vbg *x, vbg *y)            
{
    r->units = ((x->units^y->units) | (x->sign^y->sign));
    r->sign = (((x->units^y->units)^x->sign)&(y->units ^ x->sign));
    
    
    
}



vbg subgf3r(vbg x, vbg y)             

{
    vbg r;
    r.units = ((x.units^y.units) | (x.sign^y.sign));
    r.sign = (((x.units^y.units)^x.sign)&(y.units ^ x.sign));
    
    return r;
    
}



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

void isubgf3(vbg *r,vbg *x)  
{
    vec t;
    
    r->units = x->units ^ r->units;
    t  = r->units | r->sign;
    t = t ^ x->sign;
    r->sign = r->units;
    r->sign = r->sign & t;
    r->units = t | r->units;
    
    
    
    
}



void  vbg_mul( vbg *r, vbg *x, vbg *y)            {
    r->units = y->units ^ x->units ;
    r->sign = (y->sign ^ x->sign) & (r->units);
    
}






vbg vbg_mul_i(vbg x, vbg y)
{
    
    vbg r;
    r.units = x.units & y.units;
    r.sign  = (y.sign ^ x.sign) & (r.units);
    
    return r;
    
}




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




void sub_64gf3(vbg *R, vbg *A, vbg *B)
{
    int i;
    for (i
         
         = 0; i < (sizeof(vec)); i++ )
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