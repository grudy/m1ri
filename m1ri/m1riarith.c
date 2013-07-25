
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
#include "m5d.h"
#include "m7d.h"
void add_vbg(vbg * r, vbg const * x, vbg const * y)

{
    r->units = (x->units ^ y->sign) & (x->sign ^ y->units); // ///r0 ← (x0 ⊕y->1)∧(x1 ⊕y->0);
    r->sign = (ST(x->units, y->sign, x->sign ) | ST(x->sign, y->units, y->sign)); //// r1 ← s XOR t.
    
    
    
    
}



vbg add_m3dr(vbg  x, vbg const y)
{
    vec t;
    x.sign  = y.units ^ x.sign;
    t = (x.sign & x.units) ^ y.sign;
    x.units = (y.units ^ x.units) |  t;
    x.sign = t & x.sign;
    return x;
    
    
    
}

void sub_m3d( vbg *r, vbg const *x, vbg const *y)
{
    r->units = ((x->units^y->units) | (x->sign^y->sign));
    r->sign = (((x->units^y->units)^x->sign)&(y->units ^ x->sign));
    
    
    
}

void vbg_negation(vbg *r)
{


    r->units  = r->sign = r->sign ^ r->units;

}

vbg sub_m3dr(vbg const x, vbg const y)

{
    vbg r;
    r.units = ((x.units^y.units) | (x.sign^y.sign));
    r.sign = (((x.units^y.units)^x.sign)&(y.units ^ x.sign));
    
    return r;
    
}


 void iadd_vbg(vbg *r,vbg  *x)
{
    
    vec t;
    
    t = x->units ^ r->sign;
    r->sign = x->units ^ r->units;
    r->units = x->units ^ r->units;
    r->sign = r->sign & t;
    t = t ^ x->sign;
    r->units = t | r->units;
    
    
    
    
    
}

void isub_m3d(vbg  *r,vbg  *x)
{
    vec t;
    
    r->units = x->units ^ r->units;
    t  = r->units | r->sign;
    t = t ^ x->sign;
    r->sign = r->units;
    r->sign = r->sign & t;
    r->units = t | r->units;
    
    
    
    
}



void  vbg_mul( vbg *r, vbg  *x, vbg  *y)            {
    r->units = y->units ^ x->units ;
    r->sign = (y->sign ^ x->sign) & (r->units);
    
}






vbg vbg_mul_i(vbg const x, vbg const y)
{
    
    vbg r;
    r.units = x.units & y.units;
    r.sign  = (y.sign ^ x.sign) & (r.units);
    
    return r;
    
}



m3d_t * m3d_hadamard(m3d_t const *a, m3d_t const *b)
{
    
    
    m3d_t  * c = malloc(sizeof(m3d_t));
    
    if((a->nrows == b->nrows) && ( b->ncols == a->ncols))
    {
        *c = m3d_create(c, a->nrows , b->ncols);
        int i, j;
        
        for( i = 0; i < a->nrows; i++)
        {
            for(j = 0; j < (a->width ); j++)
            {
                
                c->rows[i][j] = vbg_mul_i(a->rows[i][j], b->rows[i][j]);
            }
            
            
        }
        
    }
    
    

    return c;
    
    
}




void sub_64_m3d(vbg **R, vbg  **A, vbg  **B)
{
    int i;
    for (i= 0; i < m1ri_word; i++ )
    {
        R[i][0] = sub_m3dr(A[i][0], B[i][0]);
    }
}



void add_64_m3d(vbg **R, vbg   **A, vbg  **B)
{
    int i;
    for (i = 0; i < m1ri_word; i++ )
    {
        R[i][0] = add_m3dr(A[i][0], B[i][0]);
    }
    
    
    
    
}
void m3d_sub( m3d_t *r, m3d_t const *x, m3d_t const *y)
{
    int n , i;
    for(i = 0; i < x->nrows; i++)
    {
        for(n = 0; n < x->width; n++)
        {
    sub_m3d(&r->rows[i][n], &x->rows[i][n], &y->rows[i][n]);
        }

    }

}

m3d_t m3d_add(m3d_t  *a, m3d_t  *b)
{
    
    
    m3d_t  c;
    if((a->nrows == b->nrows) && ( b->ncols == a->ncols))
    {
        c = m3d_create(&c, a->nrows , b->ncols);
        int i, j;
        
        for( i = 0; i < a->nrows; i++)
        {
            for(j = 0; j < (a->width ); j++)
            {
              
            add_vbg(&c.rows[i][j], &a->rows[i][j], &b->rows[i][j]);
            
            }
            
            
        }
        
    }
    
    return c;
    
    
}


