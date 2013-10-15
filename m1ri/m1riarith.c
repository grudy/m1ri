
/** 
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



static inline void add_vbg(vbg * restrict  r, vbg const * restrict  x, vbg const * y)

{
    r->units = (x->units ^ y->sign) & (x->sign ^ y->units); // ///r0 ← (x0 ⊕y->1)∧(x1 ⊕y->0);
    r->sign = (M1RI_ST(x->units, y->sign, x->sign ) | M1RI_ST(x->sign, y->units, y->sign)); //// r1 ← s XOR t.
}

/** 

*/
static inline vbg add_m3dr(vbg  x, vbg const y)
{ 
    vec t;
    x.sign  = y.units ^ x.sign;
    t = (x.sign & x.units) ^ y.sign;
    x.units = (y.units ^ x.units) |  t;
    x.sign = t & x.sign;
    return x; 
}

inline void sub_m3d( vbg *r, vbg const *x, vbg const *y)
{
    r->units = ((x->units^y->units) | (x->sign^y->sign));
    r->sign = (((x->units^y->units)^x->sign)&(y->units ^ x->sign));
}

void vbg_negation(vbg *r)
{
     r->sign = r->sign ^ r->units;
}

vbg sub_m3dr(vbg const x, vbg const y)

{
    vbg r;
    r.units = ((x.units^y.units) | (x.sign^y.sign));
    r.sign = (((x.units^y.units)^x.sign)&(y.units ^ x.sign));
    
    return r;
}

static inline void iadd_vbg(vbg *r,vbg  const * restrict x)
{
    vec t;
    t = x->units ^ r->sign;
    r->sign = x->units ^ r->units;
    r->units = x->units ^ r->units;
    r->sign = r->sign & t;
    t = t ^ x->sign;
    r->units = t | r->units;
}

static inline void isub_m3d(vbg  *r,vbg  *x)
{
    vec t;
    r->units = x->units ^ r->units;
    t  = r->units | r->sign;
    t = t ^ x->sign;
    r->sign = x->units ^ r->sign;
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

m3d_t * m3d_hadamard(m3d_t const * restrict a, m3d_t const * restrict b)
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


void m3d_sub_64(vbg **R, vbg  **A, vbg  **B)
{
    int i;
    for (i= 0; i < M1RI_RADIX; i++ )
    {
        R[i][0] = sub_m3dr(A[i][0], B[i][0]);
    }
    
}
/**
	Adds two 64 by 64 m3d_t matrices
*/
void m3d_add_64(vbg **R, vbg   **A, vbg  **B)
{
    int i;
    for (i = 0; i < M1RI_RADIX; i++ )
    {
        R[i][0] = add_m3dr(A[i][0], B[i][0]);
    }

}
/**
	Subtracts m3d_t a  of arbitrary 
	length and height.
	r = (x * y) 
*/

void m3d_sub( m3d_t *r, const  m3d_t  *x, const m3d_t  *y)
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

/**
	Adds two m3d_t's
	If this isn't possible it does nothing
*/
void m3d_add_r(m3d_t * c, m3d_t  * restrict a, m3d_t  * restrict b)
{
    if((a->nrows == b->nrows) && ( b->ncols == a->ncols))
    {
    	int i, j;
        for( i = 0; i < a->nrows; i++)
        {
            for(j = 0; j < (a->width ); j++)
            {
                add_vbg(&c->rows[i][j], &a->rows[i][j], &b->rows[i][j]);        
            }
        }
    }  
}

void *  m3d_combine3(vbg *table, vbg *input )
{
    vbg t, a, b, c;
    t.sign = t.units = 0;
    a = input[0];
    b = input[1];
    c = input[2];
    table[0] = t;
    table[1] = a;
    table[2] = b;
    table[4] = c;

    add_vbg(&t, &a, &b);
    table[3] = t;
    iadd_vbg(&t, &c);
    table[7] = t;
    isub_m3d(&t, &a);
    table[6] = t;
    add_vbg((table + 5), &a , &b);

    return 0;
    
}


void m3d_combine4(vbg *table, vbg *input )
{
    vbg t, a, b, c , d;
    t.sign = t.units = 0;
    a = input[0];
    b = input[1];
    c = input[2];
    d = input[3];
    
    table[0] = t;
    table[1] = a;
    table[2] = b;
    table[4] = c;
    table[8] = d;
    
    add_vbg(&t, &c, &d);
    
    table[12] = t;
    
    add_vbg(&t,&b,&c);
    table[6] = t;
    iadd_vbg(&t,&d);
    table[14] = t;
    isub_m3d(&t,&c);
    table[10] = t;
    
    add_vbg(&t,&b,&c);
    table[3] = t;
    iadd_vbg(&t, &d);
    
    
    
    table[11] = t;
    iadd_vbg(&t, &c);
    table[15] = t;
    isub_m3d(&t, &d);
    table[7] = t;
    isub_m3d(&t, &b);
    table[5] = t;
    iadd_vbg(&t, &d);
    table[13] = t;
    isub_m3d(&t, &c);
    table[9] = t;
    
    
}


void m3d_combine5(vbg *table, vbg *input )
{
	int i;
    vbg e, *t4;

    m3d_combine4(table, input);
    e = input[4];
    t4 = table + 16;
    table[16] = e;
    
    for ( i = 1; i < 16 ; i ++ ) {
        add_vbg(t4 + i, table + i, &e);
    }
}


void m3d_combine6(vbg *table, vbg *input )

{
    vbg f, *t5;
    int i;
    m3d_combine5(table, input);
    f = input[5];
    t5 = (table + 32);
    table [32] = f;
    for (i = 1; i < 32; i++)
        add_vbg((t5 + i), (table + i), &f);
    
}

 void m3d_combine7(vbg *table, vbg *input )

{
    
    vbg g, *t6;
    int i;
    m3d_combine6(table, input);
    g = input[6];
    t6 = (table+64);
    table[64] = g;
    for (i = 1; i < 64; i = i +1) {
        add_vbg((t6 + i), (table + i), &g );
    }
 
}


void m3d_combine8(vbg *table, vbg *input)

{
    vbg h, *t7;
    int i;
    
    m3d_combine7(table, input);
    h = input[7];
    t7 = (table+128);
    table[128] = h;
    
    for (i = 1; i < 128; i++)
        add_vbg((t7 + i), (table+i), &h);
}


void m3d_mul_64(vbg **R, vbg ** restrict A, vbg ** restrict B)
{
    int i;
    vbg t1, t2, r1, r2, a;
    vec v1, v2;
    
    vbg  tables6[9][64];
    vbg tables5[2][32];
    
	for (i = 0; i < 9; i ++)
    {
        m3d_combine6(&tables6[i][0], &(B [6*i][0]));
    }
   
    for (i = 0; i < 2; i ++)
    {
        m3d_combine5(&tables5[i][0], &(B[54 + (5 * i)][0]));
    }

    for (i = 0; i < 64; i ++  )//i from 0 <= i < 64
    {
        a = A[i][0];
        v2 = a.sign;
    
        v1 = (a.units ^ v2);		
        r1 = tables6[0][v1&63];
        v1 >>= 6;
        r2 = tables6[0][v2&63];
        v2 >>= 6;
        t1 = tables6[1][v1&63]; iadd_vbg(&r1, &t1);v1 >>= 6;
        t2 = tables6[1][v2&63]; iadd_vbg(&r2, &t2); v2 >>= 6;
        t1 = tables6[2][v1&63]; iadd_vbg(&r1, &t1); v1 >>= 6;
        t2 = tables6[2][v2&63]; iadd_vbg(&r2, &t2); v2 >>= 6;
        t1 = tables6[3][v1&63]; iadd_vbg(&r1, &t1); v1 >>= 6;
        t2 = tables6[3][v2&63]; iadd_vbg(&r2, &t2); v2 >>= 6;
        t1 = tables6[4][v1&63]; iadd_vbg(&r1, &t1); v1 >>= 6;
        t2 = tables6[4][v2&63]; iadd_vbg(&r2, &t2); v2 >>= 6;
        t1 = tables6[5][v1&63]; iadd_vbg(&r1, &t1); v1 >>= 6;
        t2 = tables6[5][v2&63]; iadd_vbg(&r2, &t2); v2 >>= 6;
        t1 = tables6[6][v1&63]; iadd_vbg(&r1, &t1); v1 >>= 6;
        t2 = tables6[6][v2&63]; iadd_vbg(&r2, &t2); v2 >>= 6;
        t1 = tables6[7][v1&63]; iadd_vbg(&r1, &t1); v1 >>= 6;
        t2 = tables6[7][v2&63]; iadd_vbg(&r2, &t2); v2 >>= 6;
        t1 = tables6[8][v1&63]; iadd_vbg(&r1, &t1); v1 >>= 6;
        t2 = tables6[8][v2&63]; iadd_vbg(&r2, &t2); v2 >>= 6;
        t1 = tables5[0][v1&31]; iadd_vbg(&r1, &t1); v1 >>= 5;
        t2 = tables5[0][v2&31]; iadd_vbg(&r2, &t2); v2 >>= 5;
        t1 = tables5[1][v1&31]; iadd_vbg(&r1, &t1);
        t2 = tables5[1][v2&31]; iadd_vbg(&r2, &t2);
        
        isub_m3d(&r1, &r2);
        R[i][0] = r1;
       // */
    }
    
}

//32 * 64,2048 bit, 256 byte matrix(slice) multiplication
void mul_32_m3d(vbg *R, vbg *A, vbg *B)
{
    long i;
    vbg t1, t2, r1, r2, a;
    long v1, v2;
    
    vbg tables5[4][32];
    vbg tables4[3][16];
    for (i = 1; i < 4; i ++)
        
        m3d_combine5(tables5[i], B + 0 + 5*i);
    for (i = 0; i < 3; i++)
        m3d_combine4(tables4[i], B + 20 + 4*i);
    
    for (i = 0;i < 32; i++)
    {
        
        a = A[i];
        v2 = a.sign;
        v1 = a.units ^ v2;
        t1 = tables5[0][v1&31]; v1 >>= 5;
        t2 = tables5[0][v2&31]; v2 >>= 5;
        t1 = tables5[1][v1&31]; iadd_vbg(&r1, &t1); v1 >>= 5;
        t2 = tables5[1][v2&31]; iadd_vbg(&r2, &t2); v2 >>= 5;
        t1 = tables5[2][v1&31]; iadd_vbg(&r1, &t1); v1 >>= 5;
        t2 = tables5[2][v2&31]; iadd_vbg(&r2, &t2); v2 >>= 5;
        t1 = tables5[3][v1&31]; iadd_vbg(&r1, &t1); v1 >>= 5;
        t2 = tables5[3][v2&31]; iadd_vbg(&r2, &t2); v2 >>= 5;
        t1 = tables4[0][v1&15]; iadd_vbg(&r1, &t1); v1 >>= 4;
        t2 = tables4[0][v2&15]; iadd_vbg(&r2, &t2); v2 >>= 4;
        t1 = tables4[1][v1&15]; iadd_vbg(&r1, &t1); v1 >>= 4;
        t2 = tables4[1][v2&15]; iadd_vbg(&r2, &t2); v2 >>= 4;
        t1 = tables4[2][v1&15]; iadd_vbg(&r1, &t1);
        t2 = tables4[2][v2&15]; iadd_vbg(&r2, &t2);
        
        isub_m3d(&r1, &r2);
        R[i] = r1;
    }
    
}

//16 * 64,1024 bit, 128 byte matrix(slice) multiplication
void mul_16_m3d(vbg *R, vbg *A, vbg *B)
{
    long i;
    vbg t1, t2, r1, r2, a;
    long v1, v2;
    
    vbg tables4[4][16];
    for (i = 0; i < 4; i++)
        m3d_combine4(tables4[i], B + (4*i));
    for (i = 0;  i < 16; i++)
    {
        a = A[i];
        v2 = a.sign;
        v1 = a.units ^ v2;
        r1 = tables4[0][v1&15]; v1 >>= 4;
        r2 = tables4[0][v2&15]; v2 >>= 4;
        t1 = tables4[1][v1&15]; iadd_vbg(&r1, &t1); v1 >>= 4;
        t2 = tables4[1][v2&15]; iadd_vbg(&r2, &t2); v2 >>= 4;
        t1 = tables4[2][v1&15]; iadd_vbg(&r1, &t1); v1 >>= 4;
        t2 = tables4[2][v2&15]; iadd_vbg(&r2, &t2); v2 >>= 4;
        t1 = tables4[3][v1&15]; iadd_vbg(&r1, &t1);
        t2 = tables4[3][v2&15]; iadd_vbg(&r2, &t2);
    
        isub_m3d(&r1, &r2);
        R[i] = r1;
    }
}

//8 * 64,512 bit, m1ri_word byte matrix(slice) multiplication
void mul_8_m3d(vbg *R, vbg *A, vbg *B)

{
    int i;
    vbg t1, t2, r1, r2, a;
    vec v1, v2;
    
    vbg tables4[2][16];
    for (i = 0; i < 2; i++)
        m3d_combine4(tables4[i], B + (4*i));
    for (i = 0; i < 8; i++)
    {
        a = A[i];
    v2 = a.sign;
    v1 = a.units ^ v2;
    r1 = tables4[0][v1&15]; v1 >>= 4;
    r2 = tables4[0][v2&15]; v2 >>= 4;
    t1 = tables4[1][v1&15]; iadd_vbg(&r1, &t1);
    t2 = tables4[1][v2&15]; iadd_vbg(&r2, &t2);
    
    isub_m3d(&r1, &r2);
    R[i] = r1;
    }
}






//4 * 64,256 bit, 32 byte matrix(slice) multiplication
void mul_4_m3d(vbg *R, vbg *A, vbg *B)
{
    int i;
    vbg r1, r2, a;
    vec v1, v2;
    
    vbg table4[16];
    for (i = 0; i < 1; i++)
        m3d_combine4(table4, B + (4*i));
    for(i = 0; i < 4; i++)
    {
        a = A[i];
        v2 = a.sign;
        v1 = a.units ^ v2;
        r1 = table4[v1&15];
        r2 = table4[v2&15];
        
        isub_m3d(&r1, &r2);
        R[i] = r1;
    }
    
}





