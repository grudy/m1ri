//
// m1ri_combine.c
//
//// computions for small matrices over gf3
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
 */

#include <m1ri/m1ri_combine.h>


void *  combine3(vbg *table, vbg *input )
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


void combine4(vbg *table, vbg *input )
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


void combine5(vbg *table, vbg *input )
{
	int i;
    vbg e, *t4;
    
    
    
    combine4(table, input);
    e = input[4];
    t4 = table + 16;
    table[16] = e;
    
    for ( i = 1; i < 16 ; i ++ ) {
        add_vbg(t4 + i, table + i, &e);
    }
    
    
    
}


void combine6(vbg *table, vbg *input )

{
    vbg f, *t5;
    int i;
    combine5(table, input);
    f = input[5];
    t5 = (table + 32);
    table [32] = f;
    
    for (i = 1; i < 32; i++)
        add_vbg((t5 + i), (table + i), &f);
    
    
}

void combine7(vbg *table, vbg *input )

{
    
    vbg g, *t6;
    int i;
    
    combine6(table, input);
    g = input[6];
    t6 = (table+64);
    table[64] = g;
    
    for (i = 1; i < 64; i = i +1) {
        add_vbg((t6 + i), (table + i), &g );
    }
    
    
    
    
    
}


void combine8(vbg *table, vbg *input)

{
    vbg h, *t7;
    int i;
    
    combine7(table, input);
    h = input[7];
    t7 = (table+128);
    table[128] = h;
    
    for (i = 1; i < 128; i++)
        add_vbg((t7 + i), (table+i), &h);
}





void *  m5_combine3(vfd *table, vfd *input )
{
    vfd t, a, b, c;
    t.sign =  t.middle = t.units = 0;
    a = input[0]; 
    b = input[1];
    c = input[2];
    
    table[0] = t;
    table[1] = a;
    table[2] = b;
    table[4] = c;
    
    addgf5(&t, &a, &b);
    table[3] = t;
    iaddgf5(&t, &c);
    table[7] = t;
    isubgf5(&t, &a);
    table[6] = t;
    
    addgf5((table + 5), &a , &b);
    
    
    
    return 0;
    
}


void m5_combine4(vfd *table, vfd *input )
{
    vfd t, a, b, c , d;
    t.sign =  t.middle = t.units = 0;
    a = input[0];
    b = input[1];
    c = input[2];
    d = input[3];
    
    table[0] = t;
    table[1] = a;
    table[2] = b;
    table[4] = c;
    table[8] = d;
    
    addgf5(&t, &c, &d);
    
    table[12] = t;
    
    addgf5(&t,&b,&c);
    table[6] = t;
    iaddgf5(&t,&d);
    table[14] = t;
    isubgf5(&t,&c);
    table[10] = t;
    
    addgf5(&t,&b,&c);
    table[3] = t;
    iaddgf5(&t, &d);
    
    
    
    table[11] = t;
    iaddgf5(&t, &c);
    table[15] = t;
    isubgf5(&t, &d);
    table[7] = t;
    isubgf5(&t, &b);
    table[5] = t;
    iaddgf5(&t, &d);
    table[13] = t;
    isubgf5(&t, &c);
    table[9] = t;
    
    
}


void m5_combine5(vfd *table, vfd *input )
{
	int i;
    vfd e, *t4;
    
    
    
    m5_combine4(table, input);
    e = input[4];
    t4 = table + 16;
    table[16] = e;
    
    for ( i = 1; i < 16 ; i ++ ) {
        addgf5(t4 + i, table + i, &e);
    }
    
    
    
}


void m5_combine6(vfd *table, vfd *input )

{
    vfd f, *t5;
    int i;
    m5_combine5(table, input);
    f = input[5];
    t5 = (table + 32);
    table [32] = f;
    
    for (i = 1; i < 32; i++)
        addgf5((t5 + i), (table + i), &f);
    
    
}

void m5_combine7(vfd *table, vfd *input )

{
    
    vfd g, *t6;
    int i;
    
    m5_combine6(table, input);
    g = input[6];
    t6 = (table+64);
    table[64] = g;
    
    for (i = 1; i < 64; i = i +1) {
        addgf5((t6 + i), (table + i), &g );
    }
    
    
    
    
    
}


void m5_combine8(vfd *table, vfd *input)

{
    vfd h, *t7;
    int i;
    
    m5_combine7(table, input);
    h = input[7];
    t7 = (table+128);
    table[128] = h;
    
    for (i = 1; i < 128; i++)
        addgf5((t7 + i), (table+i), &h);
}


 void *   m7_combine3(vtri  *table, vtri  *input )
{
    vtri  t, a, b, c;
    t.sign =  t.middle = t.units = 0;
    a = input[0]; 
    b = input[1];
    c = input[2];
    
    table[0] = t;
    table[1] = a;
    table[2] = b;
    table[4] = c;
    
    addgf7(&t, &a, &b);
    table[3] = t;
    addgf7r(&t, &c);
    table[7] = t;
    isubgf7(&t, &a);
    table[6] = t;
    
    addgf7((table + 5), &a , &b);
    
    
    
    return 0;
    
}


void m7_combine4(vtri  *table, vtri  *input )
{
    vtri  t, a, b, c , d;
    t.sign =  t.middle = t.units = 0;
    a = input[0];
    b = input[1];
    c = input[2];
    d = input[3];
    
    table[0] = t;
    table[1] = a;
    table[2] = b;
    table[4] = c;
    table[8] = d;
    
    addgf7(&t, &c, &d);
    
    table[12] = t;
    
    addgf7(&t,&b,&c);
    table[6] = t;
    addgf7r(&t,&d);
    table[14] = t;
    isubgf7(&t,&c);
    table[10] = t;
    
    addgf7(&t,&b,&c);
    table[3] = t;
    addgf7r(&t, &d);
    
    
    
    table[11] = t;
    addgf7r(&t, &c);
    table[15] = t;
    isubgf7(&t, &d);
    table[7] = t;
    isubgf7(&t, &b);
    table[5] = t;
    addgf7r(&t, &d);
    table[13] = t;
    isubgf7(&t, &c);
    table[9] = t;
    
    
}


void m7_combine5(vtri  *table, vtri  *input )
{
	int i;
    vtri  e, *t4;
    
    
    
    m7_combine4(table, input);
    e = input[4];
    t4 = table + 16;
    table[16] = e;
    
    for ( i = 1; i < 16 ; i ++ ) {
        addgf7(t4 + i, table + i, &e);
    }
    
    
    
}


void m7_combine6(vtri  *table, vtri  *input )

{
    vtri  f, *t5;
    int i;
    m7_combine5(table, input);
    f = input[5];
    t5 = (table + 32);
    table [32] = f;
    
    for (i = 1; i < 32; i++)
        addgf7((t5 + i), (table + i), &f);
    
    
}

void m7_combine7(vtri  *table, vtri  *input )

{
    
    vtri  g, *t6;
    int i;
    
    m7_combine6(table, input);
    g = input[6];
    t6 = (table+64);
    table[64] = g;
    
    for (i = 1; i < 64; i = i +1) {
        addgf7((t6 + i), (table + i), &g );
    }
    
    
    
    
    
}

    

void m7_combine8(vtri *table, vtri *input)

{
    vtri h, *t7;
    int i;
    
    m7_combine7(table, input);
    h = input[7];
    t7 = (table+128);
    table[128] = h;
    
    for (i = 1; i < 128; i++)
        addgf7((t7 + i), (table+i), &h);
}




