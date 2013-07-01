//
//  combine.h
//  m1riproject
//
//  Created by grudy on 6/25/13.
//  Copyright (c) 2013 William Alumbaugh. All rights reserved.
//

#ifndef M1RIPROJECT_COMBINE_H
#define M1RIPROJECT_COMBINE_H

#include "m1riwrappers.h"
#include "m1ri_3dt.h"
#include "m1rielarith.h"

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
    
    addgf3(&t, &a, &b);
    table[3] = t;
    iaddgf3(&t, &c);
    table[7] = t;
    isubgf3(&t, &a);
    table[6] = t;
    
    addgf3((table + 5), &a , &b);
    
    
    
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
    
    addgf3(&t, &c, &d);
    
    table[12] = t;
    
    addgf3(&t,&b,&c);
    table[6] = t;
    iaddgf3(&t,&d);
    table[14] = t;
    isubgf3(&t,&c);
    table[10] = t;
    
    addgf3(&t,&b,&c);
    table[3] = t;
    iaddgf3(&t, &d);
    
    
    
    table[11] = t;
    iaddgf3(&t, &c);
    table[15] = t;
    isubgf3(&t, &d);
    table[7] = t;
    isubgf3(&t, &b);
    table[5] = t;
    iaddgf3(&t, &d);
    table[13] = t;
    isubgf3(&t, &c);
    table[9] = t;
    
    
}


void combine5(vbg *table, vbg *input )
{
    vbg e, *t4;
    
    
    
    combine4(table, input);
    e = input[4];
    t4 = table + 16;
    table[16] = e;
    
    for (int i = 1; i < 16 ; i ++ ) {
        addgf3(t4 + i, table + i, &e);
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
        addgf3((t5 + i), (table + i), &f);
    
    
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
        addgf3((t6 + i), (table + i), &g );
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
        addgf3((t7 + i), (table+i), &h);
}




#endif
