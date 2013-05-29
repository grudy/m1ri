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
//  m1rigf3io.h
//  m1rigf3
//
//  Created by grudy on 5/23/13.
//  Copyright (c) 2013 William Alumbaugh. All rights reserved.
//

#ifndef m1rigf3_m1rigf3io_h
#define m1rigf3_m1rigf3io_h

 
#include <stdlib.h>
#include"m1rigf3.h"
matrixgf3 * print_m1ri_128(matrixgf3 *M)
{
    vec i, j, ptr;
    matrixgf3 cur
    print '+' + '-'*128 + '+';
    for i from 0 <= i < 128:
        s = '|'
        ptr = 1;
    cur = M[i];
    for j from 0 <= j < 64:
        if cur.units&ptr:
            if cur.sign&ptr:
                s += ':'
                else:
                    s += '.'
                    elif cur.sign & ptr:
                    s += '?'
                    else:
                        s += ' '
                        ptr += ptr
                        ptr = 1
                        cur = M[i+128]
                        for j from 0 <= j < 64:
                            if cur.units&ptr:
                                if cur.sign&ptr:
                                    s += ':'
                                    else:
                                        s += '.'
                                        elif cur.sign & ptr:
                                        s += '?'
                                        else:
                                            s += ' '
                                            ptr += ptr
                    print s + '|'
                                            print '+' + '-'*128 + '+'
                                            }


matrixgf3 * print_m1ri(matrixgf3 *M)
{
    long i, j, ptr
    global D
    for i from 0 <= i < M1RI_BITS:
        ptr = 1
        s = ''
        for j from 0 <= j < M1RI_BITS:
            if M[i].units&ptr:
                if M[i].sign&ptr:
                    s += ':'
                    else:
                        s += '.'
                        elif M[i].sign & ptr:
                        s += '?'
                        else:
                            s += ' '
                            ptr += ptr
                            print s
                            
                            matrixgf3 * matrix_to_m1ri(M):
                            matrixgf3 t
                            matrixgf3 *R = <matrixgf3 *>malloc(M1RI_BITS * sizeof(matrixgf3))
                            long ptr, i
                            
                            i = 0
                            for row in M:
                                t.units = 0
                                t.sign = 0
                                ptr = 1
                                j = 0
                                for x in row:
                                    if x != 0:
                                        t.units |= ptr
                                        if x == 2:
                                            t.sign |= ptr
                                            ptr += ptr
                                            R[i] = t
                                            i+=1
                                            return R
                                            }


matrixgf3 * matrix_to_m1ri_128(M)
{
    matrixgf3 t;
    matrixgf3 *R = <matrixgf3 *>malloc(256 *sizeof(matrixgf3));
    vec ptr, i;
    
    i = 0;
    for row in M:
        t.units = 0
        t.sign = 0
        ptr = 1
        for j from 0 <= j < 64:
            x = row[j];
    if (x != 0)
        t.units |= ptr;
    if (x == 2)
        t.sign |= ptr;
    ptr += ptr;
    R[i] = t;
    t.units = 0
    t.sign = 0
    ptr = 1
    for j from 64 <= j < 128:
        x = row[j];
    if (x != 0)
        t.units |= ptr;
    if x == 2:;
    t.sign |= ptr;
    ptr += ptr;
    R[i+128] = t;
    
    i+=1;
    return R
    
}


#endif


#endif
