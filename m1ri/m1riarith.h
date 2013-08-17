
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
 
 m1riarith.h
 */

#ifndef M1RIPROJECT_M1RIELARITH_H
#define M1RIPROJECT_M1RIELARITH_H

#include <m1ri/m1ri_3dt.h>




void add_vbg(vbg * r, vbg const * x, vbg const * y);

void vbg_negation(vbg * );

vbg add_m3dr(vbg  , vbg );

void sub_m3d( vbg *, vbg const *, vbg const * );               //multiply matrix x by by matrix y.   The product is matrix r.




vbg sub_m3dr(vbg , vbg );               //multiply matrix x by by matrix y.   The product is matrix r.




/********************************************
 matrix r = (direct sum matrix r + matrix x)
 ********************************************/
void iadd_vbg(vbg *,vbg  *);

void isub_m3d(vbg *,vbg  *);



void  vbg_mul( vbg *, vbg  *, vbg  *);             //multiply matrix x by y assinging the output to r


void m3d_sub( m3d_t *, m3d_t const *, m3d_t const *);

//return the value of the matrix multiplied


vbg vbg_mul_i(vbg const , vbg const);
void m3d_add_r(m3d_t * , m3d_t  *, m3d_t  *);
/*
    Hadamard multiplication
*/
m3d_t * m3d_hadamard(m3d_t const *, m3d_t const *);



/* * * * * * * * * * * * * * * * * * * *
 Subtract a 1 kilobyte Matrix from another
 1 kilobyte Matrix
 * * * * * * * * * * * * * * * * * * * * */



void sub_64_m3d(vbg **, vbg  **, vbg  **);


/* * * * * * * * * * * * * * * * * * * * * *
 Add a 1 kilobyte Matrix from another
 1 kilobyte Matrix
 * * * * * * * * * * * * * * * * * * * * * * */

void add_64_m3d(vbg **, vbg   **, vbg    ** );


m3d_t m3d_add(m3d_t  *, m3d_t  *);




#endif
