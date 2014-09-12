
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
 
 m1riarith.h
 */

#ifndef M1RIPROJECT_M1RIELARITH_H
#define M1RIPROJECT_M1RIELARITH_H

#include <m1ri/m3d.h>
#include <m5d.h>
#include <m7d.h>


/**
  Negates the  input vbg
*/
void vbg_negation(vbg * );


/**
	
*/
void sub_m3d( vbg *, vbg const *  , vbg const * );               /* multiply matrix x by by matrix y.   The product is matrix r. */


vbg sub_m3dr(vbg , vbg );               /*  */


/** *******************************************
 matrix r = (direct sum matrix r + matrix x)
 ********************************************/
/* void iadd_vbg(vbg *,vbg  *); */

/* void m3d_dec(vbg *,vbg  *); */



void  vbg_mul( vbg *, vbg  *, vbg  *);             /* multiply matrix x by y assinging the output to r */

/**
	m3d_t subtraction
	
*/
m3d_t * m3d_sub( m3d_t const *,const  m3d_t  *);

/**
	Return the value of the matrix multiplied
*/
vbg vbg_mul_i(vbg const , vbg const);

/** 
    Hadamard multiplication
*/
m3d_t * m3d_hadamard(m3d_t const * , m3d_t const * );









/** * * * * * * * * * * * * * * * * * * * *
 Subtract a 1 kilobyte Matrix from another
 1 kilobyte Matrixhgg
 * * * * * * * * * * * * * * * * * * * * */

void m3d_sub_64(vbg ** , vbg  ** , vbg  ** );

/** * * * * * * * * * * * * * * * * * * * * * *
 Add a 1 kilobyte Matrix from another
 1 kilobyte Matrix
 * * * * * * * * * * * * * * * * * * * * * * */
void m3d_add_64(vbg **, vbg   **  , vbg    **  );

m3d_t  * m3d_add(m3d_t  *, m3d_t  *);


/*
void * m3d_combine3(vbg *, vbg * );

void m3d_combine4(vbg *, vbg * );

void m3d_combine5(vbg *, vbg * );

void m3d_combine6(vbg *, vbg * );

void m3d_combine7(vbg *, vbg * );

void m3d_combine8(vbg *, vbg *);
*/
/** ***************************************************************************
								GF(3)
*****************************************************************************/
/* 64 * 64,4096 bit, 512 byte matrix(slice) multiplication */
void m3d_mul_64(vbg **, vbg ** , vbg ** );

/* 32 * 64,2048 bit, 256 byte matrix(slice) multiplication */
void mul_32_m3d(vbg *, vbg *, vbg *);

/* 16 * 64,1024 bit, 128 byte matrix(slice) multiplication */
void mul_16_m3d(vbg *, vbg *, vbg *);

/* 8 * 64,512 bit, m1ri_word byte matrix(slice) multiplication */
void mul_8_m3d(vbg *, vbg *, vbg *);

/* 4 * 64,256 bit, 32 byte matrix(slice) multiplication */
void mul_4_m3d(vbg *R, vbg *A, vbg *B);
#endif
