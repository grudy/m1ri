//Multiply Matrix slices
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
 
 
 smallm1ri.h
 */
#ifndef M1RIGF3_SMALLGF3_H
#define M1RIGF3_SMALLGF3_H
#include"m1ri_3dt.h"
#include "m1ri_combine.h"
#include "m1riarith.h"
#include "m1riwrappers.h"









//64 * 64,4096 bit, 512 byte matrix(slice) multiplication
void mul_64_m3d(vbg *, vbg *, vbg *);

//32 * 64,2048 bit, 256 byte matrix(slice) multiplication
void mul_32_m3d(vbg *, vbg *, vbg *);

//16 * 64,1024 bit, 128 byte matrix(slice) multiplication
void mul_16_m3d(vbg *, vbg *, vbg *);

//8 * 64,512 bit, 64 byte matrix(slice) multiplication
void mul_8_m3d(vbg *, vbg *, vbg *);





//4 * 64,256 bit, 32 byte matrix(slice) multiplication
void mul_4_m3d(vbg *R, vbg *A, vbg *B);



#endif