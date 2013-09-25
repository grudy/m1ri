
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
 
 m1ri_strassen.h
 */


#ifndef M1RIGF3_STRASSEN_H
#define M1RIGF3_STRASSEN_H
#include <m1ri/m1ri_3dt.h>
#include <m1ri/m1riarith.h>
#include <m1ri/m1ri_small.h>
#include <m1ri/m1ri_cubes.h>
#include <m1ri/m5d.h>
#include <m1ri/m7d.h>

void  m3d_strassen(m3d_t *, m3d_t *, m3d_t*);
void m3d_qrt_mul(m3d_t * ,m3d_t *, m3d_t *  );
void m5d_strassen(m5d_t * ,m5d_t *, m5d_t *);
void m5d_qrt_mul(m5d_t * ,m5d_t *, m5d_t *  );
void m7d_strassen(m7d_t * ,m7d_t *, m7d_t *);
void m7d_qrt_mul(m7d_t * ,m7d_t *, m7d_t *  );

#endif

