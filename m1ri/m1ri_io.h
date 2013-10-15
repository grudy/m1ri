/** * M1RI
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
 
 m1ri_io.h
*/
#ifndef M1RIPROJECT_M1RO_IO_H
#define M1RIPROJECT_M1RO_IO_H
#include <stdio.h>
#include <m1ri/m1ri_3dt.h>
#include <m1ri/m7d.h>
#include "m5d.h"

/** 
	Prints an m3d_t matrix
*/
void m3d_print(m3d_t * restrict);

/** 
	Prints an m5d_t matrix
*/
void m5d_print(m5d_t *);
/** 
	Prints an m7d_t matrix
*/
void m7d_print(m7d_t * );

void m3d_specs(m3d_t *);

void m3d_fullinfo(m3d_t *);

void m5d_specs(m5d_t *);

void m5d_fullinfo(m5d_t *);

void m7d_specs(m7d_t *);

void m7d_info(m7d_t *);


#endif
