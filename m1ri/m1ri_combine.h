
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
 
 m1ri_combine.h
 */
#ifndef M1RIPROJECT_COMBINE_H
#define M1RIPROJECT_COMBINE_H

#include <m1ri/m1riwrappers.h>
#include <m1ri/m1ri_3dt.h>
#include <m1ri/m1riarith.h>
#include <m1ri/m5d.h>
#include <m1ri/m7d.h>
void *  combine3(vbg *, vbg * );


void combine4(vbg *, vbg * );

void combine5(vbg *, vbg * );
void combine6(vbg *, vbg * );

void combine7(vbg *, vbg * );


void combine8(vbg *, vbg *);

void *  m5_combine3(vfd *, vfd * );


void m5_combine4(vfd *, vfd * );

void m5_combine5(vfd *, vfd * );
void m5_combine6(vfd *, vfd * );

void m5_combine7(vfd *, vfd * );


void m5_combine8(vfd *, vfd *);

void *  m7_combine3(vtri  *, vtri  * );


void m7_combine4(vtri  *, vtri  * );

void m7_combine5(vtri  *, vtri  * );
void m7_combine6(vtri  *, vtri  * );

void m7_combine7(vtri  *, vtri  * );


void m7_combine8(vtri  *, vtri  *);



#endif
