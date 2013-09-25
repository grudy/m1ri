/* M1RI
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
 
 m3d_arith_test.c
 */ 





#include <m1ri/m1ri.h>
#include "time.h"
int main(int argc, const char * argv[])
{
 	 m3d_t *a, *b, *c, *d, *e, *f;
    a = malloc(sizeof(m3d_t));
    b = malloc(sizeof(m3d_t));
    c =  malloc(sizeof(m3d_t));
    d = malloc(sizeof(m3d_t));
    e = malloc(sizeof(m3d_t));
    f = malloc(sizeof(m3d_t));
 	m3d_create(a, 256, 256);
 	m3d_create(b, 256, 256);
    m3d_create(d, 256, 256);
 	m3d_create(e, 256, 256);
 
    m3d_rand(a);
    m3d_rand(b);
    m3d_rand(d);
    m3d_rand(e);
	m3d_classic_mul(c, a, b);
	//m3d_strassen(f, d, e);
    m3d_print(c);
    return 0;
    
    
}

