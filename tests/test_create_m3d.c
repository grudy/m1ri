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
 
 test_create_m3d.c
 */ 





#include <m1ri/m1ri.h>

int main(int argc, const char * argv[])
{
 	
 
    m3d_t * a;
    a = m1ri_malloc(sizeof(m3d_t)); 
	m3d_create(a, 64, 64);
    m3d_rand(a);
    //m3d_print(a);
    int compression = 1;
    char * name = "m3d_testmatrix.png";
    char * comment = "comment";
    m3d_to_png( a,  name,  compression,  comment ,  0);
    //printf("d");
    //m3d_free(a);
    //#if __M1RI_HAVE_LIBPNG
    //printf("hi")
    //#endif 

     return 0;
    
    
}

