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
/ along with this program; if not, write to the Free Software
 Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 
 m1ri_hadamard.c
 */


#include "m1riwrappers.h"
#include "m1ri_3dt.h"
#include "m1riarith.h"
#include "m1ri_cubes.h"
#include "m1ri_io.h"
#include "time.h"
int main(int argc, const char * argv[])
{
    
    m3d_t * a = m1ri_malloc(sizeof(m3d_t));
    m3d_create(a, 128, 64);
    m3d_rand(a);
    m3d_print(a);
    m3d_t b;
    
    //m3d_t * b1 = m1ri_malloc(sizeof(m3d_t));
   // m3d_t *b2 = m1ri_malloc(sizeof(m3d_t));
    //m3d_t *b3 = m1ri_malloc(sizeof(m3d_t));
    //m3d_t *b4 = m1ri_malloc(sizeof(m3d_t));
  //  m3d_window_create(a, b1, 0, 0, 1, 1);
   // m3d_window_create(a, b2, 0, 1, 1, 1);
    //m3d_window_create(a, b3, 1, 0, 1, 1);
    //m3d_window_create(a, b4, 1, 1, 1, 1);
    
   m3d_transpose(a, &b);
   // m3d_identity_set(b2);
   // m3d_transposewin(b2);
    //m3d_transposewin(b3);
    m3d_print(&b);
   // int f[x];
    
    free(a);
    
    return 0;
}
