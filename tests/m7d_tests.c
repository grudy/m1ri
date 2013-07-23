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

 m7d_tests.c
 */ 

#include "m1ri/m1riwrappers.h"
#include "m1ri/m1ri/m1ri_3dt.h"
#include "m1ri/m1riarith.h"
#include "m1ri/m1ri_cubes.h"
#include "m1ri/m1ri_small.h"
#include "m1ri/m1ri_strassen.h"
#include "m1ri/m1ri_combine.h"
#include "m1ri/m1ri_classical.h"
#include "m1ri/m1ri_io.h"
#include "time.h"
int main(int argc, const char * argv[])
{
 
  
    
    m7d_t test_m7d_output  = m7d_create( &test_m7d_output, 3   ,3);
  
    m7d_rand(&test_m7d_output);
    
    m7d_print(&test_m7d_output);
    
    m7d_write_elem(&test_m7d_output, 1, 1, 1, 0, 1);
    
    
   
    
    m7d_print(&test_m7d_output);
    
    
    
    return 0;
}

