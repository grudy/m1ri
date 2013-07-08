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

 m1ri_test.c
 */ 

#include "m1ri_3dt.h"
#include "m1riwrappers.h"
#include "m1ri_3dt.h"
#include "m1riarith.h"
#include "m1ri_cubes.h"
#include "m1ri_small.h"
#include "m1ri_strassen.h"
#include "m1ri_combine.h"
#include "m1ri_classical.h"
#include "m1ri_io.h"
int main(int argc, const char * argv[])
{
 
    printf("Row Swap Test \n");
    printf("-------------  \n");
    m3d_t row_swap_test;
    row_swap_test =  m3d_create(&row_swap_test,8     , 8);
    m3d_identity_set(&row_swap_test);
    printf("Before the swap \n");
    printf("-------------  \n");    m3d_print(&row_swap_test);
    
    m3d_rowswap(&row_swap_test, 2, 4);
    printf("After the swap");
    printf("-------------  \n");    m3d_print(&row_swap_test);
      printf("Col Swap Test \n");
      printf("-------------  \n");
    m3d_t col_swap_test;
    col_swap_test =  m3d_create(&col_swap_test,6     , 6);
    m3d_identity_set(&col_swap_test);
    printf("Before the swap \n");
    printf("-------------  \n");
    m3d_print(&col_swap_test);
    
    m3d_rowswap(&col_swap_test, 2, 4);
    
    printf("After the swap");
   printf("-------------  \n");
    m3d_print(&col_swap_test);
    
    
    
    return 0;
}

