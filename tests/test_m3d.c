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

 m1ri_functions
 */ 

#include <m1ri/m1ri.h>
#include "time.h"
int main(int argc, const char * argv[])
{
 
    m3d_t id_test;
  	id_test =  m3d_create(&id_test,6     , 6);
  	m3d_identity_set(&id_test);
    
   m3d_print(&id_test);
   m3d_t rand_test;
   rand_test =  m3d_create(&rand_test,    5  , 5);
   rand_test =  m3d_rand(&rand_test);
   printf("%d", id_test.width);

	m3d_print(&rand_test);
    
    m3d_print(&rand_test);
    printf("\n\n\n\n\n");
    printf("Write Test: ");
    /** 
    write test
    */
     m3d_t write_test = m3d_create(&write_test, 10  , 10);
    printf("Before write  \n");
    m3d_print(&write_test);
    m3d_write_elem(&write_test , 4, 4, 1, 0);
     
    printf("after write \n");
    
    m3d_print(&write_test);
    
    printf("\n\n\n\n\n");
    /** 
    
        Testing m3d_equal
     
   */
    m3d_t eqla = m3d_identity(&eqla, 30);
    m3d_t eqlb = m3d_identity(&eqla, 30);
    
    
    
    
    int equal_test = m3d_equal(&eqla, &eqlb);
    if(equal_test == 0)
        printf("Equal test failed");
    if (equal_test == 1)
        printf("Equal test was successful");
    
    
    
    printf("Row Swap Test \n");
    printf("-------------  \n");
    m3d_t row_swap_test;
    row_swap_test =  m3d_create(&row_swap_test,7     , 7);
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
  
    m3d_t *a = m1ri_malloc(sizeof(m3d_t));
    
    m3d_t *b = m1ri_malloc(sizeof(m3d_t));
    m3d_create(a, 64, 64);
    m3d_create(b, 64, 64);
    m3d_rand(a);
    m3d_rand(b);
    
    
    printf("Addition Test \n");
    printf("Matrix a: \n");
    m3d_print(a);
    
    
    printf("Matrix b: \n");
    m3d_print(b);
    
    
    m3d_t c = m3d_add(a,b );
    
    m3d_print(&c);
    
    
    
    return 0;
}

