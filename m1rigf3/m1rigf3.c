//// GF3 arithmatic
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
 */
// m1rigf3.c



#include<stdlib.h>
#include "stdio.h" 
#include "m1rigf3.h"
#include "m1riwrappers.h"
#include "m1rielarith.h"
#define true 1
#define false 0
#define fn(a, b, c, d) (a^b)&(c^d) //for finding R[0]# (the first half of the value representingthe sum of vectory and vectorx, vectorr)
#define st(a, b , c) (a^b^c) //performing the (S= x[0] XOR y[1] XOR [x1]) and (T = x[1] XOR Y[0] XOR Y[1]) operations of addition









int main(int argc, const char * argv[])
{
   srand((unsigned int)time(0));
    m3d_t a;
   m3d_create(&a, 300  , 300);
    a.block =  m3d_rand(&a);
     vbg x = m3d_read_elems(&a, (a.ncols - 3)    , (a.nrows - 44) , 1);


    
    
    printf("Sign: %lld \n\n", x.sign);
    printf("Units: %lld \n\n", x.units);
    m3d_write_elem(&a, (a.ncols -3), (a.nrows - 44), 1, 0);
    
    
    x = m3d_read_elems(&a, (a.ncols - 3)    , (a.nrows - 44) , 1);
    printf("Signbit: %lld \n\n", x.sign);
    printf("Unitbit: %lld \n\n", x.units);
    
    
    return 0;
}

