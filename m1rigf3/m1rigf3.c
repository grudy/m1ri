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



#include <stdio.h>
#include<stddef.h>
#include<stdlib.h>
#include"m1rigf3.h"
#include"m1rigf3combine.h"
#define true 1
#define false 0
#define fn(a, b, c, d) (a^b)&(c^d) //for finding R[0]# (the first half of the value representingthe sum of vectory and vectorx, vectorr)
#define st(a, b , c) (a^b^c) //performing the (S= x[0] XOR y[1] XOR [x1]) and (T = x[1] XOR Y[0] XOR Y[1]) operations of addition








int main(int argc, const char * argv[])
{
    
    matrixgf3 r;
    matrixgf3 x;
    matrixgf3 y;


    
  
    //setting the values
    x.units.v = 5454535452452435;
y.sign.v = 42545454545353452;
    
    y.units.v = 245240352043592345;
    x.sign.v   =  45245234523452345;
    

    
    r = addgf3r(x, y);
    
    

    
    
    
    return 0;
}

