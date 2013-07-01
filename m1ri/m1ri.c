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

 m1ri.c
 */ 



#include<stdlib.h>
#include "stdio.h" 
#include "m1ri_3dt.h"
#include "m1riwrappers.h"
#include "m1rielarith.h"
#include "m1ristrassen.h"
#include "m1ri_small.h"
#include "m1ri_classical.h"
#include "m7d.h"
#include "m5d.h"
#include "m1ri_cubes.h"
#include "m1ri_combine.h"
#include <time.h>





int main(int argc, const char * argv[])
{
   srand((unsigned int)time(0));
    m3d_t a = m3d_identity(&a, 256);
    //a.block =  m3d_rand(&a);
    vbg    x;
    x.sign =  ibits;
    x.units = leftbit;
  
    a =  m3d_identity_set(&a);

  
  
 
    printf("Signbit: %lld \n\n", x.sign);
    printf("Unitbit: %lld \n\n", x.units);

    
    
    return 0;
}

