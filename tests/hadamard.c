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
 
 m1ri_hadamard.c
 */


#include <m1ri/m1ri.h>

#include "time.h"
int main(int argc, const char * argv[])
{
    time_t before;
    time(&before);
    int isequal;
    
    m3d_t bunches;
    m3d_t plenitude;
    m3d_t  oodles;
    
    
    
    
    plenitude = m3d_create(&bunches, 16384, 16384);
    bunches  = m3d_create(&plenitude,  16384, 16384);
    oodles  = m3d_create(&oodles, 16384, 16384);
    plenitude  = m3d_rand(&plenitude);
    bunches  = m3d_rand(&bunches);
  //  oodles =  *m3d_hadamard(&plenitude, &bunches);
    
    
    
    
    time_t after;
    time(&after);
    double time_test_m1ri = difftime( after, before);
    
    printf("Time: %9f seconds \n", time_test_m1ri );
    
    isequal = m3d_equal(&oodles, &plenitude);
    
    if(!isequal)
    {
        printf("Hadamard: passed ");
        
    }
    
    if(isequal)
    {
        printf("Hadamard: failed ");
        return 1;
        
    }
    
    return 0;
}


