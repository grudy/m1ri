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


#include "m1ri/m1riwrappers.h"
#include "m1ri/m1ri/m1ri_3dt.h"
#include "m1ri/m1riarith.h"

#include "m1ri/m1ri_io.h"
#include "time.h"
int main(int argc, const char * argv[])
{
    time_t before;
    time(&before);
    
    
    m3d_t bunches[30];
    m3d_t plenitude[30];
    m3d_t  oodles[30];
    for (int x = 0; x < 30; x ++)
    {
        
        
        
        plenitude[x] = m3d_create(&bunches[x], 16384, 16384);
        bunches[x]  = m3d_create(&plenitude[x],  16384, 16384);
        oodles[x]  = m3d_create(&oodles[x], 16384, 16384);
        plenitude[x]  = m3d_rand(&plenitude[x]);
        bunches[x]  = m3d_rand(&bunches[x]);
        plenitude[x]  = m3d_rand(&plenitude[x]);
        plenitude[x] =  m3d_hadamard(&plenitude[x], &plenitude[x]);
        
        
    }
    
    time_t after;
    time(&after);
    double time_test_m1ri = difftime( after, before);
    
    printf("%9f", time_test_m1ri );
    
    return 0;
}
