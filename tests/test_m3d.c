/** * M1RI
TOMAS J. 
 AND ROBERT W. BRADSHAW "BITSLICING AND THE METHOD OF FOUR
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
  /*
   [a][b]
   [c][d]
   
  */
  m3d_t * o = malloc(sizeof(m3d_t));
  m3d_t * ab = malloc(sizeof(m3d_t));
  
  m3d_t * cd = malloc(sizeof(m3d_t));
  m3d_t * abcd = malloc(sizeof(m3d_t));
  m3_slice * s = malloc(sizeof(m3_slice));
    
 
  m3d_create(o, 128, 128);
   m3d_rand(o);
//  


  m3d_quarter(s, o);
  m3d_concat(ab, &s->row[0][0], &s->row[1][0]);
  m3d_concat(cd, &s->row[0][1], &s->row[1][1]);

  m3d_specs(&s->row[0][0]);
  m3d_specs(o);
  if(m3d_equal(o, abcd))
  {
     printf("\nfail\n");
  
  }  
  m3d_free(o);
  m3d_free(ab);
  m3d_free(cd);
  m3d_free(abcd);
  m1ri_free(s);  
    return 0;
}

