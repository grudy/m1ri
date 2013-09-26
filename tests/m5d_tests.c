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
#include <time.h>
int main(int argc, const char * argv[])
{
    int isequal;
    m5d_t * a, * b,  *d, *e, *f, *g, *h;
    a = m1ri_malloc(sizeof(m5d_t));
    b  = m1ri_malloc(sizeof(m5d_t));
    d = m1ri_malloc(sizeof(m5d_t));
     e = m1ri_malloc(sizeof(m5d_t));
     f =  m1ri_malloc(sizeof(m5d_t));
     g  = m1ri_malloc(sizeof(m5d_t));
     h = m1ri_malloc(sizeof(m5d_t));
    m5d_rand(a);
    
    m5d_rand(b);
    isequal = m5d_equal(a, b);
    if(isequal)
    {
        printf("Equaltest: passed ");
        
    }
    
    if(!isequal)
    {
        printf("Equaltest: failed ");
        return 1;
        
    }
    
    m5d_t test_m5d_output  = m5d_create( &test_m5d_output, 3   ,3);
    
    m5d_rand(&test_m5d_output);
    
    m5d_print(&test_m5d_output);
    
    m5d_write_elem(&test_m5d_output, 1, 1, 1, 1, 1);
    
    
    m5d_identity(d,64);
    m5d_identity(e,64);
    isequal = m5d_equal(d, e);
    
    if(isequal)
    {
        printf("Equaltest: passed ");
        
    }
    
    if(!isequal)
    {
        printf("Equaltest: failed ");
        return 1;
        
    }
    m5d_print(&test_m5d_output);
    
    
    m5d_create(f, 256, 256);
    m5d_create(g, 256, 256);
    printf("specs f:");
    m5d_specs(f); 
    printf("specs g:");
    m5d_specs(g); 
    m5d_rand(f);
    m5d_rand(g);
    // m5d_classic_mul(h, f, g);
    m5d_strassen(h, f, g);
    m5d_print(h);
	m5d_specs(h);    
    return 0;
}
