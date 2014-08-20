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

 m3d_mulbenchmark.c
 */

#include <m1ri/m1ri.h>
#include "time.h"



void m3d_strassen_test(int y, int z)
{
	m3d_t * a, *b, *c;
	a = m3d_create(y, z);
    b = m3d_create(y, z);

    m3d_rand(a);
    m3d_rand(b);




	clock_t begin, end;
	double time_spent;
	begin = clock();
	c = m3d_strassen(c, a, b);
    time_t after;
    time(&after);
	end = clock();
	time_spent = (double)(end - begin) / CLOCKS_PER_SEC;

    printf("----------------------------------------------------------------------");
    printf("\nm3d_strassen on two %d by %d matrix matrices.", y, z  );
    printf(" \n------------------->Runs in%9f seconds. \n", time_spent);
    printf("----------------------------------------------------------------------");
    m3d_free(a);
    m3d_free(b);
    m3d_free(c);



}


int main(int argc, const char * argv[])
{
 	//m3d_strassen_test(64, 64);
   	m3d_strassen_test(256, 256);
    m3d_strassen_test(4196, 4196);
   // m3d_strassen_test(64000, 64000);


   return 0 ;
}
