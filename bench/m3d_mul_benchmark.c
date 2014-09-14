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



void m3d_classic_mul_test(int y, int z)
{
	m3d_t * a, *b, *c;
	a = m3d_create(y, z);
    b = m3d_create(y, z);

    m3d_rand(a);
    m3d_rand(b);




	clock_t begin, end;
	double time_spent;
	begin = clock();
	c = m3d_classic_mul(c, a, b);
    time_t after;
    time(&after);
	end = clock();
	time_spent = (double)(end - begin) / CLOCKS_PER_SEC;

    printf("----------------------------------------------------------------------");
    printf("\nm3d_classic_mul on two %d by %d matrix matrices.", y, z  );
    printf(" \n------------------->Runs in %9f seconds. \n", time_spent);
    printf("----------------------------------------------------------------------");

    m3d_free(a);
    m3d_free(b);
    m3d_free(c);

	

}





void m5d_strassen_test(int y, int z)
{
	m5d_t * a, *b, *c;
	a = m5d_create(y, z);
    b = m5d_create(y, z);

    m5d_rand(a);
    m5d_rand(b);




	clock_t begin, end;
	double time_spent;
	begin = clock();
	c = m5d_strassen(c, a, b);
    time_t after;
    time(&after);
	end = clock();
	time_spent = (double)(end - begin) / CLOCKS_PER_SEC;

    printf("----------------------------------------------------------------------");
    printf("\nm5d_strassen on two %d by %d matrix matrices.", y, z  );
    printf(" \n------------------->Runs in%9f seconds. \n", time_spent);
    printf("----------------------------------------------------------------------");
    m5d_free(a);
    m5d_free(b);
    m5d_free(c);

	

}



void m5d_classic_mul_test(int y, int z)
{
	m5d_t * a, *b, *c;
	a = m5d_create(y, z);
    b = m5d_create(y, z);

    m5d_rand(a);
    m5d_rand(b);




	clock_t begin, end;
	double time_spent;
	begin = clock();
	c = m5d_classic_mul(c, a, b);
    time_t after;
    time(&after);
	end = clock();
	time_spent = (double)(end - begin) / CLOCKS_PER_SEC;

    printf("----------------------------------------------------------------------");
    printf("\nm5d_classic_mul on two %d by %d matrix matrices.", y, z  );
    printf(" \n------------------->Runs in %9f seconds. \n", time_spent);
    printf("----------------------------------------------------------------------");

    m5d_free(a);
    m5d_free(b);
    m5d_free(c);

	

}



void m7d_strassen_test(int y, int z)
{
	m7d_t * a, *b, *c;
	a = m7d_create(y, z);
    b = m7d_create(y, z);

    m7d_rand(a);
    m7d_rand(b);




	clock_t begin, end;
	double time_spent;
	begin = clock();
	c = m7d_strassen(c, a, b);
    time_t after;
    time(&after);
	end = clock();
	time_spent = (double)(end - begin) / CLOCKS_PER_SEC;

    printf("----------------------------------------------------------------------");
    printf("\nm7d_strassen on two %d by %d matrix matrices.", y, z  );
    printf(" \n------------------->Runs in%9f seconds. \n", time_spent);
    printf("----------------------------------------------------------------------");
    m7d_free(a);
    m7d_free(b);
    m7d_free(c);

	

}



void m7d_classic_mul_test(int y, int z)
{
	m7d_t * a, *b, *c;
	a = m7d_create(y, z);
    b = m7d_create(y, z);

    m7d_rand(a);
    m7d_rand(b);




	clock_t begin, end;
	double time_spent;
	begin = clock();
	c = m7d_classic_mul(c, a, b);
    time_t after;
    time(&after);
	end = clock();
	time_spent = (double)(end - begin) / CLOCKS_PER_SEC;

    printf("----------------------------------------------------------------------");
    printf("\nm7d_classic_mul on two %d by %d matrix matrices.", y, z  );
    printf(" \n------------------->Runs in %9f seconds. \n", time_spent);
    printf("----------------------------------------------------------------------");

    m7d_free(a);
    m7d_free(b);
    m7d_free(c);

	

}


int main(int argc, const char * argv[])
{

	/*
 	m3d_strassen_test(64, 64);
   	m3d_strassen_test(256, 256);
   	m3d_strassen_test(512, 512);
	m3d_strassen_test(1024, 1024);
	m3d_strassen_test(2048, 2048);
    m3d_strassen_test(4096, 4096);
    m3d_strassen_test(8192, 8192);
	m3d_strassen_test(16384, 16384);
	
	*/
 	m3d_classic_mul_test(64, 64);
   	//m3d_classic_mul_test(256, 256);
   //	m3d_classic_mul_test(512, 512);
	//m3d_classic_mul_test(1024, 1024);
	//m3d_classic_mul_test(2048, 2048);
    m3d_classic_mul_test(4096, 4096);
   // m3d_classic_mul_test(8192, 8192);
	//m3d_classic_mul_test(16384, 16384);     


m5d_strassen_test(4096, 4096);
	
/**


		what could be m5d_mul_benchmark.c
*/
	/*
	m5d_strassen_test(64, 64);
   	m5d_strassen_test(256, 256);
   	m5d_strassen_test(512, 512);
	m5d_strassen_test(1024, 1024);
	m5d_strassen_test(2048, 2048);
		m5d_strassen_test(4096, 4096);
    
    m5d_strassen_test(8192, 8192);
	m5d_strassen_test(16384, 16384);
	
	
 	m7d_classic_mul_test(64, 64);
   	m7d_classic_mul_test(256, 256);
   	m7d_classic_mul_test(512, 512);
	m7d_classic_mul_test(1024, 1024);
	m7d_classic_mul_test(2048, 2048);
	*/
    m7d_classic_mul_test(4096, 4096);
    /*
    m7d_classic_mul_test(8192, 8192);
	m7d_classic_mul_test(16384, 16384);     

    */


			m7d_strassen_test(4096, 4096);

	


   return 0 ;
}
