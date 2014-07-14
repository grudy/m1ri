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
 
 m1ri_io.h
*/
#ifndef M1RIPROJECT_M1RO_IO_H
#define M1RIPROJECT_M1RO_IO_H
#include <stdio.h>
#include <m1ri/m3d.h>
#include <m1ri/m7d.h>
#include <m1ri/m5d.h>
#include <m1ri/m1ri_permutations.h>



/**
 \  brief prints an m3d_t matrix
 \  param  a matrix to be printed
 \
 \	Prints an m3d_t to standard output 
 */
void m3d_print(const m3d_t * );
/**
 \  brief prints an m5d_t matrix
 \  param  a matrix to be printed
 \
 \	Prints an m5d_t to standard output 
 */
void m5d_print(const m5d_t *);



/**
 \  brief prints an m7d_t matrix
 \  param  a matrix to be printed
 \
 \	Prints an m7d_t to standard output 
 */
void m7d_print(const m7d_t * );


/**
	\ brief information on matrix
	\ param a matrix that will have its output shown
	\
	\Displays:
	\1.  Number of columns
	\2.  Number of rows
	\3.  Width in vbg's .
*/
void m3d_specs(m3d_t *);
/**
	\ brief same a running m3d_print and then m3d_specs
	\ param a input matrix
*/
void m3d_fullinfo(m3d_t *);


/**
	\ brief information on matrix
	\ param a matrix that will have its output shown
	\
	\Displays:
	\1.  Number of columns
	\2.  Number of rows
	\3.  Width in vfd's .
*/
void m5d_specs(m5d_t *);


/**
	\ brief same a running m5d_print and then m5d_specs
	\ param a input matrix
*/
void m5d_fullinfo(m5d_t *);


/**
	\ brief information on matrix
	\ param a matrix that will have its output shown
	\
	\Displays:
	\1.  Number of columns
	\2.  Number of rows
	\3.  Width in vtri's .
*/
void m7d_specs(m7d_t *);



 m3d_t m3d_read_textfile(const char * fn);
 m5d_t m5d_read_textfile(const char * fn);
 m7d_t m7d_read_textfile (const char * fn);


/**
 * Print  permutation matrices
 * 
 *
 */

void m3p_print(m3p_t const *);
void m5p_print(m3p_t const *);
void m7p_print(m3p_t const *);

#if __M1RI_HAVE_LIBPNG
int m3d_to_png(const m3d_t *, const char *, int , const char *, int );
#endif 




#endif
