/** *  m1riproject
 Matrix Represenations and basic operations
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
 

 */
 
#ifndef M1RIPROJECT_M3D_DECOM_H
#define M1RIPROJECT_M3D_DECOM_h

#include <stdlib.h>
#include <m1ri/m1riwrappers.h>






/**
	\brief Solves L X = B with X and B matrices and L upper triangular
	\ X replaces B
*/
void m3d_upper_triangular(m3d_t * L, m3d_t B);
void m5d_upper_triangular(m5d_t * L, m5d_t B);
void m7d_upper_triangular(m7d_t * L, m7d_t * B);
 
 
 
 /**
   \brief Solves U X = B with X and B matrices and U upper left triangular, X replaces B\
   \param U matrix
   \param 
 */

void m3d_upper_left_triangular(m3d_t const *U, m3d_t *B  );
void m5d_upper_left_triangular(m5d_t const *U, m5d_t *B );
void m7d_upper_left_triangular(m7d_t const *U, m7d_t *B);

 
 /**
   Solves U X = B with X and B matrices and U lower triangular
   X replaces B
 */


void m3d_lower_triangular(m3d_t const *U, m3d_t *B);
void m5d_lower_triangular(m5d_t const *U, m5d_t *B);
void m7d_lower_triangular(m7d_t const *U, m7d_t *B);
 
  /**
   Solves U X = B with X and B matrices and U  lower right triangular
   X replaces B
 */
void m3d_lower_right_triangular(m3d_t const *U, m3d_t *B );
void m5d_lower_right_triangular(m5d_t const *U, m5d_t *B); 
void m7d_lower_right_triangular(m7d_t const *U, m7d_t *B);
 
 
 /**
   Invert triangular matrix a
 */
 void m3d_inverse_triangular(m3d_t * a);
 void m5d_inverse_triangular(m5d_t * a);
 void m7d_inverse_triangular(m7d_t * a);
 





#endif