/** *  m1riproject
 Matrix Represenations and basic operations
 TOMAS J. BOOTHBY AND ROBERT W. BRADSHAW "BITSLICING AND THE METHOD OF FOUR
 RUSSIANS OVER LARGER FINITE FIELDS"
 Rank-profile revealing Gaussian elimination and the CUP matrix decomposition
 Claude-Pierre Jeannerod, Clément Pernet, Arne Storjohann
 
 http://arxiv.org/abs/1112.5717
 
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
 
#ifndef M1RIPROJECT_DECOM_H
#define M1RIPROJECT_DECOM_h

#include <stdlib.h>
#include <m1ri/m1riwrappers.h>
#include <m1ri/m3d.h>
#include <m1ri/m5d.h>
#include <m1ri/m7d.h>
#include <m1ri/m1ri_io.h>




void  m3d_transpose_64( vbg **, vbg ** );
m3d_t * m3d_transpose(m3d_t *, m3d_t const  *);

void  m5d_transpose_64( vfd **, vfd ** );
m5d_t * m5d_transpose(m5d_t *, m5d_t const  *);

void  m7d_transpose_64( vtri **, vtri ** );
m7d_t * m7d_transpose(m7d_t *, m7d_t const  *);
/**
	\brief Solves L X = B with X and B matrices and L upper right triangular
	\ X replaces B
*/
void m3d_tsrm_ur(m3d_t * L, m3d_t B);


void m5d_tsrm_ur(m5d_t * L, m5d_t B);


void m7d_tsrm_ur(m7d_t * L, m7d_t * B);
 
 
 
 /**
   \brief Solves U X = B with X and B matrices and U upper left triangular, X replaces B\
   \param U matrix
   \param 
 */

void m3d_tsrm_ul(m3d_t const *U, m3d_t *B  );


void m5d_tsrm_ul(m5d_t const *U, m5d_t *B );


void m7d_tsrm_ul(m7d_t const *U, m7d_t *B);

 
 /**
   Solves U X = B with X and B matrices and U lower left triangular
   X replaces B
 */


void m3d_tsrm_ll(m3d_t const *U, m3d_t *B);


void m5d_tsrm_ll(m5d_t const *U, m5d_t *B);


void m7d_tsrm_ll(m7d_t const *U, m7d_t *B);
 
  /**
   Solves U X = B with X and B matrices and U  lower right triangular
   X replaces B
 */

void m3d_tsrm_lr(m3d_t const *U, m3d_t *B );


void m5d_lower_r_triangular(m5d_t const *U, m5d_t *B); 


void m7d_tsrm_lr(m7d_t const *U, m7d_t *B);
 
 
 /**
   \brief Invert upper triangular matrix a
 */
void m3d_trtri_upper(m3d_t * a);


void m5d_trtri_upper(m5d_t * a);




void m7d_trtri_upper(m7d_t * a);
 
 






#endif