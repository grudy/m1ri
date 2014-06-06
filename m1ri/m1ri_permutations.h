 
/** *
 
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
 
 
 Matrix Represenations and basic operations over GF(3)
 m1ri_permutations.h
 */

#ifndef M1RI_PERMUTATIONS
#define M1RI_PERMUTATIONS

#include <m1ri/m1riwrappers.h>
#include <m1ri/m1ri_3dt.h>
#include <m1ri/m5d.h>
#include <m1ri/m7d.h>




/**
    Permutations.
 */
//Structures
typedef struct m3p_t {
  /**
   *  LAPACK format.
   */
  rci_t *values;

  /**
   * The length of the swap array.
   */

  rci_t length;

} m3p_t;

typedef struct m5p_t {
  /**
   * LAPACK format.
   */
  rci_t *values;

  /**
   * The length of the swap array.
   */

  rci_t length;

} m5p_t;
typedef struct m7p_t {
  /**
   * The swap operations in LAPACK format.
   */
  rci_t *values;

  /**
   * The length of the swap array.
   */

  rci_t length;

} m7p_t;

/**
 * Construct an identity permutation.
 *
 * \param length Length of the permutation.
 */

m3p_t *m3p_init(rci_t length);

m5p_t *m5p_init(rci_t length);
m7p_t *m7p_init(rci_t length);

/**
 * Free a permutation.
 *
 * \param P Permutation to free.
 */

void m3p_free(m3p_t *);
void m5p_free(m5p_t *);
void m7p_free(m7p_t *);

/**
 * \brief Create a window/view into the permutation P.
 * \param P Permutation matrix
 * \param begin Starting index (inclusive)
 * \param end   Ending index   (exclusive)
 *
 */

m3p_t *m3p_init_window(m3p_t *, rci_t , rci_t );
m5p_t *m5p_init_window(m5p_t *, rci_t , rci_t );
m7p_t *m7p_init_window(m7p_t *, rci_t , rci_t );

/**
 * \brief Free a permutation window
 * \param condemned Permutation Matrix
 */

void m3p_free_window(m3p_t *);
void m5p_free_window(m5p_t *);
void m7p_free_window(m7p_t *);

/**
 * \brief copy permutation Q to P
 *
 * \param P Target permutation matrix (may be NULL)
 * \param Q Source permutation matrix (must not be NULL)
 */

m3p_t *m3p_copy(m3p_t *, const m3p_t *);
m5p_t *m5p_copy(m5p_t *, const m5p_t *);
m7p_t *m7p_copy(m7p_t *, const m7p_t *);


/**
 * \brief Set the permutation P to the identity permutation. The only
 * allowed value is 1.
 *
 *
 * \param P Permutation
 * \param value 1
 *

 */

void m3p_set_ui(m3p_t *, unsigned int );
void m5p_set_ui(m5p_t *, unsigned int );
void m7p_set_ui(m7p_t *, unsigned int );

/**
 * Apply the permutation P to A from the left.
 *
 * This is equivalent to row swaps walking from 0 to length-1.
 *
 * \param A Matrix.
 * \param P Permutation.
 */

void m3d_apply_p_left(m3d_t , m3p_t const *);
void m5d_apply_p_left(m5d_t , m5p_t const *);
void m7d_apply_p_left(m7d_t , m7p_t const *);

/**
 * Apply the permutation P to A from the left but transpose P before.
 *
 * This is equivalent to row swaps walking from length-1 to 0.
 *
 * \param A Matrix.
 * \param P Permutation.
 */

void m3d_apply_p_left_trans(m3d_t , m3p_t const *);
void m5d_apply_p_left_trans(m5d_t , m5p_t const *);
void m7d_apply_p_left_trans(m7d_t , m7p_t const *);

/**
 * Apply the permutation P to A from the right.
 *
 * This is equivalent to column swaps walking from length-1 to 0.
 *
 * \param A Matrix.
 * \param P Permutation.
 */

void m3d_apply_p_right(m3d_t , m3p_t const *);
void m5d_apply_p_right(m5d_t , m5p_t const *);
void m7d_apply_p_right(m7d_t , m7p_t const *);

/**
 * Apply the permutation P to A from the right but transpose P before.
 *
 * This is equivalent to column swaps walking from 0 to length-1.
 *
 * \param A Matrix.
 * \param P Permutation.
 */

void m3d_apply_p_right_trans(m3d_t , m3p_t const *);
void m5d_apply_p_right_trans(m5d_t , m5p_t const *);
void m7d_apply_p_right_trans(m7d_t , m7p_t const *);


/**
 * Apply the permutation P to A from the right starting at start_row.
 */

void m3d_apply_p_right_even_capped(m3d_t , m3p_t const *, rci_t , rci_t );
void m5d_apply_p_right_even_capped(m5d_t , m5p_t const *, rci_t , rci_t );
void m7d_apply_p_right_even_capped(m7d_t , m7p_t const *, rci_t , rci_t );

/**
 * Apply the permutation P^T to A from the right starting at start_row.
 *
 * This is equivalent to column swaps walking from 0 to length-1.
 *
 * \param A Matrix.
 * \param P Permutation.
 * \param start_row Start swapping at this row.
 * \param start_col Start swapping at this column.
 *
 * \wordoffset
 */

void m3d_apply_p_right_trans_even_capped(m3d_t , m3p_t const *, rci_t , rci_t );
void m5d_apply_p_right_trans_even_capped(m5d_t , m5p_t const *, rci_t , rci_t );
void m7d_apply_p_right_trans_even_capped(m7d_t , m7p_t const *, rci_t , rci_t );


/**
 * Apply the permutation matrix P to A from the right but transpose P before.
 *
 * This is equivalent to column swaps walking from 0 to length-1.
 *
 * \param A Matrix.
 * \param P Permutation.
 */

void m3d_apply_p_right_trans(m3d_t , m3p_t const *);
void m5d_apply_p_right_trans(m5d_t , m5p_t const *);
void m7d_apply_p_right_trans(m7d_t , m7p_t const *);




/**
 * Apply the permutation P to A from the right, but only on the upper
 * the matrix A above the main diagonal.
 *
 * This is equivalent to column swaps walking from length-1 to 0.
 *
 * \param A Matrix.
 * \param Q Permutation.
 */
void  m3d_apply_p_right_trans_tri(m3d_t , m3p_t const *);
void  m5d_apply_p_right_trans_tri(m3d_t , m3p_t const *);
void  m7d_apply_p_right_trans_tri(m3d_t , m3p_t const *);



/**
 * Print  permutation matrices
 *
 *
 */

void m3p_print(m3p_t const *);
void m5p_print(m3p_t const *)
void m7p_print(m3p_t const *)


/**
 * Compresses the matrix L in a step in blockwise-recursive PLE
 * decomposition.
 *
 * \param A Matrix.
 * \param r1 Rank of left matrix.
 * \param n1 Column cut which separates left and right matrix.
 * \param r2 Rank of right matrix.
 */

void _m3d_compress_l(m3d_t , rci_t r1, rci_t n1, rci_t r2);
void _m5d_compress_l(m3d_t , rci_t r1, rci_t n1, rci_t r2);
void _m7d_compress_l(m3d_t , rci_t r1, rci_t n1, rci_t r2);


   


 

 

 
#endif M1RI_PERMUTATIONS
