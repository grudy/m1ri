 
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
 m3d.h
 */

#ifndef M1RIGF3_M1RIGF3_H
#define M1RIGF3_M1RIGF3_H

#include <stdlib.h>
#include <m1ri/m1riwrappers.h>


/******************************
Creates  a struct of 128 bits
******************************/
typedef struct vbg{
    
    vec units;
    vec sign;
} vbg;

/**********************************
    GF(3) Matrix structure
******************************/
typedef struct {
    
    rci_t nrows; /// Number of rows
    
    rci_t ncols; /// Number of columns
    
    wi_t width; //// The number of vbg's needed to hold columns
    
    vbg * block;  /// Block containing the data contiguous in memory
    
    vbg ** rows;  ///  Pointers to rows of the matrix
    
    vec  svbg;   /// Identifies first vbg used in row
    u_int64_t a;
    u_int32_t  lblock; //  first block pointed to in a window
    u_int32_t fcol;  ///column offset of first block
    u_int8_t flags;    //IsWindowed, NotWindowed    
    
    
} m3d_t;


typedef struct
{

    m3d_t * block;
    m3d_t ** row;
    wi_t slicesize;// (slicesize ^ 2) * 64
    wi_t width;   ///width in slices horizaontally per row
    rci_t nrows;
    rci_t ncols;
    
}m3_slice;



/**
 * \brief 
 *
 * \param c 
 * \param a 
 * \param b 
 *
 * \wordoffset
 */
/**
 Read n bits from a s portion of an element
 x = rows
 y = columns
 M = Matrix read from
 */
void m3d_transpose( m3d_t   *);

/**
 \Brief Read n sign bits
 \param M = Matrix read from 
 \param x = rows
 \param y = columns
 \param n = number of bits to read 
*/
vec m3d_rs_bits(m3d_t const *, rci_t  , rci_t  , int  );

/**
 \Brief Read n unit bits
 \param M = Matrix read from 
 \param x = rows
 \param y = columns
 \param n = number of bits to read 
*/

vec m3d_ru_bits(m3d_t const  *, rci_t  , rci_t  , int  );

/**
 \Brief Read n elements
 \param M = Matrix read from 
 \param x = rows
 \param y = columns
 \param n = elements to read from
*/

vbg m3d_read_elems(m3d_t const *, rci_t  , rci_t  , int  );



m3d_t *  m3d_transposewin(const m3d_t   * );

/**
 \Brief Swap rows in m3d_t 
 \param M = Matrix to swap rows of
 \param a = first set of rows
 \param b = second set of rows to swap
*/
void * m3d_rowswap (m3d_t  * , rci_t , rci_t );

/**
 \Brief Swap columns in m3d_t 
 \param M = Matrix to swap columns of
 \param a = first set of columns
 \param b = second set of columns to swap
*/
void  m3d_colswap(m3d_t *, rci_t , rci_t );

/**
	\Brief Write an value in a matrix
	\param M matrix to write to 
	\param x row of value to change 
	\param y column of value to change
	\param s value of sign 
	\param u value of units
*/
void   m3d_write_elem( m3d_t * ,rci_t , rci_t , vec , vec  );


/**
 \Brief allocate the block pointers for a m5d_t
 \param block m3d_t->block to allocate
 \param nrows rows in matrix
 \param width vfds the matrix takes up
*/
static inline void  * m3d_block_allocate(rci_t  nrows,  wi_t  width)
{
	
    m3d_t * block = m1ri_calloc(nrows * width ,  sizeof(vbg) );
    return block;
}

/**
 \Brief allocate the row pointers for a m3d_t
 \param block m3d_t->block to allocate
 \param nrows rows in matrix
*/
static inline vbg ** m3d_row_alloc(vbg * block, wi_t width, rci_t nrows)
{
	int i;
	vbg ** rows;
    rows = m1ri_malloc( nrows * width * sizeof(vbg *));
    for ( i = 0; i <  nrows;  i++ )
    {
        rows[i]  = block + i * width;
    };
    return rows;
}

/** 
  \brief make a m3d_t
  \nrows rows in m3d_t
  \ncols columns in m3d_t 
  \
*/
m3d_t *  m3d_create(  rci_t , rci_t );


m3d_t m3d_rand(m3d_t * );

/** 
 \brief Make an Identity Matrix times a scalar 'length'
 \param a = Identity matrix 
 \param length = matrix size (row length and column width)
  
*/
void m3d_set_ui(m3d_t *A,unsigned int );



/** 
 \brief identity matrix of size n * n
 \param matrix 
 \param n size of rows and column of matrix
 \
 \Returns an n * n identity matrix
*/
m3d_t  *  m3d_identity(m3d_t  *, rci_t );

/** 
 \brief windows in m1ri_word rows * m1ri_word column incriments
 \param stvbg = the vbg/width offset from the base matrix
 \param strow = row offset in increments of 64
 \param sizecol  = cols * 64
 \param sizerow  = rows * 64
*/
 
m3d_t *    m3d_init_window(const m3d_t  *, rci_t , rci_t , rci_t , rci_t );




/** 
 Concat b on the end of a, the result is c
   [a] [b] ----->  [a b]   ===  C
*/


/** 
    
    \brief concat matrix a and matrix b, write result to c
    \param c = matrix stacked 
    \param a = left matrix
    \param b = right matrix
    \[a]
    \     ===  C
    \[b]
 	\
*/
m3d_t *  m3d_concat( m3d_t * , m3d_t * );

/** 
    
    \brief stack matrix a on matrix b, write result to c
    \param c = matrix stacked 
    \param a = top matrix
    \param b = bottom matrix
    \[a]
    \     ===  C
    \[b]
 	\
*/
m3d_t  * m3d_stack(m3d_t * ,const   m3d_t * , const m3d_t * );
 
/** 
	\brief  Checks if two m3d_t's is equal to another
	\param a = first matrix
	\param b = second matrix
	\return 1 if equal, 0 if false
 */

int m3d_equal(m3d_t const  *, m3d_t const  *);

/**
 * \brief copy matrix b to a
 * \param a matrix to hold copy
 * \param b matrix to copy
 */
m3d_t *  m3d_copy(m3d_t  * , m3d_t const * );
m3d_t * m3d_copy_cutoff(m3d_t  * , m3d_t const * );


/**
 * \brief Releases a m3d_t into the wilderness.  
 * \param a GF(3) matrix
 *
 * \Frees allocated memory in matrix
 */
void m3d_free( m3d_t *  );





/**
 * \brief Data structure  for holding m3d_t matrix windows  
 * \param c Previously malloced structure for holding windows   
 * \param a GF(3) matrix
 * \param slicesize n*n size of slices(matrix windows), where n is a multiple of 64
 * \
 * \
 */
void  m3d_slices(m3_slice *  ,const m3d_t * , wi_t );



/**
 * \brief Creates 4 equally sized windows  
 * \param a Matrix over GF(3) 
 * \
 * \Returns a structure holding windows to four quadrants of  matrix a
 *
 * \[0][1]
 * \[2][2] 
 */
 
m3_slice *  m3d_quarter( const m3d_t * );


/**
 * \brief Creates 4 equally sized windows  
 * \param a Matrix over GF(3) 
 * \
 * \ 
 */
m3d_t  * m3d_transpose_sliced(m3d_t * );

/**
 \Brief Swap columns in m3d_t after a certain row
 \param M = Matrix to swap columns of
 \param a = first set of columns
 \param b = second set of columns to swap
 \param start_row  starting row
*/
void  m3d_colswap_capped_row(m3d_t *, rci_t , rci_t, rci_t );
int m3d_cmp(m3d_t *A, m3d_t *B);

/**
  \brief Negates the  input vbg
  \param r vector to negate
*/
void vbg_negation(vbg * );


/**
	
*/
void sub_m3d( vbg *, vbg const *  , vbg const * );       


/**
 * \brief subtract a by by vector b.   The difference is vector r.
 * \param x = minuend
 * \param y = subtrahend
 * \
 * \returns r with difference
 */

vbg sub_m3dr(vbg , vbg );               


/** *****************************
 matrix r = (direct sum matrix r + matrix x)
 ******************************/
 
/* void iadd_vbg(vbg *,vbg  *); */

/* void isub_m3d(vbg *,vbg  *); */



void  vbg_mul( vbg *, vbg  *, vbg  *);             /* multiply matrix x by y assinging the output to r */

/**
 * \brief subtract a by by vector b.   The difference is vector r.
 * \param r = matrix, may be null 
 * \param x = minuend
 * \param y = subtrahend
 * \
 * \returns r with difference
 */

m3d_t *  m3d_sub(m3d_t *,   const  m3d_t  *, const m3d_t  *);

/**
	Return the value of the matrix multiplied
*/
vbg vbg_mul_i(vbg const , vbg const);

/** 

 * \brief Hadamard Multiplication
 * \param c product of matrix a and b 
 * \param a GF(3) matrix
 * \param slicesize n*n size of slices(matrix windows), where n is a multiple of 64
 * \
 * \
*/
m3d_t * m3d_hadamard(m3d_t * , m3d_t const * , m3d_t const * );




/** * * * * * * * * * * * * * * * * * * * *
 Subtract a 1 kilobyte Matrix from another
 1 kilobyte Matrixhgg
 * * * * * * * * * * * * * * * * * * * * */

void m3d_sub_64(vbg ** , vbg  ** , vbg  ** );

/** * * * * * * * * * * * * * * * * * * * * * *
 Add a 1 kilobyte Matrix from another
 1 kilobyte Matrix
 * * * * * * * * * * * * * * * * * * * * * * */
 /**
 \Brief Add a 64 by 64 m3d_t matrix where the 
 \param R = Where sum is written, may be Null 
 \param A = Matrix to Sum
 \param B = Matrix to Sum
 \Assumes there are 64 values, doesn't check
*/
void m3d_add_64(vbg **, vbg   **  , vbg    **  );

 /**
 \Brief Add matrix a + b = c
 \param c = Where sum is written, may be Null 
 \param a = Matrix to Sum
 \param b = Matrix to Sum
 */
m3d_t  * m3d_add(m3d_t *, const m3d_t  *,const  m3d_t  *);


 /**
 \Brief Return submatrix S from matrix M 
 \param S = Submatrix to be, must be null
 \param M = Matrix to gain a submatrix
 \param b = Matrix to Sum
 */
m3d_t * m3d_submatrix(m3d_t *S, const m3d_t *M, const rci_t lowr, const rci_t lowc, const rci_t highr, const rci_t highc);


/** **************************************************
								GF(3)
****************************************************/
/* 64 * 64,4096 bit, 512 byte matrix(slice) multiplication */
void m3d_mul_64(vbg **, vbg ** , vbg ** );

/* 32 * 64,2048 bit, 256 byte matrix(slice) multiplication */
void mul_32_m3d(vbg *, vbg *, vbg *);

/* 16 * 64,1024 bit, 128 byte matrix(slice) multiplication */
void mul_16_m3d(vbg *, vbg *, vbg *);

/* 8 * 64,512 bit, m1ri_word byte matrix(slice) multiplication */
void mul_8_m3d(vbg *, vbg *, vbg *);

/* 4 * 64,256 bit, 32 byte matrix(slice) multiplication */
void mul_4_m3d(vbg *R, vbg *A, vbg *B);
#endif
