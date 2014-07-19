/** 
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
 
 m5d.h
 */

#ifndef M1RIPORJECT_M5D_H
#define M1RIPROJECT_M5D_H

#include <stdlib.h>
#include <m1ri/m1riwrappers.h>

#ifndef _VFDM5_
#define _VFDM5_
/** 
The Representation for GF5
(Boothby's Additive Method) 
000 0
001 1
011 2
101 3
111 4

*/

/** *****************************
 Creates  a structure for GF(5) Matrices
 ******************************/

typedef struct {
    
    vec units;
    
    vec middle;
    
    vec sign;
    
} vfd;



/** 

This is unused because of the current representation 
static inline  vfd fold5( vec s3, vec s2, vec s1, vec s0)
{
    vfd r;
    vec t = s2 | s1;
    r.units = s0 ^ t;
    r.middle = (r.units & s0) ^ (s3 ^ s1);
    r.sign = (t ^ s2 ) | (r.middle & s3  );
    return r;
}
*/

/** 
 GF(5) Matrix structure
 
 */

typedef struct m5d_t{
    
    rci_t nrows; /* < number of rows */
    
    rci_t ncols; /* < number of columns */
    
    wi_t width; /* < the number vfd's needed to hold columns */

    vfd * block;  /* < block containing the data contiguous in memory */
    
    vfd ** rows;  /*  < pointers to rows of the matrix */
    u_int32_t  fblock; /*   first block pointed to in a window */
    u_int32_t fcol;  /* column offset of first block */
 
    u_int8_t flags;
	vec  svfd;   /* Identifies first vfd used in row */
 
    u_int32_t  lblock; /*   first block pointed to in a window */
    
    
    
} m5d_t;

 
 /**
 * \brief Data structure  for holding m5d_t matrix windows  
 * \param c Previously malloced structure for holding windows   
 * \param a GF(5) matrix
 * \param slicesize n*n size of slices(matrix windows), where n is a multiple of 64
 * \
 * \
 */
typedef struct
{
	/** 
	
	 */
    m5d_t * block;
    /** 
    (slicesize ^ 2) * 64
    */
    m5d_t ** row;
    /*  (slicesize ^ 2) * 64 */
    wi_t slicesize;
    /* width in slices horizaontally per row */
    wi_t width;   
    rci_t nrows;
    rci_t ncols;
    
}m5_slice;


/**
 \Brief allocate the block pointers for a m5d_t
 \param block m5d_t->block to allocate
 \param nrows rows in matrix
 \param width vfds the matrix takes up
*/
static inline vfd  * m5d_block_allocate(vfd * block, rci_t  nrows,  wi_t  width)
{
    
    block  = m1ri_malloc(nrows * width * sizeof(vfd) );
    return block;
     
}

#endif

/**
 * \brief copy matrix b to a
 * \param a matrix to hold copy
 * \param b matrix to copy
 */
m5d_t *  m5d_copy(m5d_t  * , const m5d_t  * );
void m5d_copy_cutoff(m5d_t  * , m5d_t const * );

/** 
 

/* 
 Matrix Windows
 ______________
 
 
 |   [A0 | A1]
 A = | --------
 |   |A2 | A3]
 
 
 
 */

/**
 \Brief Read n middle bits
 \param M = Matrix read from 
 \param x = rows
 \param y = columns
 \param n = number of bits to read 
*/

vec m5d_rm_bits(m5d_t *, rci_t  , rci_t  , int  ) ;


/**
 \Brief Read n sign bits
 \param M = Matrix read from 
 \param x = rows
 \param y = columns
 \param n = number of bits to read 
*/

vec m5d_rs_bits(m5d_t *, rci_t  , rci_t  , int  );


/**
 \Brief Read n sign bits
 \param M = Matrix read from 
 \param x = rows
 \param y = columns
 \param n = number of bits to read 
*/
vec m5d_ru_bits(m5d_t *, rci_t  , rci_t  , int  );

/**
 \Brief Read n elements
 \param x = rows
 \param y = columns
 \param M = Matrix read from 
 \param n = elements to read from

*/
vfd m5d_read_elems(m5d_t *M, rci_t  x, rci_t  y, int  );

/**
 \Brief Swap rows in m5d_t 
 \param M = Matrix to swap rows of
 \param a = first set of rows
 \param b = second set of rows to swap
*/

void * m5d_rowswap (m5d_t * , rci_t , rci_t  );


/**
	\Brief Write an value in a matrix
	\param M matrix to write to 
	\param x row of value to change 
	\param y column of value to change
	\param s value of sign 
	\param m value of middle
	\param u value of units
*/

void *  m5d_write_elem( m5d_t * ,rci_t , rci_t , vec , vec , vec);




/**
 \Brief allocate the row pointers for a m5d_t
 \param block m5d_t->block to allocate
 \param  rows in matrix
 
*/
vfd ** m5d_row_alloc(vfd * , vfd ** , wi_t , rci_t );



/** 
  \brief make a m5d_t
  \nrows rows in m5d_t
  \ncols columns in m5d_t 
  \
*/
m5d_t * m5d_create( rci_t nrows, rci_t ncols);




vfd * m5d_rand(m5d_t * );

/** 
 \brief Make an Identity Matrix
 \a = Identity matrix
 \n = matrix size (row length and column width)
*/
void   m5d_set_ui(m5d_t *, rci_t );


/** 
 \brief identity matrix of size n * n
 \param matrix 
 \param n size of rows and column of matrix
 \
 \Returns an n * n identity matrix
*/
m5d_t  *  m5d_identity(rci_t );


/**
 * \brief Releases a m5d_t into the wilderness.  
 * \param a GF(5) matrix
 *
 * \Frees allocated memory in matrix
 */
void m5d_free( m5d_t *  );


/**
 * \brief subtract a by by vector b.   The difference is vector r.
 * \param r vector to right result to
 * \param a = minuend
 * \param b = subtrahend
 * \
 */
void vfd_sub( vfd *, vfd *, vfd *);               
/** *****************************
 matrix r = (direct sum matrix r + matrix x)
 ******************************/


m5d_t * m5d_add(m5d_t * , const m5d_t *,const  m5d_t *);

/** 
 * \brief subtract a by by vector b.   The difference is vector c.
 * \param c = matrix, may be null 
 * \param a = minuend
 * \param b = subtrahend
 * \
 * \returns c with difference
 */

m5d_t * m5d_sub(m5d_t * , const m5d_t *,const  m5d_t *);


void vfd_sub( vfd *, vfd *, vfd *);               /* subtract vector x by by vector y.   The product is vector r. */



/** 
 * \brief subtract a by by vector b.   The difference is written to a
 * \param a minuend initially
 * \param b subtrahend
 \
 \ a becomes difference
 */
 
void m5d_sub_d(m5d_t  *  , m5d_t * );

/** 
 * \brief subtract a by by vector b.   The difference is vector r.
 * \param r = matrix, may be null 
 * \param a = minuend
 * \param b = subtrahend
 * \
 * \returns r with difference
 */
void add_vfd(vfd *, vfd * , vfd *);



/*
 * \brief like add_vfd but sum is multiplied by 2 
 */
void m5d_add2(vfd * , vfd * , vfd * );


/*
 * \brief like add_vfd but sum is multiplied by 2 
 
 */
void m5d_add2_i(vfd * , vfd * );


/** 
	\brief  Checks if two m5d_t's is equal to another
	\param a = first matrix
	\param b = second matrix
	\return 1 if equal, 0 if false
 */
int m5d_equal(m5d_t const *, m5d_t const *);


void m5d_add_64(vfd **, vfd **  , vfd ** );


/** 
 \brief windows in M1RI_RADIX rows * M1RI_RADIX column incriments
 \param stvbg = the vbg or width offset from the base matrix
 \param strow = row offset in increments of 64
 \param sizecol  = cols * 64
 \param sizerow  = rows * 64
 */

m5d_t   * m5d_init_window( m5d_t const *, const rci_t , const  rci_t ,const rci_t , const rci_t );

void m5d_sub_64(vfd **c ,vfd  ** a , vfd **b);


  /**
 * \brief subtract a by by vector b.   The difference is vector r.
 * \param x = minuend
 * \param y = subtrahend
 * \
 * \returns r with difference
 */      
void m5d_sub_i(vfd *,vfd *);

/**
 * \brief scalar multiplication of vfd vector
 * \param a vfd to be multiplied
 * \returns each value in the vector multiplied by 2
 */
vfd m5d_mul2(vfd);

/**
 * \brief scalar multiplication of vfd vector
 * \param a vfd to be multiplied
 * \returns each value in the vector multiplied by 3
 */
vfd m5d_mul3(vfd);

/**
 * \brief scalar multiplication of vfd vector
 * \param a vfd to be multiplied
 * \returns each value in the vector multiplied by r
 */
vfd m5d_mul4(vfd);




/* These Functions are expanded upon Toms 
 implementation of the method of four RUSSIANS
 */
void *  m5d_combine3(vfd *, vfd * );
void m5d_combine4(vfd *, vfd * );
void m5d_combine5(vfd *, vfd * );
void m5d_combine6(vfd *, vfd * );
void m5d_combine7(vfd *, vfd * );
void m5d_combine8(vfd *, vfd *);

/* ***************************************************
				GF(5) base cases for multiplication
*************************************************** */

void m5d_mul_64(vfd **, vfd **, vfd **);
void m5d_mul_32(vfd *, vfd *, vfd *);
void m5d_mul_16(vfd *, vfd *, vfd *);
void m5d_mul_8(vfd *, vfd *, vfd *);
void m5d_mul_4(vfd *R, vfd *A, vfd *B);
/* vfd * m5d_transpose_vfd(vfd  **, vfd **); */
void  m5d_transpose(m5d_t   * a);



/**
 * \brief Data structure  for holding m5d_t matrix windows  
 * \param c Previously malloced structure for holding windows   
 * \param a GF(5) matrix
 * \param slicesize n*n size of slices(matrix windows), where n is a multiple of 64
 * \
 * \
 */
void  m5d_slices(m5_slice *  , m5d_t * , wi_t );


/**
 * \brief Creates 4 equally sized windows  
 * \param a Matrix over GF(5) 
 * \
 * \Returns a structure holding windows to four quadrants of  matrix a
 *
 * \[0][1]
 * \[2][2] 
 */

m5_slice *  m5d_quarter(const m5d_t * );

m5d_t *  m5d_transpose_sliced(m5d_t * );

m5d_t  * m5_blockslice_allocate( rci_t  ,  wi_t  );

m5d_t ** m5_rowslice_allocate(m5d_t * ,  wi_t , rci_t );


/** 

 * \brief Hadamard Multiplication
 * \param c product of matrix a and b 
 * \param a GF(5) matrix
 * \param slicesize n*n size of slices(matrix windows), where n is a multiple of 64
 * \
 * \
*/
m5d_t * m5d_hadamard(m5d_t *, m5d_t const * , m5d_t const *  );


/**
 \Brief Swap columns in m5d_t 
 \param M = Matrix to swap columns of
 \param a = first set of columns
 \param b = second set of columns to swap
*/
void  m5d_colswap(m5d_t *, rci_t , rci_t );

/**
 \Brief Swap columns in m5d_t after a certain row
 \param M = Matrix to swap columns of
 \param a = first set of columns
 \param b = second set of columns to swap
 \param start_row  starting row
*/
void m5d_colswap_capped_row(m5d_t *, rci_t , rci_t, rci_t );

int m5d_cmp(m5d_t *A, m5d_t *B);

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
m5d_t *  m5d_concat( m5d_t * , m5d_t * );

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
m5d_t  * m5d_stack(m5d_t * ,const   m5d_t * , const m5d_t * );
 
int m5d_is_zero(const m5d_t *);
 /**
 \Brief Return a scalar product of the input Matrix B, return result to C
 \param C Matrix to return, can be null
 \param a Scalar
 \param B Matrix to Multiply
 */
m5d_t *m5d_mul_scalar(m5d_t *, const long , const m5d_t *);


 
 

#endif
