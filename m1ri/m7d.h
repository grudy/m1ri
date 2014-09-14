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
 
 m7d.h
 */
#ifndef M1RIPROJECT_M7D_H
#define M1RIPROJECT_M7D_H
#include <stdlib.h>
#include <m1ri/m1riwrappers.h>


/** *****************************
 Creates  a structure of 192 bits
 ******************************/

typedef struct {
    
    vec units;
    vec middle;
    vec sign;
    
} vtri;




/** 
 GF(7) Matrix structure
 GF7

000 = 0    
100  = 1
010 = 2
110 = 3
001 = 4
101 = 5
011 = 6

 */

typedef struct {
    
    rci_t nrows; /* < number of rows */
    
    rci_t ncols; /* < number of columns */
    
    wi_t width; /* < the number of vbg's needed to hold columns */
    
    vtri * block;  /* < block containing the data contiguous in memory */
    
    vtri ** rows;  /*  < pointers to rows of the matrix */
    
    vec  svtri;   /* Identifies first vbg used in row */
    
    /*  wi_t rowstride;  //vbg's in block to traverse to  get to first */
    
    u_int32_t  lblock; /*   first block pointed to in a window */
    u_int32_t fcol;  /* column offset of first block */
    u_int8_t flags;    /* IsWindowed, NotWindowed */

    
    
} m7d_t;


static inline void m7d_inc(vtri  *x, vtri *y)
{
     vtri  r;
    vec s;
    vec t;

    r.units = x->units ^ y->units;
    s = (x->units & y->units);
    r.middle = s^ x->middle ^ y->middle;
    t = (((s) & (x->middle | y->middle)) | (x->middle & y->middle) );
    r.sign = x->sign ^ y->sign ^ t;
    s =  ((t) & (x->sign | y->sign)) | (x->sign & y->sign);
    
    t = (r.units & s );
    x->units  = s ^ r.units;
    s = t & r.middle;
    x->middle  = r.middle ^ t ;
    x->sign  = r.sign ^ s;

}

typedef struct
{
     
    m7d_t * block;
    m7d_t ** row;
    wi_t slicesize;/*  (slicesize ^ 2) * 64 */
    wi_t width;   /* /width in slices horizaontally per row */
    rci_t nrows;
    rci_t ncols;
    
}m7_slice;




/** 
 Matrix Windows
 ______________
  
 | [A0 | A1]
 A = | --------
 | |A2 | A3]

 */


/**
 \Brief Read n middle bits
 \param x = rows
 \param y = columns
 \param M = Matrix read from 
 \param n = elements to read from
*/

vec m7d_rm_bits(m7d_t *M, rci_t  x, rci_t  y, int  n) ;




/**
 \Brief Read n sign bits
 \param M = Matrix read from 
 \param x = rows
 \param y = columns
 \param n = number of bits to read 
*/

vec m7d_rs_bits(m7d_t *M, rci_t  x, rci_t  y, int  n);

/**
 \Brief Read n units bits
 \param M = Matrix read from 
 \param x = rows
 \param y = columns
 \param n = number of bits to read 
*/
vec m7d_ru_bits(m7d_t *M, rci_t  x, rci_t  y, int  n);



/**
 \Brief Read n sign bits
 \param x = rows
 \param y = columns
 \param M = Matrix read from 
 \param n = elements to read from
*/


vtri m7d_read_elems(m7d_t *M, rci_t  x, rci_t  y, int  n);


/**
 \Brief Swap rows in m7d_t 
 \param M = Matrix to swap rows of
 \param a = first set of rows
 \param b = second set of rows to swap
*/
void  m7d_rowswap (m7d_t * M, rci_t row_a, rci_t  row_b);

/**
 * \brief copy matrix b to a
 * \param a matrix to hold copy
 * \param b matrix to copy
 */
m7d_t *  m7d_copy(m7d_t  * , m7d_t const * );


m7d_t * m7d_copy_cutoff(m7d_t  * , m7d_t const * );





/**
	\Brief Write an value in a matrix
	\param M matrix to write to 
	\param x row of value to change 
	\param y column of value to change
	\param s value of sign 
	\param m value of middle
	\param u value of units
*/

void *  m7d_write_elem( m7d_t * M,rci_t x, rci_t y, vec s,  vec m , vec u );
/** 
 
 */

vtri  * m7d_block_allocate(vtri * block, rci_t  nrows,  wi_t  width);


/**
 \Brief allocate the row pointers for a m7d_t
 \param rows pointed to 
 \param block m7d_t->block to allocate
 \param nrows rows in matrix
*/
vtri ** m7d_row_alloc(vtri * block, vtri ** rows, wi_t width, rci_t nrows);
/** 
 
 */

m7d_t  * m7d_create(rci_t nrows, rci_t ncols);
/** 
 
 */
void reduce_vtri( vtri * );
void m7d_rand(m7d_t * a);
/** 
 \brief Make an Identity Matrix
 \a = Identity matrix
 \n = matrix size (row length and column width)
 
 
 */

void  m7d_set_ui(m7d_t *, rci_t );


/** 
 \brief identity matrix of size n * n
 \param matrix 
 \param n size of rows and column of matrix
 \
 \Returns an n * n identity matrix
*/
m7d_t  *  m7d_identity(m7d_t * , rci_t );


/** 
 \brief random matrix of size n * n
 \param a null m7d_t 
 \param n size of rows and column of matrix
 \
 \Returns an n * n identity matrix
*/
static inline m7d_t * m7d_create_rand(m7d_t * a, rci_t n)
{
	 
	 a = m7d_create( n, n);
	 m7d_rand(a);
	 return a;
	

}


/** 
	\brief  Checks if two m7d_t's is equal to another
	\param a = first matrix
	\param b = second matrix
	\return 1 if equal, 0 if false
 */
 
int m7d_equal(m7d_t const *, m7d_t const *);


/**
 * \brief Releases a m7d_t into the wilderness.  
 * \param a GF(7) matrix
 *
 * \Frees allocated memory in matrix
 */
void m7d_free( m7d_t *  );


/**
 * \brief Releases a m7_slice  allocated by a m7d_quarter function 
 * \param a GF(7) matrix slice structure
 *
 * \Frees allocated and memory in slices, including windows
 */
static inline void m7d_quarter_free(m7_slice *a)
{
	
	m7d_free(a->row[0]);
	m7d_free(a->row[1]);
	m7d_free(a->row[2]);
	m7d_free(a->row[3]);
	
	m1ri_free(a->row);
	m1ri_free(a);
	
}
void vtri_negate( vtri * );

/** *****************************
 matrix r = (direct sum matrix r + matrix x)
 ******************************/

static inline void add_vtri(vtri * r, vtri  const * x, vtri const * y)

{
    
    
    /*
    s3s2s1s0 ← add(a2a1a0, b2b1b0)
    r2r1r0 ← add(s2s1s0, s3)

    */
    vec s;
    vec t;
    
    r->units = x->units ^ y->units;

    s = (x->units & y->units);
    r->middle = s^ x->middle ^ y->middle;
    t = ((s) & (x->middle | y->middle)) | (x->middle & y->middle);
    r->sign = x->sign ^ y->sign ^ t;
    /* to here I know is right */
    
    
    s = ((t) & (x->sign | y->sign)) | (x->sign & y->sign);
    
    t = s & r->units;
    r->units = s ^ r->units;
    s= t & r->middle;
    r->middle = r->middle ^ t;
    r->sign = r->sign | s;

    
    
}



/**

 \Brief Add a 64 by 64 m3d_t matrix where the 
 \param R = Where sum is written, may be Null 
 \param A = augend
 \param B = addend
 \Assumes there are 64 values, doesn't check
*/



static inline void m7d_add_64(vtri **R, vtri   * const *A, vtri  * const *B)
{
    int i;
    for (i = 0; i < M1RI_RADIX; i++ )
    {
    	add_vtri(&R[i][0], &A[i][0], &B[i][0]);
       /*  R[i][0] = add_m7dr(A[i][0], B[i][0]); */
    }

}

static inline void m7d_dec(vtri  * x, vtri * y)
{
	vtri temp;
   	temp.units = ~y->units;
   	temp.middle = ~y->middle;
   	temp.sign = ~y->sign;
   	
   	m7d_inc( x, &temp);

}



/**
 * \brief subtract a by by vector b.   The difference is vector r.
 * \param r = matrix, may be null 
 * \param x = minuend
 * \param y = subtrahend
 * \
 * \returns r with difference
 */

m7d_t *  m7d_sub(m7d_t *,   const  m7d_t  *, const m7d_t  *);

void m7d_vtri_sub(vtri *, vtri *, vtri *  );
void m7d_sub_64(vtri **, vtri   **, vtri  **);




/**
 * \brief Data structure  for holding m7d_t matrix windows  
 * \param c Previously malloced structure for holding windows   
 * \param a GF(7) matrix
 * \param slicesize n*n size of slices(matrix windows), where n is a multiple of 64
 * \
 * \
 */


 
 
vtri sub_m7dr(vtri const x, vtri const y);
/** 
	GF(7) Addition on a single M1RI word.
*/
m7d_t * m7d_add( m7d_t *, const m7d_t *,const  m7d_t *);

static inline void m7d_add_2r(vtri *, vtri *);

static void m7d_add_4r( vtri *, vtri *);



/** 
 \brief windows in M1RI_RADIX rows * M1RI_RADIX column incriments
 \param stvbg = the vbg or width offset from the base matrix
 \param strow = row offset in increments of 64
 \param sizecol  = cols * 64
 \param sizerow  = rows * 64
 */
m7d_t   * m7d_init_window(const m7d_t *, const rci_t , const  rci_t ,const rci_t , const rci_t );


void  m7d_slices(m7_slice *  , m7d_t * , wi_t );

/**
 * \brief Creates 4 equally sized windows  
 * \param a Matrix over GF(7) 
 * \
 * \Returns a structure holding windows to four quadrants of  matrix a
 *
 * \[0][1]
 * \[2][2] 
 */
 
m7_slice * m7d_quarter(const m7d_t * );


m7d_t  * m7_blockslice_allocate( rci_t  ,  wi_t  );
m7d_t ** m7_rowslice_allocate(m7d_t * ,  wi_t , rci_t );




void m7d_mul_64(vtri **, vtri **, vtri **);




m7d_t * m7d_hadamard(m7d_t * ,const m7d_t  * ,const m7d_t  *  );


/**
 \Brief Swap columns in m7d_t 
 \param M = Matrix to swap columns of
 \param a = first set of columns
 \param b = second set of columns to swap
*/
void m7d_colswap(m7d_t *, rci_t , rci_t);
/**
 \Brief Swap columns in m5d_t after a certain row
 \param M = Matrix to swap columns of
 \param a = first set of columns
 \param b = second set of columns to swap
 \param start_row  starting row
 \param end_row end row 
*/
static inline void m7d_col_swap_in_rows(m7d_t *M, rci_t col_a, rci_t col_b, rci_t start_row, rci_t end_row)
{
      if((M->ncols >= (col_a ) && (M->ncols >= col_b)))
    {
        int i;
        vec block_a, block_b, dif_a, dif_b, a_place, b_place; 
        vtri tempa, tempb;
         block_a = (col_a-1)/M1RI_RADIX;
         block_b = (col_b-1)/M1RI_RADIX;
         dif_a = col_a%M1RI_RADIX;
         dif_b = col_b%M1RI_RADIX;
         a_place =  rightbit <<  dif_a ;
         b_place =  rightbit << dif_b ;
        if(block_a == block_b)
        { 

              
          for( i = start_row; i < end_row; i++)
          {
		     
		  
		    tempa.units  = (b_place  & M->rows[i][block_b].units) ? (a_place  ): 0;
		    tempb.units  = (a_place  & M->rows[i][block_a].units) ? (b_place  ): 0; 
		     
		    tempa.middle  = (b_place  & M->rows[i][block_b].middle) ? (a_place  ): 0;
		    tempb.middle  = (a_place  & M->rows[i][block_a].middle) ? (b_place  ): 0; 
		     
		    tempa.sign  = (b_place  & M->rows[i][block_b].sign) ? (a_place  ): 0;
		    tempb.sign  = (a_place  & M->rows[i][block_a].sign) ? (b_place  ): 0; 
		    M->rows[i][block_a].units  = (tempa.units)  ? M->rows[i][block_a].units  | tempa.units :   M->rows[i][block_a].units  & ~a_place; 
		    M->rows[i][block_b].units  = (tempb.units)  ? M->rows[i][block_a].units  | tempb.units :   M->rows[i][block_b].units  & ~b_place;  
		     
		    M->rows[i][block_a].sign  = (tempa.sign)  ? (M->rows[i][block_a].sign  | tempa.sign) :   M->rows[i][block_a].sign  & ~a_place; 
		    M->rows[i][block_b].sign  = (tempb.sign)  ? (M->rows[i][block_a].sign  | tempb.sign) :   M->rows[i][block_b].sign  & ~b_place; 
		      
		    M->rows[i][block_a].middle  = (tempa.middle)  ? (M->rows[i][block_a].middle  | tempa.middle) :   M->rows[i][block_a].middle  & ~a_place; 
		    M->rows[i][block_b].middle  = (tempb.middle)  ? (M->rows[i][block_a].middle  | tempb.middle) :   M->rows[i][block_b].middle  & ~b_place; 
		     
		       
		       

		       
          }
    
        }
        
      
        
        
    } 

}

/**
 \Brief Swap columns in m7d_t after a certain row
 \param M = Matrix to swap columns of
 \param a = first set of columns
 \param b = second set of columns to swap
 \param start_row  starting row
*/
void m7d_colswap_capped_row(m7d_t *, rci_t , rci_t, rci_t );



int m7d_cmp(m7d_t *A, m7d_t *B);


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
m7d_t *  m7d_concat(m7d_t * ,const   m7d_t * , const m7d_t * );

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
m7d_t  * m7d_stack(m7d_t * ,const   m7d_t * , const m7d_t * );
 
 
/** 
    
    \brief find if the input Matrix is 0
    \param a input matrix, must NOT be NULL
 	\
 	\returns 1 if zero, else 0
 	
*/ 
int m7d_is_zero(const m7d_t *);

void m7d_add_i(m7d_t * , m7d_t *); 

void m7d_sub_i(m7d_t * , m7d_t *); 




/*

	\brief incremental subtraction, but where the subtrahend is changed
	
*/
void m7d_sub_r(m7d_t   *,  m7d_t const *);
/**
	Return the value of the matrix multiplied
*/

 /**
 \Brief Return a scalar product of the input Matrix B, return result to C
 \param C Matrix to return, can be null
 \param a Scalar
 \param B Matrix to Multiply
 */
 
m7d_t *m7d_mul_scalar(m7d_t *, const long , const m7d_t *);


 /**
 \Brief Row A[ar] = A[ar] + B[ar] after the startcol.
 \param A Matrix to return, can NOT be null
 \param ar Row to sum of Matrix A
 \param B Matrix remains unchanged
 \param br Row to add to sum of Matrix A
 \pram  startcol column to start adding the row
 \
 */
void m7d_add_row(m7d_t *, rci_t , const m7d_t *, rci_t , rci_t );

 

#endif
