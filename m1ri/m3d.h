 
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
typedef struct vbg
{
    vec sign;
    vec units;

} vbg;

/**********************************
    \brief GF(3) Matrix structure
    \
    \00 = 0
	\01 = 1
	\11 = 2  
	\
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
 Read n bits from a s portion of an element
 x = rows
 y = columns
 M = Matrix read from
 *
void m3d_transpose( m3d_t   *);
*/

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



//m3d_t *  m3d_transposewin(const m3d_t   * );

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
	
    //m3d_t * block = m1ri_calloc(nrows * width ,  sizeof(vbg) );
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
    rows = m1ri_calloc( nrows ,  sizeof(vbg *));
    for ( i = 0; i <  nrows;  i++ )
    {
		rows[i]  = block + (i * width);
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

/**
 * \brief Fill matrix with random values
 * \
 * \
 * \
 * \wordoffset
 */
void m3d_rand(m3d_t * );

/** 
 \brief random matrix of size n * n
 \param a null m3d_t 
 \param n size of rows and column of matrix
 \
 \Returns an n * n identity matrix
*/
static inline m3d_t * m3d_create_rand(m3d_t * a, rci_t n)
{
	 
	 a = m3d_create( n, n);
	 m3d_rand(a);
	 return a;
	

}
/** 
 \brief Make an Identity Matrix times a scalar 'length'
 \param a = Identity matrix 
 \param length = matrix size (row length and column width)
  
*/
void m3d_set_ui(m3d_t *A,unsigned int );



/** 
 \brief identity matrix of size n * n
 \param a matrix 
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
 \Brief Swap columns in m3d_t after a certain row
 \param M = Matrix to swap columns of
 \param a = first set of columns
 \param b = second set of columns to swap
 \param start_row  starting row
 \param end_row end row 
*/
static inline void m3d_col_swap_in_rows(m3d_t *M, rci_t col_a, rci_t col_b, rci_t start_row, rci_t end_row)
{
	
		 if (col_a == col_b)
    		return;
        int i;
        vec block_a, block_b, dif_a, dif_b, a_place, b_place; 
        vbg tempa, tempb;
         block_a = (col_a-1)/M1RI_RADIX;
         block_b = (col_b-1)/M1RI_RADIX;
         dif_a = col_a%M1RI_RADIX;
         dif_b = col_b%M1RI_RADIX;
         a_place =  leftbit >>  dif_a ;
         b_place =  leftbit >> dif_b ;
    

              
          for( i = start_row; i <  end_row; i++)
          {
		     
		  
		           tempa.units  = (b_place  & M->rows[i][block_b].units) ? (a_place  ): 0;
		     tempb.units  = (a_place  & M->rows[i][block_a].units) ? (b_place  ): 0; 
		       M->rows[i][block_a].units  = (tempa.units)  ? M->rows[i][block_a].units  | tempa.units :   M->rows[i][block_a].units  & ~a_place; 
		       M->rows[i][block_b].units  = (tempb.units)  ? M->rows[i][block_a].units  | tempb.units :   M->rows[i][block_b].units  & ~b_place;  
		         tempa.sign  = (b_place  & M->rows[i][block_b].sign) ? (a_place  ): 0;
		     tempb.sign  = (a_place  & M->rows[i][block_a].sign) ? (b_place  ): 0; 
		       M->rows[i][block_a].sign  = (tempa.sign)  ? (M->rows[i][block_a].sign  | tempa.sign) :   M->rows[i][block_a].sign  & ~a_place; 
		       M->rows[i][block_b].sign  = (tempb.sign)  ? (M->rows[i][block_a].sign  | tempb.sign) :   M->rows[i][block_b].sign  & ~b_place; 
		     
		       

		       
          }
    
   


}



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
m3d_t *  m3d_concat(m3d_t * ,const   m3d_t * , const m3d_t * );

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
 * \brief Releases a m3_slice  allocated by a m3d_quarter function 
 * \param a GF(3) matrix slice structure
 *
 * \Frees allocated and memory in slices, including windows
 */
static inline void m3d_quarter_free(m3_slice *a)
{
	
	m3d_free(a->row[0]);
	m3d_free(a->row[1]);
	m3d_free(a->row[2]);
	m3d_free(a->row[3]);
	
	m1ri_free(a->row);
	m1ri_free(a);
	
}


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
 
m3_slice *  m3d_quarter( const m3d_t * a );





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


/** 
	\brief  Compares two m3d_t matrices
	\param a = first matrix
	\param b = second matrix
	\return 1 if equal, 0 
 */
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

static inline void sub_m3dr(vbg  * r,vbg const * x, vbg const * y)

{
    r->units = ((x->units^y->units) | (x->sign^y->sign));
    r->sign = (((x->units^y->units)^x->sign)&(y->units ^ x->sign));
    
 
}
  
         
/** 
 \brief Optimized subtract for  64 * 64 matrix
 * \param R = minuend
 * \param A = subtrahend
 * \param B = difference
*/
static inline void m3d_sub_64(vbg **R, vbg  **A, vbg  **B)
{
    int i;
    for (i= 0; i < M1RI_RADIX; i++ )
    {
      sub_m3dr(R[i], A[i], B[i]);
    }
    
}


/** *****************************
 matrix r = (direct sum matrix r + matrix x)
 ******************************/
 


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


static inline void  m3d_sub_unshackled(m3d_t * r,   const  m3d_t  *x, const m3d_t  *y)
{
	
	int n , i;
	for(i = 0; i < x->nrows; i++)
    {
    	for(n = 0; n < x->width; n++)
        {
		  sub_m3d(r->rows[i] + n, x->rows[i] + n , y->rows[i] + n );
        }
    }
   



}


static inline void isub_m3d(vbg  *r,vbg  *x)
{
    vec t;
    r->units = x->units ^ r->units;
    t  = r->units | r->sign;
    t = t ^ x->sign;
    r->sign = x->units ^ r->sign;
    r->sign = r->sign & t;
    r->units = t | r->units;     
}


static inline void add_vbg(vbg *   r, vbg const *   x, vbg const * y)

{
    /* r->units = (x->units ^ y->sign) & (x->sign ^ y->units); // ///r0 ← (x0 ⊕y->1)∧(x1 ⊕y->0); */
    /* r->sign = (M1RI_ST(x->units, y->sign, x->sign ) | M1RI_ST(x->sign, y->units, y->sign)); //// r1 ← s XOR t. */



 	r->units = y->sign ^ x->units;
    r->sign = y->units ^ x->sign;
    r->sign = r->sign & r->units;
    r->units = r->units ^ x->sign;
    r->units = (y->units ^ x->units) | r->units;




}



static inline void   m3d_add_unshackled(m3d_t *c, const m3d_t  *a,const m3d_t  *b)
{

   int i, j;
        for( i = 0; i < a->nrows; i++)
        {
            for(j = 0; j < (a->width ); j++)
            {
            	add_vbg(&c->rows[i][j], &a->rows[i][j], &b->rows[i][j]);
        	}   
        }
}

/** 
	\brief, addition base case
	\param x augend
	param  y addend
	\return x as sum
*/
static inline void m3d_inc(vbg  * x, vbg const *  y)
{ 
    vec t;
    x->sign  = y->units ^ x->sign;
    t = (x->sign & x->units) ^ y->sign;
    x->units = (y->units ^ x->units) |  t;
    x->sign = t & x->sign;
    
}





/**
	Return the value of the matrix multiplied
*/
vbg vbg_mul_elementwise(vbg const , vbg const);

/** 

 * \brief Hadamard Multiplication
 * \param c product of matrix a and b 
 * \param a GF(3) matrix
 * \param slicesize n*n size of slices(matrix windows), where n is a multiple of 64
 * \
 * \
*/
m3d_t * m3d_hadamard(m3d_t * , m3d_t const * , m3d_t const * );







/** * * * * * * * * * * * * * * * * * * * * * *
 Add a 1 kilobyte Matrix from another
 1 kilobyte Matrix
 * * * * * * * * * * * * * * * * * * * * * * */
 /**
 \Brief Add a 64 by 64 m3d_t matrix where the 
 \param R = Where sum is written, may be Null 
 \param A = augend
 \param B = addend
 \Assumes there are 64 values, doesn't check
*/



static inline void m3d_add_64(vbg **R, vbg    ** const A, vbg   ** const B)
{
	
	int i;
    for (i = 0; i < M1RI_RADIX; i++ )
    {
         add_vbg(R[i], A[i], B[i]);
    }
    
    
    
    
    

}
 /**
 \Brief Add matrix a + b = c
 \param c  Where sum is written, may be Null 
 \param a  augend
 \param b  addend
 */
m3d_t  * m3d_add(m3d_t *, const m3d_t  *,const  m3d_t  *);

//void   m3d_add_unshackled(m3d_t *, const m3d_t  *,const  m3d_t  *);

 /**
 \Brief Return submatrix S from matrix M 
 \param S Submatrix to be, must be null
 \param M  Matrix to gain a submatrix
 \param lowr lower row
 \param lowc lower column
 \param highr high row
 \param highc high column 
 */
m3d_t * m3d_submatrix(m3d_t *, const m3d_t *, const rci_t , const rci_t , const rci_t , const rci_t);

 /**
 \Brief Return a scalar product of the input Matrix B, return result to C
 \param C Matrix to return, can be null
 \param a Scalar
 \param B Matrix to Multiply
 */
m3d_t *m3d_mul_scalar(m3d_t *, const long , const m3d_t *);



 /**
 \Brief Row A[ar] = A[ar] + B[ar] after the startcol.
 \param A Matrix to return, can NOT be null
 \param ar Row to sum of Matrix A
 \param B Matrix remains unchanged
 \param br Row to add to sum of Matrix A
 \pram  startcol column to start adding the row
 \
 */
void m3d_add_row(m3d_t *A, rci_t ar, const m3d_t *B, rci_t br, rci_t start_col);

 /** 
    
    \brief find if the input Matrix is 0
    \param a input matrix, must NOT be NULL
 	\
 	\returns 1 if zero, else 0
 	
*/
int m3d_is_zero(const m3d_t *);

static inline void m3d_set_zero(m3d_t * a)
{
  if(a == NULL)
  {
  	  m1ri_die("m3d_set_zero: bad argument to m3d_set_zero\n");

  
  }
  for(int i = 0; i < a->nrows; i++)
  {
    for(int j = 0; j < a->width; j++)
    {
       a->rows[i][j].units = 0;
       a->rows[i][j].sign = 0;
    }
  
  }

}





/** **************************************************
								GF(3)
****************************************************/
/* 64 * 64,4096 bit, 512 byte matrix(slice) multiplication */
void m3d_mul_64(vbg **, vbg  **  const, vbg  **  const);


#endif
