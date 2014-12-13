
/** *





 TOMAS J. BOOTHBY AND ROBERT W. BRADSHAW "BITSLICING AND THE METHOD OF FOUR
 RUSSIANS OVER LARGER FINITE FIELDS"

 Copyright (C) 2008 Martin Albrecht <malb@informatik.uni-bremen.de>
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
 m1ri_permutations.c
 */




#include "m1ri_permutations.h"


/**
 * Construct an identity permutation.
 *
 * \param length Length of the permutation.
 */

m3p_t *m3p_init(rci_t length)
{
  m3p_t * a = m1ri_malloc(sizeof(m3p_t));
  a->values = m1ri_malloc(sizeof(rci_t) * length);
  a->length = length;
  for(rci_t i = 0; i < length; i++)
  {
    a->values[i] = i;
  }


  return a;
}

m5p_t *m5p_init(rci_t length)
{
	m5p_t * a = m1ri_malloc(sizeof(m3p_t));
  	a->values = m1ri_malloc(sizeof(rci_t) * length);
  	a->length = length;
  	for(rci_t i = 0; i < length; i++)
  	{
   	 a->values[i] = i;
  	}


  return a;

}
m7p_t *m7p_init(rci_t length)
{

  m7p_t * a = m1ri_malloc(sizeof(m3p_t));
  a->values = m1ri_malloc(sizeof(rci_t) * length);
  a->length = length;
  for(rci_t i = 0; i < length; i++)
  	{
   	 a->values[i] = i;
  	}


  return a;

}

/**
 * Free a permutation.
 *
 * \param P Permutation to free.
 */

void m3p_free(m3p_t *P)
{
  m1ri_free(P->values);
  m1ri_free(P);
}
void m5p_free(m5p_t *P)
{
  m1ri_free(P->values);
  m1ri_free(P);
}
void m7p_free(m7p_t *P)
{
  m1ri_free(P->values);
  m1ri_free(P);
}

/**
 * \brief Create a window/view into the permutation P.
 * \param P Permutation matrix
 * \param begin Starting index (inclusive)
 * \param end   Ending index   (exclusive)
 *
 */

m3p_t * m3p_init_window(m3p_t *P, rci_t begin, rci_t end)
{
  m3p_t* window = ( m1ri_malloc(sizeof(m3p_t)));
  window->values = P->values + begin;
  window->length = end - begin;
  return window;
}

m5p_t * m5p_init_window(m5p_t *P, rci_t begin, rci_t end)
{
  m5p_t* window = ( m1ri_malloc(sizeof(m5p_t)));
  window->values = P->values + begin;
  window->length = end - begin;
  return window;
}
m7p_t * m7p_init_window(m7p_t *P, rci_t begin, rci_t end)
{
  m7p_t* window = ( m1ri_malloc(sizeof(m7p_t)));
  window->values = P->values + begin;
  window->length = end - begin;
  return window;
}

/**
 * \brief Free a permutation window
 * \param condemned Permutation Matrix
 */

void m3p_free_window(m3p_t *condemned)
{
  m1ri_free(condemned->values);
}
void m5p_free_window(m5p_t *condemned)
{
  m1ri_free(condemned->values);
}
void m7p_free_window(m7p_t *condemned)
{
  m1ri_free(condemned->values);
}



static inline void m3d_write_col_to_rows_blockd(m3d_t *A, m3d_t const *B, rci_t const *permutation, vec const *write_mask, rci_t const start_row, rci_t const stop_row, rci_t length)
{
	for(rci_t i = 0; i < length; i += M1RI_RADIX) {
    /* optimisation for identity permutations */
    if (write_mask[i / M1RI_RADIX] == allbits)
    	continue;
    int const todo = MIN(M1RI_RADIX, length - i);
    wi_t const a_vec = i / M1RI_RADIX;
    wi_t vecs[M1RI_RADIX];
    int bits[M1RI_RADIX];
    vec bitmasks[M1RI_RADIX];

    /* we pre-compute bit access in advance */
    for(int k = 0; k < todo; ++k)
    {
        rci_t const colb = permutation[i + k];
        vecs[k] = colb / M1RI_RADIX;
        bits[k] = colb % M1RI_RADIX;
        bitmasks[k] = rightbit << bits[k];
    }

    for (rci_t r = start_row; r < stop_row; ++r)
    {
      vbg const *Brow = B->rows[r-start_row];
      vbg *Arow = A->rows[r];
      register vec value = 0;

      switch(todo-1)
      {
    	  case 63: value |= ((Brow[vecs[63]].units & bitmasks[63]) >> bits[63]) << 63;
	      case 62: value |= ((Brow[vecs[62]].units & bitmasks[62]) >> bits[62]) << 62;
	      case 61: value |= ((Brow[vecs[61]].units & bitmasks[61]) >> bits[61]) << 61;
	      case 60: value |= ((Brow[vecs[60]].units & bitmasks[60]) >> bits[60]) << 60;
	      case 59: value |= ((Brow[vecs[59]].units & bitmasks[59]) >> bits[59]) << 59;
	      case 58: value |= ((Brow[vecs[58]].units & bitmasks[58]) >> bits[58]) << 58;
	      case 57: value |= ((Brow[vecs[57]].units & bitmasks[57]) >> bits[57]) << 57;
	      case 56: value |= ((Brow[vecs[56]].units & bitmasks[56]) >> bits[56]) << 56;
	      case 55: value |= ((Brow[vecs[55]].units & bitmasks[55]) >> bits[55]) << 55;
	      case 54: value |= ((Brow[vecs[54]].units & bitmasks[54]) >> bits[54]) << 54;
	      case 53: value |= ((Brow[vecs[53]].units & bitmasks[53]) >> bits[53]) << 53;
	      case 52: value |= ((Brow[vecs[52]].units & bitmasks[52]) >> bits[52]) << 52;
	      case 51: value |= ((Brow[vecs[51]].units & bitmasks[51]) >> bits[51]) << 51;
	      case 50: value |= ((Brow[vecs[50]].units & bitmasks[50]) >> bits[50]) << 50;
	      case 49: value |= ((Brow[vecs[49]].units & bitmasks[49]) >> bits[49]) << 49;
	      case 48: value |= ((Brow[vecs[48]].units & bitmasks[48]) >> bits[48]) << 48;
	      case 47: value |= ((Brow[vecs[47]].units & bitmasks[47]) >> bits[47]) << 47;
	      case 46: value |= ((Brow[vecs[46]].units & bitmasks[46]) >> bits[46]) << 46;
	      case 45: value |= ((Brow[vecs[45]].units & bitmasks[45]) >> bits[45]) << 45;
	      case 44: value |= ((Brow[vecs[44]].units & bitmasks[44]) >> bits[44]) << 44;
	      case 43: value |= ((Brow[vecs[43]].units & bitmasks[43]) >> bits[43]) << 43;
	      case 42: value |= ((Brow[vecs[42]].units & bitmasks[42]) >> bits[42]) << 42;
	      case 41: value |= ((Brow[vecs[41]].units & bitmasks[41]) >> bits[41]) << 41;
	      case 40: value |= ((Brow[vecs[40]].units & bitmasks[40]) >> bits[40]) << 40;
	      case 39: value |= ((Brow[vecs[39]].units & bitmasks[39]) >> bits[39]) << 39;
	      case 38: value |= ((Brow[vecs[38]].units & bitmasks[38]) >> bits[38]) << 38;
	      case 37: value |= ((Brow[vecs[37]].units & bitmasks[37]) >> bits[37]) << 37;
	      case 36: value |= ((Brow[vecs[36]].units & bitmasks[36]) >> bits[36]) << 36;
	      case 35: value |= ((Brow[vecs[35]].units & bitmasks[35]) >> bits[35]) << 35;
	      case 34: value |= ((Brow[vecs[34]].units & bitmasks[34]) >> bits[34]) << 34;
	      case 33: value |= ((Brow[vecs[33]].units & bitmasks[33]) >> bits[33]) << 33;
	      case 32: value |= ((Brow[vecs[32]].units & bitmasks[32]) >> bits[32]) << 32;
	      case 31: value |= ((Brow[vecs[31]].units & bitmasks[31]) >> bits[31]) << 31;
	      case 30: value |= ((Brow[vecs[30]].units & bitmasks[30]) >> bits[30]) << 30;
	      case 29: value |= ((Brow[vecs[29]].units & bitmasks[29]) >> bits[29]) << 29;
	      case 28: value |= ((Brow[vecs[28]].units & bitmasks[28]) >> bits[28]) << 28;
	      case 27: value |= ((Brow[vecs[27]].units & bitmasks[27]) >> bits[27]) << 27;
	      case 26: value |= ((Brow[vecs[26]].units & bitmasks[26]) >> bits[26]) << 26;
	      case 25: value |= ((Brow[vecs[25]].units & bitmasks[25]) >> bits[25]) << 25;
	      case 24: value |= ((Brow[vecs[24]].units & bitmasks[24]) >> bits[24]) << 24;
	      case 23: value |= ((Brow[vecs[23]].units & bitmasks[23]) >> bits[23]) << 23;
	      case 22: value |= ((Brow[vecs[22]].units & bitmasks[22]) >> bits[22]) << 22;
	      case 21: value |= ((Brow[vecs[21]].units & bitmasks[21]) >> bits[21]) << 21;
	      case 20: value |= ((Brow[vecs[20]].units & bitmasks[20]) >> bits[20]) << 20;
	      case 19: value |= ((Brow[vecs[19]].units & bitmasks[19]) >> bits[19]) << 19;
	      case 18: value |= ((Brow[vecs[18]].units & bitmasks[18]) >> bits[18]) << 18;
	      case 17: value |= ((Brow[vecs[17]].units & bitmasks[17]) >> bits[17]) << 17;
	      case 16: value |= ((Brow[vecs[16]].units & bitmasks[16]) >> bits[16]) << 16;
	      case 15: value |= ((Brow[vecs[15]].units & bitmasks[15]) >> bits[15]) << 15;
	      case 14: value |= ((Brow[vecs[14]].units & bitmasks[14]) >> bits[14]) << 14;
	      case 13: value |= ((Brow[vecs[13]].units & bitmasks[13]) >> bits[13]) << 13;
	      case 12: value |= ((Brow[vecs[12]].units & bitmasks[12]) >> bits[12]) << 12;
	      case 11: value |= ((Brow[vecs[11]].units & bitmasks[11]) >> bits[11]) << 11;
	      case 10: value |= ((Brow[vecs[10]].units & bitmasks[10]) >> bits[10]) << 10;
	      case  9: value |= ((Brow[vecs[ 9]].units & bitmasks[ 9]) >> bits[ 9]) <<  9;
	      case  8: value |= ((Brow[vecs[ 8]].units & bitmasks[ 8]) >> bits[ 8]) <<  8;
	      case  7: value |= ((Brow[vecs[ 7]].units & bitmasks[ 7]) >> bits[ 7]) <<  7;
	      case  6: value |= ((Brow[vecs[ 6]].units & bitmasks[ 6]) >> bits[ 6]) <<  6;
	      case  5: value |= ((Brow[vecs[ 5]].units & bitmasks[ 5]) >> bits[ 5]) <<  5;
	      case  4: value |= ((Brow[vecs[ 4]].units & bitmasks[ 4]) >> bits[ 4]) <<  4;
	      case  3: value |= ((Brow[vecs[ 3]].units & bitmasks[ 3]) >> bits[ 3]) <<  3;
	      case  2: value |= ((Brow[vecs[ 2]].units & bitmasks[ 2]) >> bits[ 2]) <<  2;
	      case  1: value |= ((Brow[vecs[ 1]].units & bitmasks[ 1]) >> bits[ 1]) <<  1;
	      case  0: value |= ((Brow[vecs[ 0]].units & bitmasks[ 0]) >> bits[ 0]) <<  0;
     	 default:
      	  break;
    }
/*       for(int k = 0; k < todo; ++k) { */
/*         value |= ((Brow[vecs[k]].units & bitmasks[k]) << bits[k]) >> k; */
/*       } */
      /* and write the vec once */
      Arow[a_vec].units |= value;

      value = 0;

      /* gathering sign bits*/
      switch(todo-1)
		{
    	case 63: value |= ((Brow[vecs[63]].sign & bitmasks[63]) >> bits[63]) << 63;
	      case 62: value |= ((Brow[vecs[62]].sign & bitmasks[62]) >> bits[62]) << 62;
	      case 61: value |= ((Brow[vecs[61]].sign & bitmasks[61]) >> bits[61]) << 61;
	      case 60: value |= ((Brow[vecs[60]].sign & bitmasks[60]) >> bits[60]) << 60;
	      case 59: value |= ((Brow[vecs[59]].sign & bitmasks[59]) >> bits[59]) << 59;
	      case 58: value |= ((Brow[vecs[58]].sign & bitmasks[58]) >> bits[58]) << 58;
	      case 57: value |= ((Brow[vecs[57]].sign & bitmasks[57]) >> bits[57]) << 57;
	      case 56: value |= ((Brow[vecs[56]].sign & bitmasks[56]) >> bits[56]) << 56;
	      case 55: value |= ((Brow[vecs[55]].sign & bitmasks[55]) >> bits[55]) << 55;
	      case 54: value |= ((Brow[vecs[54]].sign & bitmasks[54]) >> bits[54]) << 54;
	      case 53: value |= ((Brow[vecs[53]].sign & bitmasks[53]) >> bits[53]) << 53;
	      case 52: value |= ((Brow[vecs[52]].sign & bitmasks[52]) >> bits[52]) << 52;
	      case 51: value |= ((Brow[vecs[51]].sign & bitmasks[51]) >> bits[51]) << 51;
	      case 50: value |= ((Brow[vecs[50]].sign & bitmasks[50]) >> bits[50]) << 50;
	      case 49: value |= ((Brow[vecs[49]].sign & bitmasks[49]) >> bits[49]) << 49;
	      case 48: value |= ((Brow[vecs[48]].sign & bitmasks[48]) >> bits[48]) << 48;
	      case 47: value |= ((Brow[vecs[47]].sign & bitmasks[47]) >> bits[47]) << 47;
	      case 46: value |= ((Brow[vecs[46]].sign & bitmasks[46]) >> bits[46]) << 46;
	      case 45: value |= ((Brow[vecs[45]].sign & bitmasks[45]) >> bits[45]) << 45;
	      case 44: value |= ((Brow[vecs[44]].sign & bitmasks[44]) >> bits[44]) << 44;
	      case 43: value |= ((Brow[vecs[43]].sign & bitmasks[43]) >> bits[43]) << 43;
	      case 42: value |= ((Brow[vecs[42]].sign & bitmasks[42]) >> bits[42]) << 42;
	      case 41: value |= ((Brow[vecs[41]].sign & bitmasks[41]) >> bits[41]) << 41;
	      case 40: value |= ((Brow[vecs[40]].sign & bitmasks[40]) >> bits[40]) << 40;
	      case 39: value |= ((Brow[vecs[39]].sign & bitmasks[39]) >> bits[39]) << 39;
	      case 38: value |= ((Brow[vecs[38]].sign & bitmasks[38]) >> bits[38]) << 38;
	      case 37: value |= ((Brow[vecs[37]].sign & bitmasks[37]) >> bits[37]) << 37;
	      case 36: value |= ((Brow[vecs[36]].sign & bitmasks[36]) >> bits[36]) << 36;
	      case 35: value |= ((Brow[vecs[35]].sign & bitmasks[35]) >> bits[35]) << 35;
	      case 34: value |= ((Brow[vecs[34]].sign & bitmasks[34]) >> bits[34]) << 34;
	      case 33: value |= ((Brow[vecs[33]].sign & bitmasks[33]) >> bits[33]) << 33;
	      case 32: value |= ((Brow[vecs[32]].sign & bitmasks[32]) >> bits[32]) << 32;
	      case 31: value |= ((Brow[vecs[31]].sign & bitmasks[31]) >> bits[31]) << 31;
	      case 30: value |= ((Brow[vecs[30]].sign & bitmasks[30]) >> bits[30]) << 30;
	      case 29: value |= ((Brow[vecs[29]].sign & bitmasks[29]) >> bits[29]) << 29;
	      case 28: value |= ((Brow[vecs[28]].sign & bitmasks[28]) >> bits[28]) << 28;
	      case 27: value |= ((Brow[vecs[27]].sign & bitmasks[27]) >> bits[27]) << 27;
	      case 26: value |= ((Brow[vecs[26]].sign & bitmasks[26]) >> bits[26]) << 26;
	      case 25: value |= ((Brow[vecs[25]].sign & bitmasks[25]) >> bits[25]) << 25;
	      case 24: value |= ((Brow[vecs[24]].sign & bitmasks[24]) >> bits[24]) << 24;
	      case 23: value |= ((Brow[vecs[23]].sign & bitmasks[23]) >> bits[23]) << 23;
	      case 22: value |= ((Brow[vecs[22]].sign & bitmasks[22]) >> bits[22]) << 22;
	      case 21: value |= ((Brow[vecs[21]].sign & bitmasks[21]) >> bits[21]) << 21;
	      case 20: value |= ((Brow[vecs[20]].sign & bitmasks[20]) >> bits[20]) << 20;
	      case 19: value |= ((Brow[vecs[19]].sign & bitmasks[19]) >> bits[19]) << 19;
	      case 18: value |= ((Brow[vecs[18]].sign & bitmasks[18]) >> bits[18]) << 18;
	      case 17: value |= ((Brow[vecs[17]].sign & bitmasks[17]) >> bits[17]) << 17;
	      case 16: value |= ((Brow[vecs[16]].sign & bitmasks[16]) >> bits[16]) << 16;
	      case 15: value |= ((Brow[vecs[15]].sign & bitmasks[15]) >> bits[15]) << 15;
	      case 14: value |= ((Brow[vecs[14]].sign & bitmasks[14]) >> bits[14]) << 14;
	      case 13: value |= ((Brow[vecs[13]].sign & bitmasks[13]) >> bits[13]) << 13;
	      case 12: value |= ((Brow[vecs[12]].sign & bitmasks[12]) >> bits[12]) << 12;
	      case 11: value |= ((Brow[vecs[11]].sign & bitmasks[11]) >> bits[11]) << 11;
	      case 10: value |= ((Brow[vecs[10]].sign & bitmasks[10]) >> bits[10]) << 10;
	      case  9: value |= ((Brow[vecs[ 9]].sign & bitmasks[ 9]) >> bits[ 9]) <<  9;
	      case  8: value |= ((Brow[vecs[ 8]].sign & bitmasks[ 8]) >> bits[ 8]) <<  8;
	      case  7: value |= ((Brow[vecs[ 7]].sign & bitmasks[ 7]) >> bits[ 7]) <<  7;
	      case  6: value |= ((Brow[vecs[ 6]].sign & bitmasks[ 6]) >> bits[ 6]) <<  6;
	      case  5: value |= ((Brow[vecs[ 5]].sign & bitmasks[ 5]) >> bits[ 5]) <<  5;
	      case  4: value |= ((Brow[vecs[ 4]].sign & bitmasks[ 4]) >> bits[ 4]) <<  4;
	      case  3: value |= ((Brow[vecs[ 3]].sign & bitmasks[ 3]) >> bits[ 3]) <<  3;
	      case  2: value |= ((Brow[vecs[ 2]].sign & bitmasks[ 2]) >> bits[ 2]) <<  2;
	      case  1: value |= ((Brow[vecs[ 1]].sign & bitmasks[ 1]) >> bits[ 1]) <<  1;
	      case  0: value |= ((Brow[vecs[ 0]].sign & bitmasks[ 0]) >> bits[ 0]) <<  0;
      default:
        break;
      }

		Arow[a_vec].sign |= value;

    	}
	}

}

void m3d_apply_p_right_even(m3d_t *A, m3p_t const *P, rci_t start_row, rci_t start_col, int notrans)
{
  	if(A->nrows - start_row == 0)
   		return;
  	rci_t const length = MIN(P->length, A->ncols);
  	wi_t const width = A->width;
  	int step_size = MIN(A->nrows - start_row, 4096);

 	/* our temporary where we store the columns we want to swap around */
  	m3d_t *B = m3d_create(step_size, A->ncols);
  	vbg *Arow;
  	vbg *Brow;

	  /* setup mathematical permutation */
  	rci_t *permutation = (rci_t*)m1ri_calloc(A->ncols, sizeof(rci_t));
  	for(rci_t i = 0; i < A->ncols; ++i)
    	permutation[i] = i;

	if (!notrans)
	{
    	for(rci_t i = start_col; i < length; ++i)
    	{
      	rci_t t = permutation[i];
      	permutation[i] = permutation[P->values[i]];
      	permutation[P->values[i]] = t;
    	}
  	}
  	else {

    for(rci_t i = start_col; i < length; ++i)
    {
    	rci_t t = permutation[length - i - 1];
      	permutation[length - i - 1] = permutation[P->values[length - i - 1]];
     	 permutation[P->values[length - i - 1]] = t;
    }
  }

  /* we have a bitmask to encode where to write to */
  vec *write_mask = (vec*)m1ri_calloc(width, sizeof(vec));

  for(rci_t i = 0; i < A->ncols; i += M1RI_RADIX)
  {
	int const todo = MIN(M1RI_RADIX, A->ncols - i);
    for(int k = 0; k < todo; ++k)
    {
      if(permutation[i + k] == i + k)
      {
        write_mask[i / M1RI_RADIX] |= rightbit << k;
      }
    }

  }
  write_mask[width-1] |=  ~(allbits >> (	M1RI_RADIX - (A->ncols)) % M1RI_RADIX);

  	for(rci_t i = start_row; i < A->nrows; i += step_size)
  	{
   		step_size = MIN(step_size, A->nrows - i);
    	for(int k = 0; k < step_size; ++k)
    	{
	    	Arow = A->rows[i+k];
      		Brow = B->rows[k];


      /*copy row & clear those values which will be overwritten */
		for(wi_t j = 0; j < width; ++j)
		{
        	Brow[j].units = Arow[j].units;
        	Brow[j].sign  = Arow[j].sign;

        	Arow[j].units = Arow[j].units & write_mask[j];
        	Arow[j].sign  = Arow[j].sign & write_mask[j];

      	}
    }

    	/* write out the permutation */
    	m3d_write_col_to_rows_blockd(A, B, permutation, write_mask, i, i + step_size, length);
  }

  m1ri_free(permutation);
  m1ri_free(write_mask);
  m3d_free(B);

}


static inline void m5d_write_col_to_rows_blockd(m5d_t *A, m5d_t const *B, rci_t const *permutation, vec const *write_mask, rci_t const start_row, rci_t const stop_row, rci_t length)
{
	for(rci_t i = 0; i < length; i += M1RI_RADIX) {
    /* optimisation for identity permutations */
    if (write_mask[i / M1RI_RADIX] == allbits)
    	continue;
    int const todo = MIN(M1RI_RADIX, length - i);
    wi_t const a_vec = i / M1RI_RADIX;
    wi_t vecs[M1RI_RADIX];
    int bits[M1RI_RADIX];
    vec bitmasks[M1RI_RADIX];

    /* we pre-compute bit access in advance */
    for(int k = 0; k < todo; ++k)
    {
        rci_t const colb = permutation[i + k];
        vecs[k] = colb / M1RI_RADIX;
        bits[k] = colb % M1RI_RADIX;
        bitmasks[k] = rightbit << bits[k];
    }

    for (rci_t r = start_row; r < stop_row; ++r)
    {
      vfd const *Brow = B->rows[r-start_row];
      vfd *Arow = A->rows[r];
      register vec value = 0;

      switch(todo-1)
      {
    	  case 63: value |= ((Brow[vecs[63]].units & bitmasks[63]) >> bits[63]) << 63;
	      case 62: value |= ((Brow[vecs[62]].units & bitmasks[62]) >> bits[62]) << 62;
	      case 61: value |= ((Brow[vecs[61]].units & bitmasks[61]) >> bits[61]) << 61;
	      case 60: value |= ((Brow[vecs[60]].units & bitmasks[60]) >> bits[60]) << 60;
	      case 59: value |= ((Brow[vecs[59]].units & bitmasks[59]) >> bits[59]) << 59;
	      case 58: value |= ((Brow[vecs[58]].units & bitmasks[58]) >> bits[58]) << 58;
	      case 57: value |= ((Brow[vecs[57]].units & bitmasks[57]) >> bits[57]) << 57;
	      case 56: value |= ((Brow[vecs[56]].units & bitmasks[56]) >> bits[56]) << 56;
	      case 55: value |= ((Brow[vecs[55]].units & bitmasks[55]) >> bits[55]) << 55;
	      case 54: value |= ((Brow[vecs[54]].units & bitmasks[54]) >> bits[54]) << 54;
	      case 53: value |= ((Brow[vecs[53]].units & bitmasks[53]) >> bits[53]) << 53;
	      case 52: value |= ((Brow[vecs[52]].units & bitmasks[52]) >> bits[52]) << 52;
	      case 51: value |= ((Brow[vecs[51]].units & bitmasks[51]) >> bits[51]) << 51;
	      case 50: value |= ((Brow[vecs[50]].units & bitmasks[50]) >> bits[50]) << 50;
	      case 49: value |= ((Brow[vecs[49]].units & bitmasks[49]) >> bits[49]) << 49;
	      case 48: value |= ((Brow[vecs[48]].units & bitmasks[48]) >> bits[48]) << 48;
	      case 47: value |= ((Brow[vecs[47]].units & bitmasks[47]) >> bits[47]) << 47;
	      case 46: value |= ((Brow[vecs[46]].units & bitmasks[46]) >> bits[46]) << 46;
	      case 45: value |= ((Brow[vecs[45]].units & bitmasks[45]) >> bits[45]) << 45;
	      case 44: value |= ((Brow[vecs[44]].units & bitmasks[44]) >> bits[44]) << 44;
	      case 43: value |= ((Brow[vecs[43]].units & bitmasks[43]) >> bits[43]) << 43;
	      case 42: value |= ((Brow[vecs[42]].units & bitmasks[42]) >> bits[42]) << 42;
	      case 41: value |= ((Brow[vecs[41]].units & bitmasks[41]) >> bits[41]) << 41;
	      case 40: value |= ((Brow[vecs[40]].units & bitmasks[40]) >> bits[40]) << 40;
	      case 39: value |= ((Brow[vecs[39]].units & bitmasks[39]) >> bits[39]) << 39;
	      case 38: value |= ((Brow[vecs[38]].units & bitmasks[38]) >> bits[38]) << 38;
	      case 37: value |= ((Brow[vecs[37]].units & bitmasks[37]) >> bits[37]) << 37;
	      case 36: value |= ((Brow[vecs[36]].units & bitmasks[36]) >> bits[36]) << 36;
	      case 35: value |= ((Brow[vecs[35]].units & bitmasks[35]) >> bits[35]) << 35;
	      case 34: value |= ((Brow[vecs[34]].units & bitmasks[34]) >> bits[34]) << 34;
	      case 33: value |= ((Brow[vecs[33]].units & bitmasks[33]) >> bits[33]) << 33;
	      case 32: value |= ((Brow[vecs[32]].units & bitmasks[32]) >> bits[32]) << 32;
	      case 31: value |= ((Brow[vecs[31]].units & bitmasks[31]) >> bits[31]) << 31;
	      case 30: value |= ((Brow[vecs[30]].units & bitmasks[30]) >> bits[30]) << 30;
	      case 29: value |= ((Brow[vecs[29]].units & bitmasks[29]) >> bits[29]) << 29;
	      case 28: value |= ((Brow[vecs[28]].units & bitmasks[28]) >> bits[28]) << 28;
	      case 27: value |= ((Brow[vecs[27]].units & bitmasks[27]) >> bits[27]) << 27;
	      case 26: value |= ((Brow[vecs[26]].units & bitmasks[26]) >> bits[26]) << 26;
	      case 25: value |= ((Brow[vecs[25]].units & bitmasks[25]) >> bits[25]) << 25;
	      case 24: value |= ((Brow[vecs[24]].units & bitmasks[24]) >> bits[24]) << 24;
	      case 23: value |= ((Brow[vecs[23]].units & bitmasks[23]) >> bits[23]) << 23;
	      case 22: value |= ((Brow[vecs[22]].units & bitmasks[22]) >> bits[22]) << 22;
	      case 21: value |= ((Brow[vecs[21]].units & bitmasks[21]) >> bits[21]) << 21;
	      case 20: value |= ((Brow[vecs[20]].units & bitmasks[20]) >> bits[20]) << 20;
	      case 19: value |= ((Brow[vecs[19]].units & bitmasks[19]) >> bits[19]) << 19;
	      case 18: value |= ((Brow[vecs[18]].units & bitmasks[18]) >> bits[18]) << 18;
	      case 17: value |= ((Brow[vecs[17]].units & bitmasks[17]) >> bits[17]) << 17;
	      case 16: value |= ((Brow[vecs[16]].units & bitmasks[16]) >> bits[16]) << 16;
	      case 15: value |= ((Brow[vecs[15]].units & bitmasks[15]) >> bits[15]) << 15;
	      case 14: value |= ((Brow[vecs[14]].units & bitmasks[14]) >> bits[14]) << 14;
	      case 13: value |= ((Brow[vecs[13]].units & bitmasks[13]) >> bits[13]) << 13;
	      case 12: value |= ((Brow[vecs[12]].units & bitmasks[12]) >> bits[12]) << 12;
	      case 11: value |= ((Brow[vecs[11]].units & bitmasks[11]) >> bits[11]) << 11;
	      case 10: value |= ((Brow[vecs[10]].units & bitmasks[10]) >> bits[10]) << 10;
	      case  9: value |= ((Brow[vecs[ 9]].units & bitmasks[ 9]) >> bits[ 9]) <<  9;
	      case  8: value |= ((Brow[vecs[ 8]].units & bitmasks[ 8]) >> bits[ 8]) <<  8;
	      case  7: value |= ((Brow[vecs[ 7]].units & bitmasks[ 7]) >> bits[ 7]) <<  7;
	      case  6: value |= ((Brow[vecs[ 6]].units & bitmasks[ 6]) >> bits[ 6]) <<  6;
	      case  5: value |= ((Brow[vecs[ 5]].units & bitmasks[ 5]) >> bits[ 5]) <<  5;
	      case  4: value |= ((Brow[vecs[ 4]].units & bitmasks[ 4]) >> bits[ 4]) <<  4;
	      case  3: value |= ((Brow[vecs[ 3]].units & bitmasks[ 3]) >> bits[ 3]) <<  3;
	      case  2: value |= ((Brow[vecs[ 2]].units & bitmasks[ 2]) >> bits[ 2]) <<  2;
	      case  1: value |= ((Brow[vecs[ 1]].units & bitmasks[ 1]) >> bits[ 1]) <<  1;
	      case  0: value |= ((Brow[vecs[ 0]].units & bitmasks[ 0]) >> bits[ 0]) <<  0;
     	 default:
      	  break;
    }
/*       for(int k = 0; k < todo; ++k) { */
/*         value |= ((Brow[vecs[k]].units & bitmasks[k]) << bits[k]) >> k; */
/*       } */
      /* and write the vec once */
      Arow[a_vec].units |= value;

      value = 0;
      switch(todo-1)
       {
		  case 63: value |= ((Brow[vecs[63]].middle& bitmasks[63]) >> bits[63]) << 63;
	      case 62: value |= ((Brow[vecs[62]].middle& bitmasks[62]) >> bits[62]) << 62;
	      case 61: value |= ((Brow[vecs[61]].middle& bitmasks[61]) >> bits[61]) << 61;
	      case 60: value |= ((Brow[vecs[60]].middle& bitmasks[60]) >> bits[60]) << 60;
	      case 59: value |= ((Brow[vecs[59]].middle& bitmasks[59]) >> bits[59]) << 59;
	      case 58: value |= ((Brow[vecs[58]].middle& bitmasks[58]) >> bits[58]) << 58;
	      case 57: value |= ((Brow[vecs[57]].middle& bitmasks[57]) >> bits[57]) << 57;
	      case 56: value |= ((Brow[vecs[56]].middle& bitmasks[56]) >> bits[56]) << 56;
	      case 55: value |= ((Brow[vecs[55]].middle& bitmasks[55]) >> bits[55]) << 55;
	      case 54: value |= ((Brow[vecs[54]].middle& bitmasks[54]) >> bits[54]) << 54;
	      case 53: value |= ((Brow[vecs[53]].middle& bitmasks[53]) >> bits[53]) << 53;
	      case 52: value |= ((Brow[vecs[52]].middle& bitmasks[52]) >> bits[52]) << 52;
	      case 51: value |= ((Brow[vecs[51]].middle& bitmasks[51]) >> bits[51]) << 51;
	      case 50: value |= ((Brow[vecs[50]].middle& bitmasks[50]) >> bits[50]) << 50;
	      case 49: value |= ((Brow[vecs[49]].middle& bitmasks[49]) >> bits[49]) << 49;
	      case 48: value |= ((Brow[vecs[48]].middle& bitmasks[48]) >> bits[48]) << 48;
	      case 47: value |= ((Brow[vecs[47]].middle& bitmasks[47]) >> bits[47]) << 47;
	      case 46: value |= ((Brow[vecs[46]].middle& bitmasks[46]) >> bits[46]) << 46;
	      case 45: value |= ((Brow[vecs[45]].middle& bitmasks[45]) >> bits[45]) << 45;
	      case 44: value |= ((Brow[vecs[44]].middle& bitmasks[44]) >> bits[44]) << 44;
	      case 43: value |= ((Brow[vecs[43]].middle& bitmasks[43]) >> bits[43]) << 43;
	      case 42: value |= ((Brow[vecs[42]].middle& bitmasks[42]) >> bits[42]) << 42;
	      case 41: value |= ((Brow[vecs[41]].middle& bitmasks[41]) >> bits[41]) << 41;
	      case 40: value |= ((Brow[vecs[40]].middle& bitmasks[40]) >> bits[40]) << 40;
	      case 39: value |= ((Brow[vecs[39]].middle& bitmasks[39]) >> bits[39]) << 39;
	      case 38: value |= ((Brow[vecs[38]].middle& bitmasks[38]) >> bits[38]) << 38;
	      case 37: value |= ((Brow[vecs[37]].middle& bitmasks[37]) >> bits[37]) << 37;
	      case 36: value |= ((Brow[vecs[36]].middle& bitmasks[36]) >> bits[36]) << 36;
	      case 35: value |= ((Brow[vecs[35]].middle& bitmasks[35]) >> bits[35]) << 35;
	      case 34: value |= ((Brow[vecs[34]].middle& bitmasks[34]) >> bits[34]) << 34;
	      case 33: value |= ((Brow[vecs[33]].middle& bitmasks[33]) >> bits[33]) << 33;
	      case 32: value |= ((Brow[vecs[32]].middle& bitmasks[32]) >> bits[32]) << 32;
	      case 31: value |= ((Brow[vecs[31]].middle& bitmasks[31]) >> bits[31]) << 31;
	      case 30: value |= ((Brow[vecs[30]].middle& bitmasks[30]) >> bits[30]) << 30;
	      case 29: value |= ((Brow[vecs[29]].middle& bitmasks[29]) >> bits[29]) << 29;
	      case 28: value |= ((Brow[vecs[28]].middle& bitmasks[28]) >> bits[28]) << 28;
	      case 27: value |= ((Brow[vecs[27]].middle& bitmasks[27]) >> bits[27]) << 27;
	      case 26: value |= ((Brow[vecs[26]].middle& bitmasks[26]) >> bits[26]) << 26;
	      case 25: value |= ((Brow[vecs[25]].middle& bitmasks[25]) >> bits[25]) << 25;
	      case 24: value |= ((Brow[vecs[24]].middle& bitmasks[24]) >> bits[24]) << 24;
	      case 23: value |= ((Brow[vecs[23]].middle& bitmasks[23]) >> bits[23]) << 23;
	      case 22: value |= ((Brow[vecs[22]].middle& bitmasks[22]) >> bits[22]) << 22;
	      case 21: value |= ((Brow[vecs[21]].middle& bitmasks[21]) >> bits[21]) << 21;
	      case 20: value |= ((Brow[vecs[20]].middle& bitmasks[20]) >> bits[20]) << 20;
	      case 19: value |= ((Brow[vecs[19]].middle& bitmasks[19]) >> bits[19]) << 19;
	      case 18: value |= ((Brow[vecs[18]].middle& bitmasks[18]) >> bits[18]) << 18;
	      case 17: value |= ((Brow[vecs[17]].middle& bitmasks[17]) >> bits[17]) << 17;
	      case 16: value |= ((Brow[vecs[16]].middle& bitmasks[16]) >> bits[16]) << 16;
	      case 15: value |= ((Brow[vecs[15]].middle& bitmasks[15]) >> bits[15]) << 15;
	      case 14: value |= ((Brow[vecs[14]].middle& bitmasks[14]) >> bits[14]) << 14;
	      case 13: value |= ((Brow[vecs[13]].middle& bitmasks[13]) >> bits[13]) << 13;
	      case 12: value |= ((Brow[vecs[12]].middle& bitmasks[12]) >> bits[12]) << 12;
	      case 11: value |= ((Brow[vecs[11]].middle& bitmasks[11]) >> bits[11]) << 11;
	      case 10: value |= ((Brow[vecs[10]].middle& bitmasks[10]) >> bits[10]) << 10;
	      case  9: value |= ((Brow[vecs[ 9]].middle& bitmasks[ 9]) >> bits[ 9]) <<  9;
	      case  8: value |= ((Brow[vecs[ 8]].middle& bitmasks[ 8]) >> bits[ 8]) <<  8;
	      case  7: value |= ((Brow[vecs[ 7]].middle& bitmasks[ 7]) >> bits[ 7]) <<  7;
	      case  6: value |= ((Brow[vecs[ 6]].middle& bitmasks[ 6]) >> bits[ 6]) <<  6;
	      case  5: value |= ((Brow[vecs[ 5]].middle& bitmasks[ 5]) >> bits[ 5]) <<  5;
	      case  4: value |= ((Brow[vecs[ 4]].middle& bitmasks[ 4]) >> bits[ 4]) <<  4;
	      case  3: value |= ((Brow[vecs[ 3]].middle& bitmasks[ 3]) >> bits[ 3]) <<  3;
	      case  2: value |= ((Brow[vecs[ 2]].middle& bitmasks[ 2]) >> bits[ 2]) <<  2;
	      case  1: value |= ((Brow[vecs[ 1]].middle& bitmasks[ 1]) >> bits[ 1]) <<  1;
	      case  0: value |= ((Brow[vecs[ 0]].middle& bitmasks[ 0]) >> bits[ 0]) <<  0;
     	 default:
      	  break;
      }
      /* gathering sign bits*/

      Arow[a_vec].middle |= value;

      value = 0;
      switch(todo-1)
		{
    	  case 63: value |= ((Brow[vecs[63]].sign & bitmasks[63]) >> bits[63]) << 63;
	      case 62: value |= ((Brow[vecs[62]].sign & bitmasks[62]) >> bits[62]) << 62;
	      case 61: value |= ((Brow[vecs[61]].sign & bitmasks[61]) >> bits[61]) << 61;
	      case 60: value |= ((Brow[vecs[60]].sign & bitmasks[60]) >> bits[60]) << 60;
	      case 59: value |= ((Brow[vecs[59]].sign & bitmasks[59]) >> bits[59]) << 59;
	      case 58: value |= ((Brow[vecs[58]].sign & bitmasks[58]) >> bits[58]) << 58;
	      case 57: value |= ((Brow[vecs[57]].sign & bitmasks[57]) >> bits[57]) << 57;
	      case 56: value |= ((Brow[vecs[56]].sign & bitmasks[56]) >> bits[56]) << 56;
	      case 55: value |= ((Brow[vecs[55]].sign & bitmasks[55]) >> bits[55]) << 55;
	      case 54: value |= ((Brow[vecs[54]].sign & bitmasks[54]) >> bits[54]) << 54;
	      case 53: value |= ((Brow[vecs[53]].sign & bitmasks[53]) >> bits[53]) << 53;
	      case 52: value |= ((Brow[vecs[52]].sign & bitmasks[52]) >> bits[52]) << 52;
	      case 51: value |= ((Brow[vecs[51]].sign & bitmasks[51]) >> bits[51]) << 51;
	      case 50: value |= ((Brow[vecs[50]].sign & bitmasks[50]) >> bits[50]) << 50;
	      case 49: value |= ((Brow[vecs[49]].sign & bitmasks[49]) >> bits[49]) << 49;
	      case 48: value |= ((Brow[vecs[48]].sign & bitmasks[48]) >> bits[48]) << 48;
	      case 47: value |= ((Brow[vecs[47]].sign & bitmasks[47]) >> bits[47]) << 47;
	      case 46: value |= ((Brow[vecs[46]].sign & bitmasks[46]) >> bits[46]) << 46;
	      case 45: value |= ((Brow[vecs[45]].sign & bitmasks[45]) >> bits[45]) << 45;
	      case 44: value |= ((Brow[vecs[44]].sign & bitmasks[44]) >> bits[44]) << 44;
	      case 43: value |= ((Brow[vecs[43]].sign & bitmasks[43]) >> bits[43]) << 43;
	      case 42: value |= ((Brow[vecs[42]].sign & bitmasks[42]) >> bits[42]) << 42;
	      case 41: value |= ((Brow[vecs[41]].sign & bitmasks[41]) >> bits[41]) << 41;
	      case 40: value |= ((Brow[vecs[40]].sign & bitmasks[40]) >> bits[40]) << 40;
	      case 39: value |= ((Brow[vecs[39]].sign & bitmasks[39]) >> bits[39]) << 39;
	      case 38: value |= ((Brow[vecs[38]].sign & bitmasks[38]) >> bits[38]) << 38;
	      case 37: value |= ((Brow[vecs[37]].sign & bitmasks[37]) >> bits[37]) << 37;
	      case 36: value |= ((Brow[vecs[36]].sign & bitmasks[36]) >> bits[36]) << 36;
	      case 35: value |= ((Brow[vecs[35]].sign & bitmasks[35]) >> bits[35]) << 35;
	      case 34: value |= ((Brow[vecs[34]].sign & bitmasks[34]) >> bits[34]) << 34;
	      case 33: value |= ((Brow[vecs[33]].sign & bitmasks[33]) >> bits[33]) << 33;
	      case 32: value |= ((Brow[vecs[32]].sign & bitmasks[32]) >> bits[32]) << 32;
	      case 31: value |= ((Brow[vecs[31]].sign & bitmasks[31]) >> bits[31]) << 31;
	      case 30: value |= ((Brow[vecs[30]].sign & bitmasks[30]) >> bits[30]) << 30;
	      case 29: value |= ((Brow[vecs[29]].sign & bitmasks[29]) >> bits[29]) << 29;
	      case 28: value |= ((Brow[vecs[28]].sign & bitmasks[28]) >> bits[28]) << 28;
	      case 27: value |= ((Brow[vecs[27]].sign & bitmasks[27]) >> bits[27]) << 27;
	      case 26: value |= ((Brow[vecs[26]].sign & bitmasks[26]) >> bits[26]) << 26;
	      case 25: value |= ((Brow[vecs[25]].sign & bitmasks[25]) >> bits[25]) << 25;
	      case 24: value |= ((Brow[vecs[24]].sign & bitmasks[24]) >> bits[24]) << 24;
	      case 23: value |= ((Brow[vecs[23]].sign & bitmasks[23]) >> bits[23]) << 23;
	      case 22: value |= ((Brow[vecs[22]].sign & bitmasks[22]) >> bits[22]) << 22;
	      case 21: value |= ((Brow[vecs[21]].sign & bitmasks[21]) >> bits[21]) << 21;
	      case 20: value |= ((Brow[vecs[20]].sign & bitmasks[20]) >> bits[20]) << 20;
	      case 19: value |= ((Brow[vecs[19]].sign & bitmasks[19]) >> bits[19]) << 19;
	      case 18: value |= ((Brow[vecs[18]].sign & bitmasks[18]) >> bits[18]) << 18;
	      case 17: value |= ((Brow[vecs[17]].sign & bitmasks[17]) >> bits[17]) << 17;
	      case 16: value |= ((Brow[vecs[16]].sign & bitmasks[16]) >> bits[16]) << 16;
	      case 15: value |= ((Brow[vecs[15]].sign & bitmasks[15]) >> bits[15]) << 15;
	      case 14: value |= ((Brow[vecs[14]].sign & bitmasks[14]) >> bits[14]) << 14;
	      case 13: value |= ((Brow[vecs[13]].sign & bitmasks[13]) >> bits[13]) << 13;
	      case 12: value |= ((Brow[vecs[12]].sign & bitmasks[12]) >> bits[12]) << 12;
	      case 11: value |= ((Brow[vecs[11]].sign & bitmasks[11]) >> bits[11]) << 11;
	      case 10: value |= ((Brow[vecs[10]].sign & bitmasks[10]) >> bits[10]) << 10;
	      case  9: value |= ((Brow[vecs[ 9]].sign & bitmasks[ 9]) >> bits[ 9]) <<  9;
	      case  8: value |= ((Brow[vecs[ 8]].sign & bitmasks[ 8]) >> bits[ 8]) <<  8;
	      case  7: value |= ((Brow[vecs[ 7]].sign & bitmasks[ 7]) >> bits[ 7]) <<  7;
	      case  6: value |= ((Brow[vecs[ 6]].sign & bitmasks[ 6]) >> bits[ 6]) <<  6;
	      case  5: value |= ((Brow[vecs[ 5]].sign & bitmasks[ 5]) >> bits[ 5]) <<  5;
	      case  4: value |= ((Brow[vecs[ 4]].sign & bitmasks[ 4]) >> bits[ 4]) <<  4;
	      case  3: value |= ((Brow[vecs[ 3]].sign & bitmasks[ 3]) >> bits[ 3]) <<  3;
	      case  2: value |= ((Brow[vecs[ 2]].sign & bitmasks[ 2]) >> bits[ 2]) <<  2;
	      case  1: value |= ((Brow[vecs[ 1]].sign & bitmasks[ 1]) >> bits[ 1]) <<  1;
	      case  0: value |= ((Brow[vecs[ 0]].sign & bitmasks[ 0]) >> bits[ 0]) <<  0;
      default:
        break;
      }

		Arow[a_vec].sign |= value;

/*       for(int k = 0; k < todo; ++k) { */
/*         value |= ((Brow[vecs[k]].units & bitmasks[k]) << bits[k]) >> k; */
/*       } */
      /* and write the vec once */


    }
  }


}




void m5d_apply_p_right_even(m5d_t *A, m5p_t const *P, rci_t start_row, rci_t start_col, int notrans)
{
  	if(A->nrows - start_row == 0)
   		return;
  	rci_t const length = MIN(P->length, A->ncols);
  	wi_t const width = A->width;
  	int step_size = MIN(A->nrows - start_row, 4096);

 	/* our temporary where we store the columns we want to swap around */
  	m5d_t *B = m5d_create(step_size, A->ncols);
  	vfd *Arow;
  	vfd *Brow;

	  /* setup mathematical permutation */
  	rci_t *permutation = (rci_t*)m1ri_calloc(A->ncols, sizeof(rci_t));
  	for(rci_t i = 0; i < A->ncols; ++i)
    	permutation[i] = i;

	if (!notrans)
	{
    	for(rci_t i = start_col; i < length; ++i)
    	{
      	rci_t t = permutation[i];
      	permutation[i] = permutation[P->values[i]];
      	permutation[P->values[i]] = t;
    	}
  	}
  	else {

    for(rci_t i = start_col; i < length; ++i)
    {
    	rci_t t = permutation[length - i - 1];
      	permutation[length - i - 1] = permutation[P->values[length - i - 1]];
     	 permutation[P->values[length - i - 1]] = t;
    }
  }

  /* we have a bitmask to encode where to write to */
  vec *write_mask = (vec*)m1ri_calloc(width, sizeof(vec));

  for(rci_t i = 0; i < A->ncols; i += M1RI_RADIX)
  {
	int const todo = MIN(M1RI_RADIX, A->ncols - i);
    for(int k = 0; k < todo; ++k)
    {
      if(permutation[i + k] == i + k)
      {
        write_mask[i / M1RI_RADIX] |= rightbit << k;
      }
    }

  }
  write_mask[width-1] |=  ~(allbits >> (	M1RI_RADIX - (A->ncols)) % M1RI_RADIX);

  	for(rci_t i = start_row; i < A->nrows; i += step_size)
  	{
   		step_size = MIN(step_size, A->nrows - i);
    	for(int k = 0; k < step_size; ++k)
    	{
	    	Arow = A->rows[i+k];
      		Brow = B->rows[k];


      /*copy row & clear those values which will be overwritten */
		for(wi_t j = 0; j < width; ++j)
		{
        	Brow[j].units = Arow[j].units;
        	Brow[j].sign  = Arow[j].sign;

        	Arow[j].units = Arow[j].units & write_mask[j];
        	Arow[j].sign  = Arow[j].sign & write_mask[j];

      	}
    }

    	/* write out the permutation */
    	m5d_write_col_to_rows_blockd(A, B, permutation, write_mask, i, i + step_size, length);
  }

  m1ri_free(permutation);
  m1ri_free(write_mask);
  m5d_free(B);

}



static inline void m7d_write_col_to_rows_blockd(m7d_t *A, m7d_t const *B, rci_t const *permutation, vec const *write_mask, rci_t const start_row, rci_t const stop_row, rci_t length)
{
	for(rci_t i = 0; i < length; i += M1RI_RADIX) {
    /* optimisation for identity permutations */
    if (write_mask[i / M1RI_RADIX] == allbits)
    	continue;
    int const todo = MIN(M1RI_RADIX, length - i);
    wi_t const a_vec = i / M1RI_RADIX;
    wi_t vecs[M1RI_RADIX];
    int bits[M1RI_RADIX];
    vec bitmasks[M1RI_RADIX];

    /* we pre-compute bit access in advance */
    for(int k = 0; k < todo; ++k)
    {
        rci_t const colb = permutation[i + k];
        vecs[k] = colb / M1RI_RADIX;
        bits[k] = colb % M1RI_RADIX;
        bitmasks[k] = rightbit << bits[k];
    }

    for (rci_t r = start_row; r < stop_row; ++r)
    {
      vtri const *Brow = B->rows[r-start_row];
      vtri *Arow = A->rows[r];
      register vec value = 0;

      switch(todo-1)
      {
    	  case 63: value |= ((Brow[vecs[63]].units & bitmasks[63]) >> bits[63]) << 63;
	      case 62: value |= ((Brow[vecs[62]].units & bitmasks[62]) >> bits[62]) << 62;
	      case 61: value |= ((Brow[vecs[61]].units & bitmasks[61]) >> bits[61]) << 61;
	      case 60: value |= ((Brow[vecs[60]].units & bitmasks[60]) >> bits[60]) << 60;
	      case 59: value |= ((Brow[vecs[59]].units & bitmasks[59]) >> bits[59]) << 59;
	      case 58: value |= ((Brow[vecs[58]].units & bitmasks[58]) >> bits[58]) << 58;
	      case 57: value |= ((Brow[vecs[57]].units & bitmasks[57]) >> bits[57]) << 57;
	      case 56: value |= ((Brow[vecs[56]].units & bitmasks[56]) >> bits[56]) << 56;
	      case 55: value |= ((Brow[vecs[55]].units & bitmasks[55]) >> bits[55]) << 55;
	      case 54: value |= ((Brow[vecs[54]].units & bitmasks[54]) >> bits[54]) << 54;
	      case 53: value |= ((Brow[vecs[53]].units & bitmasks[53]) >> bits[53]) << 53;
	      case 52: value |= ((Brow[vecs[52]].units & bitmasks[52]) >> bits[52]) << 52;
	      case 51: value |= ((Brow[vecs[51]].units & bitmasks[51]) >> bits[51]) << 51;
	      case 50: value |= ((Brow[vecs[50]].units & bitmasks[50]) >> bits[50]) << 50;
	      case 49: value |= ((Brow[vecs[49]].units & bitmasks[49]) >> bits[49]) << 49;
	      case 48: value |= ((Brow[vecs[48]].units & bitmasks[48]) >> bits[48]) << 48;
	      case 47: value |= ((Brow[vecs[47]].units & bitmasks[47]) >> bits[47]) << 47;
	      case 46: value |= ((Brow[vecs[46]].units & bitmasks[46]) >> bits[46]) << 46;
	      case 45: value |= ((Brow[vecs[45]].units & bitmasks[45]) >> bits[45]) << 45;
	      case 44: value |= ((Brow[vecs[44]].units & bitmasks[44]) >> bits[44]) << 44;
	      case 43: value |= ((Brow[vecs[43]].units & bitmasks[43]) >> bits[43]) << 43;
	      case 42: value |= ((Brow[vecs[42]].units & bitmasks[42]) >> bits[42]) << 42;
	      case 41: value |= ((Brow[vecs[41]].units & bitmasks[41]) >> bits[41]) << 41;
	      case 40: value |= ((Brow[vecs[40]].units & bitmasks[40]) >> bits[40]) << 40;
	      case 39: value |= ((Brow[vecs[39]].units & bitmasks[39]) >> bits[39]) << 39;
	      case 38: value |= ((Brow[vecs[38]].units & bitmasks[38]) >> bits[38]) << 38;
	      case 37: value |= ((Brow[vecs[37]].units & bitmasks[37]) >> bits[37]) << 37;
	      case 36: value |= ((Brow[vecs[36]].units & bitmasks[36]) >> bits[36]) << 36;
	      case 35: value |= ((Brow[vecs[35]].units & bitmasks[35]) >> bits[35]) << 35;
	      case 34: value |= ((Brow[vecs[34]].units & bitmasks[34]) >> bits[34]) << 34;
	      case 33: value |= ((Brow[vecs[33]].units & bitmasks[33]) >> bits[33]) << 33;
	      case 32: value |= ((Brow[vecs[32]].units & bitmasks[32]) >> bits[32]) << 32;
	      case 31: value |= ((Brow[vecs[31]].units & bitmasks[31]) >> bits[31]) << 31;
	      case 30: value |= ((Brow[vecs[30]].units & bitmasks[30]) >> bits[30]) << 30;
	      case 29: value |= ((Brow[vecs[29]].units & bitmasks[29]) >> bits[29]) << 29;
	      case 28: value |= ((Brow[vecs[28]].units & bitmasks[28]) >> bits[28]) << 28;
	      case 27: value |= ((Brow[vecs[27]].units & bitmasks[27]) >> bits[27]) << 27;
	      case 26: value |= ((Brow[vecs[26]].units & bitmasks[26]) >> bits[26]) << 26;
	      case 25: value |= ((Brow[vecs[25]].units & bitmasks[25]) >> bits[25]) << 25;
	      case 24: value |= ((Brow[vecs[24]].units & bitmasks[24]) >> bits[24]) << 24;
	      case 23: value |= ((Brow[vecs[23]].units & bitmasks[23]) >> bits[23]) << 23;
	      case 22: value |= ((Brow[vecs[22]].units & bitmasks[22]) >> bits[22]) << 22;
	      case 21: value |= ((Brow[vecs[21]].units & bitmasks[21]) >> bits[21]) << 21;
	      case 20: value |= ((Brow[vecs[20]].units & bitmasks[20]) >> bits[20]) << 20;
	      case 19: value |= ((Brow[vecs[19]].units & bitmasks[19]) >> bits[19]) << 19;
	      case 18: value |= ((Brow[vecs[18]].units & bitmasks[18]) >> bits[18]) << 18;
	      case 17: value |= ((Brow[vecs[17]].units & bitmasks[17]) >> bits[17]) << 17;
	      case 16: value |= ((Brow[vecs[16]].units & bitmasks[16]) >> bits[16]) << 16;
	      case 15: value |= ((Brow[vecs[15]].units & bitmasks[15]) >> bits[15]) << 15;
	      case 14: value |= ((Brow[vecs[14]].units & bitmasks[14]) >> bits[14]) << 14;
	      case 13: value |= ((Brow[vecs[13]].units & bitmasks[13]) >> bits[13]) << 13;
	      case 12: value |= ((Brow[vecs[12]].units & bitmasks[12]) >> bits[12]) << 12;
	      case 11: value |= ((Brow[vecs[11]].units & bitmasks[11]) >> bits[11]) << 11;
	      case 10: value |= ((Brow[vecs[10]].units & bitmasks[10]) >> bits[10]) << 10;
	      case  9: value |= ((Brow[vecs[ 9]].units & bitmasks[ 9]) >> bits[ 9]) <<  9;
	      case  8: value |= ((Brow[vecs[ 8]].units & bitmasks[ 8]) >> bits[ 8]) <<  8;
	      case  7: value |= ((Brow[vecs[ 7]].units & bitmasks[ 7]) >> bits[ 7]) <<  7;
	      case  6: value |= ((Brow[vecs[ 6]].units & bitmasks[ 6]) >> bits[ 6]) <<  6;
	      case  5: value |= ((Brow[vecs[ 5]].units & bitmasks[ 5]) >> bits[ 5]) <<  5;
	      case  4: value |= ((Brow[vecs[ 4]].units & bitmasks[ 4]) >> bits[ 4]) <<  4;
	      case  3: value |= ((Brow[vecs[ 3]].units & bitmasks[ 3]) >> bits[ 3]) <<  3;
	      case  2: value |= ((Brow[vecs[ 2]].units & bitmasks[ 2]) >> bits[ 2]) <<  2;
	      case  1: value |= ((Brow[vecs[ 1]].units & bitmasks[ 1]) >> bits[ 1]) <<  1;
	      case  0: value |= ((Brow[vecs[ 0]].units & bitmasks[ 0]) >> bits[ 0]) <<  0;
     	 default:
      	  break;
    }
/*       for(int k = 0; k < todo; ++k) { */
/*         value |= ((Brow[vecs[k]].units & bitmasks[k]) << bits[k]) >> k; */
/*       } */
      /* and write the vec once */
      Arow[a_vec].units |= value;

      value = 0;
      switch(todo-1)
       {
		  case 63: value |= ((Brow[vecs[63]].middle& bitmasks[63]) >> bits[63]) << 63;
	      case 62: value |= ((Brow[vecs[62]].middle& bitmasks[62]) >> bits[62]) << 62;
	      case 61: value |= ((Brow[vecs[61]].middle& bitmasks[61]) >> bits[61]) << 61;
	      case 60: value |= ((Brow[vecs[60]].middle& bitmasks[60]) >> bits[60]) << 60;
	      case 59: value |= ((Brow[vecs[59]].middle& bitmasks[59]) >> bits[59]) << 59;
	      case 58: value |= ((Brow[vecs[58]].middle& bitmasks[58]) >> bits[58]) << 58;
	      case 57: value |= ((Brow[vecs[57]].middle& bitmasks[57]) >> bits[57]) << 57;
	      case 56: value |= ((Brow[vecs[56]].middle& bitmasks[56]) >> bits[56]) << 56;
	      case 55: value |= ((Brow[vecs[55]].middle& bitmasks[55]) >> bits[55]) << 55;
	      case 54: value |= ((Brow[vecs[54]].middle& bitmasks[54]) >> bits[54]) << 54;
	      case 53: value |= ((Brow[vecs[53]].middle& bitmasks[53]) >> bits[53]) << 53;
	      case 52: value |= ((Brow[vecs[52]].middle& bitmasks[52]) >> bits[52]) << 52;
	      case 51: value |= ((Brow[vecs[51]].middle& bitmasks[51]) >> bits[51]) << 51;
	      case 50: value |= ((Brow[vecs[50]].middle& bitmasks[50]) >> bits[50]) << 50;
	      case 49: value |= ((Brow[vecs[49]].middle& bitmasks[49]) >> bits[49]) << 49;
	      case 48: value |= ((Brow[vecs[48]].middle& bitmasks[48]) >> bits[48]) << 48;
	      case 47: value |= ((Brow[vecs[47]].middle& bitmasks[47]) >> bits[47]) << 47;
	      case 46: value |= ((Brow[vecs[46]].middle& bitmasks[46]) >> bits[46]) << 46;
	      case 45: value |= ((Brow[vecs[45]].middle& bitmasks[45]) >> bits[45]) << 45;
	      case 44: value |= ((Brow[vecs[44]].middle& bitmasks[44]) >> bits[44]) << 44;
	      case 43: value |= ((Brow[vecs[43]].middle& bitmasks[43]) >> bits[43]) << 43;
	      case 42: value |= ((Brow[vecs[42]].middle& bitmasks[42]) >> bits[42]) << 42;
	      case 41: value |= ((Brow[vecs[41]].middle& bitmasks[41]) >> bits[41]) << 41;
	      case 40: value |= ((Brow[vecs[40]].middle& bitmasks[40]) >> bits[40]) << 40;
	      case 39: value |= ((Brow[vecs[39]].middle& bitmasks[39]) >> bits[39]) << 39;
	      case 38: value |= ((Brow[vecs[38]].middle& bitmasks[38]) >> bits[38]) << 38;
	      case 37: value |= ((Brow[vecs[37]].middle& bitmasks[37]) >> bits[37]) << 37;
	      case 36: value |= ((Brow[vecs[36]].middle& bitmasks[36]) >> bits[36]) << 36;
	      case 35: value |= ((Brow[vecs[35]].middle& bitmasks[35]) >> bits[35]) << 35;
	      case 34: value |= ((Brow[vecs[34]].middle& bitmasks[34]) >> bits[34]) << 34;
	      case 33: value |= ((Brow[vecs[33]].middle& bitmasks[33]) >> bits[33]) << 33;
	      case 32: value |= ((Brow[vecs[32]].middle& bitmasks[32]) >> bits[32]) << 32;
	      case 31: value |= ((Brow[vecs[31]].middle& bitmasks[31]) >> bits[31]) << 31;
	      case 30: value |= ((Brow[vecs[30]].middle& bitmasks[30]) >> bits[30]) << 30;
	      case 29: value |= ((Brow[vecs[29]].middle& bitmasks[29]) >> bits[29]) << 29;
	      case 28: value |= ((Brow[vecs[28]].middle& bitmasks[28]) >> bits[28]) << 28;
	      case 27: value |= ((Brow[vecs[27]].middle& bitmasks[27]) >> bits[27]) << 27;
	      case 26: value |= ((Brow[vecs[26]].middle& bitmasks[26]) >> bits[26]) << 26;
	      case 25: value |= ((Brow[vecs[25]].middle& bitmasks[25]) >> bits[25]) << 25;
	      case 24: value |= ((Brow[vecs[24]].middle& bitmasks[24]) >> bits[24]) << 24;
	      case 23: value |= ((Brow[vecs[23]].middle& bitmasks[23]) >> bits[23]) << 23;
	      case 22: value |= ((Brow[vecs[22]].middle& bitmasks[22]) >> bits[22]) << 22;
	      case 21: value |= ((Brow[vecs[21]].middle& bitmasks[21]) >> bits[21]) << 21;
	      case 20: value |= ((Brow[vecs[20]].middle& bitmasks[20]) >> bits[20]) << 20;
	      case 19: value |= ((Brow[vecs[19]].middle& bitmasks[19]) >> bits[19]) << 19;
	      case 18: value |= ((Brow[vecs[18]].middle& bitmasks[18]) >> bits[18]) << 18;
	      case 17: value |= ((Brow[vecs[17]].middle& bitmasks[17]) >> bits[17]) << 17;
	      case 16: value |= ((Brow[vecs[16]].middle& bitmasks[16]) >> bits[16]) << 16;
	      case 15: value |= ((Brow[vecs[15]].middle& bitmasks[15]) >> bits[15]) << 15;
	      case 14: value |= ((Brow[vecs[14]].middle& bitmasks[14]) >> bits[14]) << 14;
	      case 13: value |= ((Brow[vecs[13]].middle& bitmasks[13]) >> bits[13]) << 13;
	      case 12: value |= ((Brow[vecs[12]].middle& bitmasks[12]) >> bits[12]) << 12;
	      case 11: value |= ((Brow[vecs[11]].middle& bitmasks[11]) >> bits[11]) << 11;
	      case 10: value |= ((Brow[vecs[10]].middle& bitmasks[10]) >> bits[10]) << 10;
	      case  9: value |= ((Brow[vecs[ 9]].middle& bitmasks[ 9]) >> bits[ 9]) <<  9;
	      case  8: value |= ((Brow[vecs[ 8]].middle& bitmasks[ 8]) >> bits[ 8]) <<  8;
	      case  7: value |= ((Brow[vecs[ 7]].middle& bitmasks[ 7]) >> bits[ 7]) <<  7;
	      case  6: value |= ((Brow[vecs[ 6]].middle& bitmasks[ 6]) >> bits[ 6]) <<  6;
	      case  5: value |= ((Brow[vecs[ 5]].middle& bitmasks[ 5]) >> bits[ 5]) <<  5;
	      case  4: value |= ((Brow[vecs[ 4]].middle& bitmasks[ 4]) >> bits[ 4]) <<  4;
	      case  3: value |= ((Brow[vecs[ 3]].middle& bitmasks[ 3]) >> bits[ 3]) <<  3;
	      case  2: value |= ((Brow[vecs[ 2]].middle& bitmasks[ 2]) >> bits[ 2]) <<  2;
	      case  1: value |= ((Brow[vecs[ 1]].middle& bitmasks[ 1]) >> bits[ 1]) <<  1;
	      case  0: value |= ((Brow[vecs[ 0]].middle& bitmasks[ 0]) >> bits[ 0]) <<  0;
     	 default:
      	  break;
      }
      /* gathering sign bits*/

      Arow[a_vec].middle |= value;

      value = 0;
      switch(todo-1)
		{
    	  case 63: value |= ((Brow[vecs[63]].sign & bitmasks[63]) >> bits[63]) << 63;
	      case 62: value |= ((Brow[vecs[62]].sign & bitmasks[62]) >> bits[62]) << 62;
	      case 61: value |= ((Brow[vecs[61]].sign & bitmasks[61]) >> bits[61]) << 61;
	      case 60: value |= ((Brow[vecs[60]].sign & bitmasks[60]) >> bits[60]) << 60;
	      case 59: value |= ((Brow[vecs[59]].sign & bitmasks[59]) >> bits[59]) << 59;
	      case 58: value |= ((Brow[vecs[58]].sign & bitmasks[58]) >> bits[58]) << 58;
	      case 57: value |= ((Brow[vecs[57]].sign & bitmasks[57]) >> bits[57]) << 57;
	      case 56: value |= ((Brow[vecs[56]].sign & bitmasks[56]) >> bits[56]) << 56;
	      case 55: value |= ((Brow[vecs[55]].sign & bitmasks[55]) >> bits[55]) << 55;
	      case 54: value |= ((Brow[vecs[54]].sign & bitmasks[54]) >> bits[54]) << 54;
	      case 53: value |= ((Brow[vecs[53]].sign & bitmasks[53]) >> bits[53]) << 53;
	      case 52: value |= ((Brow[vecs[52]].sign & bitmasks[52]) >> bits[52]) << 52;
	      case 51: value |= ((Brow[vecs[51]].sign & bitmasks[51]) >> bits[51]) << 51;
	      case 50: value |= ((Brow[vecs[50]].sign & bitmasks[50]) >> bits[50]) << 50;
	      case 49: value |= ((Brow[vecs[49]].sign & bitmasks[49]) >> bits[49]) << 49;
	      case 48: value |= ((Brow[vecs[48]].sign & bitmasks[48]) >> bits[48]) << 48;
	      case 47: value |= ((Brow[vecs[47]].sign & bitmasks[47]) >> bits[47]) << 47;
	      case 46: value |= ((Brow[vecs[46]].sign & bitmasks[46]) >> bits[46]) << 46;
	      case 45: value |= ((Brow[vecs[45]].sign & bitmasks[45]) >> bits[45]) << 45;
	      case 44: value |= ((Brow[vecs[44]].sign & bitmasks[44]) >> bits[44]) << 44;
	      case 43: value |= ((Brow[vecs[43]].sign & bitmasks[43]) >> bits[43]) << 43;
	      case 42: value |= ((Brow[vecs[42]].sign & bitmasks[42]) >> bits[42]) << 42;
	      case 41: value |= ((Brow[vecs[41]].sign & bitmasks[41]) >> bits[41]) << 41;
	      case 40: value |= ((Brow[vecs[40]].sign & bitmasks[40]) >> bits[40]) << 40;
	      case 39: value |= ((Brow[vecs[39]].sign & bitmasks[39]) >> bits[39]) << 39;
	      case 38: value |= ((Brow[vecs[38]].sign & bitmasks[38]) >> bits[38]) << 38;
	      case 37: value |= ((Brow[vecs[37]].sign & bitmasks[37]) >> bits[37]) << 37;
	      case 36: value |= ((Brow[vecs[36]].sign & bitmasks[36]) >> bits[36]) << 36;
	      case 35: value |= ((Brow[vecs[35]].sign & bitmasks[35]) >> bits[35]) << 35;
	      case 34: value |= ((Brow[vecs[34]].sign & bitmasks[34]) >> bits[34]) << 34;
	      case 33: value |= ((Brow[vecs[33]].sign & bitmasks[33]) >> bits[33]) << 33;
	      case 32: value |= ((Brow[vecs[32]].sign & bitmasks[32]) >> bits[32]) << 32;
	      case 31: value |= ((Brow[vecs[31]].sign & bitmasks[31]) >> bits[31]) << 31;
	      case 30: value |= ((Brow[vecs[30]].sign & bitmasks[30]) >> bits[30]) << 30;
	      case 29: value |= ((Brow[vecs[29]].sign & bitmasks[29]) >> bits[29]) << 29;
	      case 28: value |= ((Brow[vecs[28]].sign & bitmasks[28]) >> bits[28]) << 28;
	      case 27: value |= ((Brow[vecs[27]].sign & bitmasks[27]) >> bits[27]) << 27;
	      case 26: value |= ((Brow[vecs[26]].sign & bitmasks[26]) >> bits[26]) << 26;
	      case 25: value |= ((Brow[vecs[25]].sign & bitmasks[25]) >> bits[25]) << 25;
	      case 24: value |= ((Brow[vecs[24]].sign & bitmasks[24]) >> bits[24]) << 24;
	      case 23: value |= ((Brow[vecs[23]].sign & bitmasks[23]) >> bits[23]) << 23;
	      case 22: value |= ((Brow[vecs[22]].sign & bitmasks[22]) >> bits[22]) << 22;
	      case 21: value |= ((Brow[vecs[21]].sign & bitmasks[21]) >> bits[21]) << 21;
	      case 20: value |= ((Brow[vecs[20]].sign & bitmasks[20]) >> bits[20]) << 20;
	      case 19: value |= ((Brow[vecs[19]].sign & bitmasks[19]) >> bits[19]) << 19;
	      case 18: value |= ((Brow[vecs[18]].sign & bitmasks[18]) >> bits[18]) << 18;
	      case 17: value |= ((Brow[vecs[17]].sign & bitmasks[17]) >> bits[17]) << 17;
	      case 16: value |= ((Brow[vecs[16]].sign & bitmasks[16]) >> bits[16]) << 16;
	      case 15: value |= ((Brow[vecs[15]].sign & bitmasks[15]) >> bits[15]) << 15;
	      case 14: value |= ((Brow[vecs[14]].sign & bitmasks[14]) >> bits[14]) << 14;
	      case 13: value |= ((Brow[vecs[13]].sign & bitmasks[13]) >> bits[13]) << 13;
	      case 12: value |= ((Brow[vecs[12]].sign & bitmasks[12]) >> bits[12]) << 12;
	      case 11: value |= ((Brow[vecs[11]].sign & bitmasks[11]) >> bits[11]) << 11;
	      case 10: value |= ((Brow[vecs[10]].sign & bitmasks[10]) >> bits[10]) << 10;
	      case  9: value |= ((Brow[vecs[ 9]].sign & bitmasks[ 9]) >> bits[ 9]) <<  9;
	      case  8: value |= ((Brow[vecs[ 8]].sign & bitmasks[ 8]) >> bits[ 8]) <<  8;
	      case  7: value |= ((Brow[vecs[ 7]].sign & bitmasks[ 7]) >> bits[ 7]) <<  7;
	      case  6: value |= ((Brow[vecs[ 6]].sign & bitmasks[ 6]) >> bits[ 6]) <<  6;
	      case  5: value |= ((Brow[vecs[ 5]].sign & bitmasks[ 5]) >> bits[ 5]) <<  5;
	      case  4: value |= ((Brow[vecs[ 4]].sign & bitmasks[ 4]) >> bits[ 4]) <<  4;
	      case  3: value |= ((Brow[vecs[ 3]].sign & bitmasks[ 3]) >> bits[ 3]) <<  3;
	      case  2: value |= ((Brow[vecs[ 2]].sign & bitmasks[ 2]) >> bits[ 2]) <<  2;
	      case  1: value |= ((Brow[vecs[ 1]].sign & bitmasks[ 1]) >> bits[ 1]) <<  1;
	      case  0: value |= ((Brow[vecs[ 0]].sign & bitmasks[ 0]) >> bits[ 0]) <<  0;
      default:
        break;
      }

		Arow[a_vec].sign |= value;

/*       for(int k = 0; k < todo; ++k) { */
/*         value |= ((Brow[vecs[k]].units & bitmasks[k]) << bits[k]) >> k; */
/*       } */
      /* and write the vec once */


    }
  }


}

void m7d_apply_p_right_even(m7d_t *A, m7p_t const *P, rci_t start_row, rci_t start_col, int notrans)
{
  	if(A->nrows - start_row == 0)
   		return;
  	rci_t const length = MIN(P->length, A->ncols);
  	wi_t const width = A->width;
  	int step_size = MIN(A->nrows - start_row, 4096);

 	/* our temporary where we store the columns we want to swap around */
  	m7d_t *B = m7d_create(step_size, A->ncols);
  	vtri *Arow;
  	vtri *Brow;

	  /* setup mathematical permutation */
  	rci_t *permutation = (rci_t*)m1ri_calloc(A->ncols, sizeof(rci_t));
  	for(rci_t i = 0; i < A->ncols; ++i)
    	permutation[i] = i;

	if (!notrans)
	{
    	for(rci_t i = start_col; i < length; ++i)
    	{
      	rci_t t = permutation[i];
      	permutation[i] = permutation[P->values[i]];
      	permutation[P->values[i]] = t;
    	}
  	}
  	else {

    for(rci_t i = start_col; i < length; ++i)
    {
    	rci_t t = permutation[length - i - 1];
      	permutation[length - i - 1] = permutation[P->values[length - i - 1]];
     	 permutation[P->values[length - i - 1]] = t;
    }
  }

  /* we have a bitmask to encode where to write to */
  vec *write_mask = (vec*)m1ri_calloc(width, sizeof(vec));

  for(rci_t i = 0; i < A->ncols; i += M1RI_RADIX)
  {
	int const todo = MIN(M1RI_RADIX, A->ncols - i);
    for(int k = 0; k < todo; ++k)
    {
      if(permutation[i + k] == i + k)
      {
        write_mask[i / M1RI_RADIX] |= rightbit << k;
      }
    }

  }
  write_mask[width-1] |=  ~(allbits >> (	M1RI_RADIX - (A->ncols)) % M1RI_RADIX);

  	for(rci_t i = start_row; i < A->nrows; i += step_size)
  	{
   		step_size = MIN(step_size, A->nrows - i);
    	for(int k = 0; k < step_size; ++k)
    	{
	    	Arow = A->rows[i+k];
      		Brow = B->rows[k];


      /*copy row & clear those values which will be overwritten */
		for(wi_t j = 0; j < width; ++j)
		{
        	Brow[j].units = Arow[j].units;
        	Brow[j].sign  = Arow[j].sign;

        	Arow[j].units = Arow[j].units & write_mask[j];
        	Arow[j].sign  = Arow[j].sign & write_mask[j];

      	}
    }

    	/* write out the permutation */
    	m7d_write_col_to_rows_blockd(A, B, permutation, write_mask, i, i + step_size, length);
  }

  m1ri_free(permutation);
  m1ri_free(write_mask);
  m7d_free(B);

}


/**
 * \brief copy permutation Q to P
 *
 * \param P Target permutation matrix (may be NULL)
 * \param Q Source permutation matrix (must not be NULL)
 */

m3p_t * m3p_copy(m3p_t *P, const m3p_t *Q)
{
  if(Q == NULL)
  {
    m1ri_die("Matrix to be copied is NULL");

  }

  if(P->values != NULL)
  {
    m1ri_free(P->values);
  }

  P->values = m1ri_malloc(sizeof(rci_t) * Q->length);
  P->length = Q->length;

  for(int i = 0; i < Q->length; i++)
  {
    P->values[i] = Q->values[i];

  }


}
m5p_t * m5p_copy(m5p_t *P, const m5p_t *Q)
{
    if(Q == NULL)
  {
    m1ri_die("Matrix to be copied is NULL");

  }
    if(P->values != NULL)
  {
    m1ri_free(P->values);
  }

  P->values = m1ri_malloc(sizeof(rci_t) * Q->length);
  P->length = Q->length;

  for(int i = 0; i < Q->length; i++)
  {
    P->values[i] = Q->values[i];

  }

}
m7p_t * m7p_copy(m7p_t *P, const m7p_t *Q)
{
    if(Q->values == NULL)
  {
    m1ri_die("Matrix to be copied is NULL");

  }
    if(P->values != NULL)
  {
    m1ri_free(P->values);
  }

  P->values = m1ri_malloc(sizeof(rci_t) * Q->length);
  P->length = Q->length;

  for(int i = 0; i < Q->length; i++)
  {
    P->values[i] = Q->values[i];

  }

}


/**
 * \brief Set the permutation P to the identity permutation. The only
 * allowed value is 1.
 *
 *
 * \param P Permutation
 * \param value 1
 *

 */

void m3p_set_ui(m3p_t *P, unsigned int value)
{
  assert(value == 1);
  for(int i = 0; i < P->length; i++)
  {
     P->values[i] = i;
  }
}
void m5p_set_ui(m5p_t *P, unsigned int value)
{
   assert(value == 1);
  for(int i = 0; i < P->length; i++)
  {
     P->values[i] = i;
  }
}
void m7p_set_ui(m7p_t *P, unsigned int value)
{
  assert(value == 1);
  for(int i = 0; i < P->length; i++)
  {
     P->values[i] = i;
  }


}

/**
 * Apply the permutation P to A from the left.
 *
 * This is equivalent to row swaps walking from 0 to length-1.
 *
 * \param A Matrix.
 * \param P Permutation.
 */

void m3d_apply_p_left(m3d_t *A, m3p_t const *P)
{

 	if(A->ncols == 0)
  	  return;
  	rci_t const length = MIN(P->length, A->nrows);
  	for (rci_t i = 0; i < length; ++i)
  	 {
   		assert(P->values[i] >= i);
    	m3d_rowswap(A, i, P->values[i]);
  	}
}


void m5d_apply_p_left(m5d_t *A, m5p_t const *P)
{
 	if(A->ncols == 0)
  	  return;
  	rci_t const length = MIN(P->length, A->nrows);
  	for (rci_t i = 0; i < length; ++i)
  	{
   		assert(P->values[i] >= i);
   		m5d_rowswap(A, i, P->values[i]);
  	}
}


void m7d_apply_p_left(m7d_t *A, m7p_t const *P)
{

	if(A->ncols == 0)
  	  return;
  	rci_t const length = MIN(P->length, A->nrows);
  	for (rci_t i = 0; i < length; ++i)
  	{
   		assert(P->values[i] >= i);
    	m7d_rowswap(A, i, P->values[i]);
  	}

}

/**
 * Apply the permutation P to A from the left but transpose P before.
 *
 * This is equivalent to row swaps walking from length-1 to 0.
 *
 * \param A Matrix.
 * \param P Permutation.
 */

void m3d_apply_p_left_trans(m3d_t *A, m3p_t const *P)
{
	if(A->ncols == 0)
    	return;
	rci_t const length = MIN(P->length, A->nrows);
  	for (rci_t i = length - 1; i >= 0; --i)
  	{
    	assert(P->values[i] >= i);
    	m3d_rowswap(A, i, P->values[i]);
  	}


}
void m5d_apply_p_left_trans(m5d_t *A, m5p_t const *P)
{

	if(A->ncols == 0)
    	return;
	rci_t const length = MIN(P->length, A->nrows);
  	for (rci_t i = length - 1; i >= 0; --i)
  	{
    	assert(P->values[i] >= i);
    	m5d_rowswap(A, i, P->values[i]);
  	}



}
void m7d_apply_p_left_trans(m7d_t *A, m7p_t const *P)
{

	if(A->ncols == 0)
    	return;
	rci_t const length = MIN(P->length, A->nrows);
  	for (rci_t i = length - 1; i >= 0; --i)
  	{
    	assert(P->values[i] >= i);
    	m7d_rowswap(A, i, P->values[i]);
  	}


}

/**
 * Apply the permutation P to A from the right.
 *
 * This is equivalent to column swaps walking from length-1 to 0.
 *
 * \param A Matrix.
 * \param P Permutation.
 */

void m3d_apply_p_right(m3d_t *A, m3p_t const *P)
{


  for(rci_t i = (P->length - 1); i >= 0; i--)
  {

      m3d_colswap( A, i, P->values[i]);


  }
}
void m5d_apply_p_right(m5d_t *A, m5p_t const *P)
{
    for(rci_t i = (P->length - 1); i >= 0; i--)
  {

       m5d_colswap( A, i, P->values[i]);


  }
}
void m7d_apply_p_right(m7d_t *A, m7p_t const *P)
{

  for(rci_t i = (P->length - 1); i >= 0; i--)
  {

      m7d_colswap( A, i, P->values[i]);


  }

}




/**
 * Apply the permutation P to A from the right starting at start_row.
 */

void m3d_apply_p_right_even_capped(m3d_t *A, m3p_t const *P, rci_t start_row, rci_t start_col)
{

  	if(!A->nrows)
    	return;
   	m3d_apply_p_right_even(A, P, start_row, start_col, 1);


}
void m5d_apply_p_right_even_capped(m5d_t *A, m5p_t const *P, rci_t start_row, rci_t start_col)
{

	if(!A->nrows)
    	return;
   	m5d_apply_p_right_even(A, P, start_row, start_col, 1);

}
void m7d_apply_p_right_even_capped(m7d_t *A, m7p_t const *P, rci_t start_row, rci_t start_col)
{

	if(!A->nrows)
    	return;
   	m7d_apply_p_right_even(A, P, start_row, start_col, 1);


}


void m3d_apply_p_right_trans_even_capped(m3d_t *A, m3p_t const *P, rci_t start_row, rci_t start_col)
{
	if(!A->nrows)
    	return;
   	m3d_apply_p_right_even(A, P, start_row, start_col, 0);
}
void m5d_apply_p_right_trans_even_capped(m5d_t *A, m5p_t const *P, rci_t start_row, rci_t start_col)
{
	if(!A->nrows)
    	return;
   	m5d_apply_p_right_even(A, P, start_row, start_col, 0);
}
void m7d_apply_p_right_trans_even_capped(m7d_t *A, m7p_t const *P, rci_t start_row, rci_t start_col)
{
	if(!A->nrows)
    	return;
   	m7d_apply_p_right_even(A, P, start_row, start_col, 0);
}


/**
 * Apply the permutation matrix P to A from the right but transpose P before.
 *
 * This is equivalent to column swaps walking from 0 to length-1.
 *
 * \param A Matrix.
 * \param P Permutation.
 */

void m3d_apply_p_right_trans(m3d_t *A, m3p_t const *P)
{
  	if(A->nrows == 0)
    return;
  	rci_t const length = MIN(P->length, A->ncols);
  	int const step_size = MAX(  4096 / A->width, 1);
  	for(rci_t j = 0; j < A->nrows; j += step_size)
  	{
   		rci_t stop_row = MIN(j + step_size, A->nrows);

    	for (rci_t i = 0; i < length; ++i)
    	{
      		assert(P->values[i] >= i);
      		m3d_col_swap_in_rows(A, i, P->values[i], j, stop_row);
    	}
  }
}


void m5d_apply_p_right_trans(m5d_t *A, m5p_t const *P)
{
  	if(A->nrows == 0)
    	return;
  	rci_t const length = MIN(P->length, A->ncols);

  	int const step_size = MAX(  4096 / A->width, 1);

  	for(rci_t j = 0; j < A->nrows; j += step_size)
  	{
   		rci_t stop_row = MIN(j + step_size, A->nrows);

    	for (rci_t i = 0; i < length; ++i)
    	{
      		assert(P->values[i] >= i);
      		m5d_col_swap_in_rows(A, i, P->values[i], j, stop_row);
    	}
  }
}


void m7d_apply_p_right_trans(m7d_t *A, m7p_t const *P)
{
  if(A->nrows == 0)
    return;
  	rci_t const length = MIN(P->length, A->ncols);
  	int const step_size = MAX(  4096 / A->width, 1);
  	for(rci_t j = 0; j < A->nrows; j += step_size)
  	{
   		rci_t stop_row = MIN(j + step_size, A->nrows);

    	for (rci_t i = 0; i < length; ++i)
    	{
      		assert(P->values[i] >= i);
      		m7d_col_swap_in_rows(A, i, P->values[i], j, stop_row);
    	}
  }
}



void  m3d_apply_p_right_trans_tri(m3d_t *A, m3p_t const *Q)
{
 	if(A->nrows == 0)
    return;
  	rci_t const length = MIN(Q->length, A->ncols);
  	int const step_size = MAX(  4096 / A->width, 1);
  	for(rci_t j = 0; j < A->nrows; j += step_size)
  	{
   		rci_t stop_row = MIN(j + step_size, A->nrows);

    	for (rci_t i = 0; i < length; ++i)
    	{
      		assert(Q->values[i] >= i);
      		m3d_col_swap_in_rows(A, i, Q->values[i], j, stop_row);
    	}
  }

}



void  m5d_apply_p_right_trans_tri(m5d_t *A, m5p_t const *Q)
{
	if(A->nrows == 0)
	{
    	return;
  	}

  	rci_t const length = MIN(Q->length, A->ncols);
  	int const step_size = MAX(  4096 / A->width, 1);
  	for(rci_t j = 0; j < A->nrows; j += step_size)
  	{
   		rci_t stop_row = MIN(j + step_size, A->nrows);

    	for (rci_t i = 0; i < length; ++i)
    	{
      		assert(Q->values[i] >= i);
      		m5d_col_swap_in_rows(A, i, Q->values[i], j, stop_row);
    	}
  }
}

void  m7d_apply_p_right_trans_tri(m7d_t *A, m7p_t const *Q)
{
	if(A->nrows == 0)
    {
        	return;
    }
  	rci_t const length = MIN(Q->length, A->ncols);
  	int const step_size = MAX(  4096 / A->width, 1);
  	for(rci_t j = 0; j < A->nrows; j += step_size)
  	{
   		rci_t stop_row = MIN(j + step_size, A->nrows);

    	for (rci_t i = 0; i < length; ++i)
    	{
      		assert(Q->values[i] >= i);
      		m7d_col_swap_in_rows(A, i, Q->values[i], j, stop_row);
    	}
  }
}



