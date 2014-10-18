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
 
 m1ri_io.c
 
 */
#include "m1ri_io.h"
#if __M1RI_HAVE_LIBPNG
#include <png.h>
#endif /* __M1RI_HAVE_LIBPNG */
/** 
 
 Print a block of an m3d
 a = the unit bits
 b =  the sign bits
 l_unused = space to the left of block unused in the matrix
 r_unused =  space to the right of block unused in the matrix
 */
 
 



/* /print_m3d_block_buffered original title */

static inline void print_m3d_block(vec a, vec b, u_int32_t l_unused, u_int32_t r_unused)
{
	int i = 0;
    bool out;
    char buffer[512];
    for(int x = (0  + l_unused); x < (64 - r_unused); x = x + 1)
    {
    
    	buffer[i++] = '[';
    	buffer[i++] = ' ';
        out = (( a & (rightbit <<  x)) == (b & (rightbit  << x))) ? 0:  1;
        
        if((out == 0) && (b & (rightbit  << x)))
        {
        	
        	buffer[i++] = '2'; 
              /*  printf("[ 1 ]"); */
        }
        
        else if(out == 1) /* && (b & (rightbit  << x))) */
        {
        	buffer[i++] = '1'; 
            /* printf("[ 2 ]"); */
        }
        
		else
        {   
           buffer[i++] = '0';
        }
        
        buffer[i++] = ' ';
        buffer[i++] = ']';
    }
    buffer[i] = '\0';
    printf("%s", buffer);
    
}




/**  
	The function is  to print m3d_t
*/
 void m3d_print(const m3d_t *  a)
{
    int i, m;
    printf("\n \n"); 
    
    /* testing for io improvement  */
    if(a->width >=  16)
    {
    
    char buffer[8192];
    for( i  = 0; i < a->nrows ; i++)
    {
        
      
         print_m3d_block(a->rows[i][0].units, a->rows[i][0].sign, a->fcol, 0 );
             
             m = 1;
             while(m < (a->width -1))
             {
                 print_m3d_block(a->rows[i][m].units, a->rows[i][m].sign, 0, 0);
                 ++m;
                 
             }
             
             if(a->ncols%64 == 0)
             {
                 print_m3d_block(a->rows[i][m].units, a->rows[i][m].sign, 0, 0 );
                 
             }
             
             if(a->ncols%64 != 0)
             {
             
              print_m3d_block(a->rows[i][m].units, a->rows[i][m].sign, 0, (64 - a->ncols%64) );
              printf("got here \n");
             }
         
    
      
        printf("\n");
          }
      }
    
    else
    {   
    for( i  = 0; i < a->nrows ; i++)
    {
        
         if(a->width > 1)
         {
         print_m3d_block(a->rows[i][0].units, a->rows[i][0].sign, a->fcol, 0);
             
             m = 1;
             while(m < (a->width -1))
             {
                 print_m3d_block(a->rows[i][m].units, a->rows[i][m].sign, 0, 0);
                 ++m;
             }
             
             if(a->ncols%64 == 0)
             {
                 print_m3d_block(a->rows[i][m].units, a->rows[i][m].sign, 0, 0 );
                 
             }
             
             if(a->ncols%64 != 0)
             {
              print_m3d_block(a->rows[i][m].units, a->rows[i][m].sign, 0, (64 - a->ncols%64) );
             }
         }
    
        if(a->width  ==  1)
        {
            
            if(a->ncols%64 != 0)
            {
             
         print_m3d_block(a->rows[i][0].units, a->rows[i][0].sign, a->fcol, (64 - (a->ncols + a->fcol)%64) );
             
             
         }
          if(a->ncols%64 == 0)
          {
             print_m3d_block(a->rows[i][0].units, a->rows[i][0].sign, a->fcol, 0 );
           }
      
          }
        printf("\n");
          }
     } 


    printf("\n \n \n");
    
}
void print_m5d_block(vec a, vec b, vec c,  u_int32_t l_unused, u_int32_t r_unused)
{
    bool out[3];
    short value, x ;
    for( x = (0  + l_unused); x < (64 - r_unused); x = x + 1)
    {
        value = 0;
		out[0] = ( a & (rightbit <<  x));
        out[1] =  ( b & (rightbit <<  x));
        out[2]  = ( c & (rightbit <<  x));
        
        if (out[2] > 0) {
            value =  value + 1;
        
        }
        
        if (out[1] > 0) {
            value = value + 1;
        }
        
        if(out[0] > 0)
        {
            value = value + 2;
        }
        printf("[ %d ]", value);
    }
}

void m5d_print(const m5d_t *a)
{
    int i, m;
    printf("\n \n");
  	for( i  = 0; i < a->nrows ; i++)
        {
            
            if(a->width > 1)
            {
                print_m5d_block(a->rows[i][0].units, a->rows[i][0].middle,  a->rows[i][0].sign, a->fcol, 0);
                
                m = 1;
                while(m < (a->width -1))
                {
                    print_m5d_block(a->rows[i][m].units,a->rows[i][m].middle, a->rows[i][m].sign, 0, 0);
                    ++m;
                }
                
                if(a->ncols%64 == 0)
                {
                    print_m5d_block(a->rows[i][m].units,a->rows[i][m].middle, a->rows[i][m].sign, 0, 0 );
                    
                }
                
                if(a->ncols%64 != 0)
                {
                    print_m5d_block(a->rows[i][m].units,a->rows[i][m].middle, a->rows[i][m].sign, 0, (64 - a->ncols%64) );
                }
            }
            
            if(a->width  ==  1)
            {
                if(a->ncols%64 != 0)
                {
                    print_m5d_block(a->rows[i][0].units,a->rows[i][0].middle,  a->rows[i][0].sign, a->fcol, (64 - (a->ncols + a->fcol)%64) );
                }
                if(a->ncols%64 == 0)
                {
                    print_m5d_block(a->rows[i][0].units,a->rows[i][0].middle, a->rows[i][0].sign, a->fcol, 0 );
                }   
            }
            printf("\n");
    }
    
    printf("\n \n \n");
}

static inline void print_m7d_block(vec a, vec b, vec c,  u_int32_t l_unused, u_int32_t r_unused)
{
	
    u_int64_t out[3];
    short value;/* , x ; */
    for( int x = (0  + l_unused); x < (64 - r_unused); x = x + 1)
    {
        out[0] =  a & (rightbit <<  x);
        out[1] =   b & (rightbit <<  x);
        out[2]  =  c & (rightbit <<  x);
        
    	if((out[2]  ==  out[1]) & (out[1] == out[0]))
       	{
       	printf("[ 0 ]");	
       	}
       
        else
        {
        value = 0;
        if (out[0]  > 0) {
            value = value + 1;
        }
        
        if (out[1] > 0) {
            value = value + 2;
        }
        
        if(out[2] > 0 )
        {
            value = value + 4;           
        }
        printf("[ %d ]", value);
		}
     
    
    }
    
}
/** 
	Prints an m7d_t matrix
*/

void m7d_print(const m7d_t *a)
{
    int i, m;
    printf("\n \n");
  	for( i  = 0; i < a->nrows ; i++)
        {
            
            if(a->width > 1)
            {
                print_m7d_block(a->rows[i][0].units, a->rows[i][0].middle,  a->rows[i][0].sign, a->fcol, 0);
                
                m = 1;
                while(m < (a->width -1))
                {
                    print_m7d_block(a->rows[i][m].units,a->rows[i][m].middle, a->rows[i][m].sign, 0, 0);
                    ++m;
                }
                
                if(a->ncols%64 == 0)
                {
                	
                    print_m7d_block(a->rows[i][m].units,a->rows[i][m].middle, a->rows[i][m].sign, 0, 0 );
                    
                }
                
                if(a->ncols%64 != 0)
                {
                	
                    print_m7d_block(a->rows[i][m].units,a->rows[i][m].middle, a->rows[i][m].sign, 0, (64 - a->ncols%64) );
                }
            }
            
            if(a->width  ==  1)
            {
                if(a->ncols%64 != 0)
                {
                    print_m7d_block(a->rows[i][0].units,a->rows[i][0].middle,  a->rows[i][0].sign, a->fcol, (64 - (a->ncols + a->fcol)%64) );
                }
                if(a->ncols%64 == 0)
                {
                    print_m7d_block(a->rows[i][0].units,a->rows[i][0].middle, a->rows[i][0].sign, a->fcol, 0 );
                }   
            }
            printf("\n");
    }
    
    printf("\n \n \n");
}

void m3d_specs(const m3d_t * a)
{
    
    
    if (a->flags & iswindowed) {
        printf("Is Windowed   \n");
    }
    else if (a->flags & iswindowed == 0)
    {
        printf("Is not windowed   \n");
    }
    
    printf("Number of rows   : %d \n", a->nrows );
    printf("Number of columns: %d \n", a->ncols );

    printf("Width------------: %d \n", a->width );
    
    
}

void m3d_fullinfo(m3d_t * a)
{
    m3d_print(a);
    m3d_specs(a);
  
}


void m5d_specs(m5d_t * a)
{
    
    
    if (a->flags & iswindowed) {
        printf("Is Windowed   \n");
    }
    else if ((a->flags & iswindowed) == 0)
    {
        printf("Is not windowed   \n");
    }
   
    printf("Number of rows   : %d \n", a->nrows );
    printf("Number of columns: %d \n", a->ncols );
    printf("Width------------: %d \n", a->width );
    
    
}
/**
	Prints an m5d_t to standard output and then displays: 
	1.  Number of columns
	2.  Number of rows
	3.  Width in vfd's .
*/
void m5d_fullinfo(m5d_t * a)
{
    m5d_print(a);
    m5d_specs(a);
  
}

void m7d_specs(m7d_t * a)
{
    if ((a->flags == iswindowed)) {
        printf("Is Windowed   \n");
    }
    else if (a->flags  == notwindowed)
    {
        printf("Is not windowed   \n");
    }
    printf("Number of rows   : %d \n", a->nrows );
    printf("Number of columns: %d \n", a->ncols );
   	printf("Width------------: %d \n", a->width );
    
}

/**
	Prints an m7d_t to standard output and then displays: 
	1.  Number of columns
	2.  Number of rows
	3.  Width in vtri's .
*/
void m7d_fullinfo(m7d_t * a)
{
    m7d_print(a);
    m7d_specs(a);
  
}







void m3p_print(m3p_t const * P )
{
 	printf("[ ");
  	for(rci_t i = 0; i < P->length; ++i) 
  	{
    	printf("%zd ", (size_t)P->values[i]);
  	}
  	printf("]");
}

void m5p_print(m5p_t const * P)
{
 	printf("[ ");
  	for(rci_t i = 0; i < P->length; ++i) 
  	{
    	printf("%zd ", (size_t)P->values[i]);
  	}
  	printf("]");
}



void m7p_print(m7p_t const * P)
{
 	printf("[ ");
  	for(rci_t i = 0; i < P->length; ++i) 
  	{
    	printf("%zd ", (size_t)P->values[i]);
  	}
  	printf("]");
}

   

 

 


 m3d_t m3d_read_textfile(const char * fn)
 {
             FILE *in_file  = fopen("name_of_file", "r"); /*  read only  */
             
 
 }
 m5d_t m5d_read_textfile(const char * fn)
 {
             FILE *in_file  = fopen("name_of_file", "r"); /*  read only  */
             
 }
 
 m7d_t m7d_read_textfile (const char * fn)
{
            FILE *in_file  = fopen("name_of_file", "r"); /*  read only  */
            

}


#if __M1RI_HAVE_LIBPNG
#define PNGSIGSIZE 8
int m3d_to_png(const m3d_t *A, const char *fn, int compression_level, const char *comment, int verbose) {
  FILE *fh = fopen(fn, "wb");

  if (!fh) {
    if(verbose)
      printf("Could not open file '%s' for writing\n",fn);
    return 1;
  }

  png_structp png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);

  if (!png_ptr) {
    if(verbose)
      printf("failed to initialise PNG write struct.\n");
    fclose(fh);
    return 3;
  }
  png_set_user_limits(png_ptr, 0x7fffffffL,  0x7fffffffL);

  png_infop info_ptr = png_create_info_struct(png_ptr);

  if (!info_ptr) {
    if (verbose)
      printf("failed to initialise PNG info struct\n");
    png_destroy_write_struct(&png_ptr, &info_ptr);
    fclose(fh);
    return 3;
  }

  if (setjmp(png_jmpbuf(png_ptr))) {
    if (verbose)
      printf("error writing PNG file\n");
    png_destroy_write_struct(&png_ptr, &info_ptr);
    fclose(fh);
    return 1;
  }

  png_init_io(png_ptr, fh);
  png_set_compression_level(png_ptr, compression_level);

  png_set_IHDR(png_ptr, info_ptr, A->ncols, A->nrows, 2, \
               PNG_COLOR_TYPE_PALETTE, \
               PNG_INTERLACE_NONE, \
               PNG_COMPRESSION_TYPE_DEFAULT, \
               PNG_FILTER_TYPE_DEFAULT);

  png_text txt_ptr[3];

  char pdate[21];
  time_t ptime=time(NULL);
  struct tm *ltime=localtime(&ptime);
  sprintf(pdate,"%04d/%02d/%02d %02d:%02d:%02d",ltime->tm_year+1900,ltime->tm_mon+1,ltime->tm_mday,ltime->tm_hour,ltime->tm_min,ltime->tm_sec);

  txt_ptr[0].key="Software";
  txt_ptr[0].text="m1ri";
  txt_ptr[0].compression=PNG_TEXT_COMPRESSION_NONE;
  txt_ptr[1].key="Date";
  txt_ptr[1].text=pdate;
  txt_ptr[1].compression=PNG_TEXT_COMPRESSION_NONE;
  txt_ptr[2].key="Comment";
  txt_ptr[2].text=(char*)comment;
  txt_ptr[2].compression=PNG_TEXT_COMPRESSION_NONE;

  png_set_text(png_ptr, info_ptr, txt_ptr, 3);

  png_write_info(png_ptr, info_ptr);

  png_set_packswap(png_ptr);
  /* png_set_invert_mono(png_ptr); */
  
		  png_bytep row = m1ri_calloc(sizeof(u_int8_t),A->ncols/4+16);

		  wi_t j=0;
		  vbg  tmp;
		  u_int64_t p_row[2];/*   Packed Row  */
		  tmp.units = 0;
		  tmp.sign = 0;  /*  = 0; */
  		  u_int64_t ep_bit[64];
  
  
		  for(int q  = 0; q < 64 ; q++)
		  {
			ep_bit[q] = rightbit  << q;
  
		  } 
  
		  
		   for(rci_t i=0; i<A->nrows; i++) 
			{
			  vbg *rowptr = A->rows[i];
			  for(j=0; j<A->width-1; j++) {
	
	
			  tmp.units = rowptr[j].units;
			  tmp.sign = rowptr[j].sign;
	  

			  p_row[0] = p_row[0] | (ep_bit[63] & tmp.units);
			  p_row[0] = p_row[0] | ((ep_bit[63] & tmp.sign) >> 1);
			  
			  
			  
	  
	  
	  
	  
	  
	  
	  
	  
	  
			  row[8*j+0] = (png_byte)((tmp.units>> 0) & 0xff);
			  row[8*j+1] = (png_byte)((tmp.units>> 8) & 0xff);
			  row[8*j+2] = (png_byte)((tmp.units>>16) & 0xff);
			  row[8*j+3] = (png_byte)((tmp.units>>24) & 0xff);
			  row[8*j+4] = (png_byte)((tmp.units>>32) & 0xff);
			  row[8*j+5] = (png_byte)((tmp.units>>40) & 0xff);
			  row[8*j+6] = (png_byte)((tmp.units>>48) & 0xff);
			  row[8*j+7] = (png_byte)((tmp.units>>56) & 0xff);
	  
	  
	  
			  /* ////////////////////////////////////////// */
			}
  
	
			tmp.units = rowptr[j].units;
			switch( (A->ncols/8 + ((A->ncols%8) ? 1 : 0)) %8 ) {
    case 0: row[8*j+7] = (png_byte)((tmp.units>>56) & 0xff);
    case 7: row[8*j+6] = (png_byte)((tmp.units>>48) & 0xff);
    case 6: row[8*j+5] = (png_byte)((tmp.units>>40) & 0xff); 
    case 5: row[8*j+4] = (png_byte)((tmp.units>>32) & 0xff); 
    case 4: row[8*j+3] = (png_byte)((tmp.units>>24) & 0xff);
    case 3: row[8*j+2] = (png_byte)((tmp.units>>16) & 0xff);
    case 2: row[8*j+1] = (png_byte)((tmp.units>> 8) & 0xff);
    case 1: row[8*j+0] = (png_byte)((tmp.units>> 0) & 0xff);
    
    
    
    
       
       /*
           row[8*j+0] = (png_byte)((tmp.units>> 0) & 0xff);
      row[8*j+1] = (png_byte)((tmp.units>> 8) & 0xff);
      row[8*j+2] = (png_byte)((tmp.units>>16) & 0xff);
      row[8*j+3] = (png_byte)((tmp.units>>24) & 0xff);
      row[8*j+4] = (png_byte)((tmp.units>>32) & 0xff);
      row[8*j+5] = (png_byte)((tmp.units>>40) & 0xff);
      row[8*j+6] = (png_byte)((tmp.units>>48) & 0xff);
      row[8*j+7] = (png_byte)((tmp.units>>56) & 0xff);
      
      */
          /* ////////////////////////////////////////// */

    };
    png_write_row(png_ptr, row);
  }
  m1ri_free(row);

  png_write_end(png_ptr, info_ptr);
  png_destroy_write_struct(&png_ptr, &info_ptr);
  fclose(fh);
  return 0;
}




#endif 