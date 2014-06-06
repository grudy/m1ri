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
/** 
 
 Print a block of an m3d
 a = the unit bits
 b =  the sign bits
 l_unused = space to the left of block unused in the matrix
 r_unused =  space to the right of block unused in the matrix
 */
 
 

 /*
static inline void print_m3d_block(vec a, vec b, u_int32_t l_unused, u_int32_t r_unused)
{
	int x;
    bool out;
    for( x = (0  + l_unused); x < (64 - r_unused); x = x + 1)
    {
        out = (( a & (leftbit >>  x)) == (b & (leftbit  >> x))) ? 0:  1;
        
        if((out == 0) && (b & (leftbit  >> x)))
        {
               printf("[ 1 ]");
        }
        
        else if((out == 1) && (b & (leftbit  >> x)))
        {
            printf("[ 2 ]");
        }
        
		else
        {   
           printf("[ %d ]", out);
        }
    }
  
    
}

*/

///print_m3d_block_buffered original title

static inline void print_m3d_block(vec a, vec b, u_int32_t l_unused, u_int32_t r_unused)
{
	int i = 0;
    bool out;
    char buffer[512];
    for(int x = (0  + l_unused); x < (64 - r_unused); x = x + 1)
    {
    
    	buffer[i++] = '[';
    	buffer[i++] = ' ';
        out = (( a & (leftbit >>  x)) == (b & (leftbit  >> x))) ? 0:  1;
        
        if((out == 0) && (b & (leftbit  >> x)))
        {
        	
        	buffer[i++] = '2'; 
              // printf("[ 1 ]");
        }
        
        else if(out == 1) //&& (b & (leftbit  >> x)))
        {
        	buffer[i++] = '1'; 
            //printf("[ 2 ]");
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

static inline void print_m3d_block_more_buffered(vec a, vec b, u_int32_t l_unused, u_int32_t r_unused, char buffer[], int * i)
{
	 
    bool out;
    
    for(int x = (0  + l_unused); x < (64 - r_unused); x = x + 1)
    {
    
    	buffer[*i++] = '[';
    	buffer[*i++] = ' ';
        out = (( a & (leftbit >>  x)) == (b & (leftbit  >> x))) ? 0:  1;
        
        if((out == 0) && (b & (leftbit  >> x)))
        {
        	
        	buffer[*i++] = '2'; 
              // printf("[ 1 ]");
        }
        
        else if(out == 1) //&& (a & (leftbit  >> x)))
        {
        	buffer[*i++] = '1'; 
            //printf("[ 2 ]");
        }
        
		else
        {   
           buffer[*i++] = '0';
        }
        
        buffer[*i++] = ' ';
        buffer[*i++] = ']';
    }
   
    
}





/**  
	The function is  to print m3d_t
*/
 void m3d_print(m3d_t *  a)
{
    int i, m;
    printf("\n \n"); 
    
    //testing for io improvement 
    if(a->width >=  16)
    {
    
    int  * j = m1ri_malloc(sizeof(int));
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


    printf("\n \n \n ");
    
}
void print_m5d_block(vec a, vec b, vec c,  u_int32_t l_unused, u_int32_t r_unused)
{
    bool out[3];
    short value, x ;
    for( x = (0  + l_unused); x < (64 - r_unused); x = x + 1)
    {
        value = 0;
		out[0] = ( a & (leftbit >>  x));
        out[1] =  ( b & (leftbit >>  x));
        out[2]  = ( c & (leftbit >>  x));
        
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
        printf("[%d]", value);
    }
}

void m5d_print(m5d_t *a)
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
    
    printf("\n \n \n ");
}

static inline void print_m7d_block(vec a, vec b, vec c,  u_int32_t l_unused, u_int32_t r_unused)
{
	
    u_int64_t out[3];
    short value;//, x ;
    for( int x = (0  + l_unused); x < (64 - r_unused); x = x + 1)
    {
        out[0] =  a & (leftbit >>  x);
        out[1] =   b & (leftbit >>  x);
        out[2]  =  c & (leftbit >>  x);
        
    	if((out[2]  ==  out[1]) & (out[1] == out[0]))
       	{
       	printf("[0]");	
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
        printf("[%d]", value);
		}
     
    
    }
    
}
/** 
	Prints an m7d_t matrix
*/
void m7d_print(m7d_t *a)
{
    int i, m;
    printf("\n \n");
    for( i  = 0; i < a->nrows ; i++)
    {
      if(a->width > 1)
        {
            print_m7d_block(a->rows[i][0].units, a->rows[i][0].middle, a->rows[i][0].sign, a->fcol, 0);
            
            m = 1;
            while(m < (a->width -1))
            {
                print_m7d_block(a->rows[i][m].units,a->rows[i][m].middle, a->rows[i][m].sign, 0, 0);
                ++m;
            }
            
            if(a->ncols%64 == 0)
            {
                print_m7d_block(a->rows[i][m].units,a->rows[i][m].middle, a->rows[i][m].sign, a->fcol, 0 );
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
                
                print_m7d_block(a->rows[i][0].units,a->rows[i][0].middle, a->rows[i][0].sign, 0, (64 - a->ncols%64) );
            }
            if(a->ncols%64 == 0)
            {
                print_m7d_block(a->rows[i][0].units, a->rows[i][0].middle,  a->rows[i][0].sign, 0 , 0 );
            }
            
        }
        
        printf("\n");
        
    }
    
    printf("\n \n \n ");
    
    
    
}
/**
	Prints an m3d_t to standard output and then displays: 
	1.  Number of columns
	2.  Number of rows
	3.  Width in vbg's .
*/
void m3d_specs(m3d_t * a)
{
    
    
    if (a->flags & iswindowed) {
        printf("Is Windowed   \n");
    }
    else if (a->flags & iswindowed == 0)
    {
        printf("Is not windowed   \n");
    }
    printf("Number of columns: %d \n", a->ncols );
    printf("Number of rows   : %d \n", a->nrows );
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
    printf("Number of columns: %d \n", a->ncols );
    printf("Number of rows   : %d \n", a->nrows );
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
    
    printf("Number of columns: %d \n", a->ncols );
    printf("Number of rows   : %d \n", a->nrows );
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





/**
 * Print  permutation matrices
 *
 *
 */

void m3p_print(m3p_t const *P)
{

  

}
void m5p_print(m3p_t const *P)
{

}

void m7p_print(m3p_t const *P)
{

}


