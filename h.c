//
//  main.c
//  simpleaddingproofofconcept
// This is intended to represent  the addition logic  used 
//    in  the paper TOMAS J. BOOTHBY AND ROBERT W. BRADSHAW,  BITSLICING AND THE METHOD OF FOUR RUSSIANS OVER LARGER FINITE FIELDS
//   
//



#include <stdio.h>

#include <stdbool.h>   //For a declaration of bools
                        //This representation isn't bitpacked properly 

int main(int argc, const char * argv[])
{
    
    
    
    
    
    const int n = 4;
    
    const int nn = n * n;
    
    bool x[nn];  //boolean array, not a  bit vector.   Values actually take up an entire integer
    
    
    bool y[nn];   //second boolean array 
    
    bool r[nn];   //will represent the sum of x and y 
    
    
    int  q;   //
    
    
    
    
    for (q = 0; q < nn; q = q + 2) {  //This shows the logic of addition over GF(3)   
        
                                            
        
        bool s = (x[q]^y[q + 1]^x[q + 1]);      //  For  represnting the S=  x[0]  XOR  y[1] XOR [x1] 
        bool t = (x[q+1]^y[q]^y[q+1]);          //          
        r[q] = ((x[q]^y[q+1])&&(x[q+1]^y[q]));
        r[q + 1] = (s || t);
        
        
        
        
        
        
        
        
        return 0;
        
        
        
        
        
        
    };
    
    
    
    
    
    return 0;
}

