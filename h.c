//
//  main.c
//  simpleaddingproofofconcept
//
//   2013 William Alumbaugh.
//



#include <stdio.h>

#include <stdbool.h>


int main(int argc, const char * argv[])
{
    
    
    
    
    
    const int n = 4;
    
    const int nn = n * n;
    
    bool x[nn];
    
    
    
    bool y[nn];
    
    bool r[nn];
    
    
    int  q;
    
    
    
    
    for (q = 0; q < nn; q = q + 2) {
        
        
        
        bool s = (x[q]^y[q + 1]^x[q + 1]);
        bool t = (x[q+1]^y[q]^y[q+1]);
        r[q] = ((x[q]^y[q+1])&&(x[q+1]^y[q]));
        r[q + 1] = (s || t);
        
        
        
        
        
        
        
        
        return 0;
        
        
        
        
        
        
    };
    
    
    
    
    
    return 0;
}

