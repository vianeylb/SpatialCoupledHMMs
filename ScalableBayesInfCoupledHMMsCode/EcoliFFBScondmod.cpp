#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <math.h>
#include <stdlib.h>

extern "C" {
    
SEXP EcoliFFBScondmod(SEXP Yt1, SEXP Yt, SEXP nrow, SEXP I, SEXP I1, SEXP a, SEXP b, SEXP pr){
long double a1, b1, g1, p1, p2;
int l, nrow1, A, C;
SEXP Prob;

Yt1 = coerceVector(Yt1, INTSXP);
Yt = coerceVector(Yt, INTSXP);
nrow = coerceVector(nrow, INTSXP);
I = coerceVector(I, REALSXP);
I1 = coerceVector(I1, REALSXP);
a = coerceVector(a, REALSXP);
b = coerceVector(b, REALSXP);
pr = coerceVector(pr, REALSXP);

nrow1 = INTEGER(nrow)[0];
a1 = REAL(a)[0];
b1 = REAL(b)[0];
g1 = REAL(pr)[0];

PROTECT(Prob = allocVector(REALSXP, 2));

p1 = 1;
p2 = 1;
for (l = 0; l < nrow1; l++) {
              A = INTEGER(Yt1)[l];
              C = INTEGER(Yt)[l];
              if ((C == 0) && (A == 0))
    	      {p1 = p1 * exp(-a1 - b1 * REAL(I)[0]);
    	       p2 = p2 * exp(-a1 - b1 * REAL(I1)[0]);
              } else if ((C == 0) && (A == 1))
              { p1 = p1*(1 - exp(-a1 - b1* REAL(I)[0]));
                p2 = p2*(1 - exp(-a1 - b1* REAL(I1)[0]));
              } else if ( (C == 1) && (A == 0) )
              { p1 = p1*g1;
                p2 = p2*g1;
              } else { p1 = p1*(1-g1); 
                       p2 = p2*(1-g1);} 
}       

REAL(Prob)[0] = p1;
REAL(Prob)[1] = p2;

UNPROTECT(1);
return Prob;
}        
}      
