#ifndef BLASMATH_H_
#define BLASMATH_H_


/*#ifdef ROWS
 #include <cblas.h>
 #else*/

#ifdef _USE_ATLAS
#include <cblas.h>
#endif

#ifdef _USE_ACML
#include <acml.h>
#endif

#ifdef _USE_MKL
#include <mkl.h>
#endif

extern double matMult(int n3, int n4, int n1, double *CAprime, double *b);
extern void makeCAprime(int n3, int n1, int n2, char *a, double *c, double *CAprime, int *corder);
extern void permuteCAprime(int n3, int n1, int n2, char *a, double *c, double *CAprime, long *idum, int *iordre, int *corder);

#endif /*BLASMATH_H_*/
