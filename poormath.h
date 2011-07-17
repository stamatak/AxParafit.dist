#ifndef POORMATH_H_
#define POORMATH_H_

extern double matMult(int n3, int n4, int n1, double **CAprime, double **b);
extern void makeCAprime(int n3, int n1, int n2, char **a, double **c, double **CAprime, int *corder);
extern void permuteCAprime(int n3, int n1, int n2, char **a, double **c, double **CAprime, long *idum, int *iordre, int *corder);

#endif /*POORMATH_H_*/
