#ifndef COMMON_H_
#define COMMON_H_

#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0 / IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1 + IMM1 / NTAB)
#define EPS 1.2e-7
#define RNMX (1.0 - EPS)
#define NTURN 10

extern float ran2(long *idum);
extern double gettime();
extern int getNonZero(int n, char *a, int *corder);
extern int permute(long *idum, int n, int *iordre, char *a, int *corder);
extern int filexists(char *filename);

#endif /*COMMON_H_*/
