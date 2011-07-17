/*  Axelerated Parafit and Axelerated Parallel Parafit version 1.0              */
/*  Copyright  October 2006 by Alexandros Stamatakis                            */
/*  All rights reserved.                                                        */
/*                                                                              */
/*  This is a highly optimized and parallelized C porting of the original       */
/*  Parafit Fortran code by Pierre Legendre available at                        */
/*                                                                              */
/*  http://www.bio.umontreal.ca/Casgrain/en/labo/parafit.html                   */
/*                                                                              */
/*  and published as:                                                           */
/*                                                                              */ 
/*  Legendre, P., Y. Desdevises and E. Bazin. 2002.                             */ 
/*  A statistical test for host-parasite coevolution.                           */
/*  Systematic Biology 51(2): 217-234.                                          */
/*                                                                              */
/*                                                                              */
/* This code is made available under GNU GPL version 3 or higher                */
/*                                                                              */
/* When using this code please cite:                                            */
/* A. Stamatakis, A. Auch, J. Meier-Kolthoff, M. Göker: “AxPcoords & Parallel   */
/* AxParafit: Statistical Co-Phylogenetic Analyses on Thousands of Taxa”.       */
/* In BMC Bioinformatics, 8:405, 2007.                                          */
/*                                                                              */
/* Dr. Alexandros Stamatakis                                                    */
/*                                                                              */
/* Group Leader: Scientific Computing Group (Exelixis Lab & HPC Infrastructure) */
/* Heidelberg Institute for Theoretical Studies (HITS gGmbH)                    */
/*                                                                              */
/* Schloss-Wolfsbrunnenweg 35                                                   */
/* D-69118 Heidelberg                                                           */
/* Germany                                                                      */
/*                                                                              */
/* Tel:   +49 1511 7496080 (Mobile)                                             */
/*       +49 6221 533240 (Office)                                               */
/* Fax:   +49 6221 533298                                                       */
/* Skype: stamatak                                                              */
/* Email: Alexandros.Stamatakis@h-its.org                                       */
/* WWW:   http://www.exelixis-lab.org                                           */


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
