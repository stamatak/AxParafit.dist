#include "blasmath.h"
#include <stdlib.h> 



double matMult(int n3, int n4, int n1, double *CAprime, double *b)
{  
  int i;
  double dACC, *C;  

  C = (double *)calloc(n3 * n4, sizeof(double));

  /*#ifdef ROWS
  printf("ATLAS CBLAS \n");
  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, n3, n4, n1, 1.0, CAprime, n1, b, n1, 1.0, C, n4);
  #else*/
  
  {
    char transa='n';
    char transb='t';
    double alpha=1.0;
    double beta=1.0;
  
#ifdef _USE_ATLAS   
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans, n3, n4, n1, 1.0, CAprime, n3, b, n4, 1.0, C, n3);
#endif
    
#ifdef _USE_ACML
    dgemm(transa, transb, n3, n4, n1, alpha, CAprime, n3, b, n4, beta, C, n3); // ACML
    // void dgemm(char transa, char transb, int m, int n, int k, double alpha, double *a, int lda, double *b, int ldb, double beta, double *c, int ldc);
#endif
          
#ifdef _USE_MKL  	       
    dgemm(&transa, &transb, &n3, &n4, &n1, &alpha, CAprime, &n3, b, &n4, &beta, C, &n3); // MKL
    // void dgemm(char *transa,char *transb,int *m,int *n,int *k,double *alpha,double *a,int *lda,double *b,int *ldb,double *beta,double *c,int *ldc);
#endif
  }
         
  dACC = 0.0;
  
  for(i = 0; i < n3 * n4; i++)     
    dACC += C[i] * C[i];
  
  free(C);
     
  return dACC;
}

void makeCAprime(int n3, int n1, int n2, char *a, double *c, double *CAprime, int *corder)
{
  register int i, j, k, m;
  register double temp;
 

  /*#ifdef ROWS 
   for(i = 0; i < n1; i++)    
     {
       m = getNonZero(n2, &a[i * n2], corder);
       for(j = 0; j < n3; j++)
	 {
	   temp = 0.0;
	   for(k = 0; k < m; k++)	 
	     temp += c[j * n2 + corder[k]]; 	  
	   CAprime[j * n1 + i] = temp;	
	 }
     }  
     #else  */
  for(i = 0; i < n1; i++)   
    {
      m = getNonZero(n2, &a[i * n2], corder);
      for(j = 0; j < n3; j++)
	{
	  temp = 0.0;
	  for(k = 0; k < m; k++)      
	    temp += c[j * n2 + corder[k]]; 	      
	  CAprime[i * n3 + j] = temp;	
	}
    }
  
  /*#endif  */

  return;
}

void permuteCAprime(int n3, int n1, int n2, char *a, double *c, double *CAprime, long *idum, int *iordre, int *corder)
{
  register int i, j, k;
  register double temp;
  int m;
  /*#ifdef ROWS  
  for(i = 0; i < n1; i++)	 	 
    {	   
      m = permute(idum, n2, iordre, &a[i * n2], corder);
      for(j = 0; j < n3; j++)
	{
	  temp = 0.0;	       
	  for(k = 0; k < m; k++)
	    temp += c[j * n2 + corder[k]];	      
	  CAprime[j * n1 + i] = temp;	      	      
	}	 	   
    }
    #else*/
      
   for(i = 0; i < n1; i++)	 	 
    {	         
      m = permute(idum, n2, iordre, &a[i * n2], corder);
      for(j = 0; j < n3; j++)
	{
	  temp = 0.0;	       
	  for(k = 0; k < m; k++)
	    temp += c[j * n2 + corder[k]];	      
	  CAprime[i * n3 + j] = temp;	      	      
	}	 	   
    }
   /*#endif*/
  return;
}
