#include "poormath.h"

double matMult(int n3, int n4, int n1, double **CAprime, double **b)
{
#ifdef _HAND_OPTIMIZED_MM
  register int i, j, k, N3, N4;
  register double dACC, t1, t2, t3, t4;
  register double *cPtr, *cPtr2, *bPtr, *bPtr2;
  int n3b, n4b, mod; 

  dACC = 0.0;
      
  n3b = 2;
  n4b = 2;

  mod = (n3 % n3b);
  if(mod != 0)      
    N3 = n3 - mod;            
  else      
    N3 = n3;    
  
  mod = (n4 % n4b);
  if(mod != 0)      
    N4 = n4 - mod;            
  else      
    N4 = n4;      
  
  for(i = 0; i < N3; i += n3b)    
    {
      
      cPtr  = CAprime[i];
      cPtr2 = CAprime[i + 1];
      
      for(j = 0; j < N4; j += n4b)	 
	{
	  t1 = t2 = t3 = t4 = 0.0;
	  bPtr = b[j];
	  bPtr2 = b[j + 1];	 
	
	  for(k = 0; k < n1; k++)	  
	    {
	      t1 += cPtr[k]  * bPtr[k];  
	      t2 += cPtr2[k] * bPtr[k]; 
	      t3 += cPtr[k]  * bPtr2[k]; 
	      t4 += cPtr2[k] * bPtr2[k]; 	       
	    }
	  
	  dACC += (t1 * t1) + (t2 * t2) + (t3 * t3) + (t4 * t4);	  	 
	}
      
      for(;j < n4; j++)
	{
	  bPtr = b[j];
	  
	  t1 = t2 = 0.0;
	  for(k = 0; k < n1; k++)	  	    
	    {
	      t1 += cPtr[k]  * bPtr[k];
	      t2 += cPtr2[k] * bPtr[k];	     
	    }	  
	  
	  dACC += (t1 * t1) + (t2 * t2);
	}      
    }

  for(; i < n3; i++)    
    {
      cPtr  = CAprime[i];     
     
      for(j = 0; j < N4; j += 2)	 
	{
	  t1 = t2 = 0.0;	 
	  bPtr = b[j];
	  bPtr2 = b[j + 1];
	 
	  for(k = 0; k < n1; k++)	  
	    {
	      t1 += cPtr[k]  * bPtr[k];  	 
	      t2 += cPtr[k]  * bPtr2[k]; 	      	     	      
	    }
	  
	  dACC += (t1 * t1) + (t2 * t2);	  
	}
	
      for(;j < n4; j++)
	{
	  bPtr = b[j];
	  
	  t1 = 0.0;
	  for(k = 0; k < n1; k++)	  	    
	    {
	      t1 += cPtr[k]  * bPtr[k];	        
	    }	  
	  
	  dACC += (t1 * t1);
	}      
    }
#else
  
  register int i, j, k;
  register double dACC, t;

  dACC = 0.0;

  for(i = 0; i < n3; i++)
    for(j = 0; j < n4; j++)
      {
	t = 0.0;
	for(k = 0; k < n1; k++)
	  t += CAprime[i][k] * b[j][k];

	dACC += (t * t);
      }

#endif

  /*  printf("TIME %f, dacc %f\n", gettime() - time, dACC);*/
  
  return dACC;
}

void makeCAprime(int n3, int n1, int n2, char **a, double **c, double **CAprime, int *corder)
{
  register int i, j, k, m;
  register double temp;
 
  for(i = 0; i < n1; i++)    
    {
	      m = getNonZero(n2, a[i], corder);
	      for(j = 0; j < n3; j++)
		{
		  temp = 0.0;
		  for(k = 0; k < m; k++)	 
		    temp += c[j][corder[k]]; 	  
		  CAprime[j][i] = temp;	
		}
    }  
}

void permuteCAprime(int n3, int n1, int n2, char **a, double **c, double **CAprime, long *idum, int *iordre, int *corder)
{
  register int i, j, k;
  register double temp;
  int m;
  
  for(i = 0; i < n1; i++)	 	 
    {	   
      m = permute(idum, n2, iordre, a[i], corder);
      for(j = 0; j < n3; j++)
	{
	  temp = 0.0;	       
	  for(k = 0; k < m; k++)
	    temp += c[j][corder[k]];	      
	  CAprime[j][i] = temp;	      	      
	}	 	   
    }
  return;
}
