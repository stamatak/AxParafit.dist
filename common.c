#include "common.h"

#include <sys/time.h>
#include <stdio.h>

#ifdef _WIN32
#include <time.h>
#endif


float ran2(long *idum)
{
  int j;
  long k;
  static long idum2 = 123456789;
  static long iy = 0;
  static long iv[NTAB];
  float temp;
 
  if(*idum <= 0)
    {
      if(-(*idum) < 1) 
	*idum = 1;
      else
	*idum = -(*idum);
      idum2 = (*idum);
      for(j = NTAB + 7; j >= 0; j--)
	{
	  k = (*idum) / IQ1;
	  *idum = IA1 * (*idum - k * IQ1) - k * IR1;
	  if(*idum < 0)
	    *idum += IM1;
	  if(j < NTAB)
	    iv[j] = (*idum);
	}
      iy = iv[0];
    }
  
  k = (*idum) / IQ1;
  *idum = IA1 * (*idum -k * IQ1) - k * IR1;
  if(*idum < 0)
    *idum += IM1;

  k = idum2 / IQ2;
  idum2 = IA2 * (idum2 - k * IQ2) - k * IR2;
  if(idum2 < 0)
    idum2 += IM2;

  j = iy / NDIV;
  iy = iv[j] - idum2;
  iv[j] = *idum;
  if(iy < 1)
    iy += IMM1;

  if((temp = AM * iy) > RNMX) 
    return RNMX;   
  else    
    return temp;
}

double gettime()
{
  #ifndef _WIN32
  struct timeval ttime;
  gettimeofday(&ttime , NULL);
  return ttime.tv_sec + ttime.tv_usec * 0.000001;
  #else
  clock_t time=clock();
  return ((double)time) / ((double)CLOCKS_PER_SEC);
  #endif
}

int getNonZero(int n, char *a, int *corder)
{
  int m = 0, i;
  
  m = 0;
  for(i = 0; i < n; i++)
    {     
      if(a[i] != 0)
	{
	  corder[m] = i;	  
	  m++;
	}
    }
 
  return m;
}

int permute(long *idum, int n, int *iordre, char *a, int *corder)
{
  int m, km1, i, j, itemp;  

  m=n;
  km1 = n-1;
  for(i=1; i <= km1; i++)
    {
    RAND:
      j = 1 + ran2(idum) * m;
      if(j > m) 
	goto RAND;
      itemp = iordre[m - 1];
      iordre[m - 1] = iordre[j - 1];
      iordre[j - 1] = itemp;     
      m = m - 1;
    }  
  
  m = 0;
  for(i = 0; i < n; i++)
    {     
      if(a[iordre[i]] != 0)
	{
	  corder[m] = i;	  
	  m++;
	}
    }
 
  return m;
}

int filexists(char *filename)
{
  FILE *fp;
  int res;
  fp = fopen(filename,"r");
  
  if(fp) 
    {
      res = 1;
      fclose(fp);
    }
  else 
    res = 0;
       
  return res;
} 
