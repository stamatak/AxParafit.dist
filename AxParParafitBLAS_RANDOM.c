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

#include <stdio.h>
#include <stdlib.h> 
#include <string.h>
#include <mpi.h>

#include "common.h"
#include "blasmath.h"

/* ACML include */

#define COMPUTE 0
#define FINALIZE 1
#define JOB_REQUEST 2
#define RESULT 3

typedef struct {
  int i;
  int j;
} jobQueue;

typedef struct {
  int ii;
  int jj;
  double F1;
  double prob1;
  double F2;
  double prob2;
} resultVector;

int processID;
int numOfWorkers;


static int gettimeSrand()
{
  struct timeval ttime;
  gettimeofday(&ttime , NULL);
  return ttime.tv_sec + ttime.tv_usec;
}

static int randomInt(int n)
{
  return rand() %n;
}


/* Interface Stuff ***************************************************************************************************************************/

typedef struct {
  int n1;
  int n2;
  int n3;
  int n4;
  int permutations;
  int computeTracefileOnly;
  int randomizedSelection;
  int useExternalTraceFile;
  char fileNameA[1024];
  char fileNameB[1024]; 
  char fileNameC[1024];
  /*char traceFileName[1024];*/
  char externalTraceFileName[1024];
  char outFileName[1024]; 
} parameters;



static void printHelp()
{
  printf("\n\n\n\n");
  
  printf("AxParParafitBLAS_RANDOM version 1.0, February 2007 by Alexandros Stamatakis\n");
  printf("Please also consult the Manual\n");
  printf("To report bugs send an email to Alexandros.Stamatakis@epfl.ch\n\n\n");
  
  printf("AxParParafitBLAS_RANDOM -p numberOfPermutations -n1 N1 -n2 N2 -n3 N3 -n4 N4\n");
  printf("                        -A associationMatrix -B parasiteMatrix -C hostMatrix\n");
  printf("                        -n runID -t traceFileName [-h] [-r numberOfRandomLinks]\n"); 
  printf("\n");
  printf("       -p     Specify the number of permutations you want to execute.\n");  
  printf("\n");
  printf("       -n1    Specify number of Rows in associationMatrix\n");
  printf("       -n2    Specify number of Columns in associationMatrix\n");
  printf("       -n3    Specify number of Rows in hostMatrix\n");    
  printf("       -n4    Specify number of Columns in parasiteMatrix\n");
  printf("\n");
  printf("       -A     Specify file name of association matrix\n");
  printf("       -B     Specify file name of parasite matrix\n");
  printf("       -C     Specify file name of host matrix \n");
  printf("\n");
  printf("       -n     Specify a run Name/ID for this run which will be appended to the output files\n");
  printf("\n");
  printf("       -t     Specify the name of a binary trace-file written by a global cospeciation test with sequential AxParafit[BLAS].\n");
  printf("\n");
  printf("       -h     Display this help message\n");  
  printf("\n"); 
  printf("       -r     Specify the number of randomly selected individual cospeciation tests\n");
  printf("              DEFAULT: number of non-zero entries in associationMatrix (equivalent to full analysis)\n");
  printf("\n\n\n\n");
}


static void get_args(int argc, char *argv[], parameters *params, int processID)
{
  int i, k, res;
  int number;
#define NUM_OPT 12
  char *options[NUM_OPT]= {"-n1", "-n2", "-n3", "-n4", "-A", "-B", "-C", "-p", "-n", "-t", "-h", "-r"};
  char fileName[1024];
  char runName[1024];
  int set[NUM_OPT];

  /* init */

  params->n1 = 0;
  params->n2 = 0;
  params->n3 = 0;
  params->n4 = 0;
  params->useExternalTraceFile = 0;
  params->permutations = 0;
  params->computeTracefileOnly = 0;
  params->randomizedSelection = 0;
  strcpy(params->outFileName, "outfile");

  /***************/

  for(i = 0; i < NUM_OPT; i++)
    set[i] = 0;

  for(i = 1; i < argc; i++)
    {
      int found = 0;
      
      for(k = 0; k < NUM_OPT && !found; k++)	
	if(strcmp(options[k], argv[i]) == 0)	    	    	    
	  found = 1;	      	      	         
      k--;
      
      if(found)
	{	 
	  switch(k)
	    {
	    case 0:
	      if(i == argc - 1)
		{
		  printf("Error, argument expected after option %s\n", options[k]);
		  exit(-1);
		}
	      else
		{
		  res = sscanf(argv[++i],"%d", &number);
		  if(res == 0 || res == EOF)
		    {
		      printf("argument of %s must be an integer number\n", options[k]);
		      exit(-1);
		    }
		}
	      params->n1 = number;
	      set[k] = 1;		
	      break;
	    case 1:
	       if(i == argc - 1)
		{
		  printf("Error, argument expected after option %s\n", options[k]);
		  exit(-1);
		}
	      else
		{
		  res = sscanf(argv[++i],"%d", &number);
		  if(res == 0 || res == EOF)
		    {
		      printf("argument of %s must be an integer number\n", options[k]);
		      exit(-1);
		    }		 
		}
	       params->n2 = number;
	       set[k] = 1;
	       break;	    
	    case 2:
	       if(i == argc - 1)
		{
		  printf("Error, argument expected after option %s\n", options[k]);
		  exit(-1);
		}
	      else
		{
		  res = sscanf(argv[++i],"%d", &number);
		  if(res == 0 || res == EOF)
		    {
		      printf("argument of %s must be an integer number\n", options[k]);
		      exit(-1);
		    }		 
		}	
	       params->n3 = number;
	       set[k] = 1;
	       break;
	    case 3: 
	      if(i == argc - 1)
		{
		  printf("Error, argument expected after option %s\n", options[k]);
		  exit(-1);
		}
	      else
		{
		  res = sscanf(argv[++i],"%d", &number);
		  if(res == 0 || res == EOF)
		    {
		      printf("argument of %s must be an integer number\n", options[k]);
		      exit(-1);
		    }		  
		}	
	      params->n4 = number;
	      set[k] = 1;
	      break; 
	    case 4:
	       if(i == argc - 1)
		{
		  printf("Error, argument expected after option %s\n", options[k]);
		  exit(-1);
		}
	       else
		 {
		   strcpy(fileName, argv[++i]);		   
		   if(!filexists(fileName))
		     {
		       printf("file %s does not exist\n", fileName);
		       exit(-1);
		     }
		 }
	       strcpy(params->fileNameA, fileName);
	       set[k] = 1;
	       break;
	    case 5:
	       if(i == argc - 1)
		{
		  printf("Error, argument expected after option %s\n", options[k]);
		  exit(-1);
		}
	       else
		 {
		   strcpy(fileName, argv[++i]);		
		   if(!filexists(fileName))
		     {
		       printf("file %s does not exist\n", fileName);
		       exit(-1);
		     }
		 }
	       strcpy(params->fileNameB, fileName);
	       set[k] = 1;
	       break;
	    case 6:
	       if(i == argc - 1)
		{
		  printf("Error, argument expected after option %s\n", options[k]);
		  exit(-1);
		}
	       else
		 {
		   strcpy(fileName, argv[++i]);		 
		   if(!filexists(fileName))
		     {
		       printf("file %s does not exist\n", fileName);
		       exit(-1);
		     }
		 }	 
	       strcpy(params->fileNameC, fileName);
	       set[k] = 1;
	       break;
	    case 7:
	      if(i == argc - 1)
		{
		  printf("Error, argument expected after option %s\n", options[k]);
		  exit(-1);
		}
	      else
		{
		  res = sscanf(argv[++i],"%d", &number);
		  if(res == 0 || res == EOF)
		    {
		      printf("argument of %s must be an integer number\n", options[k]);
		      exit(-1);
		    }	  
		}	     	      	     
	      params->permutations = number;
	      set[k] = 1;
	      break; 	      	      
	    case 8:
	      if(i == argc - 1)
		{
		  printf("Error, argument expected after option %s\n", options[k]);
		  exit(-1);
		}
	       else
		 {
		   strcpy(runName, argv[++i]);
		  		   
		   strcat(params->outFileName,   ".");		  
		   strcat(params->outFileName,   runName);		  		   		   
		 }
	      set[k] = 1;
	      break;	      	  	        
	    case 9:
	      if(i == argc - 1)
		{
		  printf("Error, argument expected after option %s\n", options[k]);
		  exit(-1);
		}
	      else
		{
		  strcpy(fileName, argv[++i]);		
		  if(!filexists(fileName))
		    {
		       printf("file %s does not exist\n", fileName);
		       exit(-1);
		    }
		}	 
	       strcpy(params->externalTraceFileName, fileName);
	       params->useExternalTraceFile = 1;	   
	       set[k] = 1;
	       break;
	    case 10:	      
	      if(processID == 0)
		printHelp();
	      exit(0);
	      break;	   
	    case 11:	      
	      if(i == argc - 1)
		{
		  printf("Error, argument expected after option %s\n", options[k]);
		  exit(-1);
		}
	      else
		{
		  res = sscanf(argv[++i],"%d", &number);
		  if(res == 0 || res == EOF)
		    {
		      printf("argument of %s must be an integer number\n", options[k]);
		      exit(-1);
		    }	   
		}	
	       if(number <= 0)
		 {
		   printf("Randomized Selection must be a positive Integer Value\n");
		   exit(-1);
		 }
	       params->randomizedSelection = number;
	       break;
	    default:
	      printf("unknown option %s\n", options[k]);
	      exit(-1);
	    }
	}
      else
	{
	  printf("unknown option %s\n", argv[i]);
	  exit(-1);
	}     
    }


  for(i = 0; i < 10; i++)
    {
      if(set[i] == 0)
	{
	  printf("Error option %s must be specified!\n", options[i]);
	  exit(-1);
	}
    }


  if(params->useExternalTraceFile == 1 && params->computeTracefileOnly == 1)
    {
      printf("Conflicting Options, you are trying to read in a tracfile with \"-t\" while\n");
      printf("you only want to compute the global significance with \"-g\" \n");
      exit(-1);
    }
    
  if(processID == 0 && filexists(params->outFileName))
    {
      printf("Output File %s already exists\n", params->outFileName);
      exit(-1);
    }    

}


/************************* Interface end ***************************************************************************************************/

int main (int argc, char *argv[])
{
  int n1, n2, n3, n4;
  int numberOfPermutations;
  int i, j, k, l, p, ii, jj, m; 
  char *a; 
  int *iordre;
  int *corder;
  double *b;
  double *c;
  double *CAprime;
  double ssqeigvB;
  double ssqeigvC;
  double tracetot;
  double temp, prob;
  double *traceper;
  double dACC;
  FILE *f, *outf; 
  long idum; 
  int iGE, iGEF1, iGEF2, dummyINT;
  char fileNameA[1024], 
    fileNameB[1024], 
    fileNameC[1024],
    outFileName[1024]; 
  double F1, F2, trace0, prob1, prob2, F1per, F2per;
  double time;
  MPI_Status msgStatus; 
  int jobMsg[3];
  double resultMsg[7];
  int jobID;
  int nonZero;
  parameters *params;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &processID);
  MPI_Comm_size(MPI_COMM_WORLD, &numOfWorkers);
  
  params = (parameters *)malloc(sizeof(parameters));
  
  nonZero = 0;
  get_args(argc, argv, params, processID);

  numberOfPermutations = params->permutations;
  n1 = params->n1;
  n2 = params->n2;
  n3 = params->n3;
  n4 = params->n4;
  strcpy(fileNameA, params->fileNameA);
  strcpy(fileNameB, params->fileNameB);
  strcpy(fileNameC, params->fileNameC);
  strcpy(outFileName, params->outFileName);

  if(processID == 0)
    {
      time = gettime();  
      outf = fopen(outFileName, "w");   
      fprintf(outf, "Permutations: %d N1 %d N2 %d, N3 %d N4 %d\n", numberOfPermutations - 1, n1, n2, n3, n4);
    }   
  
  traceper = (double *)malloc(sizeof(double) * numberOfPermutations);

  iordre = (int *) malloc (sizeof(int) * n2);
  corder = (int *) malloc (sizeof(int) * n2);
  
  a = (char *)malloc(sizeof(char) * n1 * n2); 

  b = (double *)malloc(sizeof(double) * n4 * n1);
  
  CAprime = (double *)malloc(sizeof(double) * n3 * n1);  

  c = (double *)malloc(sizeof(double) * n3 * n2);

  /******** READ DATA ***********************************************/

  f = fopen(fileNameA, "r");
    
  for(i = 0; i < n1; i++)    
    for(j = 0; j < n2; j++)	      
      {
	int d, v;
	v = fscanf(f, "%d", &d);
	if(v == 0)
	  {
	    printf("Format Conversion Error while reading Matrix A(%s) at position A[%d][%d]\n", fileNameA, i, j);
	    exit(-1);
	  }
	if(v == EOF)
	  {
	    printf("End of File reached while reading Matrix A(%s) at position A[%d][%d]\n", fileNameA, i, j);
	    exit(-1);
	  }
	a[i * n2 + j] = (char)d;
	if(a[i * n2 + j] != 0)
	  nonZero++;	
      }
  
  fclose(f);   
 
  f = fopen(fileNameB, "r");
  
  for(i = 0; i < n1; i++)   
    for(j = 0; j < n4; j++)      
      {
	int v;
#ifdef ROWS
	v = fscanf(f, "%lf",&b[j * n1 + i]);
#else
	v = fscanf(f, "%lf",&b[i * n4 + j]);
#endif
	if(v == 0)
	  {
	    printf("Format Conversion Error while reading Matrix B(%s) at position B[%d][%d]\n", fileNameB, i, j);
	    exit(-1);
	  }
	if(v == EOF)
	  {
	    printf("End of File reached while reading Matrix B(%s) at position B[%d][%d]\n", fileNameB, i, j);
	    exit(-1);
	  }
      }
  
  

  fclose(f);

  f = fopen(fileNameC, "r");
  
  for(i = 0; i < n3; i++)        
    for(j = 0; j < n2; j++)	
      {
	int v;
	v = fscanf(f, "%lf",&c[i * n2 + j]);    
	if(v == 0)
	  {
	    printf("Format Conversion Error while reading Matrix C(%s) at position C[%d][%d]\n", fileNameC, i, j);
	    exit(-1);
	  }
	if(v == EOF)
	  {
	    printf("End of File reached while reading Matrix C(%s) at position C[%d][%d]\n", fileNameC, i, j);
	    exit(-1);
	  }
      }
	
  fclose(f);       
  ssqeigvB = 0.0; 

#ifdef ROWS
  for(i = 0; i < n4; i++)
    {
      temp = 0.0;
      for(j = 0; j < n1; j++)     
	temp += b[i * n1 + j] * b[i * n1 + j];      
      ssqeigvB += temp * temp;
    }
#else
  for(i = 0; i < n4; i++)
    {
      temp = 0.0;
      for(j = 0; j < n1; j++)     
	temp += b[j * n4 + i] * b[j * n4 + i];      
      ssqeigvB += temp * temp;
    }
#endif

  ssqeigvC = 0.0;

  for(i = 0; i < n3; i++)
    {
      temp = 0.0;
      for(j = 0; j < n2; j++)
	temp += c[i * n2 + j] * c[i * n2 + j];
      ssqeigvC += temp * temp;
    }    
 
  if(processID == 0)
    {
      fprintf(outf, "Sum of squared PCoA eigenvalues of B = %1.5f\n\n", ssqeigvB);
      fprintf(outf, "Sum of squared PCoA eigenvalues of C = %1.5f\n\n", ssqeigvC);
    }

  if(ssqeigvC > ssqeigvB) 
    tracetot = ssqeigvC;
  else
    tracetot = ssqeigvB;

  if(processID == 0)
    fprintf(outf, "TraceTot = %1.5f\n\n", tracetot);

  {
    FILE *t;
    int readCount;

    if(processID == 0)
      printf("READING trace file %s\n", params->externalTraceFileName);

    t = fopen(params->externalTraceFileName, "r");      
    readCount = fread(((void *)traceper), sizeof(double), numberOfPermutations, t);       
    fclose(t);           
      
    iGE = 1;

    for(p = 1; p < numberOfPermutations; p++)
      {                            	 	  	  
	if(traceper[p] >= traceper[0])
	  iGE++;          	 
      }
      
    prob = (double)(iGE) / (double)(numberOfPermutations);
    if(processID == 0)
      {
	fprintf(outf, " Global test of cospeciation:                     ParaFitGlobal = %1.5f   Prob  = %1.5f\n\n", traceper[0], prob);       
	fprintf(outf, " Test of individual host-parasite links:\n\n");  
	fprintf(outf, "                        F1 = ParaFitLink1                  F2 = ParaFitLink2\n\n\n");     
	printf("Global test of cospeciation:                     ParaFitGlobal = %1.5f   Prob  = %1.5f\n\n", traceper[0], prob);
      }
  }
  
  if(processID == 0)
    {
      int count;
      jobQueue *jobs;     
      int jobsSent, jobsReceived;         
      int resultCounter = 0;         

      jobs = (jobQueue *)malloc(nonZero * sizeof(jobQueue));
            
      count = 0;
      for(i = 0; i < n1; i++)
	for(j = 0; j < n2; j++)
	  {
	    if(a[i * n2 + j] != 0)
	      {
		jobs[count].i = i;
		jobs[count].j = j;
		count++;
	      }
	  }
     
      /* permute the job vector a bit */

      srand((unsigned int) gettimeSrand());
	
      for(i = 0; i < 10 * count; i++)
	{
	  int e_i, e_j, p, q;
	  
	  p = randomInt(count);
	  q = randomInt(count);
	  
	  e_i = jobs[q].i;
	  e_j = jobs[q].j;
	  jobs[q].i = jobs[p].i;
	  jobs[q].j = jobs[p].j;
	  jobs[p].i = e_i;
	  jobs[p].j = e_j;
	}       

      /*for(i = 0; i < count; i++)
	{
	  printf("%d %d %d\n", i, jobs[i].i, jobs[i].j);
	  }*/           
      

      fclose(outf);

      jobsReceived = nonZero;

      if(params->randomizedSelection > 0)
	{
	  if(params->randomizedSelection > jobsReceived)
	    {
	      printf("Number %d of randomly selected individual cospeciation tests higher than \n", params->randomizedSelection);
	      printf("number %d of non-zero entries in association matrix, resetting it to %d\n", nonZero, nonZero);
	    }
	  else	    
	    jobsReceived = params->randomizedSelection;
	}

      jobsSent = 0;
      while(jobsReceived > 0)
	{
	  MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &msgStatus);
	  switch(msgStatus.MPI_TAG)
	    {
	    case JOB_REQUEST: 
	      MPI_Recv(&dummyINT, 1, MPI_INT, msgStatus.MPI_SOURCE, JOB_REQUEST, MPI_COMM_WORLD, &msgStatus);	     
	      if(jobsSent < nonZero)
		 {
		   jobMsg[0] = jobsSent;
		   jobMsg[1] = jobs[jobsSent].i;
		   jobMsg[2] = jobs[jobsSent].j;		   
		   MPI_Send(jobMsg, 3, MPI_INT, msgStatus.MPI_SOURCE, COMPUTE, MPI_COMM_WORLD);
		   jobsSent++;
		 }
	       break;
	    case RESULT:
	      MPI_Recv(resultMsg, 7, MPI_DOUBLE, msgStatus.MPI_SOURCE, RESULT, MPI_COMM_WORLD, &msgStatus);	     
	      jobsReceived--;
	      jobID                = (int)resultMsg[0];
	      {
		int ii, jj;
		double F1, prob1, F2, prob2;
		ii    = (int)resultMsg[1];
		jj    = (int)resultMsg[2];
		F1    = resultMsg[3];
		prob1 = resultMsg[4];
		F2    = resultMsg[5];
		prob2 = resultMsg[6];	     	      	      	     
	      	     
		printf("Parasite  %d  Host %d   F1 =   %1.5f   Prob1 =  %1.5f   F2 =    %1.5f   Prob2 =  %1.5f\n", 
		       ii + 1,  jj + 1,  F1, prob1, F2, prob2);

		outf = fopen(outFileName, "a");
		fprintf(outf, "Parasite  %d  Host %d   F1 =   %1.5f   Prob1 =  %1.5f   F2 =    %1.5f   Prob2 =  %1.5f\n", 
		       ii + 1,  jj + 1,  F1, prob1, F2, prob2);
		fclose(outf);
		      
	      }

	      if(jobsSent < nonZero)
		{
		  jobMsg[0] = jobsSent;
		  jobMsg[1] = jobs[jobsSent].i;
		  jobMsg[2] = jobs[jobsSent].j;
		  MPI_Send(jobMsg, 3, MPI_INT, msgStatus.MPI_SOURCE, COMPUTE, MPI_COMM_WORLD);
		  jobsSent++;
		}
	      break;	      
	    }
	}     

      printf("There are %d host-parasite links in matrix A\n", nonZero);

      outf = fopen(outFileName, "a");	       
      fprintf(outf, "There are %d host-parasite links in matrix A\n", nonZero);    
      fclose(outf);

      for(i = 1; i < numOfWorkers; i++)
	{
	  MPI_Send(&dummyINT, 1, MPI_INT, i, FINALIZE, MPI_COMM_WORLD);
	}      

      printf("TIME %f\n", gettime() - time);
      goto FINISH;
    }
  else
    {
      MPI_Send(&dummyINT, 1, MPI_INT, 0, JOB_REQUEST, MPI_COMM_WORLD);

      while(1)
	{		  
	  MPI_Probe(0, MPI_ANY_TAG, MPI_COMM_WORLD, &msgStatus);	  
	  switch(msgStatus.MPI_TAG)
	    {
	    case COMPUTE:
	      MPI_Recv(jobMsg, 3, MPI_INT, 0, COMPUTE, MPI_COMM_WORLD, &msgStatus);	     
	      jobID = jobMsg[0];
	      ii = jobMsg[1];
	      jj = jobMsg[2];

	      a[ii * n2 + jj] = 0;
	      	     
	      makeCAprime(n3, n1, n2, a, c, CAprime, corder);
	      dACC = matMult(n3, n4, n1, CAprime, b);	    

	      F1 = (traceper[0] - dACC);
	      F2 = (traceper[0] - dACC)/(tracetot - traceper[0]);	  

	      for(i = 0; i < n2; i++)
		iordre[i] = i;
	      idum = -1;   
	      for(i = 0; i < NTURN; i++)
		ran2(&idum);

	      iGEF1 = 1;
	      iGEF2 = 1;
	      
	      for(p = 1; p < numberOfPermutations; p++)
		{		 
		  permuteCAprime(n3, n1, n2, a, c, CAprime, &idum, iordre, corder);
		  dACC = matMult(n3, n4, n1, CAprime, b);		  		 

		  F1per = traceper[p] - dACC;
		  F2per = (traceper[p] - dACC) / (tracetot - traceper[p]);		 

		  if(F1per >= F1)
		    iGEF1++;
		  if(F2per >= F2)
		    iGEF2++;		 				 
		}
	      prob1 = (double)(iGEF1) / (double)(numberOfPermutations);
	      prob2 = (double)(iGEF2) / (double)(numberOfPermutations);	     	      

	      a[ii * n2 + jj]    = 1; 
	      resultMsg[0] = (double)jobID;
	      resultMsg[1] = (double)ii;
	      resultMsg[2] = (double)jj;
	      resultMsg[3] = F1;
	      resultMsg[4] = prob1;
	      resultMsg[5] = F2;
	      resultMsg[6] = prob2;
	     	    
	      MPI_Send(resultMsg, 7, MPI_DOUBLE, 0, RESULT, MPI_COMM_WORLD);
	      break;
	    case FINALIZE:
	      MPI_Recv(&dummyINT, 1, MPI_INT, 0, FINALIZE, MPI_COMM_WORLD, &msgStatus);
	      goto FINISH;
	    }
	}
    }
 FINISH:
  MPI_Finalize();   	        
}


   
