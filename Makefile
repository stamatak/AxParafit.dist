# Axelerated Parafit and Axelerated Parallel Parafit Makefile 
# Copyright January 2007 by Alexandros Stamatakis

CC = gcc 

MPICC = mpicc

#added -D_HAND_OPTIMIZED_MM to compile poormath.c with the hand-tuned version of
#dense matric mult by default

CFLAGS = -O3 -fomit-frame-pointer -funroll-loops -D_HAND_OPTIMIZED_MM

######################################################################################################
# SPECIFY 32 or 64 bit architecture

#ARCH=64
ARCH=32
ARCHFLAGS=-m$(ARCH)

######################################################################################################
# STATIC or DYNAMIC linking of libraries

#LINKERFLAGS=-static
LINKERFLAGS=

######################################################################################################
# AMD ACML BLAS-specific FLAGS
# 
# ATTENTION ! ATTENTION ! adapt BLAS and BLAS_CFLAGS to your local installation/path

# LINUX/UNIX Alexi's settings on the cluster
#BLAS_CFLAGS = -O3 -fomit-frame-pointer -funroll-loops -D_USE_ACML -I/home/stamatak/BLAS/gnu$(ARCH)/include/
#BLAS = -L/home/stamatak/BLAS/gnu$(ARCH)/lib/ -lacml -lm -lg2c

# LINUX/UNIX Tuebingen settings
#BLAS_CFLAGS = -O3 -fomit-frame-pointer -funroll-loops -D_USE_ACML -I../acml3.6.0/gnu$(ARCH)/include/
#BLAS = -L../acml3.6.0/gnu$(ARCH)/lib/ -lacml -lm -lg2c

# WINDOWS SETTINGS
#BLAS_CFLAGS = -O3 -fomit-frame-pointer -funroll-loops -D_USE_ACML -Ic:/msys/1.0/mingw/acml3.6.0/gnu32/include
#BLAS = -Lc:/msys/1.0/mingw/acml3.6.0/gnu32/lib -lacml -lm -lg2c -static


######################################################################################################
# Intel MKL BLAS-specific FLAGS
#
# ATTENTION ! ATTENTION ! adapt BLAS and BLAS_CFLAGS to your local installation/path

#BLAS_CFLAGS = -O3 -fomit-frame-pointer -funroll-loops -D_USE_MKL -I../intel-mkl/include/
#BLAS = -L../intel-mkl/lib/$(ARCH) -lmkl -lmkl_lapack -lm -lguide -shared


######################################################################################################
# ATLAS BLAS-specific FLAGS
# ATTENTION ! ATTENTION ! adapt BLAS and BLAS_CFLAGS to your local installation/path

BLAS_CFLAGS = -O3 -fomit-frame-pointer -funroll-loops -D_USE_ATLAS -I/root/ATLAS-LIB/Linux_P4SSE2/include/
BLAS = -L/root/ATLAS-LIB/Linux_P4SSE2/lib/ -lm   -lcblas -latlas 

######################################################################################################

RM = rm -f

STRIP = strip

all : AxParafit  AxParafitBLAS AxParParafit AxParParafitBLAS AxParParafitBLAS_RANDOM

AxParafit : common.o AxParafit.o poormath.o
	$(CC) $(ARCHFLAGS) $(LINKERFLAGS) -o AxParafit AxParafit.o common.o poormath.o

AxParParafit : AxParParafit.o common.o poormath.o
	$(MPICC) -o AxParParafit AxParParafit.o common.o poormath.o $(LINKERFLAGS) $(ARCHFLAGS)

AxParParafit.o : AxParParafit.c
	$(MPICC) $(ARCHFLAGS) -c AxParParafit.c 

######################################################################################################

blasmath.o:	blasmath.c
	$(CC) $(ARCHFLAGS) $(BLAS_CFLAGS) -c blasmath.c

AxParParafitBLAS : AxParParafitBLAS.o common.o blasmath.o
	$(MPICC) -o AxParParafitBLAS AxParParafitBLAS.o common.o blasmath.o  $(BLAS) $(ARCHFLAGS) $(LINKERFLAGS)

AxParParafitBLAS.o : AxParParafitBLAS.c
	$(MPICC) $(ARCHFLAGS) $(BLAS_CFLAGS) -c  AxParParafitBLAS.c 


AxParParafitBLAS_RANDOM : AxParParafitBLAS_RANDOM.o common.o blasmath.o
	$(MPICC) -o AxParParafitBLAS_RANDOM AxParParafitBLAS_RANDOM.o common.o blasmath.o $(BLAS) $(ARCHFLAGS) $(LINKERFLAGS)

AxParParafitBLAS_RANDOM.o : AxParParafitBLAS_RANDOM.c
	$(MPICC) $(ARCHFLAGS) $(BLAS_CFLAGS) -c  AxParParafitBLAS_RANDOM.c 

AxParafitBLAS : AxParafitBLAS.o common.o blasmath.o
	$(CC) -o AxParafitBLAS AxParafitBLAS.o common.o blasmath.o $(BLAS) $(ARCHFLAGS) $(LINKERFLAGS)

AxParafitBLAS.o : AxParafitBLAS.c
	$(CC) $(ARCHFLAGS) $(BLAS_CFLAGS) -c  AxParafitBLAS.c

%.o: %.c
	$(CC) $(ARCHFLAGS) $(CFLAGS) -o $@ -c $<

clean : 
	$(RM) *.o
	$(RM) AxParafit AxParParafit AxParafitBLAS AxParParafitBLAS AxParParafitBLAS_RANDOM
	$(RM) AxParafit.exe AxParParafit.exe AxParafitBLAS.exe AxParParafitBLAS.exe AxParParafitBLAS_RANDOM.exe

strip :
	$(STRIP) AxParafit AxParParafit AxParafitBLAS AxParParafitBLAS AxParParafitBLAS_RANDOM
	$(STRIP) AxParafit.exe AxParParafit.exe AxParafitBLAS.exe AxParParafitBLAS.exe AxParParafitBLAS_RANDOM.exe
