/*
 * svdDynamic.h
 * Copyright (c) 2000
 * Thomas F. El-Maraghi
 *
 * Singular value decomposition (SVD) routines.
 *
 */

// Modified for stand-alone use, FEG, Jul 18, 2006

#include<math.h>
#include<stdio.h>
#include<stdlib.h>

#ifndef __SVD_dynamic
#define __SVD_dynamic

#define MAX_SVD_ITERATIONS 100
#define max(A,B) ((A)<(B)?(B):(A))

#define signof(A,B)    (((B)>=0)? (fabs(A)) : (-fabs(A)))

int SVDHelper( const int m, const int n,
	       double *U, double *w, double *V,
	       double *rv1 );
static double SVD_PYTHAG( const double a, const double b );
int SVD( const double *A, const int m, const int n,
	 double **U, double **w, double **V, double **rv1 );
void SortSV( int *svPerm, double *w, const int n );
int SolveLinearSystem( const double *A, const double *b,
		       const int m, const int n,
		       double **x, double **w );
void InvertMatrix( const double *U, const double *w, const double *V,
		   const int n, double *I );
#endif

