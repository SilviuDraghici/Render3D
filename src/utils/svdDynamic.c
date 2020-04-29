/*
 * svdDynamic.c
 * Copyright (c) 2000
 * Thomas F. El-Maraghi
 *
 * Singular value decomposition (SVD) routines.
 *
 */

// Modified for stand-alone use, FEG, Jul 18, 2006

#include "svdDynamic.h"

#define signof(A,B)    (((B)>=0)? (fabs(A)) : (-fabs(A)))

/*
 * Returns the Singular Value Decomposition U*w*V^T of the matrix A
 * (passed to the function in U).
 *
 * U orthogonal mxn matrix (mxn, m >= n)
 * w n-vector of signular values (nx1)
 * V orthogonal nxn matrix (nxn)
 * rv1  superdiagonal of singular value matrix
 *
 * This subroutine is a translation of the algol procedure svd, 
 * num. math. 14, 403-420(1970) by golub and reinsch. 
 * handbook for auto. comp., vol ii-linear algebra, 134-151(1971). 
 * See http://www.netlib.org/  for Eispack svd.f Fortran77 version 
 */
int SVDHelper( const int m, const int n,
	       double *U, double *w, double *V, 
	       double *rv1 )
{
  int flag, i, its, j, jj, k, l = 0, nm = 0;
  double c, f, h, s, x, y, z;
  double anorm = (double)0;
  double g = (double)0;
  double scale = (double)0;
  double tst;

  if( m < n )
  {
//    FATAL( "SVD: You must augment A with extra zero rows" );
   fprintf(stderr,"SVD: Thou shouldst augment A with extra zero rows, silly!\n");
   exit(0);
  }
  
  for( i = 0; i < n; i++ ) {
    l = i + 1;
    rv1[i] = (double)(scale * g);
    g = s = scale = (double)0;
    if( i < m ) {
      for( k = i; k < m; k++ ) 
	scale += fabs( U[k*n+i] );
      if( scale ) {
	for( k = i; k < m; k++ ) {
	  U[k*n + i] /= scale;
	  s += (double)(U[k*n+i] * U[k*n+i]);
	}
	f = U[i*n + i];
	g = (double)-signof(sqrt(s),f);
	h = (double)(f * g - s);
	U[i*n+i] = (double)(f - g);
	if( i != n - 1 ) {
	  for( j = l; j < n; j++ ) {
	    for( s = (double)0, k = i; k < m; k++ ) 
	      s += (double)(U[k*n+i] * U[k*n+j]);
	    f = (double)(s / h);
	    for( k = i; k < m; k++ )
	      U[k*n+j] += (double)(f * U[k*n+i]);
	  }
	}
	for( k = i; k < m; k++ ) 
	  U[k*n+i] *= scale;
      }
    }
    w[i] = (double)(scale * g);
    g = s = scale = (double)0;
    if( i < m && i != n - 1) {
      for( k = l; k < n; k++ )
	scale += fabs(U[i*n+k]);
      if( scale ) {
	for( k = l; k < n; k++ ) {
	  U[i*n+k] /= scale;
	  s += (double)(U[i*n+k] * U[i*n+k]);
	}
	f = U[i*n+l];
	g = (double)-signof(sqrt(s),f);
	h = (double)(f * g - s);
	U[i*n+l] = f - g;
	for( k = l; k < n; k++ ) 
	  rv1[k] = (double)(U[i*n+k] / h);
	if( i != m - 1) {
	  for( j = l; j < m; j++ ) {
	    for( s = (double)0, k = l; k < n; k++ )
	      s += (double)(U[j*n+k] * U[i*n+k]);
	    for( k = l; k < n; k++ )
	      U[j*n+k] += (double)(s * rv1[k]);
	  }
	}
	for( k = l; k < n; k++ )
	  U[i*n+k] *= scale;
      }
    }
    anorm = max( anorm, (fabs(w[i]) + fabs(rv1[i])) );
  }
  for( i = n - 1; i >= 0; i-- ) {
    if( i < n - 1 ) {
      if( g ) {
	for( j = l; j < n; j++ )
	  V[j*n+i] = (double)((U[i*n+j] / U[i*n+l]) / g);
	for( j = l; j < n; j++ ) {
	  for( s = (double)0, k = l; k < n; k++ )
	    s += (double)(U[i*n+k] * V[k*n+j]);
	  for( k = l; k < n; k++ )
	    V[k*n+j] += (double)(s * V[k*n+i]);
	}
      }
      for( j = l; j < n; j++ )
	V[i*n+j] = V[j*n+i] = (double)0;
    }
    V[i*n+i] = (double)1;
    g = rv1[i];
    l = i;
  }
  for( i = n - 1; i >= 0; i-- ) {
    l = i + 1;
    g = w[i];
    if( i < n - 1 )
      for( j = l; j < n; j++ )
	U[i*n+j] = (double)0;
    if( g ) {
      g = (double)((double)1 / g);
      if( i != n - 1 ) {
	for( j = l; j < n; j++ ) {
	  for( s = (double)0, k = l; k < m; k++ )
	    s += (double)(U[k*n+i] * U[k*n+j]);
	  f = (double)((s / U[i*n+i]) * g);
	  for( k = i; k < m; k++ )
	    U[k*n+j] += (double)(f * U[k*n+i]);
	}
      }
      for( j = i; j < m; j++ )
	U[j*n+i] *= g;
    } else {
      for( j = i; j < m; j++ )
	U[j*n+i] = (double)0;
    }
    ++U[i*n+i];
  }

  for( k = n - 1; k >= 0; k-- ) {
    for( its = 1; its <= MAX_SVD_ITERATIONS; its++ ) {
      flag = 1;
      for( l = k; l >= 0; l-- ) {
	nm = l - 1;
	tst = fabs(rv1[l]) + anorm;
	if( tst == anorm ) {
	  flag = 0;
	  break;
	}
	tst = fabs(w[nm]) + anorm;
	if( tst == anorm )
	  break;
      }
      /* Found a zero diagonal element w[nm] */
      if( flag ) {
	c = (double)0;
	s = (double)1;
	for( i = l; i <= k; i++ ) {
	  f = (double)(s * rv1[i]);
	  rv1[i] = (double)(c * rv1[i]);
	  tst = fabs(f) + anorm;
	  if( tst == anorm ) 
	    break;
	  else {
	    g = w[i];
	    h = SVD_PYTHAG(f,g);
	    w[i] = h;
	    h = (double)((double)1 / h);
	    c = (double)(g * h);
	    s = (double)(-f * h);
	    for( j = 0; j < m; j++ ) {
	      y = U[j*n+nm];
	      z = U[j*n+i];
	      U[j*n+nm] = (double)(y * c + z * s);
	      U[j*n+i] = (double)(z * c - y * s);
	    }
	  }
	}
      }
      z = w[k];
      if( l == k ) {
     	if( z < (double)0 ) {
	  w[k] = (double)-z;
	  for( j = 0; j < n; j++ )
	    V[j*n+k] = (double)(-V[j*n+k]);
	}
	break;
      }
      if( its >= MAX_SVD_ITERATIONS ) {
	return( -1 );
      }
      
      x = w[l];
      nm = k - 1;
      y = w[nm];
      g = rv1[nm];
      h = rv1[k];
      f = (double)(0.5 * (((g + z) / h) * ((g - z) / y) + y / h - h / y));
      g = SVD_PYTHAG(f,(double)1);
      f = (double)(x - (z / x) * z + (h / x) * (y / (f + signof(g,f)) - h));
      c = s = (double)1;
      for( j = l; j <= nm; j++ ) {
	i = j + 1;
	g = rv1[i];
	y = w[i];
	h = (double)(s * g);
	g = (double)(c * g);
	z = SVD_PYTHAG(f,h);
	rv1[j] = z;
	c = (double)(f / z);
	s = (double)(h / z);
	f = (double)(x * c + g * s);
	g = (double)(g * c - x * s);
	h = (double)(y * s);
	y = (double)(y * c);
	for( jj = 0; jj < n; jj++ ) {
	  x = V[jj*n+j];
	  z = V[jj*n+i];
	  V[jj*n+j] = (double)(x * c + z * s);
	  V[jj*n+i] = (double)(z * c - x * s);
	}
	z = SVD_PYTHAG(f,h);
	w[j] = z;
	if( z != (double)0 ) {
	  c = (double)(f / z);
	  s = (double)(h / z);
	} 
	f = (double)((c * g) + (s * y));
	x = (double)((c * y) - (s * g));
	for( jj = 0; jj < m; jj++ ) {
	  y = U[jj*n+j];
	  z = U[jj*n+i];
	  U[jj*n+j] = (double)(y * c + z * s);
	  U[jj*n+i] = (double)(z * c - y * s);
	}
      }
      rv1[l] = (double)0;
      rv1[k] = f;
      w[k] = x;
    }
  }
  return( 0 );
}

static double SVD_PYTHAG( const double a, const double b )
{
  double at = fabs(a);
  double bt = fabs(b);
  double ct;
  if( at > bt ) {
    ct = bt / at;
    return (double)( at * sqrt( 1.0 + ct * ct ) );
  }
  if( bt <= 0.0 )
    return 0.0;
  ct = at / bt;
  return (double)( bt * sqrt( 1.0 + ct * ct ) );
}


/*
 * Returns the Singular Value Decomposition A = U*w*V^T of the matrix A
 *
 * A input mxn matrix
 * U orthogonal mxn matrix (mxn)
 * w n-vector of signular values (nx1)
 * V orthogonal nxn matrix (nxn)
 * rv1  superdiagonal of singular value matrix
 */
int SVD( const double *A, const int m, const int n,
	 double **U, double **w, double **V, double **rv1 )
{
  double *tmp;
  int r, c, k, svdReturn;
  int N = m * n;
  int ownRv1 = 0;

  if( m >= n ) {

    /* allocate memory if necessary */
    if( *U == NULL ) *U = (double *)calloc(N,sizeof(double));
    if( *w == NULL ) *w = (double *)calloc(n,sizeof(double));
    if( *V == NULL ) *V = (double *)calloc(n*n,sizeof(double));

    if( *rv1 == NULL ) {
      *rv1 = (double *)calloc(n,sizeof(double));
      ownRv1 = 1;
    }

    /* copy A to U if necessary */
    if( *U != A )
      for( k = 0; k < N; k++ )
	(*U)[k] = A[k];

    svdReturn = SVDHelper( m, n, *U, *w, *V, *rv1 );
    
  } else {
    
    /* allocate memory if necessary */
    if( *U == NULL ) *U = (double *)calloc(N,sizeof(double));
    if( *w == NULL ) *w = (double *)calloc(m,sizeof(double));
    if( *V == NULL ) *V = (double *)calloc(m*m,sizeof(double));

    if( *rv1 == NULL ) {
      *rv1 = (double *)calloc(m,sizeof(double));
      ownRv1 = 1;
    }

    /* transpose A and store in U */
    tmp = *U;
    for( c = 0; c < n; c++ )
      for( r = 0; r < m; r++, tmp++ ) 
	*tmp = A[r*n + c];

    svdReturn = SVDHelper( n, m, *U, *w, *V, *rv1 );
    
    /* swap U and V */
    tmp = *U;
    *U = *V;
    *V = tmp;
  }   

  if( ownRv1 ) free(*rv1);

  return svdReturn;
}

/*
 * Sort singular values into decreasing order,
 * return permutation array svPerm *... ith sorted
 * singular value is then w[svPerm[i]] 
 */
void SortSV( int *svPerm, double *w, const int n )
{
  int i, j, iTmp;
  double tmp;

  for( i = 0; i < n; i++ )
    svPerm[i] = i;

  for( i = 0; i < n; i++ ) {
    /* Find max in remaining set i..numSamp */
    tmp = w[svPerm[i]];
    iTmp = i;
    for( j = i + 1; j < n; j++ )
      if( w[svPerm[j]] > tmp ) {
	tmp = w[svPerm[j]];
	iTmp = j;
      }
    /* Switch */
    j = svPerm[i];
    svPerm[i] = svPerm[iTmp];
    svPerm[iTmp] = j;
  }
}


/*
 * Solve the linear system Ax=b
 */
int SolveLinearSystem( const double *A, const double *b,
		       const int m, const int n,
		       double **x, double **w )
{ 
  double *U, *V, *s;
  int i, j, svdReturn;

  if( *x == NULL ) *x = (double *)calloc(n,sizeof(double));

  svdReturn = SVD( A, m, n, &U, w, &V, &s );

  for( i = 0; i < n; i++ ) {      
    s[i] = 0.0;
    for( j = 0; j < n; j++ )          
      s[i] += U[j*n+i] * b[j-1];                                   
    s[i] /= (*w)[i];                        
  }                                                    
                                   
  for( i = 0; i < n; i++ ) {        
    (*x)[i-1] = 0.0;      
    for( j = 0; j < n; j++ )                                            
      (*x)[i-1] += V[i*n+j] * s[j];
  }
 
  free(U);
  free(V);
  free(s);

  return( svdReturn );                  
}
      

/*
 * Given svd decomposition U w V^T of a matrix, compute
 * the inverse (or pseudoinverse) I = V w^-1 U^T 
 */
void InvertMatrix( const double *U, const double *w, const double *V, 
		   const int n, double *I )
{ 
  int i, j, k;
  double *scr;

  scr = (double *)calloc(n,sizeof(double));
  
  for( k = 0; k < n; k++ ) {
    /* Compute scr = kth column of (w^-1 U^T)  */
    for( i = 0; i < n; i++ )
      scr[i] = (double)(U[k*n+i] / w[i]);
      
    /* Compute kth col of Ainv = V * scr */ 
    for( i = 0; i < n; i++ ) {
      I[i*n+k] = (double)0;      
      for( j = 0; j < n; j++ )
	I[i*n+k] += (double)(V[i*n+j] * scr[j]);
    }
  }
  
  free(scr);
}
