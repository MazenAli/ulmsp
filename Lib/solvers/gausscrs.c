#ifndef GAUSSCRS_C
#define GAUSSCRS_C

#include <stdlib.h>

#include "gausscrs.h"

void
gausscrs(pccrs A,
         prealvector b,
         pcindexvector fixed)
{
/*
 *  this subroutine solves the linear system of equations
 *                 A x = b
 *  by reducing the coefficient matrix a to a diagonal matrix.
 */

  index i,j,k,n, ptr;
  real recip, *T;

/**********************************************/

  n    = b->length;
  /* Allocate memory for matrix A in full format */
  T = (real *) calloc(n*n,sizeof(crs));
  /* Copy crs matrix to full format, store A row-wise */
  for (i = 0 ; i < n ; i++) {
    for(j = A->rowptr[i]; j < A->rowptr[i+1]; j++) {
      T[i*n+A->colind[j]-INDEX_BASE] = A->vals[j];
    }
  }
  /* Set fixed rows to canonical unit vector */
  for (i = 0 ; i<fixed->length ; i++) {
    ptr = fixed->vals[i]-INDEX_BASE;
    for (j = 0 ; j<n ; j++) {
      T[n*ptr+j] = 0.0;
    }
    T[(n+1)*ptr] = 1.0;
  }

  /* consider k-th column */
  for (k = 0; k < n-1 ; k++){
    /* modify entries in A, i.e. create LU decomp*/
    recip = 1.0/T[k*n+k];
    for (j = k+1; j<n; j++){ 
      /* compute  l_jk, store in A */
      T[j*n+k] *= recip;  
      for (i = k+1; i<n; i++){
        /* mofify (n-k-1)*(n-k-1) submatrix */
        T[j*n+i] -= T[j*n+k] * T[k*n+i];
      }
      /* compute L y = b and store in b */
      b->vals[j] -= T[j*n+k] * b->vals[k];
    }
  }

  /* solve R x = y and store in y */
  k = n;
  do {
    k--;
    for (j = k+1; j < n ; j++){
      b->vals[k]-=T[k*n+j]*b->vals[j];
    }
    b->vals[k]/=T[k*n+k];
  } while (k);

  free(T);
}

#endif
