#ifndef REALMATRIX_C
#define REALMATRIX_C

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>

/* Conflicting typedefs for "index" */
#define index string_index
    #include <string.h>
#undef index

#include <math.h>

#include "realmatrix.h"

prealmatrix
new_realmatrix(index rows, index cols)
{
    prealmatrix A;
    A = (prealmatrix) malloc(sizeof(realmatrix));
    init_realmatrix(A, rows, cols);

    return A;
}


void
init_realmatrix(prealmatrix A, index rows, index cols)
{
    index nz;

    assert(A);

    nz = rows * cols;
    if (nz>0) {
        A->vals = (real*) malloc(nz*sizeof(real));
    } else {
        A->vals = NULL;
    }
    A->rows = rows;
    A->cols = cols;
}


void
resize_realmatrix(prealmatrix A, index rows, index cols)
{
    assert(A);

    del_realmatrix(A);
    init_realmatrix(A, rows, cols);
}


void
copy_realmatrix(prealmatrix dest, pcrealmatrix src)
{
    assert(dest);
    assert(src);

    if ((dest->rows!=src->rows)||(dest->cols!=src->cols)) {
        del_realmatrix(dest);
        init_realmatrix(dest, src->rows, src->cols);
    }

    memcpy(dest->vals, src->vals, (dest->rows*dest->cols)*sizeof(real));
}


void
swap_realmatrix(prealmatrix A, prealmatrix B)
{
    index tmp_dim;
    real *tmp_vals;

    assert(A);
    assert(B);

    tmp_dim = A->rows;
    A->rows = B->rows;
    B->rows = tmp_dim;

    tmp_dim = A->cols;
    A->cols = B->cols;
    B->cols = tmp_dim;

    tmp_vals = A->vals;
    A->vals  = B->vals;
    B->vals  = tmp_vals;
}

void
del_realmatrix(prealmatrix A)
{
    assert(A);

    if (A->vals) free(A->vals);
    A->rows = 0;
    A->cols = 0;
}


real
getentry_realmatrix(pcrealmatrix A, index i, index j)
{
    assert(A);
    assert(i < A->rows+INDEX_BASE);
    assert(j < A->cols+INDEX_BASE);

    return A->vals[ (j-INDEX_BASE) * A->rows + (i-INDEX_BASE)];
}


void
setentry_realmatrix(prealmatrix A, index i, index j, real entry)
{
    assert(A);
    assert(i < A->rows+INDEX_BASE);
    assert(j < A->cols+INDEX_BASE);

    A->vals[ (j-INDEX_BASE) * A->rows + (i-INDEX_BASE)] = entry;
}


void
addentry_realmatrix(prealmatrix A, index i, index j, real entry)
{
    assert(A);
    assert(i<A->rows+INDEX_BASE);
    assert(j<A->cols+INDEX_BASE);

    A->vals[ (j-INDEX_BASE) * A->cols + (i-INDEX_BASE)] += entry;
}


void
print_realmatrix(pcrealmatrix A)
{
    index i, j;

    assert(A);

    for ( i = INDEX_BASE; i < A->rows+INDEX_BASE; ++i) {
      for ( j = INDEX_BASE; j < A->cols+INDEX_BASE; ++j) {
        (void) printf("%8.4f ", getentry_realmatrix(A, i, j));
      }
      (void) printf("\n");
    }
}


/*
 * Load INDEXMATRIX with cols columns A
 * Return A   if transpose is 0
 *     or A^T if transpose is 1
 */
prealmatrix
load_realmatrix(char *fname, index rows, transpose t)
{
  FILE *file;
  index cnt, cols, i, j;
  real a;
  prealmatrix x;
  int suppr_out;

  file = fopen(fname,"r");

  if (file == NULL){
    printf("\n load_realmatrix() Error!!! (a)\n\n");
    return NULL;
  }

  cnt = 0;
  while (fscanf(file,"%lf",&a) != EOF) cnt++;
  if (  cnt - rows * (cnt / rows) ){
    fclose(file);
    printf("\n load_realmatrix() Error!!! (b), cnt = %lu\n\n", cnt);
    return NULL;
  }

  x = (prealmatrix) malloc(sizeof(realmatrix));

  cols = cnt/rows;
  fseek(file,0L,SEEK_SET);
  if (t==notrans){
   init_realmatrix(x, cols,rows);
   for (j=0; j<cols; j++) {
     for (i=0; i<rows; i++) {
        suppr_out = fscanf(file,"%lf",&(x->vals[i*cols+j]));
        (void) suppr_out;
      }
    }
  } else {
    init_realmatrix(x,rows,cols);
    for (j=0; j<rows*cols; j++) {
      suppr_out = fscanf(file,"%lf",&(x->vals[j]));
      (void) suppr_out;
    }
  }
  fclose(file);
  return x;
}


int
write_realmatrix(char *fname, prealmatrix A, transpose t)
{
  FILE *file;
  index i, j, m, n;
  real *Ax;

  file = fopen(fname,"w");

  if (file == NULL){
    printf("\n write_realmatrix() Error!!!\n\n");
    return 0;
  }

  Ax = A->vals; m = A->rows; n = A->cols;
  if (t==notrans){
    for ( i = 0; i < m; i++) {
      for ( j = 0; j < n; j++) {
        fprintf(file,"%22.15g ", Ax[i+j*m]);
      }
      fprintf(file,"\n");
    }
  } else {
    for ( j = 0; j < n; j++) {
      for ( i = 0; i < m; i++) {
        fprintf(file,"%22.15g ", Ax[i+j*m]);
      }
      fprintf(file,"\n");
    }
  }
  fclose(file);
  return 1;
}

#endif
