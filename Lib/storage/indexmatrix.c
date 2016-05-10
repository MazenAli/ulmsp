#ifndef INDEXMATRIX_C
#define INDEXMATRIX_C

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>

/* Conflicting typedefs for "index" */
#define index string_index
    #include <string.h>
#undef index

#include <math.h>

#include "indexmatrix.h"

pindexmatrix
new_indexmatrix(index rows, index cols)
{
    pindexmatrix A;
    A = (pindexmatrix) malloc(sizeof(indexmatrix));
    init_indexmatrix(A, rows, cols);

    return A;
}


void
init_indexmatrix(pindexmatrix A, index rows, index cols)
{
    index nz;

    assert(A);

    nz = rows * cols;
    if (nz>0) {
        A->vals = (index*) malloc(nz*sizeof(index));
    } else {
        A->vals = NULL;
    }
    A->rows = rows;
    A->cols = cols;
}


void
resize_indexmatrix(pindexmatrix A, index rows, index cols)
{
    assert(A);

    del_indexmatrix(A);
    init_indexmatrix(A, rows, cols);
}


void
copy_indexmatrix(pindexmatrix dest, pcindexmatrix src)
{
    assert(dest);
    assert(src);

    if ((dest->rows!=src->rows)||(dest->cols!=src->cols)) {
        del_indexmatrix(dest);
        init_indexmatrix(dest, src->rows, src->cols);
    }

    memcpy(dest->vals, src->vals, (dest->rows*dest->cols)*sizeof(index));
}

void
swap_indexmatrix(pindexmatrix A, pindexmatrix B)
{
    index tmp_dim, *tmp_vals;

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
del_indexmatrix(pindexmatrix A)
{
    assert(A);

    if (A->vals) free(A->vals);
    A->rows = 0;
    A->cols = 0;
}


index
getentry_indexmatrix(pcindexmatrix A, index i, index j)
{
    assert(A);
    assert(i < A->rows+INDEX_BASE);
    assert(j < A->cols+INDEX_BASE);

    return A->vals[ (j-INDEX_BASE) * A->rows + (i-INDEX_BASE)];
}


void
setentry_indexmatrix(pindexmatrix A, index i, index j, index entry)
{
    assert(A);
    assert(i < A->rows+INDEX_BASE);
    assert(j < A->cols+INDEX_BASE);

    A->vals[ (j-INDEX_BASE) * A->rows + (i-INDEX_BASE)] = entry;
}


void
addentry_indexmatrix(pindexmatrix A, index i, index j, index entry)
{
    assert(A);
    assert(i<A->rows+INDEX_BASE);
    assert(j<A->cols+INDEX_BASE);

    A->vals[ (j-INDEX_BASE) * A->rows + (i-INDEX_BASE)] += entry;
}


void
print_indexmatrix(pcindexmatrix A)
{
    index i, j;

    assert(A);

    for ( i = INDEX_BASE; i < A->rows+INDEX_BASE; ++i) {
      for ( j = INDEX_BASE; j < A->cols+INDEX_BASE; ++j) {
        printf("%4lu ", getentry_indexmatrix(A, i, j));
      }
      printf("\n");
    }
}

/*
 * Load INDEXMATRIX with cols columns A
 * Return A   if transpose is 0
 *     or A^T if transpose is 1
 */
pindexmatrix load_indexmatrix(char *fname, index rows, transpose t)
{
  FILE *file;
  index cnt, cols, a, i, j;
  pindexmatrix x;

  file = fopen(fname,"r");

  if (file == NULL){
    printf("\n load_indexmatrix() Error!!! (a)\n\n");
    return NULL;
  }

  cnt = 0;
  while (fscanf(file,"%lu",&a) != EOF) cnt++;
  if (  cnt - rows * (cnt / rows) ){
    fclose(file);
    printf("\n load_indexmatrix() Error!!! (b), cnt = %lu\n\n", cnt);
    return NULL;
  }

  x = (pindexmatrix) malloc(sizeof(indexmatrix));

  cols = cnt/rows;
  fseek(file,0L,SEEK_SET);
  if (t==notrans){
   init_indexmatrix(x, cols,rows);
   for (j=0; j<cols; j++) {
     for (i=0; i<rows; i++) {
        fscanf(file,"%lu",&(x->vals[i*cols+j]));
      }
    }
  } else {
    init_indexmatrix(x,rows,cols);
    for (j=0; j<rows*cols; j++) {
      fscanf(file,"%lu",&(x->vals[j]));
    }
  }
  fclose(file);
  return x;
}

int write_indexmatrix(char *fname, pindexmatrix A, transpose t)
{
  FILE *file;
  index i, j, m, n, *Ax;

  file = fopen(fname,"w");

  if (file == NULL){
    printf("\n write_indexmatrix() Error!!!\n\n");
    return 0;
  }

  Ax = A->vals; m = A->rows; n = A->cols;
  if (t==notrans){
    for ( i = 0; i < m; i++) {
      for ( j = 0; j < n; j++) {
        fprintf(file,"%4lu ", Ax[i+j*m]);
      }
      fprintf(file,"\n");
    }
  } else {
    for ( j = 0; j < n; j++) {
      for ( i = 0; i < m; i++) {
        fprintf(file,"%4lu ", Ax[i+j*m]);
      }
      fprintf(file,"\n");
    }
  }
  fclose(file);
  return 1;
}

#endif
