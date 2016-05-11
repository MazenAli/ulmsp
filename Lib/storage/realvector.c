#ifndef REALVECTOR_C
#define REALVECTOR_C

#ifdef USE_BLAS
    #include <blas.h>
#endif

#ifdef USE_CBLAS
    #include <cblas.h>
#endif

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>

/* Conflicting typedefs for "index" */
#define index string_index
    #include <string.h>
#undef index

#include <math.h>

#include "realvector.h"

prealvector
new_realvector(index length)
{
    prealvector x;
    x = (prealvector) malloc(sizeof(realvector));
    init_realvector(x, length);

    return x;
}


void
init_realvector(prealvector x, index length)
{
    assert(x);

    if (length>0) {
        x->vals = (real*) malloc(length*sizeof(real));
    } else {
        x->vals = NULL;
    }
    x->length = length;
}


void
resize_realvector(prealvector x, index length)
{
    assert(x);

    del_realvector(x);
    init_realvector(x, length);
}


void
copy_realvector(prealvector dest, pcrealvector src)
{
    assert(dest);
    assert(src);

    if (dest->length!=src->length) {
        del_realvector(dest);
        init_realvector(dest, src->length);
    }

    memcpy(dest->vals, src->vals, dest->length*sizeof(real));
}


void
swap_realvector(prealvector x, prealvector y)
{
    index tmp_length;
    real *tmp_vals;

    assert(x);
    assert(y);

    tmp_length = x->length;
    x->length  = y->length;
    y->length  = tmp_length;

    tmp_vals = x->vals;
    x->vals  = y->vals;
    y->vals  = tmp_vals;
}


void
del_realvector(prealvector x)
{
    assert(x);

    if (x->vals) free(x->vals);
    x->length = 0;
}


real
getentry_realvector(pcrealvector x, index i)
{
    assert(x);
    check_base(i);
    assert(i<x->length+INDEX_BASE);

    return x->vals[i-INDEX_BASE];
}


void
setentry_realvector(prealvector x, index i, real entry)
{
    assert(x);
    check_base(i);
    assert(i<x->length+INDEX_BASE);

    x->vals[i-INDEX_BASE] = entry;
}


void
addentry_realvector(prealvector x, index i, real entry)
{
    assert(x);
    check_base(i);
    assert(i<x->length+INDEX_BASE);

    x->vals[i-INDEX_BASE] += entry;
}

void
fill_realvector(prealvector x, real val)
{
    index i;

    assert(x);

    for (i=INDEX_BASE; i<x->length+INDEX_BASE; ++i) {
        setentry_realvector(x, i, val);
    }
}


void
print_realvector(pcrealvector x)
{
    index i;

    assert(x);

    for (i=INDEX_BASE; i<x->length+INDEX_BASE; ++i) {
        (void) printf("%8.4f\n", getentry_realvector(x, i));
    }
}


prealvector
load_realvector(char *fname)
{
  prealvector x;
  FILE *file;
  index cnt;
  real a;

  file = fopen(fname,"r");

  if (file == NULL){
    printf("\n fopen() Error!!!\n\n");
    return NULL;
  }

  cnt = 0;
  while (fscanf(file,"%lf",&a) != EOF) cnt++;

  x = (prealvector) malloc(sizeof(realvector));
  init_realvector(x, cnt);

  fseek(file,0L,SEEK_SET);
  cnt = 0;
  while (fscanf(file,"%lf",&(x->vals[cnt])) != EOF) cnt++;
  fclose(file);
  return x;
}


int
write_realvector(char *fname, prealvector v)
{
  FILE *file;
  index i;
  real *vx;

  file = fopen(fname,"w");

  if (file == NULL){
    printf("\n write_realvector() Error!!!\n\n");
    return 0;
  }

  vx = v->vals; 
  for ( i = 0 ; i < v->length ; i++ ) {
    fprintf(file,"%22.15f\n", vx[i]);
  }
  fclose(file);
  return 1;
}


void
scal_realvector(real alpha, prealvector x)
{
#ifdef USE_BLAS
    int N   = x->length;
    int one = 1;

    assert(x);

    dscal_(&N, &alpha, x->vals, &one);

#elif defined(USE_CBLAS)
    int N   = x->length;

    assert(x);

    cblas_dscal(N, alpha, x->vals, 1);

#else
    index i;

    assert(x);

    if (alpha==(real) 0) {
        memset(x->vals, 0, x->length*sizeof(real));
    } else {
        for (i=INDEX_BASE; i<x->length+INDEX_BASE; ++i) {
            setentry_realvector(x, i,
                                alpha*getentry_realvector(x, i));
        }
    }
#endif
}


real
dot_realvector(pcrealvector x, pcrealvector y)
{
#ifdef USE_BLAS
    int N   = x->length;
    int one = 1;

    assert(x);
    assert(y);
    assert(x->length==y->length);

    return ddot_(&N, x->vals, &one, y->vals, &one);

#elif defined(USE_CBLAS)
    int N = x->length;

    assert(x);
    assert(y);
    assert(x->length==y->length);


    return cblas_ddot(N, x->vals, 1, y->vals, 1);

#else
    index i;
    real ret = 0.;

    assert(x);
    assert(y);
    assert(x->length==y->length);

    for (i=INDEX_BASE; i<x->length+INDEX_BASE; ++i) {
        ret += getentry_realvector(x, i)*getentry_realvector(y, i);
    }

    return ret;
#endif
}


real
nrm2_realvector(pcrealvector x)
{
#ifdef USE_BLAS
    int N   = x->length;
    int one = 1;

    assert(x);

    return dnrm2_(&N, x->vals, &one);

#elif defined(USE_CBLAS)
    int N = x->length;

    assert(x);

    return cblas_dnrm2(N, x->vals, 1);

#else
    index i;
    real ret;

    assert(x);

    ret = 0.;
    for (i=INDEX_BASE; i<x->length+INDEX_BASE; ++i) {
        ret += getentry_realvector(x, i)*getentry_realvector(x, i);
    }

    return sqrt(ret);
#endif
}


void
axpy_realvector(real alpha, pcrealvector x, prealvector y)
{
#ifdef USE_BLAS
    int N   = x->length;
    int one = 1;

    assert(x);
    assert(y);
    assert(x->length==y->length);

    daxpy_(&N, &alpha, x->vals, &one, y->vals, &one);

#elif defined(USE_CBLAS)
    int N = x->length;

    assert(x);
    assert(y);
    assert(x->length==y->length);

    cblas_daxpy(N, alpha, x->vals, 1, y->vals, 1);

#else
    index i;

    assert(x);
    assert(y);
    assert(x->length==y->length);

    if (alpha!=(real) 0) {
        for (i=INDEX_BASE; i<y->length+INDEX_BASE; ++i) {
            setentry_realvector(y, i,
                                alpha*getentry_realvector(x, i)+
                                getentry_realvector(y, i));
        }
    }
#endif
}

#endif
