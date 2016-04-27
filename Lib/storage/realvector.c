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

    index i;

    if (dest->length!=src->length) {
        del_realvector(dest);
        init_realvector(dest, src->length);
    }

    for (i=INDEX_BASE; i<src->length+INDEX_BASE; ++i) {
        setentry_realvector(dest, i, getentry_realvector(src, i));
    }
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
    assert(i<x->length+INDEX_BASE);

    return x->vals[i-INDEX_BASE];
}


void
setentry_realvector(prealvector x, index i, real entry)
{
    assert(x);
    assert(i<x->length+INDEX_BASE);

    x->vals[i-INDEX_BASE] = entry;
}


void
addentry_realvector(prealvector x, index i, real entry)
{
    assert(x);
    assert(i<x->length+INDEX_BASE);

    x->vals[i-INDEX_BASE] += entry;
}


void
print_realvector(pcrealvector x)
{
    assert(x);

    index i;

    for (i=INDEX_BASE; i<x->length+INDEX_BASE; ++i) {
        (void) printf("%8.4f\n", getentry_realvector(x, i));
    }
}


void
scal_realvector(real alpha, prealvector x)
{
    assert(x);

    #ifdef USE_BLAS
        int N   = x->length;
        int one = 1;
        dscal_(&N, &alpha, x->vals, &one);

    #elif defined(USE_CBLAS)
        int N   = x->length;
        cblas_dscal(N, alpha, x->vals, 1);

    #else
        index i;
        for (i=INDEX_BASE; i<x->length+INDEX_BASE; ++i) {
            setentry_realvector(x, i,
                                alpha*getentry_realvector(x, i));
        }
    #endif
}


real
dot_realvector(pcrealvector x, pcrealvector y)
{
    assert(x);
    assert(y);
    assert(x->length==y->length);

    #ifdef USE_BLAS
        int N   = x->length;
        int one = 1;

        return ddot_(&N, x->vals, &one, y->vals, &one);

    #elif defined(USE_CBLAS)
        int N = x->length;

        return cblas_ddot(N, x->vals, 1, y->vals, 1);

    #else
        index i;
        real ret = 0.;

        for (i=INDEX_BASE; i<x->length+INDEX_BASE; ++i) {
            ret += getentry_realvector(x, i)*getentry_realvector(y, i);
        }

        return ret;
    #endif
}


real
nrm2_realvector(pcrealvector x)
{
    assert(x);

    #ifdef USE_BLAS
        int N   = x->length;
        int one = 1;

        return dnrm2_(&N, x->vals, &one);

    #elif defined(USE_CBLAS)
        int N = x->length;

        return cblas_dnrm2(N, x->vals, 1);

    #else
        index i;
        real ret = 0.;
        for (i=INDEX_BASE; i<x->length+INDEX_BASE; ++i) {
            ret += getentry_realvector(x, i)*getentry_realvector(x, i);
        }

        return sqrt(ret);
    #endif
}


void
axpy_realvector(real alpha, pcrealvector x, prealvector y)
{
    assert(x);
    assert(y);
    assert(x->length==y->length);

    #ifdef USE_BLAS
        int N   = x->length;
        int one = 1;
        daxpy_(&N, &alpha, x->vals, &one, y->vals, &one);

    #elif defined(USE_CBLAS)
        int N = x->length;
        cblas_daxpy(N, alpha, x->vals, 1, y->vals, 1);

    #else
        index i;

        for (i=INDEX_BASE; i<y->length+INDEX_BASE; ++i) {
            setentry_realvector(y, i,
                                alpha*getentry_realvector(x, i)+
                                getentry_realvector(y, i));
        }
    #endif
}

#endif
