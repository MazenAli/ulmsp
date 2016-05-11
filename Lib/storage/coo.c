#ifndef COO_C
#define COO_C

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

#include "coo.h"

pcoo
new_coo(index nonz, index numr, index numc)
{
    pcoo A;

    A = (pcoo) malloc(sizeof(coo));
    init_coo(A, nonz, numr, numc);

    return A;
}


void
init_coo(pcoo A, index nonz, index numr, index numc)
{
    assert(A);

    if (nonz>0) {
        assert(numr>0);
        assert(numc>0);
        A->vals = (real*) malloc(nonz*sizeof(real));
        A->rows = (index*) malloc(nonz*sizeof(index));
        A->cols = (index*) malloc(nonz*sizeof(index));
    } else {
        assert(!numr);
        assert(!numc);
        A->vals = NULL;
        A->rows = NULL;
        A->cols = NULL;
    }
    A->nonz = nonz;
    A->numr = numr;
    A->numc = numc;
}


void
resize_coo(pcoo A, index nonz, index numr, index numc)
{
    assert(A);

    if (nonz>0) {
      assert(numr>0);
      assert(numc>0);
      A->vals = (real*)  realloc(A->vals,nonz*sizeof(real));
      A->rows = (index*) realloc(A->rows,nonz*sizeof(index));
      A->cols = (index*) realloc(A->cols,nonz*sizeof(index));
      A->nonz = nonz;
      A->numr = numr;
      A->numc = numc;
    } else {
        del_coo(A);
        init_coo(A, nonz, numr, numc);
    }
}


void
del_coo(pcoo A)
{
    assert(A);

    if (A->vals) free(A->vals);
    if (A->rows) free(A->rows);
    if (A->cols) free(A->cols);
    A->nonz = 0;
    A->numr = 0;
    A->numc = 0;
}


real
getentry_coo(pccoo A, index i, index j)
{
    index    ii;
    real entry;

    check_base(i);
    check_base(j);
    assert(i<A->numr+INDEX_BASE);
    assert(j<A->numc+INDEX_BASE);

    entry = 0.;

    for (ii=0; ii<A->nonz; ++ii) {
        if (A->rows[ii]==i && A->cols[ii]==j) {
            entry = A->vals[ii];
        }
    }

    return entry;
}


void
setentry_coo(pcoo A, index i, index j, real entry)
{
    index    ii;

    check_base(i);
    check_base(j);
    assert(i<A->numr+INDEX_BASE);
    assert(j<A->numc+INDEX_BASE);

    for (ii=0; ii<A->nonz; ++ii) {
        if (A->rows[ii]==i && A->cols[ii]==j) {
            A->vals[ii] = entry;
            return;
        }
    }

    fprintf(stderr, "Index not in range of sparse matrix\n");
}


void
addentry_coo(pcoo A, index i, index j, real entry)
{
    index    ii;

    check_base(i);
    check_base(j);
    assert(i<A->numr+INDEX_BASE);
    assert(j<A->numc+INDEX_BASE);

    for (ii=0; ii<A->nonz; ++ii) {
        if (A->rows[ii]==i && A->cols[ii]==j) {
            A->vals[ii] += entry;
            return;
        }
    }

    fprintf(stderr, "Index not in range of sparse matrix\n");
}


void
print_coo(pccoo A)
{
    index i;

    assert(A);

    for (i=0; i<A->nonz; ++i) {
        (void) printf("(%lu, %lu) : %8.4f\n",
                      A->rows[i], A->cols[i], A->vals[i]);
    }
}


void
printdense_coo(pccoo A)
{
    index i, j;

    assert(A);

    for (i=INDEX_BASE; i<A->numr+INDEX_BASE; ++i) {
        for (j=INDEX_BASE; j<A->numc+INDEX_BASE; ++j) {
            (void) printf("%8.4f", getentry_coo(A, i, j));
        }
        (void) printf("\n");
    }
}

#endif
