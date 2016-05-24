#ifndef GSCRS_C
#define GSCRS_C

#include <assert.h>

#include "gscrs.h"

void
gscrs(pccrs A,
      prealvector x,
      pcrealvector b,
      index it)
{
    index i, j, k;
    real sum;

    /* Input check */
    assert(A->numc==x->length);
    assert(A->numr==b->length);
    assert(x->length==b->length);

    for (k=0; k<it; ++k) {
        for (i=0; i<A->numr; ++i) {
            sum = (real) 0;
            for (j=A->rowptr[i]+1; j<A->rowptr[i+1]; ++j) {
                sum += A->vals[j]*x->vals[A->colind[j]-INDEX_BASE];
            }
            x->vals[i] = (b->vals[i]-sum)/A->vals[A->rowptr[i]];
        }
    }
}


void
gscrs_constrains(pccrs A,
                 prealvector x,
                 pcrealvector b,
                 pcindexvector fixedMask,
                 index it)
{
    index i, j, k;
    real sum;

    /* Input check */
    assert(A->numc==x->length);
    assert(A->numr==b->length);
    assert(x->length==b->length);
    assert(fixedMask->length==x->length);

    for (k=0; k<it; ++k) {
        for (i=0; i<A->numr; ++i) {
            if (!fixedMask->vals[i]) {
                sum = (real) 0;
                for (j=A->rowptr[i]+1; j<A->rowptr[i+1]; ++j) {
                    if (!fixedMask->vals[A->colind[j]-INDEX_BASE]) {
                        sum += A->vals[j]*x->vals[A->colind[j]-INDEX_BASE];
                    }
                }
                x->vals[i] = (b->vals[i]-sum)/A->vals[A->rowptr[i]];
            }
        }
    }
}

#endif
