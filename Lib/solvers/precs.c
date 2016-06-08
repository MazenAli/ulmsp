#ifndef PRECS_C
#define PRECS_C

#include <assert.h>

#include "precs.h"

void
diagcrs(pccrs        A,
        prealvector  x,
        pcrealvector b)
{
    index i;

    /* Input check */
    assert(A->numc==x->length);
    assert(A->numr==b->length);
    assert(x->length==b->length);

    for (i=0; i<x->length; ++i) {
        x->vals[i] = b->vals[i]/A->vals[A->rowptr[i]];
    }
}


void
diagcrs_constrains(pccrs         A,
                   prealvector   x,
                   pcrealvector  b,
                   pcindexvector fixedMask)
{
    index i;

    /* Input check */
    assert(A->numc==x->length);
    assert(A->numr==b->length);
    assert(x->length==b->length);
    assert(fixedMask->length==x->length);

    for (i=0; i<x->length; ++i) {
        if (!fixedMask->vals[i]) {
            x->vals[i] = b->vals[i]/A->vals[A->rowptr[i]];
        }
    }
}


void
frwgscrs_constrains(pccrs         A,
                    prealvector   x,
                    pcrealvector  b,
                    pcindexvector fixedMask)
{
    index i, j;
    real sum;

    /* Input check */
    assert(A->numc==x->length);
    assert(A->numr==b->length);
    assert(x->length==b->length);
    assert(fixedMask->length==x->length);

    for (i=0; i<A->numr; ++i) {
        if (!fixedMask->vals[i]) {
            sum = (real) 0;
            for (j=A->rowptr[i]+1; j<A->rowptr[i+1]; ++j) {
                if (!fixedMask->vals[A->colind[j]-INDEX_BASE]
                    && i>A->colind[j]-INDEX_BASE) {
                    sum += A->vals[j]*x->vals[A->colind[j]-INDEX_BASE];
                }
            }
            x->vals[i] = (b->vals[i]-sum)/A->vals[A->rowptr[i]];
        }
    }
}


void
bkwgscrs_constrains(pccrs         A,
                    prealvector   x,
                    pcrealvector  b,
                    pcindexvector fixedMask)
{
    index j;
    long  i;
    real sum;

    /* Input check */
    assert(A->numc==x->length);
    assert(A->numr==b->length);
    assert(x->length==b->length);
    assert(fixedMask->length==x->length);

    for (i=A->numr-1; i>=0; --i) {
        if (!fixedMask->vals[i]) {
            sum = (real) 0;
            for (j=A->rowptr[i]+1; j<A->rowptr[i+1]; ++j) {
                if (!fixedMask->vals[A->colind[j]-INDEX_BASE]
                    && (unsigned) i<A->colind[j]-INDEX_BASE) {
                    sum += A->vals[j]*x->vals[A->colind[j]-INDEX_BASE];
                }
            }
            x->vals[i] = (b->vals[i]-sum)/A->vals[A->rowptr[i]];
        }
    }
}


void
symgscrs_constrains(pccrs         A,
                    prealvector   x,
                    pcrealvector  b,
                    pcindexvector fixedMask)
{
    index i;
    prealvector tmp = new_realvector(x->length);

    /* Input check */
    assert(A->numc==x->length);
    assert(A->numr==b->length);
    assert(x->length==b->length);
    assert(fixedMask->length==x->length);

    fill_realvector(tmp, 0.);

    /* Forward sub */
    frwgscrs_constrains(A, tmp, b, fixedMask);

    /* Diagonal scaling */
    for (i=0; i<tmp->length; ++i) {
        tmp->vals[i] *= A->vals[A->rowptr[i]];
    }

    /* Backward sub */
    bkwgscrs_constrains(A, x, tmp, fixedMask);

    del_realvector(tmp);
}

#endif
