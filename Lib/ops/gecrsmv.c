#ifndef GECRSMV_C
#define GECRSMV_C

#include <stdio.h>
#include <assert.h>

#include "gecrsmv.h"

void
gecrsmv(transpose t,
        real alpha,
        pccrs A,
        pcrealvector x,
        real beta,
        prealvector y)
{
    index i, j;

    /* Input parsing */
    assert(A);
    assert(x);
    assert(y);

    if (t==notrans) {
        assert(A->numc==x->length);
        assert(A->numr==y->length);
    } else if (t==trans) {
        assert(A->numr==x->length);
        assert(A->numc==y->length);
    } else {
        fprintf(stderr, "Transpose type unknown.\n");
    }

    /* Scale y */
    if (beta!=(real) 1) scal_realvector(beta, y);

    /* Perform y <- alpha*A*x + y */
    if (t==notrans) {
        for (i=0; i<A->numr; ++i) {
            for(j=A->rowptr[i]; j<A->rowptr[i+1]; ++j) {
                    addentry_realvector(y,
                                        i+INDEX_BASE,
                                        alpha*
                                        getentry_realvector(x, A->colind[j])*
                                        A->vals[j]);
            }
        }
    } else {
        for (i=0; i<A->numr; ++i) {
            for(j=A->rowptr[i]; j<A->rowptr[i+1]; ++j) {
                    addentry_realvector(y,
                                        A->colind[j],
                                        alpha*
                                        getentry_realvector(x, i+INDEX_BASE)*
                                        A->vals[j]);

            }
        }
    }
}

#endif
