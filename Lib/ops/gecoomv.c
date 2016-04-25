#ifndef GECOOMV_C
#define GECOOMV_C

#include <stdio.h>
#include <assert.h>

#include "gecoomv.h"

void
gecoomv(transpose t,
        real alpha,
        pccoo A,
        pcrealvector x,
        real beta,
        prealvector y)
{
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

    index i;

    for (i=INDEX_BASE; i<y->length+INDEX_BASE; ++i) {
        setentry_realvector(y, i,
                            beta*getentry_realvector(y, i));
    }

    for (i=0; i<A->nonz; ++i) {
        if (t==notrans) {
            addentry_realvector(y,
                                A->rows[i],
                                alpha*
                                getentry_realvector(x, A->cols[i])*
                                A->vals[i]);
        } else {
            addentry_realvector(y,
                                A->cols[i],
                                alpha*
                                getentry_realvector(x, A->rows[i])*
                                A->vals[i]);

        }
    }
}

#endif
