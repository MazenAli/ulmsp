#ifndef CGCOO_C
#define CGCOO_C

#include <assert.h>
#include <stdio.h>

#include "../ops/gecoomv.h"
#include "cgcoo.h"

index
cgcoo(pccoo A,
      prealvector x,
      pcrealvector b,
      real tol,
      index maxit)
{
    index i;
    real residual;
    prealvector r, p, temp, rold;
    transpose t;

    /* Input parsing */
    assert(A);
    assert(x);
    assert(b);
    assert(A->numc==x->length);
    assert(x->length==b->length);
    assert(A->numr==A->numc);

    r    = new_realvector(x->length);
    p    = new_realvector(x->length);
    temp = new_realvector(x->length);
    rold = new_realvector(x->length);

    /* Initial residual */
    copy_realvector(r, b);
    t = notrans;
    gecoomv(t, -1., A, x, 1., r);
    copy_realvector(p, r);

    for (i=0; i<=maxit; ++i) {
        real ak, bk;

        residual = nrm2_realvector(r);

        if (residual<=tol) {
            #ifdef VERBOSE
                (void) printf("cg: r = %f\n", residual);
            #endif

            del_realvector(r);
            del_realvector(p);
            del_realvector(temp);
            del_realvector(rold);

            return i;
        }

        /* alpha_k */
        ak = dot_realvector(r, r);
        gecoomv(t, 1., A, p, 0., temp);
        ak /= dot_realvector(p, temp);

        /* Update x_k and r_k */
        copy_realvector(rold, r);
        axpy_realvector(ak, p, x);
        axpy_realvector(-ak, temp, r);

        /* Update p_k */
        bk = dot_realvector(r, r)/dot_realvector(rold, rold);
        scal_realvector(bk, p);
        axpy_realvector(1., r, p);
    }

    fprintf(stderr,
            "Max iterations reached: maxit = %ld, r = %f\n", maxit, residual);

    del_realvector(r);
    del_realvector(p);
    del_realvector(temp);
    del_realvector(rold);

    return maxit;
}

#endif
