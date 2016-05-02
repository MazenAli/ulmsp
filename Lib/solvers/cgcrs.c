#ifndef CGCRS_C
#define CGCRS_C

#include <assert.h>
#include <stdio.h>

#include "../ops/gecrsmv.h"
#include "cgcrs.h"

index
cgcrs(pccrs A,
      prealvector x,
      pcrealvector b,
      real tol,
      index maxit)
{
    index i;
    real residual;
    prealvector r, p, Ap;
    transpose t;

    /* Input parsing */
    assert(A);
    assert(x);
    assert(b);
    assert(A->numc==x->length);
    assert(x->length==b->length);
    assert(A->numr==A->numc);

    r  = new_realvector(x->length);
    p  = new_realvector(x->length);
    Ap = new_realvector(x->length);

    /* Initial residual */
    copy_realvector(r, b);
    t = notrans;
    gecrsmv(t, -1., A, x, 1., r);
    copy_realvector(p, r);

    residual = nrm2_realvector(r);

    for (i=0; i<=maxit; ++i) {
        real ak, bk, residual_old;

        if (residual<=tol) {
            #ifdef VERBOSE
                (void) printf("cg: r = %f\n", residual);
            #endif

            del_realvector(r);
            del_realvector(p);
            del_realvector(Ap);

            return i;
        }

        /* alpha_k */
        ak = residual*residual;
        gecrsmv(t, 1., A, p, 0., Ap);
        ak /= dot_realvector(p, Ap);

        /* Update x_k and r_k */
        axpy_realvector(ak, p, x);
        axpy_realvector(-ak, Ap, r);

        /* Update p_k */
        residual_old = residual;
        residual     = nrm2_realvector(r);
        bk           = residual*residual/(residual_old*residual_old);
        scal_realvector(bk, p);
        axpy_realvector(1., r, p);
    }

    fprintf(stderr,
            "Max iterations reached: maxit = %ld, r = %f\n", maxit, residual);

    del_realvector(r);
    del_realvector(p);
    del_realvector(Ap);

    return maxit;
}

#endif
