#ifndef CGCRS_C
#define CGCRS_C

#include <assert.h>
#include <stdio.h>
#include <math.h>

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
    real rho_old, rho_new, lambda;
    prealvector r = new_realvector(x->length);
    prealvector p = new_realvector(x->length);
    prealvector a = new_realvector(x->length);
    transpose t = notrans;

    /* Initial residual, r = b-A*x */
    copy_realvector(r, b);
    gecrsmv(t, -1., A, x, 1., r);
    copy_realvector(p, r);

    rho_new = dot_realvector(r,r);

    for (i=0; i<=maxit; ++i) {

        if (rho_new<=tol*tol) {
            #ifdef VERBOSE
                printf("cg: r = %.16f\n", sqrt(rho_new));
            #endif

            del_realvector(r);
            del_realvector(p);
            del_realvector(a);

            return i;
        }

        /* alpha_k */
        gecrsmv(t, 1., A, p, 0., a);
        lambda = rho_new/dot_realvector(p, a);

        /* Update x_k and r_k */
        axpy_realvector(lambda, p, x);
        axpy_realvector(-lambda, a, r);

        /* rho_new */
        rho_old = rho_new;
        rho_new = dot_realvector(r,r);

        /* Update p_k */
        scal_realvector(rho_new/rho_old, p);
        axpy_realvector(1., r, p);
    }

    fprintf(stderr,
            "Max iterations reached: maxit = %ld, r = %.16f\n", maxit,
                                                             sqrt(rho_new));

    del_realvector(r);
    del_realvector(p);
    del_realvector(a);

    return maxit;
}


index
cgcrs_constrains(pccrs A,                   /* in     */
                 prealvector x,             /* in/out */
                 pcrealvector b,            /* in     */
                 pcindexvector fixedNodes,  /* in     */
                 real tol,                  /* in     */
                 index maxit)               /* in     */
{
    index i, k;
    real rho_old, rho_new, lambda;
    prealvector r = new_realvector(x->length);
    prealvector p = new_realvector(x->length);
    prealvector a = new_realvector(x->length);
    transpose   t = notrans;

    /* Initial residual, r = b-A*x */
    copy_realvector(r, b);
    gecrsmv(t, -1., A, x, 1., r);

    /* Incorporate constrains */
    for ( k =  0 ; k < fixedNodes->length; k++){
        r->vals[fixedNodes->vals[k]-INDEX_BASE] = 0.;
    }
    copy_realvector(p, r);

    rho_new = dot_realvector(r,r);

    for (i=0; i<=maxit; ++i) {

        if (rho_new<=tol*tol) {
            #ifdef VERBOSE
                printf("cg:  r = %.16f\n", sqrt(rho_new));
            #endif

            del_realvector(r);
            del_realvector(p);
            del_realvector(a);

            return i;
        }

        /* alpha_k */
        gecrsmv(t, 1., A, p, 0., a);

        /* Incorporate constrains */
        for ( k =  0 ; k < fixedNodes->length; k++){
            a->vals[fixedNodes->vals[k]-INDEX_BASE] = 0.;
        }
        lambda = rho_new/dot_realvector(p, a);

        /* Update x_k and r_k */
        axpy_realvector(lambda, p, x);
        axpy_realvector(-lambda, a, r);

        /* rho_new */
        rho_old = rho_new;
        rho_new = dot_realvector(r,r);

        /* Update p_k */
        scal_realvector(rho_new/rho_old, p);
        axpy_realvector(1., r, p);
    }

    fprintf(stderr,
            "Max iterations reached: maxit = %ld, r = %.16f\n", maxit,
                                                             sqrt(rho_new));

    del_realvector(r);
    del_realvector(p);
    del_realvector(a);

    return maxit;
}

#endif
