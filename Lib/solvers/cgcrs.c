#ifndef CGCRS_C
#define CGCRS_C

#include <assert.h>
#include <stdio.h>
#include <math.h>

#include "../ops/gecrsmv.h"
#include "cgcrs.h"
#include "precs.h"

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
pcgdiagcrs(pccrs A,
           prealvector x,
           pcrealvector b,
           real tol,
           index maxit)
{
    index i;
    real rho, lambda, rkzk, rkzknew;
    prealvector r          = new_realvector(x->length);
    prealvector z          = new_realvector(x->length);
    prealvector p          = new_realvector(x->length);
    prealvector a          = new_realvector(x->length);
    transpose   t          = notrans;

    /* Initial residual, r = b-A*x */
    copy_realvector(r, b);
    gecrsmv(t, -1., A, x, 1., r);

    rho = dot_realvector(r, r);

    /* Initial z_0 */
    fill_realvector(z, 0.);
    diagcrs(A, z, r);

    copy_realvector(p, z);
    rkzknew = dot_realvector(r, z);

    for (i=0; i<=maxit; ++i) {

        if (rho<=tol*tol) {
            #ifdef VERBOSE
                printf("pcg:  r = %.16f\n", sqrt(rho));
            #endif

            del_realvector(r);
            del_realvector(p);
            del_realvector(a);

            return i;
        }

        /* alpha_k */
        gecrsmv(t, 1., A, p, 0., a);
        lambda = rkzknew/dot_realvector(p, a);

        /* Update x_k and r_k */
        axpy_realvector(lambda, p, x);
        axpy_realvector(-lambda, a, r);

        /* Update z_k */
        diagcrs(A, z, r);

        /* rho_new */
        rho = dot_realvector(r, r);

        /* rkzknew */
        rkzk    = rkzknew;
        rkzknew = dot_realvector(r, z);

        /* Update p_k */
        scal_realvector(rkzknew/rkzk, p);
        axpy_realvector(1., z, p);
    }

    fprintf(stderr,
            "Max iterations reached: maxit = %ld, r = %.16f\n", maxit,
                                                             sqrt(rho));

    del_realvector(r);
    del_realvector(z);
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
    for ( k =  0 ; k < fixedNodes->length; ++k){
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
        for ( k =  0 ; k < fixedNodes->length; ++k){
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


index
pcgdiagcrs_constrains(pccrs A,
                      prealvector x,
                      pcrealvector b,
                      pcindexvector fixedNodes,
                      real tol,
                      index maxit)
{
    index i, k;
    real rho, lambda, rkzk, rkzknew;
    pindexvector fixedMask = new_indexvector(x->length);
    prealvector r          = new_realvector(x->length);
    prealvector z          = new_realvector(x->length);
    prealvector p          = new_realvector(x->length);
    prealvector a          = new_realvector(x->length);
    transpose   t          = notrans;

    /* Set fixedMask for Jacobi */
    for (k=0; k<fixedMask->length; ++k) {
        fixedMask->vals[k] = 0;
        for (i=0; i<fixedNodes->length; ++i) {
            if (k==fixedNodes->vals[i]-INDEX_BASE) ++fixedMask->vals[k];
        }
    }

    /* Initial residual, r = b-A*x */
    copy_realvector(r, b);
    gecrsmv(t, -1., A, x, 1., r);

    /* Incorporate constrains */
    for ( k =  0 ; k < fixedNodes->length; ++k){
        r->vals[fixedNodes->vals[k]-INDEX_BASE] = 0.;
    }

    rho = dot_realvector(r, r);

    /* Initial z_0 */
    diagcrs_constrains(A, z, r, fixedMask);

    copy_realvector(p, z);

    rkzknew = dot_realvector(r, z);

    for (i=0; i<=maxit; ++i) {

        if (rho<=tol*tol) {
            #ifdef VERBOSE
                printf("pcg:  r = %.16f\n", sqrt(rho));
            #endif

            del_realvector(r);
            del_realvector(p);
            del_realvector(a);

            return i;
        }

        /* alpha_k */
        gecrsmv(t, 1., A, p, 0., a);

        /* Incorporate constrains */
        for ( k =  0 ; k < fixedNodes->length; ++k){
            a->vals[fixedNodes->vals[k]-INDEX_BASE] = 0.;
        }
        lambda = rkzknew/dot_realvector(p, a);

        /* Update x_k and r_k */
        axpy_realvector(lambda, p, x);
        axpy_realvector(-lambda, a, r);

        /* Update z_k */
        diagcrs_constrains(A, z, r, fixedMask);

        /* rho_new */
        rho = dot_realvector(r, r);

        /* rkzknew */
        rkzk    = rkzknew;
        rkzknew = dot_realvector(r, z);

        /* Update p_k */
        scal_realvector(rkzknew/rkzk, p);
        axpy_realvector(1., z, p);
    }

    fprintf(stderr,
            "Max iterations reached: maxit = %ld, r = %.16f\n", maxit,
                                                             sqrt(rho));

    del_realvector(r);
    del_realvector(z);
    del_realvector(p);
    del_realvector(a);

    return maxit;
}


index
pcgsymgscrs_constrains(pccrs A,
                       prealvector x,
                       pcrealvector b,
                       pcindexvector fixedNodes,
                       real tol,
                       index maxit)
{
    index i, k;
    real rho, lambda, rkzk, rkzknew;
    pindexvector fixedMask = new_indexvector(x->length);
    prealvector r          = new_realvector(x->length);
    prealvector z          = new_realvector(x->length);
    prealvector p          = new_realvector(x->length);
    prealvector a          = new_realvector(x->length);
    transpose   t          = notrans;

    /* Set fixedMask for Jacobi */
    for (k=0; k<fixedMask->length; ++k) {
        fixedMask->vals[k] = 0;
        for (i=0; i<fixedNodes->length; ++i) {
            if (k==fixedNodes->vals[i]-INDEX_BASE) ++fixedMask->vals[k];
        }
    }

    /* Initial residual, r = b-A*x */
    copy_realvector(r, b);
    gecrsmv(t, -1., A, x, 1., r);

    /* Incorporate constrains */
    for ( k =  0 ; k < fixedNodes->length; ++k){
        r->vals[fixedNodes->vals[k]-INDEX_BASE] = 0.;
    }

    rho = dot_realvector(r, r);

    /* Initial z_0 */
    symgscrs_constrains(A, z, r, fixedMask);

    copy_realvector(p, z);

    rkzknew = dot_realvector(r, z);

    for (i=0; i<=maxit; ++i) {

        if (rho<=tol*tol) {
            #ifdef VERBOSE
                printf("pcg:  r = %.16f\n", sqrt(rho));
            #endif

            del_realvector(r);
            del_realvector(p);
            del_realvector(a);

            return i;
        }

        /* alpha_k */
        gecrsmv(t, 1., A, p, 0., a);

        /* Incorporate constrains */
        for ( k =  0 ; k < fixedNodes->length; ++k){
            a->vals[fixedNodes->vals[k]-INDEX_BASE] = 0.;
        }
        lambda = rkzknew/dot_realvector(p, a);

        /* Update x_k and r_k */
        axpy_realvector(lambda, p, x);
        axpy_realvector(-lambda, a, r);

        /* Update z_k */
        symgscrs_constrains(A, z, r, fixedMask);

        /* rho_new */
        rho = dot_realvector(r, r);

        /* rkzknew */
        rkzk    = rkzknew;
        rkzknew = dot_realvector(r, z);

        /* Update p_k */
        scal_realvector(rkzknew/rkzk, p);
        axpy_realvector(1., z, p);
    }

    fprintf(stderr,
            "Max iterations reached: maxit = %ld, r = %.16f\n", maxit,
                                                             sqrt(rho));

    del_realvector(r);
    del_realvector(z);
    del_realvector(p);
    del_realvector(a);

    return maxit;
}

#endif
