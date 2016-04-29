/**
 *  Assemble and solve a simple laplace problem with the CRS storage format.
 *  The laplace problem is posed on the unit square and discretized uniformly
 *  in each dimension with linear FEM.
 */

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#include "settings.h"
#include "crs.h"
#include "gecrsmv.h"
#include "cgcrs.h"

void
assemble_laplace(pcrs A, index N);


int
main(int argc, char **argv)
{
    if (argc!=4) {
        (void) fprintf(stderr, "Usage: ./main N maxit print(1/0)\n");
        return 1;
    }

    int N     = atoi(argv[1]); /* Nodes in each dimension */
    int maxit = atoi(argv[2]); /* Maximum iteration number for CG */
    int print = atoi(argv[3]); /* Print output: set to 0 for big N */

    index i;
    real alpha = 1., beta = 0.; /* gecrsmv parameters */
    real tol = 1e-8;            /* CG tolerance */
    index it;

    /* Assemble and print stiffness matrix */
    pcrs A = new_crs(0, 0, 0);
    assemble_laplace(A, N);

    if (print) {
        (void) printf("A =\n");
        printdense_crs(A);
    }

    /* Set up, perform and print y <- alpha*A*x+beta*y */
    prealvector x = new_realvector(N*N);
    prealvector y = new_realvector(N*N);

    transpose t = notrans;

    for (i=INDEX_BASE; i<x->length+INDEX_BASE; ++i) {
        setentry_realvector(x, i, 1.);
    }

    gecrsmv(t, alpha, A, x, beta, y);

    if (print) {
        (void) printf("y = A*ones =\n");
        print_realvector(y);
    }

    /* Apply CG solver */
    it = cgcrs(A, y, x, tol, maxit);

    (void) printf("Solution after %ld iterations\n", it);
    if (print) {
        (void) printf("x =\n");
        print_realvector(y);
    }

    del_crs(A);
    del_realvector(x);
    del_realvector(y);

    return 0;
}


void
assemble_laplace(pcrs A, index N)
{
    assert(N>1);
    resize_crs(A, 5*N*N-4*N, N*N, N*N);

    index i, j;

    /* Lower left node */
    A->vals[0] = 4.;
    A->vals[1] = -1.;
    A->vals[2] = -1.;

    A->colind[0] = INDEX_BASE;
    A->colind[1] = INDEX_BASE+1;
    A->colind[2] = INDEX_BASE+N;

    A->rowptr[0] = 0;

    /* Lower right node */
    A->vals[4*N-5] = -1.;
    A->vals[4*N-4] = 4.;
    A->vals[4*N-3] = -1.;

    A->colind[4*N-5] = N-2+INDEX_BASE;
    A->colind[4*N-4] = N-1+INDEX_BASE;
    A->colind[4*N-3] = N-1+INDEX_BASE+N;

    A->rowptr[N-1] = 4*N-5;

    /* Upper left node */
    A->vals[5*N*N-8*N+2] = -1.;
    A->vals[5*N*N-8*N+3] = 4.;
    A->vals[5*N*N-8*N+4] = -1.;

    A->colind[5*N*N-8*N+2] = N*N-2*N+INDEX_BASE;
    A->colind[5*N*N-8*N+3] = N*N-N+INDEX_BASE;
    A->colind[5*N*N-8*N+4] = N*N-N+1+INDEX_BASE;

    A->rowptr[N*N-N]   = 5*N*N-8*N+2;

    /* Upper right node */
    A->vals[5*N*N-4*N-3] = -1.;
    A->vals[5*N*N-4*N-2] = -1.;
    A->vals[5*N*N-4*N-1] = 4.;

    A->colind[5*N*N-4*N-3] = N*N-1-N+INDEX_BASE;
    A->colind[5*N*N-4*N-2] = N*N-2+INDEX_BASE;
    A->colind[5*N*N-4*N-1] = N*N-1+INDEX_BASE;

    A->rowptr[N*N-1] = 5*N*N-4*N-3;
    A->rowptr[N*N]   = 5*N*N-4*N;

    /* Lower most inner nodes */
    for (i=0; i<N-2; ++i) {
        A->vals[3+i*4] = -1.;
        A->vals[4+i*4] = 4.;
        A->vals[5+i*4] = -1.;
        A->vals[6+i*4] = -1.;

        A->colind[3+i*4] = i+INDEX_BASE;
        A->colind[4+i*4] = i+1+INDEX_BASE;
        A->colind[5+i*4] = i+2+INDEX_BASE;
        A->colind[6+i*4] = i+1+N+INDEX_BASE;

        A->rowptr[i+1] = 3+i*4;
    }

    /* Upper most inner nodes */
    for (i=0; i<N-2; ++i) {
        A->vals[5*N*N-8*N+5+i*4] = -1.;
        A->vals[5*N*N-8*N+6+i*4] = -1.;
        A->vals[5*N*N-8*N+7+i*4] = 4.;
        A->vals[5*N*N-8*N+8+i*4] = -1.;

        A->colind[5*N*N-8*N+5+i*4] = N*N-2*N+i+1+INDEX_BASE;
        A->colind[5*N*N-8*N+6+i*4] = N*N-N+i+INDEX_BASE;
        A->colind[5*N*N-8*N+7+i*4] = N*N-N+i+1+INDEX_BASE;
        A->colind[5*N*N-8*N+8+i*4] = N*N-N+i+2+INDEX_BASE;

        A->rowptr[N*N-N+i+1] = 5*N*N-8*N+5+i*4;
    }

    /* Right most inner nodes */
    for (i=0; i<N-2; ++i) {
        A->vals[9*N-8+i*(5*N-2)] = -1.;
        A->vals[9*N-7+i*(5*N-2)] = -1.;
        A->vals[9*N-6+i*(5*N-2)] = 4.;
        A->vals[9*N-5+i*(5*N-2)] = -1.;

        A->colind[9*N-8+i*(5*N-2)] = N*(i+1)-1+INDEX_BASE;
        A->colind[9*N-7+i*(5*N-2)] = N*(i+2)-2+INDEX_BASE;
        A->colind[9*N-6+i*(5*N-2)] = N*(i+2)-1+INDEX_BASE;
        A->colind[9*N-5+i*(5*N-2)] = N*(i+2)-1+N+INDEX_BASE;

        A->rowptr[N*(i+2)-1] = 9*N-8+i*(5*N-2);
    }

    /* Left most inner nodes */
    for (i=0; i<N-2; ++i) {
        A->vals[4*N-2+i*(5*N-2)] = -1.;
        A->vals[4*N-1+i*(5*N-2)] = 4.;
        A->vals[4*N+i*(5*N-2)]   = -1.;
        A->vals[4*N+1+i*(5*N-2)] = -1.;

        A->colind[4*N-2+i*(5*N-2)] = N*i+INDEX_BASE;
        A->colind[4*N-1+i*(5*N-2)] = N*(i+1)+INDEX_BASE;
        A->colind[4*N+i*(5*N-2)]   = N*(i+1)+1+INDEX_BASE;
        A->colind[4*N+1+i*(5*N-2)] = N*(i+2)+INDEX_BASE;

        A->rowptr[N*(i+1)] = 4*N-2+i*(5*N-2);
    }

    /* Purely inner nodes */
    for (i=0; i<N-2; ++i) {
        for (j=0; j<N-2; ++j) {
            A->vals[4*N+2+5*j+i*(5*N-2)] = -1.;
            A->vals[4*N+3+5*j+i*(5*N-2)] = -1.;
            A->vals[4*N+4+5*j+i*(5*N-2)] = 4.;
            A->vals[4*N+5+5*j+i*(5*N-2)] = -1.;
            A->vals[4*N+6+5*j+i*(5*N-2)] = -1.;

            A->colind[4*N+2+5*j+i*(5*N-2)] = N*(i)+1+j+INDEX_BASE;
            A->colind[4*N+3+5*j+i*(5*N-2)] = N*(i+1)+j+INDEX_BASE;
            A->colind[4*N+4+5*j+i*(5*N-2)] = N*(i+1)+1+j+INDEX_BASE;
            A->colind[4*N+5+5*j+i*(5*N-2)] = N*(i+1)+2+j+INDEX_BASE;
            A->colind[4*N+6+5*j+i*(5*N-2)] = N*(i+2)+1+j+INDEX_BASE;

            A->rowptr[N*(i+1)+1+j] = 4*N+2+5*j+i*(5*N-2);
        }
    }
}
