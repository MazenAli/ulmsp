/**
 *  Setting up and working with sparse CRS matrices.
 */

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#include "settings.h"
#include "crs.h"
#include "realvector.h"
#include "gecrsmv.h"

/* Fill CRS matrix A with some values */
void
setupdense_crs(pcrs A);


/* Add to CRS matrix A some values */
void
adddense_crs(pcrs A);


int
main(int argc, char **argv)
{
    if (argc!=3) {
        (void) fprintf(stderr, "Usage: ./main numr numc\n");
        return 1;
    }

    int rows = atoi(argv[1]);
    int cols = atoi(argv[2]);
    index min = (rows<cols) ? rows : cols;
    real alpha = 2., beta = 0.5; /* gecrsmv parameters */

    /* Set up, resize and change entries of CRS matrix A */
    pcrs A = new_crs(0, 0, 0);
    resize_crs(A, 3*min-2, rows, cols);

    setupdense_crs(A);
    adddense_crs(A);
    setentry_crs(A, rows/2, 1, 3);
    setentry_crs(A, rows/2, cols-1, -0.5);
    setentry_crs(A, 2, cols/2, -0.5);

    printdense_crs(A);
    (void) printf("A =\n");
    print_crs(A);

    /* Set up, resize and change entries of vector x */
    prealvector x = new_realvector(rows);
    resize_realvector(x, rows);

    setentry_realvector(x, INDEX_BASE, 2.5);
    addentry_realvector(x, INDEX_BASE, 1.5);
    setentry_realvector(x, rows-1+INDEX_BASE, 2.5);
    setentry_realvector(x, (INDEX_BASE+rows/2), -2.5);

    (void) printf("\nx =\n");
    print_realvector(x);

    /* Set up, resize and change entries of vector x */
    prealvector y = new_realvector(cols);
    setentry_realvector(y, INDEX_BASE, 1.);
    setentry_realvector(y, INDEX_BASE+1, 0.5);
    setentry_realvector(y, INDEX_BASE+3, 2);
    setentry_realvector(y, cols-1+INDEX_BASE, -1.5);
    (void) printf("\ny =\n");
    print_realvector(y);

    /* Perform y <- alpha*A^T*x+beta*y */
    transpose t   = trans;
    gecrsmv(t, alpha, A, x, beta, y);

    (void) printf("\nalpha*Op(A)*x + beta*y =\n");
    print_realvector(y);

    del_realvector(x);
    del_realvector(y);
    del_crs(A);

    return 0;
}


void
setupdense_crs(pcrs A)
{
    assert(A);
    index min = (A->numr<A->numc) ? A->numr : A->numc;
    assert(A->nonz==3*min-2);
    assert(A->nonz>2);

    index i;

    for (i=0; i<A->numr; ++i) {
        if (i==0) {
            A->rowptr[i]   = 0;
            A->rowptr[i+1] = 2;
        } else if (i+1>=min) {
            A->rowptr[i+1] = A->nonz;
        } else {
            A->rowptr[i+1] = A->rowptr[i]+3;
        }
    }

    A->colind[0] = 0+INDEX_BASE;
    A->colind[1] = 1+INDEX_BASE;
    A->colind[2] = 0+INDEX_BASE;
    for (i=3; i<A->nonz; ++i) {
        A->colind[i] = A->colind[i-3]+1;
    }

    for (i=INDEX_BASE; i<min+INDEX_BASE; ++i) {
        setentry_crs(A, i, i, 2);
        if (i>INDEX_BASE) setentry_crs(A, i, i-1, -0.5);
        if (i+1<min+INDEX_BASE) setentry_crs(A, i, i+1, -0.5);
    }
}


void
adddense_crs(pcrs A)
{
    assert(A);
    index min = (A->numr<A->numc) ? A->numr : A->numc;
    assert(A->nonz==3*min-2);

    index i;

    for (i=INDEX_BASE; i<min+INDEX_BASE; ++i) {
        addentry_crs(A, i, i, 2);
        if (i>INDEX_BASE) addentry_crs(A, i, i-1, -0.5);
        if (i+1<min+INDEX_BASE) addentry_crs(A, i, i+1, -0.5);
    }
}

