/**
 *  Setting up and working with sparse COO matrices.
 */

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#include "settings.h"
#include "coo.h"
#include "realvector.h"
#include "gecoomv.h"

/// Fill COO matrix A with some values
void
setupdense_coo(pcoo A);


/// Add to COO matrix A some values
void
adddense_coo(pcoo A);


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
    real alpha = 2., beta = 0.5; ///< gecoomv parameters

    /// Set up, resize and change entries of COO matrix A
    pcoo A = new_coo(0, 0, 0);
    resize_coo(A, 3*min-2, rows, cols);

    setupdense_coo(A);
    adddense_coo(A);
    setentry_coo(A, rows/2, 1, 3);
    setentry_coo(A, rows/2, cols-1, -0.5);
    setentry_coo(A, 2, cols/2, -0.5);

    printdense_coo(A);
    (void) printf("A =\n");
    print_coo(A);

    /// Set up, resize and change entries of vector x
    prealvector x = new_realvector(rows);
    resize_realvector(x, rows);

    setentry_realvector(x, INDEX_BASE, 2.5);
    addentry_realvector(x, INDEX_BASE, 1.5);
    setentry_realvector(x, rows-1+INDEX_BASE, 2.5);
    setentry_realvector(x, (INDEX_BASE+rows/2), -2.5);

    (void) printf("\nx =\n");
    print_realvector(x);

    /// Set up, resize and change entries of vector x
    prealvector y = new_realvector(cols);
    setentry_realvector(y, INDEX_BASE, 1.);
    setentry_realvector(y, INDEX_BASE+1, 0.5);
    setentry_realvector(y, INDEX_BASE+3, 2);
    setentry_realvector(y, cols-1+INDEX_BASE, -1.5);
    (void) printf("\ny =\n");
    print_realvector(y);

    /// Perform y <- alpha*A^T*x+beta*y
    transpose t   = trans;
    gecoomv(t, alpha, A, x, beta, y);

    (void) printf("\nalpha*Op(A)*x + beta*y =\n");
    print_realvector(y);

    del_realvector(x);
    del_realvector(y);
    del_coo(A);

    return 0;
}


void
setupdense_coo(pcoo A)
{
    assert(A);
    index min = (A->numr<A->numc) ? A->numr : A->numc;
    assert(A->nonz==3*min-2);

    index i;

    for (i=0; i<A->nonz-4; ++i) {
        A->rows[i+2] = i/3+INDEX_BASE+1;
        A->cols[i+2] = i%3+INDEX_BASE+i/3;
    }

    A->rows[0]          = 0+INDEX_BASE;
    A->rows[1]          = 0+INDEX_BASE;
    A->rows[A->nonz-2]  = min+INDEX_BASE-1;
    A->rows[A->nonz-1]  = min+INDEX_BASE-1;

    A->cols[0]          = 0+INDEX_BASE;
    A->cols[1]          = 1+INDEX_BASE;
    A->cols[A->nonz-2]  = min+INDEX_BASE-2;
    A->cols[A->nonz-1]  = min+INDEX_BASE-1;

    for (i=INDEX_BASE; i<min+INDEX_BASE; ++i) {
        setentry_coo(A, i, i, 2);
        if (i>INDEX_BASE) setentry_coo(A, i, i-1, -0.5);
        if (i+1<min+INDEX_BASE) setentry_coo(A, i, i+1, -0.5);
    }
}


void
adddense_coo(pcoo A)
{
    assert(A);
    index min = (A->numr<A->numc) ? A->numr : A->numc;
    assert(A->nonz==3*min-2);

    index i;

    for (i=INDEX_BASE; i<min+INDEX_BASE; ++i) {
        addentry_coo(A, i, i, 2);
        if (i>INDEX_BASE) addentry_coo(A, i, i-1, -0.5);
        if (i+1<min+INDEX_BASE) addentry_coo(A, i, i+1, -0.5);
    }
}

