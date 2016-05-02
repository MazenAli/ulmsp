/**
 *  Setting up and working with realvector's.
 */

#include <stdio.h>
#include <stdlib.h>

#include "settings.h"
#include "realvector.h"

int
main(int argc, char **argv)
{
    int n;
    index i;
    real alpha, dot, norm;
    prealvector x, y;

    if (argc!=2) {
        (void) fprintf(stderr, "Usage: ./main n\n");
        return 1;
    }

    n = atoi(argv[1]);
    alpha = 0.5;

    /* Set up vectors */
    x = new_realvector(n);
    y = new_realvector(n);

    for (i=INDEX_BASE; i<x->length+INDEX_BASE; ++i) {
        setentry_realvector(x, i, (real) i);
        setentry_realvector(y, i, (real) i);
    }

    /* Perform x <- alpha*x */
    scal_realvector(alpha, x);

    (void) printf("y =\n");
    print_realvector(y);

    (void) printf("x post scaling =\n");
    print_realvector(x);

    /* Perform x^T*y */
    dot = dot_realvector(x, y);
    (void) printf("dot product = %f\n", dot);

    /* Perform |x|_2 */
    norm = nrm2_realvector(x);
    (void) printf("nrm2 x = %f\n", norm);

    /* Perform |y|_2 */
    norm = nrm2_realvector(y);
    (void) printf("nrm2 y = %f\n", norm);

    /* Perform y <- alpha*x + y */
    axpy_realvector(alpha, x, y);
    (void) printf("alpha*x+y =\n");
    print_realvector(y);

    del_realvector(x);
    del_realvector(y);

    return 0;
}
