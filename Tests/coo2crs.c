/**
 *  Assemble a sparse matrix in COO format with possible repetitions
 *  of coordinates and convert to CRS.
 */

#include <stdio.h>

#include "settings.h"
#include "coo.h"
#include "crs.h"

int
main()
{
    pcoo A = new_coo(6, 3, 4);
    pcrs B;

    /* Setup COO with repetitions */
    A->vals[0] = 1;
    A->vals[1] = 8;
    A->vals[2] = 10;
    A->vals[3] = 7;
    A->vals[4] = 5;
    A->vals[5] = 4;

    A->rows[0] = INDEX_BASE;
    A->rows[1] = INDEX_BASE+1;
    A->rows[2] = INDEX_BASE+2;
    A->rows[3] = INDEX_BASE+1;
    A->rows[4] = INDEX_BASE+1;
    A->rows[5] = INDEX_BASE;

    A->cols[0] = INDEX_BASE;
    A->cols[1] = INDEX_BASE+3;
    A->cols[2] = INDEX_BASE+2;
    A->cols[3] = INDEX_BASE;
    A->cols[4] = INDEX_BASE;
    A->cols[5] = INDEX_BASE;

    (void) printf("COO matrix A=\n");
    print_coo(A);

    B = new_coo2crs(A);
    (void) printf("CRS conversion->\n");
    print_crs(B);

    /* COO w/o repetitions */
    A->cols[3] = INDEX_BASE+1;
    A->cols[5] = INDEX_BASE+1;

    del_crs(B);
    B = new_coo2crs(A);

    (void) printf("\nW/o repetitions\n");
    (void) printf("COO matrix A=\n");
    print_coo(A);

    B = new_coo2crs(A);
    (void) printf("CRS conversion->\n");
    print_crs(B);

    /* COO with zero row and repetitions */
    A->rows[1] = INDEX_BASE;
    A->rows[3] = INDEX_BASE;
    A->rows[4] = INDEX_BASE+2;

    del_crs(B);
    B = new_coo2crs(A);

    (void) printf("\nWith repetitions and zero rows\n");
    (void) printf("COO matrix A=\n");
    print_coo(A);

    B = new_coo2crs(A);
    (void) printf("CRS conversion->\n");
    print_crs(B);


    del_coo(A);
    del_crs(B);

    return 0;
}
