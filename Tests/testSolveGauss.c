/**
 *  Assemble a sparse matrix in COO format with possible repetitions
 *  of coordinates and convert to CRS.
 */

#include <stdio.h>

#include "settings.h"
#include "coo.h"
#include "crs.h"
#include "indexvector.h"
#include "realvector.h"
#include "realmatrix.h"
#include "indexmatrix.h"
#include "gausscrs.h"

int
main()
{
    pcoo A = new_coo(6, 3, 3);
    pcrs B = new_crs(1, 1, 1);
    prealvector b = new_realvector(3);
    pindexvector fixed = new_indexvector(1);
    pindexvector d;
    prealvector dr;
    prealmatrix dm;
    pindexmatrix di;

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
    A->cols[1] = INDEX_BASE+1;
    A->cols[2] = INDEX_BASE+2;
    A->cols[3] = INDEX_BASE;
    A->cols[4] = INDEX_BASE;
    A->cols[5] = INDEX_BASE;

    /* Setup b */
    d = load_indexvector("./Tests/b.txt");
    print_indexvector(d);
    /* Setup br */
    dr = load_realvector("./Tests/b.txt");
    print_realvector(dr);
    /* Setup bm */
    dm = load_realmatrix("./Tests/b.txt",2,0);
    print_realmatrix(dm);
    /* Setup bm */
    dm = load_realmatrix("./Tests/b.txt",2,1);
    print_realmatrix(dm);
    /* Setup bi */
    di = load_indexmatrix("./Tests/b.txt",2,0);
    print_indexmatrix(di);
    /* Setup bi */
    di = load_indexmatrix("./Tests/b.txt",2,1);
    print_indexmatrix(di);

    (void) printf("COO matrix A=\n");
    print_coo(A);

    init_coo2crs(B, A);
    (void) printf("CRS conversion->\n");
    print_crs(B);

    fixed->vals[0] = INDEX_BASE;
    b->vals[0] = 1.0;

    print_realvector(b);
    gausscrs(B, b, fixed);
    print_realvector(b);

    write_realvector("./Tests/ver.txt",dr);

    del_coo(A);
    del_crs(B);

    return 0;
}
