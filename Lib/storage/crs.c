#ifndef CRS_C
#define CRS_C

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

#include "crs.h"

pcrs
new_crs(index nonz, index numr, index numc)
{
    assert(nonz<=numr*numc);

    pcrs A;
    A = (pcrs) malloc(sizeof(crs));
    init_crs(A, nonz, numr, numc);

    return A;
}


void
init_crs(pcrs A, index nonz, index numr, index numc)
{
    assert(A);
    assert(nonz<=numr*numc);

    if (nonz>0) {
        assert(numr>0);
        assert(numc>0);
        A->vals   = (real*) malloc(nonz*sizeof(real));
        A->rowptr = (index*) malloc((numr+1)*sizeof(index));
        A->colind = (index*) malloc(nonz*sizeof(index));
    } else {
        assert(!numr);
        assert(!numc);
        A->vals   = NULL;
        A->rowptr = NULL;
        A->colind = NULL;
    }
    A->nonz = nonz;
    A->numr = numr;
    A->numc = numc;
}


void
resize_crs(pcrs A, index nonz, index numr, index numc)
{
    assert(A);
    assert(nonz<=numr*numc);

    del_crs(A);
    init_crs(A, nonz, numr, numc);
}


void
del_crs(pcrs A)
{
    assert(A);

    if (A->vals)   free(A->vals);
    if (A->rowptr) free(A->rowptr);
    if (A->colind) free(A->colind);
    A->nonz = 0;
    A->numr = 0;
    A->numc = 0;
}


real
getentry_crs(pccrs A, index i, index j)
{
    assert(i<A->numr+INDEX_BASE);
    assert(j<A->numc+INDEX_BASE);

    index    ii;
    real entry = 0.;

    for (ii=A->rowptr[i-INDEX_BASE];
         ii<A->rowptr[i-INDEX_BASE+1]; ++ii) {
        if (A->colind[ii]==j) {
            entry = A->vals[ii];
        }
    }

    return entry;
}


void
setentry_crs(pcrs A, index i, index j, real entry)
{
    assert(i<A->numr+INDEX_BASE);
    assert(j<A->numc+INDEX_BASE);

    index    ii;

    for (ii=A->rowptr[i-INDEX_BASE];
         ii<A->rowptr[i-INDEX_BASE+1]; ++ii) {
        if (A->colind[ii]==j) {
            A->vals[ii] = entry;
            return;
        }
    }

    fprintf(stderr, "Index not in range of sparse matrix\n");
}


void
addentry_crs(pcrs A, index i, index j, real entry)
{
    assert(i<A->numr+INDEX_BASE);
    assert(j<A->numc+INDEX_BASE);

    index    ii;

    for (ii=A->rowptr[i-INDEX_BASE];
         ii<A->rowptr[i-INDEX_BASE+1]; ++ii) {
        if (A->colind[ii]==j) {
            A->vals[ii] += entry;
            return;
        }
    }

    fprintf(stderr, "Index not in range of sparse matrix\n");
}


void
print_crs(pccrs A)
{
    assert(A);

    index i, j;

    for (i=0; i<A->numr; ++i) {
        for (j=A->rowptr[i]; j<A->rowptr[i+1]; ++j) {
            (void) printf("(%lu, %lu) : %8.4f\n",
                          i+INDEX_BASE,
                          A->colind[j],
                          A->vals[j]);
        }
    }
}


void
printdense_crs(pccrs A)
{
    assert(A);

    index i, j;

    for (i=INDEX_BASE; i<A->numr+INDEX_BASE; ++i) {
        for (j=INDEX_BASE; j<A->numc+INDEX_BASE; ++j) {
            (void) printf("%8.4f", getentry_crs(A, i, j));
        }
        (void) printf("\n");
    }
}



#endif
