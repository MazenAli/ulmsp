#ifndef CRS_C
#define CRS_C

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

/* Conflicting typedefs for "index"*/
#define index string_index
    #include <string.h>
#undef index

#include "crs.h"

pcrs
new_crs(index nonz, index numr, index numc)
{
    pcrs A;

    assert(nonz<=numr*numc);

    A = (pcrs) malloc(sizeof(crs));
    init_crs(A, nonz, numr, numc);

    return A;
}


pcrs
new_coo2crs(pccoo C)
{
    assert(C);

    pcrs A;
    A = (pcrs) malloc(sizeof(crs));
    init_coo2crs(A, C);

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
init_coo2crs(pcrs A, pccoo C)
{
    index k, i;
    index currentrow, offset;
    index *rowptr, *colind;
    real  *vals;
    int   *flags;

    assert(A);
    assert(C);

    /* In case C empty */
    if (!C->nonz) {
        A->vals   = NULL;
        A->rowptr = NULL;
        A->colind = NULL;
        A->nonz   = 0;
        A->numr   = C->numr;
        A->numc   = C->numc;

        return;
    }


    /* Set meta data */
    A->numr = C->numr;
    A->numc = C->numc;

    vals   = (real*)  malloc(C->nonz*sizeof(real));
    rowptr = (index*) malloc((A->numr+1)*sizeof(index));
    colind = (index*) malloc(C->nonz*sizeof(index));
    flags  = (int*)  malloc(C->nonz*sizeof(int));
    memset(rowptr, 0, (A->numr+1)*sizeof(index));
    memset(flags, 0, C->nonz*sizeof(int));

    /* Count number of each row */
    for (k=0; k<C->nonz; ++k) {
        ++rowptr[C->rows[k]-INDEX_BASE+1];
    }

    /* Accumulate */
    for (k=1; k<A->numr; ++k) {
        rowptr[k+1] += rowptr[k];
    }

    /* Add and flag repetitions */
    for (k=0; k<C->nonz; ++k) {
        for (i=rowptr[C->rows[k]-INDEX_BASE];
             i<rowptr[C->rows[k]-INDEX_BASE+1]; ++i) {

            if (flags[i]) {
                if (C->cols[k]==colind[i]) {
                    vals[i] += C->vals[k];
                    break;
                }
            } else {
                ++A->nonz;
                flags[i]  = 1;
                colind[i] = C->cols[k];
                vals[i]   = C->vals[k];
                break;
            }
        }
    }


    /* No repetitions */
    if (A->nonz==C->nonz) {
        A->vals   = vals;
        A->rowptr = rowptr;
        A->colind = colind;
        free(flags);

        return;
    }

    /* Eliminate repetitions */
    A->vals      = (real*)  malloc(A->nonz*sizeof(real));
    A->colind    = (index*) malloc(A->nonz*sizeof(index));

    currentrow = 0;
    offset     = 0;
    for (k=0; k<C->nonz; ++k) {
        if (k==rowptr[currentrow+1]) {
            rowptr[currentrow+1] -= offset;
            ++currentrow;
        }

        if (!flags[k]) {
            offset += rowptr[currentrow+1]-k;
            k       = rowptr[currentrow+1]-1;
        } else {
            A->vals[k-offset]   = vals[k];
            A->colind[k-offset] = colind[k];
        }
    }
    rowptr[A->numr] -= offset;
    A->rowptr = rowptr;

    free(vals);
    free(colind);
    free(flags);
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
    index    ii;
    real entry;

    assert(i<A->numr+INDEX_BASE);
    assert(j<A->numc+INDEX_BASE);

    entry = 0.;

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
    index    ii;

    assert(i<A->numr+INDEX_BASE);
    assert(j<A->numc+INDEX_BASE);

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
    index    ii;

    assert(i<A->numr+INDEX_BASE);
    assert(j<A->numc+INDEX_BASE);

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
    index i, j;

    assert(A);

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
   index i, j;

   assert(A);

    for (i=INDEX_BASE; i<A->numr+INDEX_BASE; ++i) {
        for (j=INDEX_BASE; j<A->numc+INDEX_BASE; ++j) {
            (void) printf("%8.4f", getentry_crs(A, i, j));
        }
        (void) printf("\n");
    }
}



#endif
