#ifndef CRS_C
#define CRS_C

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include<time.h>

/* Conflicting typedefs for "index"*/
#define index string_index
    #include <string.h>
#undef index

#include "crs.h"

#define MAX(a,b) ((a) > (b) ? (a) : (b))
#define ELAPSED(t0,t1) ((int) ((t1 - t0) / (double) CLOCKS_PER_SEC * 1000))


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
    pcrs A;

    assert(C);

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
init_coo2crs(pcrs A, pccoo T)
{
    index k, i;
    index *rowptr, *colind;
    real  *vals;
    index p, q, nz;
    index *Ti, *Tj, *work;
    real  *Tx;

    /* Checks, if C and T are no NULL pointer*/
    assert(T);
    assert(A);

    /* In case C empty */
    if (!T->nonz) {
        A->vals   = NULL; A->rowptr = NULL;    A->colind = NULL;
        A->nonz   = 0;    A->numr   = T->numr; A->numc   = T->numc;

        return;
    }

    /* Set dimensioin of matrix */
    A->numr = T->numr;  A->numc = T->numc;

    /* Allocate memory for CRS matrix and some workspace */
    vals   = (real*)  malloc( T->nonz    * sizeof(real));
    colind = (index*) malloc( T->nonz    * sizeof(index));
    rowptr = (index*) malloc((T->numr+1) * sizeof(index));
    work   = (index*) calloc(MAX(T->numr, T->numc)  , sizeof(index));
                                             /* zero-initialize */

    /* Check, if out of memory*/
    if (!vals || !rowptr || !colind || !work){ /* out of memory*/
        free(vals);
        free(rowptr);
        free(colind);
        free(work);

        return;
    }

    Tx = T->vals; Ti = T->rows; Tj = T->cols;

    /*
    * Convert coordinates format to CRS format,
    * allow duplicate entries
    */

    /* Count entries per row */
    for (k = 0 ; k < T->nonz ; k++) work[Ti[k]-INDEX_BASE]++ ;
    /* Create col pointer */
    cumsum(rowptr,work,T->numr);
    /* Copy data from T to C */
    for (k = 0 ; k < T->nonz ; k++){
        colind[p = work[Ti[k]-INDEX_BASE]++] = Tj[k];
        vals[p] = Tx[k] ;
    }

    /*
    * Remove duplicate entries
    */

    for (k = 0 ; k < T->numc ; k++) work[k]=0;   /* col k yet not seen =
    */
    /* Loop over each row and check if A(i,j) is duplicate */
    nz = 0;
    for (i = 0 ; i < T->numr ; i++){
        q = nz;                           /* Row i will start at q */
        for (p = rowptr[i] ; p < rowptr[i+1] ; p++){
            k = colind[p]-INDEX_BASE ;      /* Entry A(i,j) */
            if (work[k] > q){
                vals[work[k]-1] += vals[p] ;  /* A(i,j) is a duplicate */
            } else {
                work[k] = nz+1;               /* Record where col j occurs */
                colind[nz] = k+INDEX_BASE ;   /* Keep A(i,j) */
                vals[nz++] = vals[p] ;
            }
        }

        rowptr[i] = q ;                   /* Record start of row i */
    }

    rowptr[T->numr] = nz ;              /* Finalize A */
    free(work);                         /* Free workspace */

    /*
    * Drop entries which are zero
    */

    nz = 0;
    for (i = 0 ; i < T->numr ; i++){
        p = rowptr[i];                    /* Get current location of row i =
        */
        rowptr[i] = nz;                   /* Record new location of row i =
        */
        for ( ; p < rowptr [i+1] ; p++){
            if (vals[p]!=0){
            vals[nz] = vals[p];           /* Keep A(i,j) */
            colind[nz++] = colind[p];
            }
        }
    }

    rowptr[A->numr] = nz ;              /* Finalize A */

    /*
    * Realloc memory of compressed matrix
    */

    A->vals   = (real*)  realloc(vals,   nz * sizeof(real));
    A->colind = (index*) realloc(colind, nz * sizeof(index));
    A->rowptr = rowptr;

    /* check for memory error*/
    return;
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


real
cumsum (index *p, index *c, index n)
{
    index i, nz = 0 ;
    real  nz2 = 0 ;
    if (!p || !c) return (-1) ;     /* check inputs */
    for (i = 0 ; i < n ; i++){
    p [i] = nz ;
    nz += c [i] ;
    nz2 += c [i] ;        /* also in double to avoid CS_INT overflow =
    */
    c [i] = p [i] ;       /* also copy p[0..n-1] back into c[0..n-1]*/
    }
    p [n] = nz ;

    return (nz2);           /* return sum (c [0..n-1]) */
}

#endif
