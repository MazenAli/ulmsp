#ifndef REALVECTOR_C
#define REALVECTOR_C

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>

#include "realvector.h"

prealvector
new_realvector(index length)
{
    prealvector x;
    x = (prealvector) malloc(sizeof(realvector));
    init_realvector(x, length);

    return x;
}


void
init_realvector(prealvector x, index length)
{
    assert(x);

    if (length>0) {
        x->vals = (real*) malloc(length*sizeof(real));
    } else {
        x->vals = NULL;
    }
    x->length = length;
}


void
resize_realvector(prealvector x, index length)
{
    assert(x);

    del_realvector(x);
    init_realvector(x, length);
}


void
copy_realvector(prealvector dest, pcrealvector src)
{
    assert(src);

    index i;

    if (dest) del_realvector(dest);
    init_realvector(dest, src->length);

    for (i=INDEX_BASE; i<src->length+INDEX_BASE; ++i) {
        setentry_realvector(dest, i, getentry_realvector(src, i));
    }
}


void
del_realvector(prealvector x)
{
    assert(x);

    if (x->vals) free(x->vals);
    x->length = 0;
}


real
getentry_realvector(pcrealvector x, index i)
{
    assert(x);
    assert(i<x->length+INDEX_BASE);

    return x->vals[i-INDEX_BASE];
}


void
setentry_realvector(prealvector x, index i, real entry)
{
    assert(x);
    assert(i<x->length+INDEX_BASE);

    x->vals[i-INDEX_BASE] = entry;
}


void
addentry_realvector(prealvector x, index i, real entry)
{
    assert(x);
    assert(i<x->length+INDEX_BASE);

    x->vals[i-INDEX_BASE] += entry;
}


void
print_realvector(pcrealvector x)
{
    assert(x);

    index i;

    for (i=INDEX_BASE; i<x->length+INDEX_BASE; ++i) {
        (void) printf("%8.4f\n", getentry_realvector(x, i));
    }
}

#endif
