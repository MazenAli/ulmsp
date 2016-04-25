/**
 *  Coordinate list (COO) sparse matrix storage format.
 *  Assumes row sorting of entries.
 */

#ifndef COO_H
#define COO_H

#include "../settings.h"

typedef struct _coo coo;
typedef coo         *pcoo;
typedef const coo   *pccoo;

struct _coo
{
    real    *vals;
    index    *rows;
    index    *cols;
    index    nonz;
    index    numr;
    index    numc;
};

/**
 *  Constructors and destructors.
 */

pcoo
new_coo(index nonz, index numr, index numc);

void
init_coo(pcoo A, index nonz, index numr, index numc);

void
resize_coo(pcoo A, index nonz, index numr, index numc);

void
del_coo(pcoo A);

/**
 *  Access methods.
 */

real
getentry_coo(pccoo A, index i, index j);

void
setentry_coo(pcoo A, index i, index j, real entry);

void
addentry_coo(pcoo A, index i, index j, real entry);

void
print_coo(pccoo A);

void
printdense_coo(pccoo A);

#endif
