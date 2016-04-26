/**
 *  Compressed row storage (CRS) sparse matrix storage format.
 *  Assumes row sorting of entries.
 */

#ifndef CRS_H
#define CRS_H

#include "../settings.h"

typedef struct _crs crs;
typedef crs         *pcrs;
typedef const crs   *pccrs;

struct _crs
{
    real    *vals;
    index    *rowptr;
    index    *colind;
    index    nonz;
    index    numr;
    index    numc;
};

/**
 *  Constructors and destructors.
 */

pcrs
new_crs(index nonz, index numr, index numc);

void
init_crs(pcrs A, index nonz, index numr, index numc);

void
resize_crs(pcrs A, index nonz, index numr, index numc);

void
del_crs(pcrs A);

/**
 *  Access methods.
 */

real
getentry_crs(pccrs A, index i, index j);

void
setentry_crs(pcrs A, index i, index j, real entry);

void
addentry_crs(pcrs A, index i, index j, real entry);

void
print_crs(pccrs A);

void
printdense_crs(pccrs A);

#endif
