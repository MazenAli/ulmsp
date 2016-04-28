/**
 *  Compressed row storage (CRS) sparse matrix storage format.
 *  Assumes row sorting of entries.
 */

#ifndef CRS_H
#define CRS_H

#include "../settings.h"
#include "coo.h"

typedef struct _crs crs;
typedef crs         *pcrs;
typedef const crs   *pccrs;

struct _crs
{
    real    *vals;    ///< Pointer to matrix entries
    index    *rowptr; ///< Pointer to row indices
    index    *colind; ///< Pointer to column indices
    index    nonz;    ///< Number of non zeros
    index    numr;    ///< Number of rows
    index    numc;    ///< Number of columns
};

/**
 *  Constructors and destructors.
 */

pcrs
new_crs(index nonz, index numr, index numc);

pcrs
new_coo2crs(pccoo C);

void
init_crs(pcrs A, index nonz, index numr, index numc);

void
init_coo2crs(pcrs A, pccoo C);

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
