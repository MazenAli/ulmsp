/*! \file crs.h
 *  Compressed row (CRS) sparse matrix storage format and functions.
 */

#ifndef CRS_H
#define CRS_H

#include "../settings.h"
#include "coo.h"

typedef struct _crs crs;    /*!< Struct typedef: use crs to declare. */
typedef crs         *pcrs;  /*!< Pointer to crs. */
typedef const crs   *pccrs; /*!< Const pointer to crs. */

/*! \struct _crs "crs.h"
 *  \brief A struct for storing a CRS matrix.
 */
struct _crs
{
    real    *vals;    /*!< Pointer to matrix entries */
    index    *rowptr; /*!< Pointer to row indices */
    index    *colind; /*!< Pointer to column indices */
    index    nonz;    /*!< Number of non zeros */
    index    numr;    /*!< Number of rows */
    index    numc;    /*!< Number of columns */
};

/*! \name Constructors and destructors. */
/**@{*/
/*! Initialize new crs and return pointer. */
pcrs
new_crs(index nonz, index numr, index numc);

/*! Convert COO C to a CRS and return pointer. */
pcrs
new_coo2crs(pccoo C);

/*! Initialize new crs and set pointer A. */
void
init_crs(pcrs A, index nonz, index numr, index numc);

/*! Convert COO C to CRS A. */
void
init_coo2crs(pcrs A, pccoo C);

/*! Resize crs A. */
void
resize_crs(pcrs A, index nonz, index numr, index numc);

/*! Free memory. */
void
del_crs(pcrs A);
/**@}*/


/*! \name Access methods. */
/**@{*/
/*! Get entry \f$ A_{i,j} \f$. Inefficient access. */
real
getentry_crs(pccrs A, index i, index j);

/*! Set entry \f$ A_{i,j} \f$. Inefficient access. */
void
setentry_crs(pcrs A, index i, index j, real entry);

/*! Add to entry \f$ A_{i,j} \f$. Inefficient access. */
void
addentry_crs(pcrs A, index i, index j, real entry);

/*! Print A. */
void
print_crs(pccrs A);

/*! Print A in dense format. Only for small sizes. */
void
printdense_crs(pccrs A);
/**@}*/

/**
 *  Auxiliary functions
 */

/**
*/

/*! \name Auxiliary functions. */
/**@{*/
/*!
 * \param[in]  c [0..n-1]
 * \param[in]  n size
 * \param[out] p [0..n] = cumulative sum of c [0..n-1] and
 * \param[out] c [0..n-1] copy of p [0..n-1]
 */
real
cumsum(index *p, index *c, index n);
/**@}*/

#endif
