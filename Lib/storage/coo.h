/*! \file coo.h
 *  Coordinate list (COO) sparse matrix storage format and functions.
 */

#ifndef COO_H
#define COO_H

#include "../settings.h"

typedef struct _coo coo;    /*!< Struct typedef: use coo to declare. */
typedef coo         *pcoo;  /*!< Pointer to coo. */
typedef const coo   *pccoo; /*!< Const pointer to coo. */

/*! \struct _coo "coo.h"
 *  \brief A struct for storing a COO matrix.
 */
struct _coo
{
    real    *vals;  /*!< Pointer to matrix entries. */
    index    *rows; /*!< Pointer to row indices. */
    index    *cols; /*!< Pointer to column indices. */
    index    nonz;  /*!< Number of non zeros. */
    index    numr;  /*!< Number of rows. */
    index    numc;  /*!< Number of columns. */
};

/*! \name Constructors and destructors. */
/**@{*/
/*! Initialize new coo and return pointer. */
pcoo
new_coo(index nonz, index numr, index numc);

/*! Initialize new coo and set pointer A. */
void
init_coo(pcoo A, index nonz, index numr, index numc);

/*! Resize coo A. */
void
resize_coo(pcoo A, index nonz, index numr, index numc);

/*! Free memory. */
void
del_coo(pcoo A);
/**@}*/


/*! \name Access methods. */
/**@{*/
/*! Get entry \f$ A_{i, j} \f$. Inefficient access.*/
real
getentry_coo(pccoo A, index i, index j);

/*! Set entry \f$ A_{i, j} \f$. Inefficient access.*/
void
setentry_coo(pcoo A, index i, index j, real entry);

/*! Add to entry \f$ A_{i, j} \f$. Inefficient access.*/
void
addentry_coo(pcoo A, index i, index j, real entry);

/*! Print A. */
void
print_coo(pccoo A);

/*! Print A in dense format. Only for small sizes. */
void
printdense_coo(pccoo A);
/**@}*/

#endif
