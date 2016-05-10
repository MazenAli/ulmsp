/*! \file realmatrix.h
 *  A simple real matrix data type.
 */

#ifndef REALMATRIX_H
#define REALMATRIX_H

#include "../settings.h"

typedef struct _realmatrix realmatrix;    /*!< Struct typedef: use realmatrix
                                               to declare. */
typedef realmatrix         *prealmatrix;  /*!< Pointer to realmatrix. */
typedef const realmatrix   *pcrealmatrix; /*!< Const pointer to realmatrix. */

/*! \struct _realmatrix "realmatrix.h"
 *  \brief A struct for storing a matrix of reals.
 */
struct _realmatrix
{
    real    *vals;   /*!< Pointer to matrix entries. */
    index    rows;   /*!< Number of rows. */
    index    cols;   /*!< Number of columns. */
};

/*! \name Constructors and destructors. */
/**@{*/
/*! Initialize new realmatrix and return pointer. */
prealmatrix
new_realmatrix(index rows, index cols);

/*! Initialize new realmatrix and set pointer x. */
void
init_realmatrix(prealmatrix A, index rows, index cols);

/*! Resize realmatrix x. */
void
resize_realmatrix(prealmatrix A, index rows, index cols);

/*! Copy realmatrix dest<-src. */
void
copy_realmatrix(prealmatrix dest, pcrealmatrix src);

/*! Swap data A<->B. */
void
swap_realmatrix(prealmatrix A, prealmatrix B);

/*! Free memory. */
void
del_realmatrix(prealmatrix A);
/**@}*/


/*! \name Access methods. */
/**@{*/
/*! Get entry \f$ A_{i, j} \f$. */
real
getentry_realmatrix(pcrealmatrix A, index i, index j);

/*! Set entry \f$ A_{i, j} \f$. */
void
setentry_realmatrix(prealmatrix A, index i, index j, real entry);

/*! Add to entry \f$ A_{i, j} \f$. */
void
addentry_realmatrix(prealmatrix A, index i, index j, real entry);

/*! Print matrix.*/
void
print_realmatrix(pcrealmatrix A);

/*! Load realmatrix from file, return pointer. */
prealmatrix
load_realmatrix(char *fname, index cols, transpose t);

/*! Write realmatrix to file, return status. */
int
write_realmatrix(char *fname, prealmatrix A, transpose t);
/**@}*/

#endif
