/*! \file indexmatrix.h
 *  A simple real matrix data type.
 */

#ifndef INDEXMATRIX_H
#define INDEXMATRIX_H

#include "../settings.h"

typedef struct _indexmatrix indexmatrix;    /*!< Struct typedef: use indexmatrix
                                                 to declare. */
typedef indexmatrix         *pindexmatrix;  /*!< Pointer to indexmatrix. */
typedef const indexmatrix   *pcindexmatrix; /*!< Const pointer
                                                 to indexmatrix. */

/*! \struct _indexmatrix "indexmatrix.h"
 *  \brief A struct for storing a matrix of indices.
 */
struct _indexmatrix
{
    index    *vals;  /*!< Pointer to matrix entries. */
    index    rows;   /*!< Number of rows. */
    index    cols;   /*!< Number of columns. */
};

/*! \name Constructors and destructors. */
/**@{*/
/*! Initialize new indexmatrix and return pointer. */
pindexmatrix
new_indexmatrix(index rows, index cols);

/*! Initialize new indexmatrix and set pointer x. */
void
init_indexmatrix(pindexmatrix A, index rows, index cols);

/*! Resize indexmatrix x. */
void
resize_indexmatrix(pindexmatrix A, index rows, index cols);

/*! Copy indexmatrix dest<-src. */
void
copy_indexmatrix(pindexmatrix dest, pcindexmatrix src);

/*! Swap data A<->B. */
void
swap_indexmatrix(pindexmatrix A, pindexmatrix B);

/*! Free memory. */
void
del_indexmatrix(pindexmatrix A);
/**@}*/


/*! \name Access methods. */
/**@{*/
/*! Get entry \f$ A_{i, j} \f$. */
index
getentry_indexmatrix(pcindexmatrix A, index i, index j);

/*! Set entry \f$ A_{i, j} \f$. */
void
setentry_indexmatrix(pindexmatrix A, index i, index j, index entry);

/*! Add to entry \f$ A_{i, j} \f$. */
void
addentry_indexmatrix(pindexmatrix A, index i, index j, index entry);

/*! Print matrix.*/
void
print_indexmatrix(pcindexmatrix A);

/*! Load indexmatrix from file, return pointer. */
pindexmatrix
load_indexmatrix(char *fname, index rows, transpose t);

/*! Write indexmatrix to file, return status. */
int
write_indexmatrix(char *fname, pindexmatrix A, transpose t);
/**@}*/

#endif
