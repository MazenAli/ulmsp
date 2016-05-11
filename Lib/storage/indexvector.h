/*! \file indexvector.h
 *  A simple real matrix data type.
 */

#ifndef INDEXVECTOR_H
#define INDEXVECTOR_H

#include "../settings.h"

typedef struct _indexvector indexvector;    /*!< Struct typedef:
                                                 use indexvector to declare. */
typedef indexvector         *pindexvector;  /*!< Pointer to indexvector. */
typedef const indexvector   *pcindexvector; /*!< Const pointer to
                                                 indexvector. */

/*! \struct _indexvector "indexvector.h"
 *  \brief A struct for storing a vector of indices.
 */
struct _indexvector
{
    index    *vals;  /*!< Pointer to vector entries. */
    index    length; /*!< Length of vector. */
};

/*! \name Constructors and destructors. */
/**@{*/
/*! Initialize new indexvector and return pointer. */
pindexvector
new_indexvector(index length);

/*! Initialize new indexvector and set pointer x. */
void
init_indexvector(pindexvector x, index length);

/*! Resize indexvector x. */
void
resize_indexvector(pindexvector x, index length);

/*! Copy indexvector dest<-src. */
void
copy_indexvector(pindexvector dest, pcindexvector src);

/*! Swap data x<->y. */
void
swap_indexvector(pindexvector x, pindexvector y);

/*! Free memory. */
void
del_indexvector(pindexvector x);
/**@}*/


/*! \name Access methods. */
/**@{*/
/*! Get entry \f$ x_i \f$. */
index
getentry_indexvector(pcindexvector x, index i);

/*! Set entry \f$ x_i \f$. */
void
setentry_indexvector(pindexvector x, index i, index entry);

/*! Add to entry \f$ x_i \f$. */
void
addentry_indexvector(pindexvector x, index i, index entry);

/*! Fill vector \f$ x_i=val \f$ for all \f$ i \f$. */
void
fill_indexvector(pindexvector x, index val);

/*! Print vector.*/
void
print_indexvector(pcindexvector x);

/*! Load indexvector from file, return pointer. */
pindexvector
load_indexvector(char *fname);

/*! Write indexvector to file, return status. */
int
write_indexvector(char *fname, pindexvector v);
/**@}*/

#endif
