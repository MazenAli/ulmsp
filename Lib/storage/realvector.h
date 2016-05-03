/*! \file realvector.h
 *  A simple real vector data type.
 */

#ifndef REALVECTOR_H
#define REALVECTOR_H

#include "../settings.h"

typedef struct _realvector realvector;    /*!< Struct typedef: use realvector
                                               to declare. */
typedef realvector         *prealvector;  /*!< Pointer to realvector. */
typedef const realvector   *pcrealvector; /*!< Const pointer to realvector. */

/*! \struct _realvector "realvector.h"
 *  \brief A struct for storing a vector of reals.
 */
struct _realvector
{
    real    *vals;   /*!< Pointer to vector entries */
    index    length; /*!< Length of vector */
};


/*! \name Constructors and destructors. */
/**@{*/
/*! Initialize new realvector and return pointer. */
prealvector
new_realvector(index length);

/*! Initialize new realvector and set pointer x. */
void
init_realvector(prealvector x, index length);

/*! Resize realvector x. */
void
resize_realvector(prealvector x, index length);

/*! Copy realvector dest<-src. */
void
copy_realvector(prealvector dest, pcrealvector src);

/*! Free memory. */
void
del_realvector(prealvector x);
/**@}*/


/*! \name Access methods. */
/**@{*/
/*! Get entry \f$ x_i \f$. */
real
getentry_realvector(pcrealvector x, index i);

/*! Set entry \f$ x_i \f$. */
void
setentry_realvector(prealvector x, index i, real entry);

/*! Add to entry \f$ x_i \f$. */
void
addentry_realvector(prealvector x, index i, real entry);

/*! Print vector.*/
void
print_realvector(pcrealvector x);
/**@}*/


/*! \name Basic operations. */
/**@{*/
/*! Scale \f$x \leftarrow \alpha\cdot x\f$. */
void
scal_realvector(real alpha, prealvector x);

/*! Perform \f$ x^Ty \f$.*/
real
dot_realvector(pcrealvector x, pcrealvector y);

/*! Perform \f$ \|x\|_2 \f$. */
real
nrm2_realvector(pcrealvector x);

/*! Perform \f$ y \leftarrow \alpha\cdot x + y \f$.*/
void
axpy_realvector(real alpha, pcrealvector x, prealvector y);
/**@}*/

#endif
