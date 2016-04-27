/**
 *  A simple real vector data type.
 */

#ifndef REALVECTOR_H
#define REALVECTOR_H

#include "../settings.h"

typedef struct _realvector realvector;
typedef realvector         *prealvector;
typedef const realvector   *pcrealvector;

struct _realvector
{
    real    *vals;   ///< Pointer to vector entries
    index    length; ///< Length of vector
};

/**
 *  Constructors and destructors.
 */

prealvector
new_realvector(index length);

void
init_realvector(prealvector x, index length);

void
resize_realvector(prealvector x, index length);

void
copy_realvector(prealvector dest, pcrealvector src);

void
del_realvector(prealvector x);

/**
 *  Access methods.
 */

real
getentry_realvector(pcrealvector x, index i);

void
setentry_realvector(prealvector x, index i, real entry);

void
addentry_realvector(prealvector x, index i, real entry);

void
print_realvector(pcrealvector x);

/**
 *  Basic operations.
 */

/// Scale x <- alpha*x
void
scal_realvector(real alpha, prealvector x);

/// Perform x^T*y
real
dot_realvector(pcrealvector x, pcrealvector y);

/// Perform |x|_2
real
nrm2_realvector(pcrealvector x);

/// Perform y <- alpha*x + y
void
axpy_realvector(real alpha, pcrealvector x, prealvector y);

#endif
