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
    real    *vals;
    index    length;
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

#endif
