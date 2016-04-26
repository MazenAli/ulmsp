/**
 *  General matrix vector product with CRS sparse matrices.
 *  Computes y <- alpha*Op(A)*x + beta*y,
 *  where Op(A) = A or Op(A) = A^T.
 */

#ifndef GECRSMV_H
#define GECRSMV_H

#include "../settings.h"
#include "../storage/realvector.h"
#include "../storage/crs.h"

void
gecrsmv(transpose t,
        real alpha,
        pccrs A,
        pcrealvector x,
        real beta,
        prealvector y);

#endif
