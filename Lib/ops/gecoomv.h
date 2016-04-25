/**
 *  General matrix vector product with COO sparse matrices.
 *  Computes y <- alpha*Op(A)*x + beta*y,
 *  where Op(A) = A or Op(A) = A^T.
 */

#ifndef GECOOMV_H
#define GECOOMV_H

#include "../settings.h"
#include "../storage/realvector.h"
#include "../storage/coo.h"

void
gecoomv(transpose t,
        real alpha,
        pccoo A,
        pcrealvector x,
        real beta,
        prealvector y);

#endif
