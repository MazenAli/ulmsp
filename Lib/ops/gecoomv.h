/*! \file gecoomv.h
 *  General matrix vector product with COO sparse matrices.
*/

#ifndef GECOOMV_H
#define GECOOMV_H

#include "../settings.h"
#include "../storage/realvector.h"
#include "../storage/coo.h"

/*!
 *  Computes \f$ y \leftarrow \alpha Op(A)x + \beta y \f$,
 *  where \f$ Op(A) = A \f$ or \f$ Op(A) = A^T \f$.
 */
void
gecoomv(transpose t,
        real alpha,
        pccoo A,
        pcrealvector x,
        real beta,
        prealvector y);

#endif
