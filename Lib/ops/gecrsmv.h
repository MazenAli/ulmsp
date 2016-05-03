/*! \file gecrsmv.h
 *  General matrix vector product with CRS sparse matrices.
*/

#ifndef GECRSMV_H
#define GECRSMV_H

#include "../settings.h"
#include "../storage/realvector.h"
#include "../storage/crs.h"

/*!
 *  Computes \f$ y \leftarrow \alpha Op(A)x + \beta y \f$,
 *  where \f$ Op(A) = A \f$ or \f$ Op(A) = A^T \f$.
 */
void
gecrsmv(transpose t,
        real alpha,
        pccrs A,
        pcrealvector x,
        real beta,
        prealvector y);

#endif
