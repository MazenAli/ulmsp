/*! \file precs.h
 *  Preconditioners for iterative methods.
 *  Assume diagonal entries are sorted first.
 */

#ifndef PRECS_H
#define PRECS_H

#include "../settings.h"
#include "../storage/realvector.h"
#include "../storage/crs.h"

/*!
 * Diagonal (Jacobi) preconditioner.
 * Performs fixed number of iterations on \f$ Ax=b \f$
 */
void
diagcrs(pccrs A,
        prealvector x,
        pcrealvector b,);

#endif
