/*! \file cgcrs.h
 *  Conjugate Gradient Method (CG) for a CRS matrix.
 */

#ifndef CGCRS_H
#define CGCRS_H

#include "../settings.h"
#include "../storage/realvector.h"
#include "../storage/crs.h"

/*!
 * Solve \f$ Ax=b \f$ with \f$ A \f$ in CRS format.
 * Returns number of iterations.
 */
index
cgcrs(pccrs A,
      prealvector x,
      pcrealvector b,
      real tol,
      index maxit);

#endif
